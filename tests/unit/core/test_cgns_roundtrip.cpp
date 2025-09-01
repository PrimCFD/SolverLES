#include "master/FieldCatalog.hpp"
#include "master/Views.hpp"
#include "master/io/CGNSWriter.hpp"
#include "master/io/WriterConfig.hpp"
#include <catch2/catch_test_macros.hpp>

#include <cstring>
#include <filesystem>
#include <random>
#include <vector>

#ifdef SOLVERLES_WITH_CGNS
#include <cgnslib.h>
#endif

using namespace core::master;
using namespace core::master::io;
namespace fs = std::filesystem;

#ifdef SOLVERLES_WITH_CGNS
static std::vector<double> cgns_read_field(const std::string& path, const std::string& solname,
                                           const std::string& fieldname, int& nx, int& ny, int& nz)
{
    int file = -1;
    if (cg_open(path.c_str(), CG_MODE_READ, &file))
    {
        FAIL("cg_open failed");
    }
    int nbases = 0;
    cg_nbases(file, &nbases);
    REQUIRE(nbases >= 1);
    int base = 1;
    int cell_dim = 0, phys_dim = 0;
    cg_base_read(file, base, nullptr, &cell_dim, &phys_dim);

    int nzones = 0;
    cg_nzones(file, base, &nzones);
    REQUIRE(nzones >= 1);
    int zone = 1;
    char zname[33] = {0};
    cgsize_t size[9] = {0};
    cg_zone_read(file, base, zone, zname, size);
    // size[0..2] are vertex IJK, [3..5] are cell IJK
    nx = int(size[3]);
    ny = int(size[4]);
    nz = int(size[5]);

    // find solution by name
    int nsol = 0;
    cg_nsols(file, base, zone, &nsol);
    REQUIRE(nsol >= 1);
    int sol_id = -1;
    for (int s = 1; s <= nsol; ++s)
    {
        char sname[33] = {0};
        GridLocation_t loc{};
        cg_sol_info(file, base, zone, s, sname, &loc);
        if (solname == sname)
        {
            sol_id = s;
            break;
        }
    }
    REQUIRE(sol_id != -1);

    // read field (assume RealSingle or RealDouble)
    int nf = 0;
    cg_nfields(file, base, zone, sol_id, &nf);
    REQUIRE(nf >= 1);

    CG_DataType_t dtype = RealDouble;
    int field_id = -1;
    for (int f = 1; f <= nf; ++f)
    {
        char fname[33] = {0};
        CG_DataType_t dt{};
        cg_field_info(file, base, zone, sol_id, f, &dt, fname);
        if (fieldname == fname)
        {
            field_id = f;
            dtype = dt;
            break;
        }
    }
    REQUIRE(field_id != -1);

    std::vector<double> out(std::size_t(nx) * ny * nz);
    if (dtype == RealDouble)
    {
        if (cg_field_read(file, base, zone, sol_id, field_id, RealDouble, out.data()))
            FAIL("cg_field_read(double) failed");
    }
    else if (dtype == RealSingle)
    {
        std::vector<float> tmp(out.size());
        if (cg_field_read(file, base, zone, sol_id, field_id, RealSingle, tmp.data()))
            FAIL("cg_field_read(float) failed");
        for (std::size_t i = 0; i < out.size(); ++i)
            out[i] = tmp[i];
    }
    else
    {
        FAIL("Unsupported CGNS dtype");
    }

    cg_close(file);
    return out;
}
#endif // SOLVERLES_WITH_CGNS

TEST_CASE("CGNS roundtrip: one scalar field, 2 timesteps", "[io][cgns]")
{
#ifndef SOLVERLES_WITH_CGNS
    SUCCEED("CGNS disabled at build time");
    return;
#else
    // Arrange: small synthetic field
    const int nx = 8, ny = 5, nz = 3;
    std::vector<double> rho(std::size_t(nx) * ny * nz);
    auto idx = [=](int i, int j, int k)
    { return std::size_t(k) * ny * nx + std::size_t(j) * nx + i; };
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
                rho[idx(i, j, k)] = 1.0 * i + 10.0 * j + 100.0 * k;

    FieldCatalog fc;
    fc.register_scalar("rho", rho.data(), sizeof(double), {nx, ny, nz},
                       {(std::ptrdiff_t) sizeof(double), (std::ptrdiff_t)(sizeof(double) * nx),
                        (std::ptrdiff_t)(sizeof(double) * nx * ny)});
    fc.select_for_output("rho");

    WriterConfig cfg;
    cfg.backend = WriterConfig::Backend::CGNS;
    cfg.series = WriterConfig::Series::Single;
    cfg.precision = WriterConfig::Precision::Float64;

    fs::path outdir = fs::temp_directory_path() / "solverles_cgns_rt";
    fs::create_directories(outdir);
    cfg.path = outdir.string();

    CGNSWriter W(cfg);

    W.open_case("rt_case");

    // Step 0
    WriteRequest r0;
    r0.step = 0;
    r0.time = 0.0;
    auto sel = fc.selected_for_output();
    r0.fields.assign(sel.begin(), sel.end());
    REQUIRE(W.write(r0));

    // mutate data and write step 1
    for (auto& v : rho)
        v += 0.5; // simple change
    WriteRequest r1;
    r1.step = 1;
    r1.time = 0.1;
    r1.fields.assign(sel.begin(), sel.end());
    REQUIRE(W.write(r1));

    W.close();

    // Assert: read back both steps
    const std::string cgns = (outdir / "rt_case.cgns").string();

    int rx = 0, ry = 0, rz = 0;
    auto f0 = cgns_read_field(cgns, "FlowSolution_000000", "rho", rx, ry, rz);
    REQUIRE(rx == nx);
    REQUIRE(ry == ny);
    REQUIRE(rz == nz);
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
                REQUIRE(f0[idx(i, j, k)] == Approx(1.0 * i + 10.0 * j + 100.0 * k));

    auto f1 = cgns_read_field(cgns, "FlowSolution_000001", "rho", rx, ry, rz);
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
                REQUIRE(f1[idx(i, j, k)] == Approx(1.0 * i + 10.0 * j + 100.0 * k + 0.5));
#endif
}
