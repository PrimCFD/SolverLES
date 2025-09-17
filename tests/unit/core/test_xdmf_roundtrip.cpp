#include "master/FieldCatalog.hpp"
#include "master/Views.hpp"
#include "master/io/WriterConfig.hpp"
#include "master/io/XdmfHdf5Writer.hpp"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <hdf5.h>

#include <filesystem>
#include <string>
#include <vector>

using namespace core::master;
using namespace core::master::io;
namespace fs = std::filesystem;
using Catch::Approx;

static std::vector<double> h5_read_3d(const std::string& file, const std::string& dset, int& nx,
                                      int& ny, int& nz, bool want_double)
{
    hid_t f = H5Fopen(file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    REQUIRE(f >= 0);
    hid_t d = H5Dopen2(f, dset.c_str(), H5P_DEFAULT);
    REQUIRE(d >= 0);
    hid_t sp = H5Dget_space(d);
    REQUIRE(sp >= 0);
    int nd = H5Sget_simple_extent_ndims(sp);
    REQUIRE(nd == 3);
    hsize_t dims[3];
    H5Sget_simple_extent_dims(sp, dims, nullptr);
    nz = (int) dims[0];
    ny = (int) dims[1];
    nx = (int) dims[2];

    hid_t t = H5Dget_type(d);
    bool is_f64 = H5Tequal(t, H5T_IEEE_F64LE) > 0;
    bool is_f32 = H5Tequal(t, H5T_IEEE_F32LE) > 0;
    const bool ok = (is_f64 || is_f32);
    REQUIRE(ok);

    std::vector<double> out(std::size_t(nx) * ny * nz);
    if (is_f64)
    {
        if (want_double)
        {
            REQUIRE(H5Dread(d, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, out.data()) >= 0);
        }
        else
        {
            std::vector<double> tmp(out.size());
            REQUIRE(H5Dread(d, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp.data()) >= 0);
            for (std::size_t i = 0; i < out.size(); ++i)
                out[i] = tmp[i];
        }
    }
    else
    {
        std::vector<float> tmp(out.size());
        REQUIRE(H5Dread(d, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp.data()) >= 0);
        for (std::size_t i = 0; i < out.size(); ++i)
            out[i] = tmp[i];
    }

    H5Tclose(t);
    H5Sclose(sp);
    H5Dclose(d);
    H5Fclose(f);
    return out;
}

TEST_CASE("XDMF/HDF5 roundtrip: two scalar fields, 1 timestep", "[io][xdmf]")
{
    const int nx = 6, ny = 4, nz = 2;
    std::vector<double> a(std::size_t(nx) * ny * nz), b(a.size());
    auto idx = [=](int i, int j, int k)
    { return std::size_t(k) * ny * nx + std::size_t(j) * nx + i; };
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
            {
                a[idx(i, j, k)] = i + 2 * j + 3 * k;
                b[idx(i, j, k)] = 100 + i - j + 0.25 * k;
            }

    FieldCatalog fc;
    fc.register_scalar("A", a.data(), sizeof(double), {nx, ny, nz},
                       {(std::ptrdiff_t) sizeof(double), (std::ptrdiff_t)(sizeof(double) * nx),
                        (std::ptrdiff_t)(sizeof(double) * nx * ny)});
    fc.register_scalar("B", b.data(), sizeof(double), {nx, ny, nz},
                       {(std::ptrdiff_t) sizeof(double), (std::ptrdiff_t)(sizeof(double) * nx),
                        (std::ptrdiff_t)(sizeof(double) * nx * ny)});
    fc.select_for_output("A");
    fc.select_for_output("B");

    WriterConfig cfg;
    cfg.backend = WriterConfig::Backend::XDMF;
    cfg.series = WriterConfig::Series::Single;
    cfg.precision = WriterConfig::Precision::Float32; // cast on write to shrink

    fs::path outdir = fs::temp_directory_path() / "solverles_xdmf_rt";
    fs::create_directories(outdir);
    cfg.path = outdir.string();

    XdmfHdf5Writer W(cfg);
    W.open_case("xdmf_rt");

    WriteRequest r;
    r.step = 0;
    r.time = 0.0;
    auto sel = fc.selected_for_output();
    r.fields.assign(sel.begin(), sel.end());
    W.write(r);
    W.close();

    // Read back via HDF5
    std::string h5 = (outdir / "xdmf_rt.h5").string();
    int rx = 0, ry = 0, rz = 0;
    auto A = h5_read_3d(h5, "/Step_000000/A", rx, ry, rz, /*want_double=*/false);
    REQUIRE(rx == nx);
    REQUIRE(ry == ny);
    REQUIRE(rz == nz);
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
                REQUIRE(A[idx(i, j, k)] == Approx(i + 2 * j + 3 * k).margin(1e-6));

    auto B = h5_read_3d(h5, "/Step_000000/B", rx, ry, rz, /*want_double=*/false);
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
                REQUIRE(B[idx(i, j, k)] == Approx(100 + i - j + 0.25 * k).margin(1e-6));
}
