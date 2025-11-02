#include "master/FieldCatalog.hpp"
#include "master/Views.hpp"
#include "master/io/CGNSWriter.hpp"
#include "master/io/WriterConfig.hpp"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cstring>
#include <filesystem>
#include <random>
#include <vector>

#include <cgnslib.h>

using std::string;
using std::vector;

// ----- CGNS error reporter (portable w.r.t. your header)
extern "C" void test_cgns_err(int is_error, char* msg)
{
    // Print even non-error messages to help trace the last CGNS call
    fprintf(stderr, "\n[CGNS %s] %s\n", is_error ? "error" : "note", msg ? msg : "(null)");
    fflush(stderr);
}

static struct CgnsErrHook
{
    CgnsErrHook() { cg_error_handler(test_cgns_err); } // <-- call the function, don't assign
} _ceh_;

// ----- Backtrace on segfault
#include <cstdio>
#include <cstring>
#include <execinfo.h>
#include <signal.h>
#include <unistd.h>

static void segv_handler(int sig, siginfo_t*, void*)
{
    void* frames[64];
    int n = backtrace(frames, 64);
    fprintf(stderr, "\n*** Caught SIGSEGV (%d). Backtrace:\n", sig);
    backtrace_symbols_fd(frames, n, STDERR_FILENO);
    // flush & exit with signal-style code
    fsync(STDERR_FILENO);
    _Exit(128 + sig);
}
static void install_segv_bt()
{
    struct sigaction sa;
    std::memset(&sa, 0, sizeof(sa));
    sa.sa_sigaction = segv_handler;
    sa.sa_flags = SA_SIGINFO | SA_RESETHAND;
    sigemptyset(&sa.sa_mask);
    sigaction(SIGSEGV, &sa, nullptr);
}
static struct SegvBT
{
    SegvBT() { install_segv_bt(); }
} _segvbt_;

using namespace core::master;
using namespace core::master::io;
namespace fs = std::filesystem;
using Catch::Approx;

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
    char basename[33] = {0};
    REQUIRE(cg_base_read(file, base, basename, &cell_dim, &phys_dim) == 0);

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

    CGNS_ENUMT(DataType_t) dtype = CGNS_ENUMV(RealDouble);
    int field_id = -1;
    char field_name[33] = {0};
    for (int f = 1; f <= nf; ++f)
    {
        char fname[33] = {0};
        CGNS_ENUMT(DataType_t) dt{};
        cg_field_info(file, base, zone, sol_id, f, &dt, fname);
        if (fieldname == fname)
        {
            field_id = f;
            dtype = dt;
            std::strncpy(field_name, fname, 32);
            break;
        }
    }
    REQUIRE(field_id != -1);

    std::vector<double> out(std::size_t(nx) * ny * nz);
    // rmin/rmax for full-range read (cell-centered sizes already in nx,ny,nz)
    cgsize_t rmin[3] = {1, 1, 1};
    cgsize_t rmax[3] = {(cgsize_t) nx, (cgsize_t) ny, (cgsize_t) nz};

    if (dtype == CGNS_ENUMV(RealDouble))
    {
        if (cg_field_read(file, base, zone, sol_id, field_name, CGNS_ENUMV(RealDouble), rmin, rmax,
                          out.data()))
            FAIL("cg_field_read(double) failed");
    }
    else if (dtype == CGNS_ENUMV(RealSingle))
    {
        std::vector<float> tmp(out.size());
        if (cg_field_read(file, base, zone, sol_id, field_name, CGNS_ENUMV(RealSingle), rmin, rmax,
                          tmp.data()))
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

static std::vector<std::string> cgns_list_flow_solutions(const std::string& path)
{
    int fn = -1;
    if (cg_open(path.c_str(), CG_MODE_READ, &fn))
        FAIL("cg_open failed");
    int base = 1, zone = 1;
    // If bases/zones aren’t 1, this is a unit test—keep it simple:
    {
        int nb = 0;
        cg_nbases(fn, &nb);
        REQUIRE(nb >= 1);
    }
    {
        int nz = 0;
        cg_nzones(fn, base, &nz);
        REQUIRE(nz >= 1);
    }

    std::vector<std::string> names;
    // Prefer ZoneIterativeData FlowSolutionPointers (writer fills this on close()).
    if (!cg_goto(fn, base, "Zone_t", zone, "ZoneIterativeData_t", 1, "end"))
    {
        int narr = 0;
        cg_narrays(&narr);
        for (int a = 1; a <= narr; ++a)
        {
            char aname[33] = {0};
            DataType_t adt{};
            int ndim = 0;
            cgsize_t dims[12] = {0};
            cg_array_info(a, aname, &adt, &ndim, dims);
            if (std::string(aname) == "FlowSolutionPointers" && ndim == 2 && adt == Character)
            {
                const int WIDTH = int(dims[0]);
                const int N = int(dims[1]);
                std::vector<char> buf(std::size_t(WIDTH) * N, ' ');
                cg_array_read(a, buf.data());
                names.resize(N);
                for (int s = 0; s < N; ++s)
                {
                    std::string nm(&buf[std::size_t(s) * WIDTH], std::size_t(WIDTH));
                    // rtrim spaces (Fortran-style fixed width)
                    nm.erase(nm.find_last_not_of(' ') + 1);
                    names[s] = nm;
                }
                break;
            }
        }
    }
    // Fallback: enumerate by index if pointers are missing
    if (names.empty())
    {
        int nsol = 0;
        cg_nsols(fn, base, zone, &nsol);
        REQUIRE(nsol >= 1);
        names.resize(nsol);
        for (int s = 1; s <= nsol; ++s)
        {
            char nm[33] = {0};
            GridLocation_t loc{};
            cg_sol_info(fn, base, zone, s, nm, &loc);
            names[s - 1] = nm;
        }
    }
    cg_close(fn);
    return names;
}

TEST_CASE("CGNS roundtrip: one scalar field, 2 timesteps", "[io][cgns]")
{
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
                        (std::ptrdiff_t)(sizeof(double) * nx * ny)},
                       core::master::Stagger::Cell);
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
    W.write(r0);

    // mutate data and write step 1
    for (auto& v : rho)
        v += 0.5; // simple change
    WriteRequest r1;
    r1.step = 1;
    r1.time = 0.1;
    r1.fields.assign(sel.begin(), sel.end());
    W.write(r1);

    W.close();

    // Assert: read back both steps
    const std::string cgns = (outdir / "rt_case.cgns").string();
    auto sols = cgns_list_flow_solutions(cgns);
    REQUIRE(sols.size() >= 2);

    int rx = 0, ry = 0, rz = 0;
    auto f0 = cgns_read_field(cgns, sols[0], "rho", rx, ry, rz);
    REQUIRE(rx == nx);
    REQUIRE(ry == ny);
    REQUIRE(rz == nz);
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
                REQUIRE(f0[idx(i, j, k)] == Approx(1.0 * i + 10.0 * j + 100.0 * k));

    auto f1 = cgns_read_field(cgns, sols[1], "rho", rx, ry, rz);
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
                REQUIRE(f1[idx(i, j, k)] == Approx(1.0 * i + 10.0 * j + 100.0 * k + 0.5));
}
