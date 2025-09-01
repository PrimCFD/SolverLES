#include "master/FieldCatalog.hpp"
#include "master/Views.hpp"
#include "master/io/WriterConfig.hpp"
#include "master/io/XdmfHdf5Writer.hpp"
#include <catch2/catch_test_macros.hpp>

#ifdef HAVE_CUDA
#include <cuda_runtime.h>
#endif

#include <filesystem>
#include <vector>

using namespace core::master;
using namespace core::master::io;
namespace fs = std::filesystem;

TEST_CASE("IO works under CUDA build (UM or mirrored)", "[io][gpu]")
{
#ifndef HAVE_CUDA
    SUCCEED("CUDA not enabled; test skipped");
    return;
#else
    const int nx = 8, ny = 4, nz = 2;
    std::vector<double> f(std::size_t(nx) * ny * nz, 3.14);

    FieldCatalog fc;
    fc.register_scalar("f", f.data(), sizeof(double), {nx, ny, nz},
                       {(std::ptrdiff_t) sizeof(double), (std::ptrdiff_t)(sizeof(double) * nx),
                        (std::ptrdiff_t)(sizeof(double) * nx * ny)});
    fc.select_for_output("f");

    WriterConfig cfg;
    cfg.backend = WriterConfig::Backend::XDMF; // simple HDF5 write path
    fs::path outdir = fs::temp_directory_path() / "solverles_gpu_io";
    fs::create_directories(outdir);
    cfg.path = outdir.string();

    XdmfHdf5Writer W(cfg);
    REQUIRE(W.open_case("gpu_io"));

    WriteRequest req;
    req.step = 0;
    req.time = 0.0;
    auto sel = fc.selected_for_output();
    req.fields.assign(sel.begin(), sel.end());
    REQUIRE(W.write(req));
    W.close();

    SUCCEED();
#endif
}