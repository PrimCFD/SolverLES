#include "master/FieldCatalog.hpp"
#include "master/Views.hpp"
#include "master/io/WriterConfig.hpp"
#include "master/io/XdmfHdf5Writer.hpp"
#include "simple_bench.hpp"

#include <filesystem>
#include <iostream>
#include <random>
#include <vector>

using namespace core::master;
using namespace core::master::io;
namespace fs = std::filesystem;

int main(int argc, char**)
{
    const int nx = 64, ny = 64, nz = 32;
    const int nfields = 8; // many fields path

    FieldCatalog fc;
    std::vector<std::vector<double>> stor(nfields);
    for (int f = 0; f < nfields; ++f)
    {
        auto& v = stor[f];
        v.resize(std::size_t(nx) * ny * nz);
        for (std::size_t i = 0; i < v.size(); ++i)
            v[i] = double(f) + 0.001 * double(i);
        fc.register_scalar("F" + std::to_string(f), v.data(), sizeof(double), {nx, ny, nz},
                           {(std::ptrdiff_t) sizeof(double), (std::ptrdiff_t)(sizeof(double) * nx),
                            (std::ptrdiff_t)(sizeof(double) * nx * ny)});
        fc.select_for_output("F" + std::to_string(f));
    }

    WriterConfig cfg;
    cfg.backend = WriterConfig::Backend::XDMF;
    cfg.precision = WriterConfig::Precision::Float32; // realistic
    fs::path outdir = fs::temp_directory_path() / "solverles_bench_io";
    fs::create_directories(outdir);
    cfg.path = outdir.string();

    XdmfHdf5Writer W(cfg);
    W.open_case("bench_io");

    auto sel = fc.selected_for_output();
    WriteRequest r;
    r.fields.assign(sel.begin(), sel.end());
    r.step = 0;
    r.time = 0.0;

    auto [mean, stddev] = bench::run([&] { (void) W.write(r); });
    const std::size_t step_bytes = std::size_t(nx) * ny * nz * nfields * sizeof(float);
    bench::report("io_write_" + std::to_string(nx) + "x" + std::to_string(ny) + "x" +
                      std::to_string(nz),
                  mean, stddev, step_bytes);

    W.close();
    return 0;
}