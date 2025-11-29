// tests/performance/bench_io_write.cpp
#include "master/FieldCatalog.hpp"
#include "master/Views.hpp"
#include "master/io/CGNSWriter.hpp"
#include "master/io/WriterConfig.hpp"
#include "master/io/XdmfHdf5Writer.hpp"
#include "simple_bench.hpp"

#include "memory/MpiBox.hpp"
#include "mesh/Mesh.hpp"

#include <chrono>
#include <cstring>
#include <filesystem>
#include <iostream>
#include <mpi.h>
#include <numeric>
#include <string>
#include <vector>

using namespace core::master;
using namespace core::master::io;
namespace fs = std::filesystem;

struct Args
{
    std::string backend = "xmf"; // xmf | cgns
    int nx = 64, ny = 64, nz = 32;
    int reps = 5;
    bool f32 = true;
};

static Args parse_args(int argc, char** argv)
{
    Args a;
    for (int i = 1; i < argc; ++i)
    {
        if (!std::strncmp(argv[i], "--backend=", 10))
            a.backend = argv[i] + 10;
        else if (!std::strncmp(argv[i], "--nx=", 5))
            a.nx = std::atoi(argv[i] + 5);
        else if (!std::strncmp(argv[i], "--ny=", 5))
            a.ny = std::atoi(argv[i] + 5);
        else if (!std::strncmp(argv[i], "--nz=", 5))
            a.nz = std::atoi(argv[i] + 5);
        else if (!std::strncmp(argv[i], "--reps=", 7))
            a.reps = std::atoi(argv[i] + 7);
        else if (!std::strcmp(argv[i], "--float32"))
            a.f32 = true;
        else if (!std::strcmp(argv[i], "--float64"))
            a.f32 = false;
    }
    return a;
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int world_rank = 0, world_size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    const Args args = parse_args(argc, argv);

    // Cartesian grid dims (fill zeros)
    int dims[3] = {0, 0, 0};
    int periods[3] = {0, 0, 0};
    MPI_Dims_create(world_size, 3, dims);
    MPI_Comm cart_comm;
    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, /*reorder*/ 1, &cart_comm);

    int coords[3] = {0, 0, 0};
    if (cart_comm != MPI_COMM_NULL)
    {
        MPI_Cart_coords(cart_comm, world_rank, 3, coords);
    }

    // Partition global cells as evenly as possible.
    auto split = [](int N, int p, int r) -> std::pair<int, int>
    {
        // return (local, offset)
        const int base = N / p;
        const int rem = N % p;
        const int loc = base + (r < rem ? 1 : 0);
        const int off = r * base + std::min(r, rem);
        return {loc, off};
    };

    auto [lx, ox] = split(args.nx, dims[0] ? dims[0] : 1, dims[0] ? coords[0] : 0);
    auto [ly, oy] = split(args.ny, dims[1] ? dims[1] : 1, dims[1] ? coords[1] : 0);
    auto [lz, oz] = split(args.nz, dims[2] ? dims[2] : 1, dims[2] ? coords[2] : 0);

    // Mesh (cell-centered only here)
    core::mesh::Mesh mesh{};
    mesh.local = {lx, ly, lz};
    mesh.global = {args.nx, args.ny, args.nz};
    mesh.global_lo = {ox, oy, oz};
    mesh.proc_grid = {dims[0], dims[1], dims[2]};
    mesh.ng = 0;
    mesh.periodic = {false, false, false};

    // Create synthetic tile-local fields (cell-centered)
    FieldCatalog fc;
    const int nfields = 8; // same payload for both writers
    std::vector<std::vector<double>> stor(nfields);
    for (int f = 0; f < nfields; ++f)
    {
        auto& v = stor[f];
        v.resize(std::size_t(lx) * ly * lz);
        for (std::size_t i = 0; i < v.size(); ++i)
            v[i] = double(f) + 0.001 * double(i + ox + oy * args.nx + oz * args.nx * args.ny);

        fc.register_scalar("F" + std::to_string(f), v.data(), sizeof(double), {lx, ly, lz},
                           {(std::ptrdiff_t) sizeof(double), (std::ptrdiff_t)(sizeof(double) * lx),
                            (std::ptrdiff_t)(sizeof(double) * lx * ly)},
                           Stagger::Cell);
        fc.select_for_output("F" + std::to_string(f));
    }

    // Writer config (shared)
    WriterConfig cfg;
    cfg.path = (fs::temp_directory_path() / "kolmoplas_bench_io").string();
    fs::create_directories(cfg.path);
    cfg.precision = args.f32 ? WriterConfig::Precision::Float32 : WriterConfig::Precision::Float64;

    // pass mesh + cart comm so writers can do parallel I/O / slicing
    cfg.mesh = &mesh;
    cfg.mpi_cart_comm = mpi_box(cart_comm);

    // Backend selection
    std::unique_ptr<IWriter> writer;
    if (args.backend == "xmf" || args.backend == "xdmf")
    {
        cfg.backend = WriterConfig::Backend::XDMF;
        writer = std::make_unique<XdmfHdf5Writer>(cfg);
    }
    else if (args.backend == "cgns")
    {
        cfg.backend = WriterConfig::Backend::CGNS;
        writer = std::make_unique<CGNSWriter>(cfg);
    }
    else
    {
        if (world_rank == 0)
            std::cerr << "Unknown --backend=" << args.backend << " (use xmf|cgns)\n";
        MPI_Finalize();
        return 2;
    }

    const std::string casename = "bench_io_" + args.backend;

    // Open once
    writer->open_case(casename);

    // Common write request
    auto sel = fc.selected_for_output();
    WriteRequest r;
    r.fields.assign(sel.begin(), sel.end());

    // N-step timing using bench::run_mpi_max (µs for MAX over ranks)
    const std::size_t step_bytes =
        std::size_t(lx) * ly * lz * nfields * (args.f32 ? sizeof(float) : sizeof(double));

    int step = 0;
    auto write_one = [&]()
    {
        r.step = step;
        r.time = double(step);
        writer->write(r);
        ++step;
    };
    auto [mean_us, stddev_us] = bench::run_mpi_max(cart_comm, write_one, args.reps);

    writer->close();

    if (world_rank == 0)
    {
        const std::string tag = "io_" + args.backend + "_" + std::to_string(args.nx) + "x" +
                                std::to_string(args.ny) + "x" + std::to_string(args.nz) + "_p" +
                                std::to_string(world_size);
        bench::report(tag, mean_us, stddev_us, step_bytes);
    }

    if (cart_comm != MPI_COMM_NULL)
        MPI_Comm_free(&cart_comm);
    MPI_Finalize();
    return 0;
}
