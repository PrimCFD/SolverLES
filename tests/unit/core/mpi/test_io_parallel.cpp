#include "master/FieldCatalog.hpp"
#include "master/Views.hpp"
#include "master/io/IWriter.hpp"
#include "master/io/WriterConfig.hpp"
#include <catch2/catch_test_macros.hpp>

#include "master/io/CGNSWriter.hpp"
#include "master/io/XdmfHdf5Writer.hpp"

#include <filesystem>
#include <string>
#include <vector>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

using namespace core::master;
using namespace core::master::io;
namespace fs = std::filesystem;

static inline std::array<std::ptrdiff_t, 3> strides_bytes(int nx, int ny, std::size_t elem)
{
    return {(std::ptrdiff_t) elem, (std::ptrdiff_t)(elem * nx), (std::ptrdiff_t)(elem * nx * ny)};
}

#ifndef HAVE_MPI

TEST_CASE("MPI I/O smoke skipped (HAVE_MPI=0)", "[io][mpi]")
{
    SUCCEED();
}

#else // HAVE_MPI

TEST_CASE("MPI I/O: gather on root then write once", "[io][mpi][2rank]")
{
    int rank = 0, size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Small local block per rank
    const int nx = 6, ny = 4, nz = 2;
    const std::size_t n_local = (std::size_t) nx * ny * nz;

    std::vector<double> local(n_local, double(rank));

    // Gather to root (rank 0). Use nullptr recvbuf on non-root to avoid impl edge-cases.
    std::vector<double> global; // only allocated on root
    void* recvbuf = nullptr;
    if (rank == 0)
    {
        global.resize(n_local * (std::size_t) size);
        recvbuf = static_cast<void*>(global.data());
    }

    int gerr = MPI_Gather(local.data(), (int) local.size(), MPI_DOUBLE, recvbuf, (int) local.size(),
                          MPI_DOUBLE, 0, MPI_COMM_WORLD);
    REQUIRE(gerr == MPI_SUCCESS);

    // Root writes; broadcast status so others don't hang if I/O fails.
    int io_status = 1; // 1=success, 0=failure, 2=skipped (no backend)

    if (rank == 0)
    {
        try
        {
            // Build a catalog with one big, gathered field (stack along X).
            FieldCatalog fc;
            const int NX = nx * size;
            const std::size_t elem = sizeof(double);

            fc.register_scalar("rankval", global.data(), elem, {NX, ny, nz},
                               strides_bytes(NX, ny, elem));
            fc.select_for_output("rankval");

            // Prepare writer config
            WriterConfig cfg;
            fs::path outdir = fs::temp_directory_path() / "solverles_mpi_io";
            fs::create_directories(outdir);
            cfg.path = outdir.string();

            // XDMF
            std::unique_ptr<IWriter> W;
            cfg.backend = WriterConfig::Backend::XDMF;
            cfg.precision = WriterConfig::Precision::Float32; // realistic cast path
            W = std::make_unique<XdmfHdf5Writer>(cfg);
            const fs::path xmf = outdir / "mpi_io.xmf";
            const fs::path h5 = outdir / "mpi_io.h5";

            W->open_case("mpi_io");

            WriteRequest req;
            req.step = 0;
            req.time = 0.0;
            auto sel = fc.selected_for_output();
            req.fields.assign(sel.begin(), sel.end());

            W->write(req);
            W->close();

            REQUIRE(fs::exists(outdir / "mpi_io.xmf"));
            REQUIRE(fs::exists(outdir / "mpi_io.h5"));
        }
        catch (...)
        {
            io_status = 0;
        }
    }

    // Tell everyone whether root succeeded (or skipped).
    MPI_Bcast(&io_status, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (io_status == 2)
    {
        // Skipped due to no backend; still keep ranks in sync.
        MPI_Barrier(MPI_COMM_WORLD);
        SUCCEED();
        return;
    }

    REQUIRE(io_status == 1);

    // Keep output tidy / synchronize exit
    MPI_Barrier(MPI_COMM_WORLD);
    SUCCEED();
}

#endif // HAVE_MPI