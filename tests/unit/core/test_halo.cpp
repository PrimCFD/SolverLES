#include <catch2/catch_test_macros.hpp>
#include <mpi.h>

#include "Field.hpp"
#include "HaloExchange.hpp"
#include "MemoryManager.hpp"
#include "Mesh.hpp"

using core::Field;
using core::Mesh;

TEST_CASE("Halo swap conserves known pattern (3D cart)", "[mpi]")
{
    int rank = 0, nprocs = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    // Build the same style of Cartesian topology the halo code uses.
    int dims[3] = {0, 0, 0};
    int periods[3] = {0, 0, 0}; // non-periodic by default
    MPI_Dims_create(nprocs, 3, dims);
    MPI_Comm cart = MPI_COMM_NULL;
    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, /*reorder=*/1, &cart);
    if (cart == MPI_COMM_NULL)
        cart = MPI_COMM_WORLD;

    // Neighbor ranks for each axis (neg,pos)
    int xneg, xpos, yneg, ypos, zneg, zpos;
    MPI_Cart_shift(cart, 0, 1, &xneg, &xpos);
    MPI_Cart_shift(cart, 1, 1, &yneg, &ypos);
    MPI_Cart_shift(cart, 2, 1, &zneg, &zpos);

    const int nx = 8, ny = 6, nz = 4, ng = 1;
    Mesh mesh{{nx, ny, nz}, ng};

    // Allocate raw storage (including ghosts)
    auto& mm = core::MemoryManager::instance();
    const std::array<int, 3> ext_with_ghosts = {nx + 2 * ng, ny + 2 * ng, nz + 2 * ng};
    double* raw = mm.allocate<double>(static_cast<std::size_t>(ext_with_ghosts[0]) *
                                      static_cast<std::size_t>(ext_with_ghosts[1]) *
                                      static_cast<std::size_t>(ext_with_ghosts[2]));
    Field<double> f(raw, ext_with_ghosts, ng);

    // Initialize: interior = rank, ghosts = 0
    for (int k = -ng; k < nz + ng; ++k)
        for (int j = -ng; j < ny + ng; ++j)
            for (int i = -ng; i < nx + ng; ++i)
                f(i, j, k) =
                    (i >= 0 && i < nx && j >= 0 && j < ny && k >= 0 && k < nz) ? double(rank) : 0.0;

    // Exercise: do the exchange through the API that takes a communicator
    core::exchange_ghosts(f, mesh, MPI_COMM_WORLD);

    // Helper lambdas to check a whole face is a constant value
    auto check_x_face = [&](int i0, int expected_rank)
    {
        const double expected = expected_rank >= 0 ? double(expected_rank) : 0.0;
        for (int k = 0; k < nz; ++k)
            for (int j = 0; j < ny; ++j)
                for (int g = 0; g < ng; ++g)
                    REQUIRE(f(i0 + g, j, k) == expected);
    };
    auto check_y_face = [&](int j0, int expected_rank)
    {
        const double expected = expected_rank >= 0 ? double(expected_rank) : 0.0;
        for (int k = 0; k < nz; ++k)
            for (int i = 0; i < nx; ++i)
                for (int g = 0; g < ng; ++g)
                    REQUIRE(f(i, j0 + g, k) == expected);
    };
    auto check_z_face = [&](int k0, int expected_rank)
    {
        const double expected = expected_rank >= 0 ? double(expected_rank) : 0.0;
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
                for (int g = 0; g < ng; ++g)
                    REQUIRE(f(i, j, k0 + g) == expected);
    };

    // Faces: -X, +X, -Y, +Y, -Z, +Z
    check_x_face(-ng, xneg); // left ghosts
    check_x_face(nx, xpos);  // right ghosts
    check_y_face(-ng, yneg); // front ghosts
    check_y_face(ny, ypos);  // back ghosts
    check_z_face(-ng, zneg); // bottom ghosts
    check_z_face(nz, zpos);  // top ghosts

    if (cart != MPI_COMM_WORLD)
        MPI_Comm_free(&cart);
}
