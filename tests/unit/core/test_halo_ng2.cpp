#include "Field.hpp"
#include "HaloExchange.hpp"
#include "MemoryManager.hpp"
#include "Mesh.hpp"
#include <catch2/catch_test_macros.hpp>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

using core::Field;
using core::Mesh;

TEST_CASE("Halo exchange works for ng=2 (faces)", "[mpi][halo]")
{
#ifndef HAVE_MPI
    SUCCEED("MPI not enabled, skipping halo test");
#else
    int rank = 0, nprocs = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    // Small grid with ng=2
    const int nx = 5, ny = 4, nz = 3, ng = 2;
    Mesh mesh{{nx, ny, nz}, ng};

    auto& mm = core::MemoryManager::instance();
    const auto ext = mesh.extents();
    double* raw =
        mm.allocate<double>(static_cast<std::size_t>(ext[0]) * static_cast<std::size_t>(ext[1]) *
                            static_cast<std::size_t>(ext[2]));
    Field<double> f(raw, ext, ng);

    // Interior = rank; ghosts = sentinel
    for (int k = -ng; k < nz + ng; ++k)
        for (int j = -ng; j < ny + ng; ++j)
            for (int i = -ng; i < nx + ng; ++i)
                f(i, j, k) = (i >= 0 && i < nx && j >= 0 && j < ny && k >= 0 && k < nz)
                                 ? static_cast<double>(rank)
                                 : -777.0;

    // Do the exchange on WORLD
    core::exchange_ghosts(f, mesh, MPI_COMM_WORLD);

    // Build same cartesian to discover neighbors
    int dims[3] = {0, 0, 0}, periods[3] = {0, 0, 0};
    MPI_Dims_create(nprocs, 3, dims);
    MPI_Comm cart = MPI_COMM_NULL;
    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 1, &cart);
    if (cart == MPI_COMM_NULL)
        cart = MPI_COMM_WORLD;

    int xneg, xpos, yneg, ypos, zneg, zpos;
    MPI_Cart_shift(cart, 0, 1, &xneg, &xpos);
    MPI_Cart_shift(cart, 1, 1, &yneg, &ypos);
    MPI_Cart_shift(cart, 2, 1, &zneg, &zpos);

    auto check_x = [&](int i0, int neighbor)
    {
        const double expected = (neighbor >= 0) ? double(neighbor) : -777.0;
        for (int k = 0; k < nz; ++k)
            for (int j = 0; j < ny; ++j)
                for (int g = 0; g < ng; ++g)
                    REQUIRE(f(i0 + g, j, k) == expected);
    };
    auto check_y = [&](int j0, int neighbor)
    {
        const double expected = (neighbor >= 0) ? double(neighbor) : -777.0;
        for (int k = 0; k < nz; ++k)
            for (int i = 0; i < nx; ++i)
                for (int g = 0; g < ng; ++g)
                    REQUIRE(f(i, j0 + g, k) == expected);
    };
    auto check_z = [&](int k0, int neighbor)
    {
        const double expected = (neighbor >= 0) ? double(neighbor) : -777.0;
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
                for (int g = 0; g < ng; ++g)
                    REQUIRE(f(i, j, k0 + g) == expected);
    };

    check_x(-ng, xneg);
    check_x(nx, xpos);
    check_y(-ng, yneg);
    check_y(ny, ypos);
    check_z(-ng, zneg);
    check_z(nz, zpos);

    if (cart != MPI_COMM_WORLD)
        MPI_Comm_free(&cart);
    mm.release(raw);
#endif
}
