#pragma once
#include "Field.hpp"
#include "Mesh.hpp"
#include <mpi.h>
#include <vector>

/**
* @file HaloExchange.hpp
* @ingroup memory
* @brief MPI façade to exchange ghost/halo slabs of ghost‑padded fields.
*
* Posts six non‑blocking \c MPI_Isend/\c MPI_Irecv operations into the ghost layers
* of each field and overlaps communication with interior computation. Interfaces
* operate on \c Field\<T\> views to avoid extra copies.
*
* @rst
*.. code-block:: cpp
*
*    exchange_ghosts(rho, mesh, cart_comm); // X±, Y±, Z± slabs
* @endrst
*
* @warning All participating fields must share the same ghost width and extents.
*/

namespace core
{

// Exchange the 6 axis-aligned faces (no edges/corners). Works for ng >= 1.
// Uses a Cartesian communicator derived from the provided `comm`.
template <class T> void exchange_ghosts(Field<T>& f, const Mesh& m, MPI_Comm comm)
{
    const int nx = m.local[0], ny = m.local[1], nz = m.local[2];
    const int ng = m.ng;

    // Derive a 3D Cartesian communicator from `comm`
    int size = 1;
    MPI_Comm_size(comm, &size);
    int dims[3] = {0, 0, 0};
    int periods[3] = {0, 0, 0};
    MPI_Dims_create(size, 3, dims);
    MPI_Comm cart = MPI_COMM_NULL;
    MPI_Cart_create(comm, 3, dims, periods, /*reorder=*/1, &cart);
    if (cart == MPI_COMM_NULL)
        cart = comm;

    // Neighbor ranks for faces
    int xneg, xpos, yneg, ypos, zneg, zpos;
    MPI_Cart_shift(cart, 0, 1, &xneg, &xpos);
    MPI_Cart_shift(cart, 1, 1, &yneg, &ypos);
    MPI_Cart_shift(cart, 2, 1, &zneg, &zpos);

    // Pack/unpack helpers
    auto pack_x = [&](int i_start, std::vector<T>& buf)
    {
        buf.resize(static_cast<std::size_t>(ny) * nz * ng);
        std::size_t p = 0;
        for (int k = 0; k < nz; ++k)
            for (int j = 0; j < ny; ++j)
                for (int g = 0; g < ng; ++g)
                    buf[p++] = f(i_start + g, j, k);
    };
    auto unpack_x = [&](int i_ghost_start, const std::vector<T>& buf)
    {
        std::size_t p = 0;
        for (int k = 0; k < nz; ++k)
            for (int j = 0; j < ny; ++j)
                for (int g = 0; g < ng; ++g)
                    f(i_ghost_start + g, j, k) = buf[p++];
    };

    auto pack_y = [&](int j_start, std::vector<T>& buf)
    {
        buf.resize(static_cast<std::size_t>(nx) * nz * ng);
        std::size_t p = 0;
        for (int k = 0; k < nz; ++k)
            for (int i = 0; i < nx; ++i)
                for (int g = 0; g < ng; ++g)
                    buf[p++] = f(i, j_start + g, k);
    };
    auto unpack_y = [&](int j_ghost_start, const std::vector<T>& buf)
    {
        std::size_t p = 0;
        for (int k = 0; k < nz; ++k)
            for (int i = 0; i < nx; ++i)
                for (int g = 0; g < ng; ++g)
                    f(i, j_ghost_start + g, k) = buf[p++];
    };

    auto pack_z = [&](int k_start, std::vector<T>& buf)
    {
        buf.resize(static_cast<std::size_t>(nx) * ny * ng);
        std::size_t p = 0;
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
                for (int g = 0; g < ng; ++g)
                    buf[p++] = f(i, j, k_start + g);
    };
    auto unpack_z = [&](int k_ghost_start, const std::vector<T>& buf)
    {
        std::size_t p = 0;
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
                for (int g = 0; g < ng; ++g)
                    f(i, j, k_ghost_start + g) = buf[p++];
    };

    // Buffers & requests
    std::vector<T> s_xn, r_xn, s_xp, r_xp;
    std::vector<T> s_yn, r_yn, s_yp, r_yp;
    std::vector<T> s_zn, r_zn, s_zp, r_zp;
    MPI_Request reqs[12];
    int rcount = 0;
    auto bytes = [](std::size_t n) -> int { return static_cast<int>(n); }; // small halos

    // Post recvs (MPI_BYTE for portability across T)
    if (xneg != MPI_PROC_NULL)
    {
        r_xn.resize((size_t) ny * nz * ng);
        MPI_Irecv(r_xn.data(), bytes(r_xn.size() * sizeof(T)), MPI_BYTE, xneg, 100, cart,
                  &reqs[rcount++]);
    }
    if (xpos != MPI_PROC_NULL)
    {
        r_xp.resize((size_t) ny * nz * ng);
        MPI_Irecv(r_xp.data(), bytes(r_xp.size() * sizeof(T)), MPI_BYTE, xpos, 101, cart,
                  &reqs[rcount++]);
    }
    if (yneg != MPI_PROC_NULL)
    {
        r_yn.resize((size_t) nx * nz * ng);
        MPI_Irecv(r_yn.data(), bytes(r_yn.size() * sizeof(T)), MPI_BYTE, yneg, 200, cart,
                  &reqs[rcount++]);
    }
    if (ypos != MPI_PROC_NULL)
    {
        r_yp.resize((size_t) nx * nz * ng);
        MPI_Irecv(r_yp.data(), bytes(r_yp.size() * sizeof(T)), MPI_BYTE, ypos, 201, cart,
                  &reqs[rcount++]);
    }
    if (zneg != MPI_PROC_NULL)
    {
        r_zn.resize((size_t) nx * ny * ng);
        MPI_Irecv(r_zn.data(), bytes(r_zn.size() * sizeof(T)), MPI_BYTE, zneg, 300, cart,
                  &reqs[rcount++]);
    }
    if (zpos != MPI_PROC_NULL)
    {
        r_zp.resize((size_t) nx * ny * ng);
        MPI_Irecv(r_zp.data(), bytes(r_zp.size() * sizeof(T)), MPI_BYTE, zpos, 301, cart,
                  &reqs[rcount++]);
    }

    // Pack & send
    if (xneg != MPI_PROC_NULL)
    {
        pack_x(0, s_xn);
        MPI_Isend(s_xn.data(), bytes(s_xn.size() * sizeof(T)), MPI_BYTE, xneg, 101, cart,
                  &reqs[rcount++]);
    }
    if (xpos != MPI_PROC_NULL)
    {
        pack_x(nx - ng, s_xp);
        MPI_Isend(s_xp.data(), bytes(s_xp.size() * sizeof(T)), MPI_BYTE, xpos, 100, cart,
                  &reqs[rcount++]);
    }
    if (yneg != MPI_PROC_NULL)
    {
        pack_y(0, s_yn);
        MPI_Isend(s_yn.data(), bytes(s_yn.size() * sizeof(T)), MPI_BYTE, yneg, 201, cart,
                  &reqs[rcount++]);
    }
    if (ypos != MPI_PROC_NULL)
    {
        pack_y(ny - ng, s_yp);
        MPI_Isend(s_yp.data(), bytes(s_yp.size() * sizeof(T)), MPI_BYTE, ypos, 200, cart,
                  &reqs[rcount++]);
    }
    if (zneg != MPI_PROC_NULL)
    {
        pack_z(0, s_zn);
        MPI_Isend(s_zn.data(), bytes(s_zn.size() * sizeof(T)), MPI_BYTE, zneg, 301, cart,
                  &reqs[rcount++]);
    }
    if (zpos != MPI_PROC_NULL)
    {
        pack_z(nz - ng, s_zp);
        MPI_Isend(s_zp.data(), bytes(s_zp.size() * sizeof(T)), MPI_BYTE, zpos, 300, cart,
                  &reqs[rcount++]);
    }

    if (rcount)
        MPI_Waitall(rcount, reqs, MPI_STATUSES_IGNORE);

    // Unpack to ghosts (only where a message arrived)
    if (!r_xn.empty())
        unpack_x(-ng, r_xn);
    if (!r_xp.empty())
        unpack_x(nx, r_xp);
    if (!r_yn.empty())
        unpack_y(-ng, r_yn);
    if (!r_yp.empty())
        unpack_y(ny, r_yp);
    if (!r_zn.empty())
        unpack_z(-ng, r_zn);
    if (!r_zp.empty())
        unpack_z(nz, r_zp);

    if (cart != comm)
        MPI_Comm_free(&cart);
}

} // namespace core
