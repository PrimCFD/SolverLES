#pragma once
#include <array>
#include <cstdint>
#include <type_traits>
#include <vector>

#include "Field.hpp"
#include "MemoryManager.hpp"
#include "Mesh.hpp"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

namespace core
{

// --- MPI type mapping -------------------------------------------------------
#ifdef HAVE_MPI
template <class T> inline MPI_Datatype mpi_dt();
template <> inline MPI_Datatype mpi_dt<float>()
{
    return MPI_FLOAT;
}
template <> inline MPI_Datatype mpi_dt<double>()
{
    return MPI_DOUBLE;
}
template <> inline MPI_Datatype mpi_dt<std::int32_t>()
{
    return MPI_INT;
}
template <> inline MPI_Datatype mpi_dt<std::int64_t>()
{
    return MPI_LONG_LONG;
} // portable
#endif

// --- Public API --------------------------------------------------------------
// Convenience: world communicator
template <class T> inline void exchange_ghosts(Field<T>& f, const Mesh& m)
{
#ifdef HAVE_MPI
    exchange_ghosts(f, m, MPI_COMM_WORLD);
#else
    (void) f;
    (void) m; // no-op without MPI
#endif
}

// 3D exchange (faces + edges + corners)
template <class T>
inline void exchange_ghosts(Field<T>& f, const Mesh& m,
/*in*/
#ifdef HAVE_MPI
                            MPI_Comm user_comm
#else
                            int /*dummy*/ = 0
#endif
)
{
#ifndef HAVE_MPI
    (void) f;
    (void) m;
    return; // no-MPI build: do nothing
#else
    // Local extents and halo width
    const int nx = m.local[0];
    const int ny = m.local[1];
    const int nz = m.local[2];
    const int ng = m.ng;
    if (ng <= 0)
        return;

    // Build (or adapt to) a 3D Cartesian communicator
    int world_n = 1, world_rank = 0;
    MPI_Comm_size(user_comm, &world_n);
    MPI_Comm_rank(user_comm, &world_rank);

    int dims[3] = {0, 0, 0};    // auto factorization
    int periods[3] = {0, 0, 0}; // non-periodic by default
    MPI_Dims_create(world_n, 3, dims);
    MPI_Comm cart_comm = MPI_COMM_NULL;
    MPI_Cart_create(user_comm, 3, dims, periods, /*reorder=*/1, &cart_comm);
    if (cart_comm == MPI_COMM_NULL)
        cart_comm = user_comm; // fallback (single rank)

    int coords[3] = {0, 0, 0};
    if (cart_comm != MPI_COMM_NULL)
    {
        int ndims = 0;
        MPI_Cartdim_get(cart_comm, &ndims);
        if (ndims == 3)
            MPI_Cart_coords(cart_comm, world_rank, 3, coords);
    }

    // Describe one neighbor exchange
    struct Ex
    {
        int dx, dy, dz;    // neighbor offset (-1,0,1)
        int nbr;           // neighbor rank (or MPI_PROC_NULL)
        int sx, sy, sz;    // block size to send
        int ix0, iy0, iz0; // source start (interior)
        int ox0, oy0, oz0; // destination start (ghosts)
        std::size_t count; // sx*sy*sz
        T* sbuf{nullptr};
        T* rbuf{nullptr};
        MPI_Request reqs[2]{MPI_REQUEST_NULL, MPI_REQUEST_NULL}; // {recv, send}
    };

    std::vector<Ex> work;
    work.reserve(26);

    auto tag_of = [](int dx, int dy, int dz) -> int
    {
        // Unique tag in [0..26], offset so it's >=100 (avoid collisions with other code)
        return 100 + (dx + 1) * 9 + (dy + 1) * 3 + (dz + 1);
    };

    // Helper to compute neighbor rank for (dx,dy,dz)
    auto neighbor_of = [&](int dx, int dy, int dz) -> int
    {
        int here[3];
        MPI_Cart_coords(cart_comm, world_rank, 3, here);
        int there[3] = {here[0] + dx, here[1] + dy, here[2] + dz};
        int dims_local[3], periods_local[3], coords_dummy[3];
        MPI_Cart_get(cart_comm, 3, dims_local, periods_local, coords_dummy);

        // Handle non-periodic out-of-bounds
        for (int a = 0; a < 3; ++a)
        {
            if (!periods_local[a] && (there[a] < 0 || there[a] >= dims_local[a]))
            {
                return MPI_PROC_NULL;
            }
            if (periods_local[a])
            {
                if (there[a] < 0)
                    there[a] += dims_local[a];
                if (there[a] >= dims_local[a])
                    there[a] -= dims_local[a];
            }
        }
        int r = MPI_PROC_NULL;
        MPI_Cart_rank(cart_comm, there, &r);
        return r;
    };

    // Build the 26 neighbor exchanges
    for (int dx = -1; dx <= 1; ++dx)
        for (int dy = -1; dy <= 1; ++dy)
            for (int dz = -1; dz <= 1; ++dz)
            {
                if (dx == 0 && dy == 0 && dz == 0)
                    continue;

                Ex ex{};
                ex.dx = dx;
                ex.dy = dy;
                ex.dz = dz;
                ex.nbr = neighbor_of(dx, dy, dz);
                if (ex.nbr == MPI_PROC_NULL)
                    continue; // nothing to do

                ex.sx = (dx == 0) ? nx : ng;
                ex.sy = (dy == 0) ? ny : ng;
                ex.sz = (dz == 0) ? nz : ng;

                ex.ix0 = (dx < 0) ? 0 : (dx == 0 ? 0 : nx - ng);
                ex.iy0 = (dy < 0) ? 0 : (dy == 0 ? 0 : ny - ng);
                ex.iz0 = (dz < 0) ? 0 : (dz == 0 ? 0 : nz - ng);

                ex.ox0 = (dx < 0) ? -ng : (dx == 0 ? 0 : nx);
                ex.oy0 = (dy < 0) ? -ng : (dy == 0 ? 0 : ny);
                ex.oz0 = (dz < 0) ? -ng : (dz == 0 ? 0 : nz);

                ex.count = static_cast<std::size_t>(ex.sx) * static_cast<std::size_t>(ex.sy) *
                           static_cast<std::size_t>(ex.sz);

                work.push_back(ex);
            }

    // Allocate receive buffers and post Irecv first (to avoid races)
    auto& mm = core::MemoryManager::instance();
    for (auto& ex : work)
    {
        ex.rbuf = mm.allocate<T>(ex.count);
        const int tagOpp = tag_of(-ex.dx, -ex.dy, -ex.dz);
        MPI_Irecv(ex.rbuf, static_cast<int>(ex.count), mpi_dt<T>(), ex.nbr, tagOpp, cart_comm,
                  &ex.reqs[0]);
    }

    // Pack and Isend
    for (auto& ex : work)
    {
        ex.sbuf = mm.allocate<T>(ex.count);
        // pack: (k,j,i) fastest-varying i
        std::size_t p = 0;
        for (int kz = 0; kz < ex.sz; ++kz)
            for (int jy = 0; jy < ex.sy; ++jy)
                for (int ix = 0; ix < ex.sx; ++ix)
                    ex.sbuf[p++] = f(ex.ix0 + ix, ex.iy0 + jy, ex.iz0 + kz);

        const int tag = tag_of(ex.dx, ex.dy, ex.dz);
        MPI_Isend(ex.sbuf, static_cast<int>(ex.count), mpi_dt<T>(), ex.nbr, tag, cart_comm,
                  &ex.reqs[1]);
    }

    // Wait for all comms to complete
    std::vector<MPI_Request> all;
    all.reserve(work.size() * 2);
    for (auto& ex : work)
    {
        all.push_back(ex.reqs[0]);
        all.push_back(ex.reqs[1]);
    }
    if (!all.empty())
        MPI_Waitall(static_cast<int>(all.size()), all.data(), MPI_STATUSES_IGNORE);

    // Unpack and free
    for (auto& ex : work)
    {
        std::size_t p = 0;
        for (int kz = 0; kz < ex.sz; ++kz)
            for (int jy = 0; jy < ex.sy; ++jy)
                for (int ix = 0; ix < ex.sx; ++ix)
                    f(ex.ox0 + ix, ex.oy0 + jy, ex.oz0 + kz) = ex.rbuf[p++];

        mm.release(ex.sbuf);
        mm.release(ex.rbuf);
    }

    if (cart_comm != user_comm)
        MPI_Comm_free(&cart_comm);
#endif
}

} // namespace core
