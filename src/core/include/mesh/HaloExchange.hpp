#pragma once
#include "Field.hpp"
#include "Mesh.hpp"
#include <cassert>
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

namespace core::mesh
{

// Exchange the 6 axis-aligned faces (no edges/corners). Works for ng >= 1.
// Uses the caller-provided *Cartesian* communicator as-is.
// Returns true if halos were exchanged (or nothing was needed but the path is valid),
// false if no exchange was performed (caller may choose to fall back to local wrap).
template <class T>
bool exchange_ghosts(Field<T>& f, const Mesh& m, MPI_Comm comm, int tag_base = 100)
{
    static_assert(std::is_same_v<T, float> || std::is_same_v<T, double>,
                  "exchange_ghosts currently supports float/double element types");
    const int ng = m.ng;
    const auto e = f.extents(); // totals INCLUDING ghosts
    const int nx = e[0] - 2 * ng;
    const int ny = e[1] - 2 * ng;
    const int nz = e[2] - 2 * ng;

    // Bail early if there's nothing to exchange
    if (ng <= 0 || nx <= 0 || ny <= 0 || nz <= 0)
        return true;

    // Require initialized MPI
    int initialized = 0;
    MPI_Initialized(&initialized);
    if (!initialized)
        return false; // not usable here

    // Require a valid, caller-provided *Cartesian* communicator
    if (comm == MPI_COMM_NULL)
        return false;

    // Optional but helpful: during development, return errors instead of aborting
    // (comment out for production)
    MPI_Comm_set_errhandler(comm, MPI_ERRORS_RETURN);

    int size = 1, ierr = MPI_Comm_size(comm, &size);
    if (ierr != MPI_SUCCESS)
        return false; // nothing we can do; keep it fail-safe

    // Small local helper: periodic face wrap for faces only (matches master logic)
    auto idx_lo = [](int n, int g) -> int
    {
        if (n <= 0)
            return 0;
        int r = (g % n);
        return (r == 0) ? (n - 1) : (n - r);
    };
    auto idx_hi = [](int n, int g) -> int
    {
        if (n <= 0)
            return 0;
        return (g - 1) % n;
    };

    // Per-axis wrappers
    auto wrap_x_faces_local = [&](Field<T>& ff)
    {
        if (!m.periodic[0])
            return;
        for (int k = 0; k < nz; ++k)
            for (int j = 0; j < ny; ++j)
                for (int g = 1; g <= ng; ++g)
                {
                    ff(-g, j, k) = ff(idx_lo(nx, g), j, k);
                    ff(nx - 1 + g, j, k) = ff(idx_hi(nx, g), j, k);
                }
    };
    auto wrap_y_faces_local = [&](Field<T>& ff)
    {
        if (!m.periodic[1])
            return;
        for (int k = 0; k < nz; ++k)
            for (int i = 0; i < nx; ++i)
                for (int g = 1; g <= ng; ++g)
                {
                    ff(i, -g, k) = ff(i, idx_lo(ny, g), k);
                    ff(i, ny - 1 + g, k) = ff(i, idx_hi(ny, g), k);
                }
    };
    auto wrap_z_faces_local = [&](Field<T>& ff)
    {
        if (!m.periodic[2])
            return;
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
                for (int g = 1; g <= ng; ++g)
                {
                    ff(i, j, -g) = ff(i, j, idx_lo(nz, g));
                    ff(i, j, nz - 1 + g) = ff(i, j, idx_hi(nz, g));
                }
    };

    // If single rank:
    //  - no neighbors if not periodic -> nothing to do
    //  - periodic in any axis -> do local periodic wrap and return
    if (size == 1)
    {
        if (m.periodic[0])
            wrap_x_faces_local(f);
        if (m.periodic[1])
            wrap_y_faces_local(f);
        if (m.periodic[2])
            wrap_z_faces_local(f);
        return true;
    }

    // Verify the communicator is already Cartesian and matches Mesh periodicity
    int topo = MPI_UNDEFINED;
    MPI_Topo_test(comm, &topo);
    if (topo != MPI_CART)
    {
        // Not a Cartesian communicator: do nothing (explicit contract violation)
        return false;
    }
    // Optional consistency check (non-fatal)
    int cdims[3] = {0, 0, 0}, cperiods[3] = {0, 0, 0}, coords[3] = {0, 0, 0};
    MPI_Cart_get(comm, 3, cdims, cperiods, coords);
    // If needed: ensure periodic flags match Mesh (ignore cdims; decomposition is up to the app)
    if ((cperiods[0] != (m.periodic[0] ? 1 : 0)) || (cperiods[1] != (m.periodic[1] ? 1 : 0)) ||
        (cperiods[2] != (m.periodic[2] ? 1 : 0)))
    {
        // Mismatch between Mesh.periodic and communicator topology — safer to bail
        return false;
    }

    // Neighbor ranks for faces
    int xneg = MPI_PROC_NULL, xpos = MPI_PROC_NULL;
    int yneg = MPI_PROC_NULL, ypos = MPI_PROC_NULL;
    int zneg = MPI_PROC_NULL, zpos = MPI_PROC_NULL;
    MPI_Cart_shift(comm, 0, 1, &xneg, &xpos);
    MPI_Cart_shift(comm, 1, 1, &yneg, &ypos);
    MPI_Cart_shift(comm, 2, 1, &zneg, &zpos);
    int myrank = -1;
    MPI_Comm_rank(comm, &myrank);

    // Pack/unpack helpers
    auto pack_x = [&](int i_start, std::vector<T>& buf)
    {
        buf.resize(static_cast<std::size_t>(ny) * nz * ng);
        std::size_t p = 0;
        const int n = nx;
        const int nsend = (n > 0) ? std::min(ng, n) : 0;
        for (int k = 0; k < nz; ++k)
            for (int j = 0; j < ny; ++j)
                for (int g = 0; g < ng; ++g)
                {
                    int gi = (n > 0) ? ((i_start + (g % nsend)) % n + n) % n : 0;
                    buf[p++] = f(gi, j, k);
                }
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
        const int n = ny;
        const int nsend = (n > 0) ? std::min(ng, n) : 0;
        for (int k = 0; k < nz; ++k)
            for (int i = 0; i < nx; ++i)
                for (int g = 0; g < ng; ++g)
                {
                    int gj = (n > 0) ? ((j_start + (g % nsend)) % n + n) % n : 0;
                    buf[p++] = f(i, gj, k);
                }
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
        const int n = nz;
        const int nsend = (n > 0) ? std::min(ng, n) : 0;
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
                for (int g = 0; g < ng; ++g)
                {
                    int gk = (n > 0) ? ((k_start + (g % nsend)) % n + n) % n : 0;
                    buf[p++] = f(i, j, gk);
                }
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
    // Map T -> MPI_Datatype
    MPI_Datatype Tmpi = MPI_BYTE;
    if constexpr (std::is_same_v<T, double>)
        Tmpi = MPI_DOUBLE;
    else if constexpr (std::is_same_v<T, float>)
        Tmpi = MPI_FLOAT;

    // If a periodic axis has a self-neighbor (rank==myrank), prefer local copy for that axis.
    const bool xneg_self = (xneg == myrank);
    const bool xpos_self = (xpos == myrank);
    const bool yneg_self = (yneg == myrank);
    const bool ypos_self = (ypos == myrank);
    const bool zneg_self = (zneg == myrank);
    const bool zpos_self = (zpos == myrank);

    // Post recvs (only for true remote neighbors)
    if (xneg != MPI_PROC_NULL && !xneg_self)
    {
        r_xn.resize((size_t) ny * nz * ng);
        MPI_Irecv(r_xn.data(), (int) r_xn.size(), Tmpi, xneg, tag_base + 0, comm, &reqs[rcount++]);
    }
    if (xpos != MPI_PROC_NULL && !xpos_self)
    {
        r_xp.resize((size_t) ny * nz * ng);
        MPI_Irecv(r_xp.data(), (int) r_xp.size(), Tmpi, xpos, tag_base + 1, comm, &reqs[rcount++]);
    }
    if (yneg != MPI_PROC_NULL && !yneg_self)
    {
        r_yn.resize((size_t) nx * nz * ng);
        MPI_Irecv(r_yn.data(), (int) r_yn.size(), Tmpi, yneg, tag_base + 2, comm, &reqs[rcount++]);
    }
    if (ypos != MPI_PROC_NULL && !ypos_self)
    {
        r_yp.resize((size_t) nx * nz * ng);
        MPI_Irecv(r_yp.data(), (int) r_yp.size(), Tmpi, ypos, tag_base + 3, comm, &reqs[rcount++]);
    }
    if (zneg != MPI_PROC_NULL && !zneg_self)
    {
        r_zn.resize((size_t) nx * ny * ng);
        MPI_Irecv(r_zn.data(), (int) r_zn.size(), Tmpi, zneg, tag_base + 4, comm, &reqs[rcount++]);
    }
    if (zpos != MPI_PROC_NULL && !zpos_self)
    {
        r_zp.resize((size_t) nx * ny * ng);
        MPI_Irecv(r_zp.data(), (int) r_zp.size(), Tmpi, zpos, tag_base + 5, comm, &reqs[rcount++]);
    }

    // Pack & send
    if (xneg != MPI_PROC_NULL && !xneg_self)
    {
        pack_x(0, s_xn);
        MPI_Isend(s_xn.data(), (int) s_xn.size(), Tmpi, xneg, tag_base + 1, comm, &reqs[rcount++]);
    }
    if (xpos != MPI_PROC_NULL && !xpos_self)
    {
        int nsend = std::min(ng, nx);
        pack_x(nx - nsend, s_xp);
        MPI_Isend(s_xp.data(), (int) s_xp.size(), Tmpi, xpos, tag_base + 0, comm, &reqs[rcount++]);
    }
    if (yneg != MPI_PROC_NULL && !yneg_self)
    {
        pack_y(0, s_yn);
        MPI_Isend(s_yn.data(), (int) s_yn.size(), Tmpi, yneg, tag_base + 3, comm, &reqs[rcount++]);
    }
    if (ypos != MPI_PROC_NULL && !ypos_self)
    {
        int nsend = std::min(ng, ny);
        pack_y(ny - nsend, s_yp);
        MPI_Isend(s_yp.data(), (int) s_yp.size(), Tmpi, ypos, tag_base + 2, comm, &reqs[rcount++]);
    }
    if (zneg != MPI_PROC_NULL && !zneg_self)
    {
        pack_z(0, s_zn);
        MPI_Isend(s_zn.data(), (int) s_zn.size(), Tmpi, zneg, tag_base + 5, comm, &reqs[rcount++]);
    }
    if (zpos != MPI_PROC_NULL && !zpos_self)
    {
        int nsend = std::min(ng, nz);
        pack_z(nz - nsend, s_zp);
        MPI_Isend(s_zp.data(), (int) s_zp.size(), Tmpi, zpos, tag_base + 4, comm, &reqs[rcount++]);
    }

    if (rcount)
        MPI_Waitall(rcount, reqs, MPI_STATUSES_IGNORE);

    // Unpack to ghosts (only where a remote message arrived)
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

    // For any axis with self-neighbor (periodic, cdims==1), do a local copy on that axis only.
    // This avoids clobbering halos that were filled via MPI on other periodic-but-distributed axes.
    if (xneg_self || xpos_self)
        wrap_x_faces_local(f);
    if (yneg_self || ypos_self)
        wrap_y_faces_local(f);
    if (zneg_self || zpos_self)
        wrap_z_faces_local(f);

    return true;
}

} // namespace core::mesh
