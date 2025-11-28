#pragma once
#include "master/FieldCatalog.hpp"
#include "master/Log.hpp"
#include "memory/MpiBox.hpp"
#include "mesh/Field.hpp"
#include "mesh/Mesh.hpp"
#include <cstring>
#include <string>

#include "mesh/HaloExchange.hpp"

#include <initializer_list>
#include <string>
#include <vector>

namespace core::master
{

// ---- Deterministic per-field tag base derived from field name (rank-invariant),
// ---- and guaranteed to lie within the communicator's legal tag range.
static inline int tag_base_for_name_comm(const char* name, MPI_Comm comm)
{
    // djb2-xor hash; stable across ranks
    unsigned h = 5381u;
    const unsigned char* s = reinterpret_cast<const unsigned char*>(name ? name : "");
    while (*s) { h = ((h << 5) + h) ^ *s++; }

    // Query MPI_TAG_UB to stay within the impl's legal tag range.
    int *p_ub = nullptr, flag = 0; // MPI says tags must be in [0..MPI_TAG_UB]
    MPI_Comm_get_attr(comm, MPI_TAG_UB, &p_ub, &flag);
    int tag_ub = (flag && p_ub) ? *p_ub : 32767;  // standard guarantees >= 32767

    // Reserve a window per field large enough for all 26 directions, rounded up.
    // Leave some headroom: choose 64.
    const int window = 64;
    const int guard  = 200;                 // keep away from tiny tags used elsewhere
    const int limit  = std::max(guard + window, tag_ub - window); // safety
    const int slots  = std::max(1, (limit - guard) / window);
    const int slot   = static_cast<int>(h % static_cast<unsigned>(slots));
    return guard + slot * window;
}

// -------- periodic face wrap for faces only (non-MPI fallback) --------
template <class T>
static inline void wrap_periodic_faces(core::mesh::Field<T>& f, const core::mesh::Mesh& m)
{

    const auto e = f.extents();
    const int ng = f.ng();
    const int nx = e[0] - 2 * ng;
    const int ny = e[1] - 2 * ng;
    const int nz = e[2] - 2 * ng;
    if (ng == 0)
        return;

    auto idx_lo = [](int n, int g) -> int
    {
        // lower ghost at -g pulls from interior index n-1, n-2, ... modulo n
        // handle n==1: always 0
        if (n <= 0)
            return 0;
        int r = (g % n);
        return (r == 0) ? (n - 1) : (n - r);
    };
    auto idx_hi = [](int n, int g) -> int
    {
        // upper ghost at n-1+g pulls from 0,1,2,... modulo n
        if (n <= 0)
            return 0;
        return (g - 1) % n;
    };

    // X faces
    if (m.periodic[0])
    {
        for (int k = 0; k < nz; ++k)
            for (int j = 0; j < ny; ++j)
                for (int g = 1; g <= ng; ++g)
                {
                    f(-g, j, k) = f(idx_lo(nx, g), j, k);
                    f(nx - 1 + g, j, k) = f(idx_hi(nx, g), j, k);
                }
    }
    // Y faces
    if (m.periodic[1])
    {
        for (int k = 0; k < nz; ++k)
            for (int i = 0; i < nx; ++i)
                for (int g = 1; g <= ng; ++g)
                {
                    f(i, -g, k) = f(i, idx_lo(ny, g), k);
                    f(i, ny - 1 + g, k) = f(i, idx_hi(ny, g), k);
                }
    }
    // Z faces
    if (m.periodic[2])
    {
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
                for (int g = 1; g <= ng; ++g)
                {
                    f(i, j, -g) = f(i, j, idx_lo(nz, g));
                    f(i, j, nz - 1 + g) = f(i, j, idx_hi(nz, g));
                }
    }
    // Faces are enough for 7-point Laplacians.
}

inline void exchange_all_fields(FieldCatalog& cat, const core::mesh::Mesh& m,
                                void* mpi_comm_void /* may be nullptr */)
{
    if (mpi_comm_void)
    {
        MPI_Comm comm = mpi_unbox(mpi_comm_void);
        int comm_size = 1;
        MPI_Comm_size(comm, &comm_size);
        for (const AnyFieldView& v : cat.all_views())
        {
            const auto e = v.extents;
            const int tag_base = tag_base_for_name_comm(
                v.name.empty() ? "" : v.name.c_str(), comm);
            if (v.elem_size == sizeof(double))
            {
                core::mesh::Field<double> f(static_cast<double*>(v.host_ptr), e, m.ng);
                bool ok = core::mesh::exchange_ghosts(f, m, comm, tag_base, v.stagger);
                if (!ok && comm_size > 1)
                {
                    int rank = -1;
                    MPI_Comm_rank(comm, &rank);
                    if (rank == 0)
                    {
                        LOGW(
                            "[halo] MPI halo exchange was skipped on a multi-rank run "
                            "(non-Cartesian comm or periodicity mismatch). Halos will be stale.\n");
                    }
                }
            }
            else if (v.elem_size == sizeof(float))
            {
                core::mesh::Field<float> f(static_cast<float*>(v.host_ptr), e, m.ng);
                bool ok = core::mesh::exchange_ghosts(f, m, comm, tag_base, v.stagger);
                if (!ok && comm_size > 1)
                {
                    int rank = -1;
                    MPI_Comm_rank(comm, &rank);
                    if (rank == 0)
                    {
                        LOGW(
                            "[halo] MPI halo exchange was skipped on a multi-rank run "
                            "(non-Cartesian comm or periodicity mismatch). Halos will be stale.\n");
                    }
                }
            }
        }
        return;
    }

    // Non-MPI or comm==nullptr → periodic wrap fallback (only if some axis is periodic)
    if (!(m.periodic[0] || m.periodic[1] || m.periodic[2]))
        return;

    for (const AnyFieldView& v : cat.all_views())
    {
        const auto e = v.extents;
        if (v.elem_size == sizeof(double))
        {
            core::mesh::Field<double> f(static_cast<double*>(v.host_ptr), e, m.ng);
            wrap_periodic_faces(f, m);
        }
        else if (v.elem_size == sizeof(float))
        {
            core::mesh::Field<float> f(static_cast<float*>(v.host_ptr), e, m.ng);
            wrap_periodic_faces(f, m);
        }
    }
}

// Exchange only the given field names (e.g. {"u","v","w"})
inline void exchange_named_fields(FieldCatalog& cat, const core::mesh::Mesh& m,
                                  void* mpi_comm_void /* may be nullptr */,
                                  std::initializer_list<const char*> names)
{
    if (mpi_comm_void)
    {
        MPI_Comm comm = mpi_unbox(mpi_comm_void);
        int comm_size = 1;
        MPI_Comm_size(comm, &comm_size);
        int field_idx = 0;
        for (const char* name : names)
        {
            if (!cat.contains(name))
            {
                continue;
            }
            auto v = cat.view(name);
            const auto e = v.extents;
            if (v.elem_size == sizeof(double))
            {
                core::mesh::Field<double> f(static_cast<double*>(v.host_ptr), e, m.ng);
                const int tag_base = tag_base_for_name_comm(name, comm);
                bool ok = core::mesh::exchange_ghosts(f, m, comm, tag_base, v.stagger);
                if (!ok && comm_size > 1)
                {
                    int rank = -1;
                    MPI_Comm_rank(comm, &rank);
                    if (rank == 0)
                    {
                        LOGW(
                            "[halo] MPI halo exchange was skipped on a multi-rank run for field %s "
                            "(non-Cartesian comm or periodicity mismatch). Halos will be stale.\n",
                            name);
                    }
                }
            }
            else if (v.elem_size == sizeof(float))
            {
                core::mesh::Field<float> f(static_cast<float*>(v.host_ptr), e, m.ng);
                const int tag_base = tag_base_for_name_comm(name, comm);
                bool ok = core::mesh::exchange_ghosts(f, m, comm, tag_base, v.stagger);
                if (!ok && comm_size > 1)
                {
                    int rank = -1;
                    MPI_Comm_rank(comm, &rank);
                    if (rank == 0)
                    {
                        LOGW(
                            "[halo] MPI halo exchange was skipped on a multi-rank run for field %s "
                            "(non-Cartesian comm or periodicity mismatch). Halos will be stale.\n",
                            name);
                    }
                }
            }
        }
        return;
    }
    if (!(m.periodic[0] || m.periodic[1] || m.periodic[2]))
        return;

    for (const char* name : names)
    {
        if (!cat.contains(name))
            continue;
        auto v = cat.view(name);
        const auto e = v.extents;
        if (v.elem_size == sizeof(double))
        {
            core::mesh::Field<double> f(static_cast<double*>(v.host_ptr), e, m.ng);
            wrap_periodic_faces(f, m);
        }
        else if (v.elem_size == sizeof(float))
        {
            core::mesh::Field<float> f(static_cast<float*>(v.host_ptr), e, m.ng);
            wrap_periodic_faces(f, m);
        }
    }
}

} // namespace core::master
