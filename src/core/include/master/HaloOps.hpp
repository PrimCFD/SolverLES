#pragma once
#include "master/FieldCatalog.hpp"
#include "memory/MpiBox.hpp"
#include "mesh/Field.hpp"
#include "mesh/Mesh.hpp"

#ifdef HAVE_MPI
#include "mesh/HaloExchange.hpp"
#endif

#include <initializer_list>
#include <string>
#include <vector>

namespace core::master
{

// -------- periodic face wrap for faces only (non-MPI fallback) --------
template <class T>
static inline void wrap_periodic_faces(core::mesh::Field<T>& f, const core::mesh::Mesh& m)
{
    const int nx = m.local[0], ny = m.local[1], nz = m.local[2], ng = m.ng;
    if (ng == 0)
        return;

    // X faces
    if (m.periodic[0])
    {
        for (int k = 0; k < nz; ++k)
            for (int j = 0; j < ny; ++j)
                for (int g = 1; g <= ng; ++g)
                {
                    f(-g, j, k) = f(nx - g, j, k);        // west ghosts <- east interior
                    f(nx - 1 + g, j, k) = f(g - 1, j, k); // east ghosts <- west interior
                }
    }
    // Y faces
    if (m.periodic[1])
    {
        for (int k = 0; k < nz; ++k)
            for (int i = 0; i < nx; ++i)
                for (int g = 1; g <= ng; ++g)
                {
                    f(i, -g, k) = f(i, ny - g, k);
                    f(i, ny - 1 + g, k) = f(i, g - 1, k);
                }
    }
    // Z faces
    if (m.periodic[2])
    {
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
                for (int g = 1; g <= ng; ++g)
                {
                    f(i, j, -g) = f(i, j, nz - g);
                    f(i, j, nz - 1 + g) = f(i, j, g - 1);
                }
    }
    // Faces are enough for 7-point Laplacians.
}

inline void exchange_all_fields(FieldCatalog& cat, const core::mesh::Mesh& m,
                                void* mpi_comm_void /* may be nullptr */)
{
#ifdef HAVE_MPI
    if (mpi_comm_void)
    {
        MPI_Comm comm = mpi_unbox(mpi_comm_void);
        int comm_size = 1;
        MPI_Comm_size(comm, &comm_size);
        int field_idx = 0;
        for (const AnyFieldView& v : cat.all_views())
        {
            const auto e = v.extents;
            if (v.elem_size == sizeof(double))
            {
                core::mesh::Field<double> f(static_cast<double*>(v.host_ptr), e, m.ng);
                int tag_base = 100 + 16 * field_idx;  // unique tag range per field
                bool ok = core::mesh::exchange_ghosts(f, m, comm, tag_base);
                if (!ok && comm_size > 1) {
                    int rank = -1; MPI_Comm_rank(comm, &rank);
                    if (rank == 0) {
                        std::cerr << "[halo] WARNING: MPI halo exchange was skipped on a multi-rank run (non-Cartesian comm or periodicity mismatch). Halos will be stale.\n";
                    }
                }
                if (!ok && comm_size == 1 && (m.periodic[0] || m.periodic[1] || m.periodic[2]))
                    wrap_periodic_faces(f, m);
            }
            else if (v.elem_size == sizeof(float))
            {
                core::mesh::Field<float> f(static_cast<float*>(v.host_ptr), e, m.ng);
                int tag_base = 100 + 16 * field_idx;  // unique tag range per field
                bool ok = core::mesh::exchange_ghosts(f, m, comm, tag_base);                if (!ok && comm_size == 1 && (m.periodic[0] || m.periodic[1] || m.periodic[2]))
                    wrap_periodic_faces(f, m);
            }
            ++field_idx;
        }
        return;
    }
#endif
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
#ifdef HAVE_MPI
    if (mpi_comm_void)
    {
        MPI_Comm comm = mpi_unbox(mpi_comm_void);
        int comm_size = 1;
        MPI_Comm_size(comm, &comm_size);
        int field_idx = 0;
        for (const char* name : names)
        {
            if (!cat.contains(name)) {
                ++field_idx;
                continue; 
            }
            auto v = cat.view(name);
            const auto e = v.extents;
            if (v.elem_size == sizeof(double))
            {
                core::mesh::Field<double> f(static_cast<double*>(v.host_ptr), e, m.ng);
                int tag_base = 100 + 16 * field_idx;  // unique tag range per field
                bool ok = core::mesh::exchange_ghosts(f, m, comm, tag_base);
                if (!ok && comm_size > 1) {
                    int rank = -1; MPI_Comm_rank(comm, &rank);
                    if (rank == 0) {
                        std::cerr<< "[halo] WARNING: MPI halo exchange was skipped on a multi-rank run for field" << name << "(non-Cartesian comm or periodicity mismatch). Halos will be stale.\n";
                    }
                }
                if (!ok && comm_size == 1 && (m.periodic[0] || m.periodic[1] || m.periodic[2]))
                    wrap_periodic_faces(f, m);
            }
            else if (v.elem_size == sizeof(float))
            {
                core::mesh::Field<float> f(static_cast<float*>(v.host_ptr), e, m.ng);
                int tag_base = 100 + 16 * field_idx;  // unique tag range per field
                bool ok = core::mesh::exchange_ghosts(f, m, comm, tag_base);
                if (!ok && comm_size == 1 && (m.periodic[0] || m.periodic[1] || m.periodic[2]))
                    wrap_periodic_faces(f, m);
            }
            ++field_idx;
        }
        return;
    }
#endif
    if (!(m.periodic[0] || m.periodic[1] || m.periodic[2]))
        return;

    for (const char* name : names)
    {
        if (!cat.contains(name)) continue;
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
