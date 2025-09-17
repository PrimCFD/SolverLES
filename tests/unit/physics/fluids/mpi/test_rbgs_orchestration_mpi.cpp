#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <mpi.h>
#include <vector>
#include "kernels_fluids.h"

using Catch::Approx;

static inline size_t lin(int I, int J, int K, int nx_tot, int ny_tot)
{
    return static_cast<size_t>(I + nx_tot * (J + ny_tot * K));
}

// Exchange halos across rank interfaces along X (1D slab)
static void exchange_x(double* a, int nx_tot, int ny_tot, int nz_tot, int ng, MPI_Comm comm)
{
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    const int left = (rank == 0) ? MPI_PROC_NULL : rank - 1;
    const int right = (rank == size - 1) ? MPI_PROC_NULL : rank + 1;

    const int nyz = ny_tot * nz_tot;
    const int slab = ng * nyz;
    std::vector<double> sendL(slab), sendR(slab), recvL(slab), recvR(slab);

    for (int g = 0; g < ng; ++g)
    {
        const int IL = ng + g;
        const int IR = (nx_tot - 2 * ng) + g; // last interior plane
        for (int K = 0; K < nz_tot; ++K)
            for (int J = 0; J < ny_tot; ++J)
            {
                const size_t s = g * nyz + J + ny_tot * K;
                sendL[s] = a[lin(IL, J, K, nx_tot, ny_tot)];
                sendR[s] = a[lin(IR, J, K, nx_tot, ny_tot)];
            }
    }

    MPI_Sendrecv(sendL.data(), slab, MPI_DOUBLE, left, 100, recvR.data(), slab, MPI_DOUBLE, right,
                 100, comm, MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendR.data(), slab, MPI_DOUBLE, right, 200, recvL.data(), slab, MPI_DOUBLE, left,
                 200, comm, MPI_STATUS_IGNORE);

    for (int g = 0; g < ng; ++g)
    {
        const int IHL = g;
        const int IHR = (nx_tot - ng) + g;
        for (int K = 0; K < nz_tot; ++K)
            for (int J = 0; J < ny_tot; ++J)
            {
                const size_t s = g * nyz + J + ny_tot * K;
                a[lin(IHL, J, K, nx_tot, ny_tot)] = recvL[s];
                a[lin(IHR, J, K, nx_tot, ny_tot)] = recvR[s];
            }
    }
}

// Apply "mirror" (Neumann) physical BCs on all domain faces.
// X faces are physical only for the first/last rank; Y/Z faces are physical for all ranks.
static void mirror_physical(double* a, int nx_tot, int ny_tot, int nz_tot, int ng, int rank,
                            int size)
{
    // X- faces (physical only on end ranks)
    if (rank == 0)
    {
        for (int g = 0; g < ng; ++g)
        {
            const int Ihalo = g;
            const int Iint = ng + g;
            for (int K = 0; K < nz_tot; ++K)
                for (int J = 0; J < ny_tot; ++J)
                    a[lin(Ihalo, J, K, nx_tot, ny_tot)] = a[lin(Iint, J, K, nx_tot, ny_tot)];
        }
    }
    if (rank == size - 1)
    {
        for (int g = 0; g < ng; ++g)
        {
            const int Ihalo = (nx_tot - ng) + g;
            const int Iint = (nx_tot - ng - 1) - g;
            for (int K = 0; K < nz_tot; ++K)
                for (int J = 0; J < ny_tot; ++J)
                    a[lin(Ihalo, J, K, nx_tot, ny_tot)] = a[lin(Iint, J, K, nx_tot, ny_tot)];
        }
    }
    // Y faces (always physical)
    for (int g = 0; g < ng; ++g)
    {
        const int JLo = g, JLoInt = ng + g;
        const int JHi = (ny_tot - ng) + g, JHiInt = (ny_tot - ng - 1) - g;
        for (int K = 0; K < nz_tot; ++K)
            for (int I = 0; I < nx_tot; ++I)
            {
                a[lin(I, JLo, K, nx_tot, ny_tot)] = a[lin(I, JLoInt, K, nx_tot, ny_tot)];
                a[lin(I, JHi, K, nx_tot, ny_tot)] = a[lin(I, JHiInt, K, nx_tot, ny_tot)];
            }
    }
    // Z faces (always physical)
    for (int g = 0; g < ng; ++g)
    {
        const int KLo = g, KLoInt = ng + g;
        const int KHi = (nz_tot - ng) + g, KHiInt = (nz_tot - ng - 1) - g;
        for (int J = 0; J < ny_tot; ++J)
            for (int I = 0; I < nx_tot; ++I)
            {
                a[lin(I, J, KLo, nx_tot, ny_tot)] = a[lin(I, J, KLoInt, nx_tot, ny_tot)];
                a[lin(I, J, KHi, nx_tot, ny_tot)] = a[lin(I, J, KHiInt, nx_tot, ny_tot)];
            }
    }
}

TEST_CASE(
    "RBGS diffusion: red/black sweeps with halo + physical BC exchanges reduce residual globally",
    "[mpi][fluids][orchestration][rbgs]")
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank = 0, size = 1;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    const int ng = 1, ny = 8, nz = 6;
    const int nx_global = 24;
    if (nx_global % size != 0)
    {
        SUCCEED("Skipping: nx_global not divisible by size");
        return;
    }
    const int nx = nx_global / size;
    const int nx_tot = nx + 2 * ng, ny_tot = ny + 2 * ng, nz_tot = nz + 2 * ng;

    const double dx = 1.0, dy = 1.0, dz = 1.0, dt = 0.2;

    const size_t N = static_cast<size_t>(nx_tot) * ny_tot * nz_tot;
    std::vector<double> u(N), v(N), w(N), u_rhs(N), v_rhs(N), w_rhs(N), nu_eff(N, 1.0e-3);

    // Rank-distinct pattern to ensure cross-rank dependence
    for (int K = 0; K < nz_tot; ++K)
        for (int J = 0; J < ny_tot; ++J)
            for (int I = 0; I < nx_tot; ++I)
            {
                const size_t q = lin(I, J, K, nx_tot, ny_tot);
                const double base = 0.1 * I + 0.05 * J - 0.02 * K + 0.5 * rank;
                u[q] = v[q] = w[q] = base;
                u_rhs[q] = v_rhs[q] = w_rhs[q] = base;
            }

    // Initial physical BC mirror + inter-rank exchange like the plugin does
    mirror_physical(u.data(), nx_tot, ny_tot, nz_tot, ng, rank, size);
    mirror_physical(v.data(), nx_tot, ny_tot, nz_tot, ng, rank, size);
    mirror_physical(w.data(), nx_tot, ny_tot, nz_tot, ng, rank, size);
    exchange_x(u.data(), nx_tot, ny_tot, nz_tot, ng, comm);
    exchange_x(v.data(), nx_tot, ny_tot, nz_tot, ng, comm);
    exchange_x(w.data(), nx_tot, ny_tot, nz_tot, ng, comm);
    mirror_physical(u.data(), nx_tot, ny_tot, nz_tot, ng, rank, size);
    mirror_physical(v.data(), nx_tot, ny_tot, nz_tot, ng, rank, size);
    mirror_physical(w.data(), nx_tot, ny_tot, nz_tot, ng, rank, size);

    double res2_before = 0, rhs2_before = 0;
    diffuse_velocity_be_residual_c(u_rhs.data(), v_rhs.data(), w_rhs.data(), u.data(), v.data(),
                                   w.data(), nu_eff.data(), nx_tot, ny_tot, nz_tot, ng, dx, dy, dz,
                                   dt, &res2_before, &rhs2_before);

    // RED sweep -> exchange -> mirror -> BLACK sweep -> exchange -> mirror
    diffuse_velocity_be_gs_color_c(u.data(), v.data(), w.data(), u_rhs.data(), v_rhs.data(),
                                   w_rhs.data(), nu_eff.data(), nx_tot, ny_tot, nz_tot, ng, dx, dy,
                                   dz, dt, /*color=*/0);
    exchange_x(u.data(), nx_tot, ny_tot, nz_tot, ng, comm);
    exchange_x(v.data(), nx_tot, ny_tot, nz_tot, ng, comm);
    exchange_x(w.data(), nx_tot, ny_tot, nz_tot, ng, comm);
    mirror_physical(u.data(), nx_tot, ny_tot, nz_tot, ng, rank, size);
    mirror_physical(v.data(), nx_tot, ny_tot, nz_tot, ng, rank, size);
    mirror_physical(w.data(), nx_tot, ny_tot, nz_tot, ng, rank, size);

    diffuse_velocity_be_gs_color_c(u.data(), v.data(), w.data(), u_rhs.data(), v_rhs.data(),
                                   w_rhs.data(), nu_eff.data(), nx_tot, ny_tot, nz_tot, ng, dx, dy,
                                   dz, dt, /*color=*/1);
    exchange_x(u.data(), nx_tot, ny_tot, nz_tot, ng, comm);
    exchange_x(v.data(), nx_tot, ny_tot, nz_tot, ng, comm);
    exchange_x(w.data(), nx_tot, ny_tot, nz_tot, ng, comm);
    mirror_physical(u.data(), nx_tot, ny_tot, nz_tot, ng, rank, size);
    mirror_physical(v.data(), nx_tot, ny_tot, nz_tot, ng, rank, size);
    mirror_physical(w.data(), nx_tot, ny_tot, nz_tot, ng, rank, size);

    double res2_after = 0, rhs2_after = 0;
    diffuse_velocity_be_residual_c(u_rhs.data(), v_rhs.data(), w_rhs.data(), u.data(), v.data(),
                                   w.data(), nu_eff.data(), nx_tot, ny_tot, nz_tot, ng, dx, dy, dz,
                                   dt, &res2_after, &rhs2_after);

    REQUIRE(rhs2_after == Approx(rhs2_before));
    REQUIRE(res2_after < res2_before);
}
