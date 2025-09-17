#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <mpi.h>
#include <vector>
#include "kernels_fluids.h"

static inline size_t lin(int I, int J, int K, int nx_tot, int ny_tot)
{
    return static_cast<size_t>(I + nx_tot * (J + ny_tot * K));
}

static void exchange_x(double* a, int nx_tot, int ny_tot, int nz_tot, int ng, MPI_Comm comm)
{
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    const int left = (rank == 0) ? MPI_PROC_NULL : rank - 1;
    const int right = (rank == size - 1) ? MPI_PROC_NULL : rank + 1;
    const int nyz = ny_tot * nz_tot, slab = ng * nyz;
    std::vector<double> sendL(slab), sendR(slab), recvL(slab), recvR(slab);

    for (int g = 0; g < ng; ++g)
    {
        const int IL = ng + g;
        const int IR = (nx_tot - 2 * ng) + g;
        for (int K = 0; K < nz_tot; ++K)
            for (int J = 0; J < ny_tot; ++J)
            {
                const size_t s = g * nyz + J + ny_tot * K;
                sendL[s] = a[lin(IL, J, K, nx_tot, ny_tot)];
                sendR[s] = a[lin(IR, J, K, nx_tot, ny_tot)];
            }
    }

    MPI_Sendrecv(sendL.data(), slab, MPI_DOUBLE, left, 101, recvR.data(), slab, MPI_DOUBLE, right,
                 101, comm, MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendR.data(), slab, MPI_DOUBLE, right, 201, recvL.data(), slab, MPI_DOUBLE, left,
                 201, comm, MPI_STATUS_IGNORE);

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

static void mirror_physical(double* a, int nx_tot, int ny_tot, int nz_tot, int ng, int rank,
                            int size)
{
    if (rank == 0)
    {
        for (int g = 0; g < ng; ++g)
        {
            const int Ihalo = g, Iint = ng + g;
            for (int K = 0; K < nz_tot; ++K)
                for (int J = 0; J < ny_tot; ++J)
                    a[lin(Ihalo, J, K, nx_tot, ny_tot)] = a[lin(Iint, J, K, nx_tot, ny_tot)];
        }
    }
    if (rank == size - 1)
    {
        for (int g = 0; g < ng; ++g)
        {
            const int Ihalo = (nx_tot - ng) + g, Iint = (nx_tot - ng - 1) - g;
            for (int K = 0; K < nz_tot; ++K)
                for (int J = 0; J < ny_tot; ++J)
                    a[lin(Ihalo, J, K, nx_tot, ny_tot)] = a[lin(Iint, J, K, nx_tot, ny_tot)];
        }
    }
    for (int g = 0; g < ng; ++g)
    {
        const int JLo = g, JLoInt = ng + g, JHi = (ny_tot - ng) + g, JHiInt = (ny_tot - ng - 1) - g;
        for (int K = 0; K < nz_tot; ++K)
            for (int I = 0; I < nx_tot; ++I)
            {
                a[lin(I, JLo, K, nx_tot, ny_tot)] = a[lin(I, JLoInt, K, nx_tot, ny_tot)];
                a[lin(I, JHi, K, nx_tot, ny_tot)] = a[lin(I, JHiInt, K, nx_tot, ny_tot)];
            }
    }
    for (int g = 0; g < ng; ++g)
    {
        const int KLo = g, KLoInt = ng + g, KHi = (nz_tot - ng) + g, KHiInt = (nz_tot - ng - 1) - g;
        for (int J = 0; J < ny_tot; ++J)
            for (int I = 0; I < nx_tot; ++I)
            {
                a[lin(I, J, KLo, nx_tot, ny_tot)] = a[lin(I, J, KLoInt, nx_tot, ny_tot)];
                a[lin(I, J, KHi, nx_tot, ny_tot)] = a[lin(I, J, KHiInt, nx_tot, ny_tot)];
            }
    }
}

static double l2_residual(const std::vector<double>& p, const std::vector<double>& rhs, int nx_tot,
                          int ny_tot, int nz_tot, int ng, double dx, double dy, double dz)
{
    const double ax = 1.0 / (dx * dx), ay = 1.0 / (dy * dy), az = 1.0 / (dz * dz);
    double sum = 0.0;
    for (int K = ng; K < nz_tot - ng; ++K)
        for (int J = ng; J < ny_tot - ng; ++J)
            for (int I = ng; I < nx_tot - ng; ++I)
            {
                const size_t c = lin(I, J, K, nx_tot, ny_tot);
                const size_t ip = lin(I + 1, J, K, nx_tot, ny_tot),
                             im = lin(I - 1, J, K, nx_tot, ny_tot);
                const size_t jp = lin(I, J + 1, K, nx_tot, ny_tot),
                             jm = lin(I, J - 1, K, nx_tot, ny_tot);
                const size_t kp = lin(I, J, K + 1, nx_tot, ny_tot),
                             km = lin(I, J, K - 1, nx_tot, ny_tot);
                const double lap = ax * (p[ip] - 2 * p[c] + p[im]) +
                                   ay * (p[jp] - 2 * p[c] + p[jm]) +
                                   az * (p[kp] - 2 * p[c] + p[km]);
                const double r = lap - rhs[c];
                sum += r * r;
            }
    return sum;
}

TEST_CASE("Poisson orchestration: per-sweep exchanges + physical BC mirror + global mean removal",
          "[mpi][fluids][poisson]")
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank = 0, size = 1;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    const int ng = 1, ny = 10, nz = 8;
    const int nx_global = 20;
    if (nx_global % size != 0)
    {
        SUCCEED("Skipping: nx_global not divisible by size");
        return;
    }
    const int nx = nx_global / size;
    const int nx_tot = nx + 2 * ng, ny_tot = ny + 2 * ng, nz_tot = nz + 2 * ng;

    const double dx = 1.0, dy = 1.0, dz = 1.0;

    const size_t N = static_cast<size_t>(nx_tot) * ny_tot * nz_tot;
    std::vector<double> rhs(N, 0.0), p(N, 0.0);

    // Smooth RHS with nonzero mean (forces nullspace handling)
    for (int K = 0; K < nz_tot; ++K)
        for (int J = 0; J < ny_tot; ++J)
            for (int I = 0; I < nx_tot; ++I)
                rhs[lin(I, J, K, nx_tot, ny_tot)] = 0.2 + 0.01 * (I + J + K) + 0.05 * rank;

    // --- Make RHS globally zero-mean over the interior (compatibility condition) ---
    double rsum_local = 0.0, rcnt_local = 0.0;
    for (int K = ng; K < nz_tot - ng; ++K)
        for (int J = ng; J < ny_tot - ng; ++J)
            for (int I = ng; I < nx_tot - ng; ++I)
            {
                rsum_local += rhs[lin(I, J, K, nx_tot, ny_tot)];
                rcnt_local += 1.0;
            }
    double rsum = 0.0, rcnt = 0.0;
    MPI_Allreduce(&rsum_local, &rsum, 1, MPI_DOUBLE, MPI_SUM, comm);
    MPI_Allreduce(&rcnt_local, &rcnt, 1, MPI_DOUBLE, MPI_SUM, comm);
    const double rhs_mean = (rcnt > 0.0) ? rsum / rcnt : 0.0;
    for (int K = ng; K < nz_tot - ng; ++K)
        for (int J = ng; J < ny_tot - ng; ++J)
            for (int I = ng; I < nx_tot - ng; ++I)
                rhs[lin(I, J, K, nx_tot, ny_tot)] -= rhs_mean;

    // Ensure p halos are consistent before first residual
    mirror_physical(p.data(), nx_tot, ny_tot, nz_tot, ng, rank, size);
    exchange_x(p.data(), nx_tot, ny_tot, nz_tot, ng, comm);
    mirror_physical(p.data(), nx_tot, ny_tot, nz_tot, ng, rank, size);

    const double r0_local = l2_residual(p, rhs, nx_tot, ny_tot, nz_tot, ng, dx, dy, dz);
    double r0 = 0;
    MPI_Allreduce(&r0_local, &r0, 1, MPI_DOUBLE, MPI_SUM, comm);

    const int sweeps = 80; // a few dozen sweeps is plenty on this grid
    for (int it = 0; it < sweeps; ++it)
    {
        poisson_jacobi_c(rhs.data(), nx_tot, ny_tot, nz_tot, ng, dx, dy, dz, /*iters=*/1, p.data());
        exchange_x(p.data(), nx_tot, ny_tot, nz_tot, ng, comm);
        mirror_physical(p.data(), nx_tot, ny_tot, nz_tot, ng, rank, size);
        if ((it % 10) == 9)
        {
            const double r_local = l2_residual(p, rhs, nx_tot, ny_tot, nz_tot, ng, dx, dy, dz);
            double r = 0.0;
            MPI_Allreduce(&r_local, &r, 1, MPI_DOUBLE, MPI_SUM, comm);
            INFO("it=" << (it + 1) << "  r=" << r);
        }
    }

    // Global mean removal on interior
    double lsum = 0.0, lcnt = 0.0;
    for (int K = ng; K < nz_tot - ng; ++K)
        for (int J = ng; J < ny_tot - ng; ++J)
            for (int I = ng; I < nx_tot - ng; ++I)
            {
                lsum += p[lin(I, J, K, nx_tot, ny_tot)];
                lcnt += 1.0;
            }
    double gsum = 0.0, gcnt = 0.0;
    MPI_Allreduce(&lsum, &gsum, 1, MPI_DOUBLE, MPI_SUM, comm);
    MPI_Allreduce(&lcnt, &gcnt, 1, MPI_DOUBLE, MPI_SUM, comm);
    const double mean = (gcnt > 0.0) ? gsum / gcnt : 0.0;
    for (int K = ng; K < nz_tot - ng; ++K)
        for (int J = ng; J < ny_tot - ng; ++J)
            for (int I = ng; I < nx_tot - ng; ++I)
                p[lin(I, J, K, nx_tot, ny_tot)] -= mean;
    exchange_x(p.data(), nx_tot, ny_tot, nz_tot, ng, comm);
    mirror_physical(p.data(), nx_tot, ny_tot, nz_tot, ng, rank, size);

    const double r1_local = l2_residual(p, rhs, nx_tot, ny_tot, nz_tot, ng, dx, dy, dz);
    double r1 = 0;
    MPI_Allreduce(&r1_local, &r1, 1, MPI_DOUBLE, MPI_SUM, comm);

    REQUIRE(r1 < 0.7 * r0); // with correct halo policy the drop is very consistent
    // Mean ~ 0 (tolerate tiny rounding)
    double check_sum_local = 0.0;
    for (int K = ng; K < nz_tot - ng; ++K)
        for (int J = ng; J < ny_tot - ng; ++J)
            for (int I = ng; I < nx_tot - ng; ++I)
                check_sum_local += p[lin(I, J, K, nx_tot, ny_tot)];
    double check_sum = 0.0;
    MPI_Allreduce(&check_sum_local, &check_sum, 1, MPI_DOUBLE, MPI_SUM, comm);
    REQUIRE(std::abs(check_sum) <= 1e-12 * std::max(1.0, gcnt));
}
