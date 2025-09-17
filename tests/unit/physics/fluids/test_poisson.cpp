
#include <algorithm>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <vector>
#include "kernels_fluids.h"

using Catch::Approx;

static inline size_t idx(int I, int J, int K, int nx_tot, int ny_tot)
{
    return static_cast<size_t>(I + nx_tot * (J + ny_tot * K));
}

static double l2_residual_laplacian(const std::vector<double>& p, const std::vector<double>& rhs,
                                    int nx_tot, int ny_tot, int nz_tot, int ng, double dx,
                                    double dy, double dz)
{
    const double ax = 1.0 / (dx * dx);
    const double ay = 1.0 / (dy * dy);
    const double az = 1.0 / (dz * dz);
    double sum = 0.0;
    for (int K = ng; K < nz_tot - ng; ++K)
        for (int J = ng; J < ny_tot - ng; ++J)
            for (int I = ng; I < nx_tot - ng; ++I)
            {
                size_t c = idx(I, J, K, nx_tot, ny_tot);
                size_t ip = idx(I + 1, J, K, nx_tot, ny_tot);
                size_t im = idx(I - 1, J, K, nx_tot, ny_tot);
                size_t jp = idx(I, J + 1, K, nx_tot, ny_tot);
                size_t jm = idx(I, J - 1, K, nx_tot, ny_tot);
                size_t kp = idx(I, J, K + 1, nx_tot, ny_tot);
                size_t km = idx(I, J, K - 1, nx_tot, ny_tot);
                const double lap = ax * (p[ip] - 2 * p[c] + p[im]) +
                                   ay * (p[jp] - 2 * p[c] + p[jm]) +
                                   az * (p[kp] - 2 * p[c] + p[km]);
                const double r = lap - rhs[c];
                sum += r * r;
            }
    return std::sqrt(sum);
}

TEST_CASE("Jacobi Poisson: zero RHS keeps p identically zero", "[fluids][poisson]")
{
    const int nx = 6, ny = 5, nz = 4, ng = 1;
    const int nx_tot = nx + 2 * ng, ny_tot = ny + 2 * ng, nz_tot = nz + 2 * ng;
    const double dx = 1.0, dy = 1.0, dz = 1.0;
    const int iters = 7;

    const size_t N = static_cast<size_t>(nx_tot) * ny_tot * nz_tot;
    std::vector<double> rhs(N, 0.0), p(N, 0.0);

    poisson_jacobi_c(rhs.data(), nx_tot, ny_tot, nz_tot, ng, dx, dy, dz, iters, p.data());

    for (size_t q = 0; q < N; ++q)
    {
        REQUIRE(p[q] == Approx(0.0).margin(0.0));
    }
}

TEST_CASE("Jacobi Poisson reduces residual for point source", "[fluids][poisson]")
{
    const int nx = 10, ny = 8, nz = 6, ng = 1;
    const int nx_tot = nx + 2 * ng, ny_tot = ny + 2 * ng, nz_tot = nz + 2 * ng;
    const double dx = 1.0, dy = 1.0, dz = 1.0;
    const int iters = 30;

    const size_t N = static_cast<size_t>(nx_tot) * ny_tot * nz_tot;
    std::vector<double> rhs(N, 0.0), p(N, 0.0);

    // Place a unit source at a central interior cell
    const int I0 = nx / 2 + ng, J0 = ny / 2 + ng, K0 = nz / 2 + ng;
    rhs[idx(I0, J0, K0, nx_tot, ny_tot)] = 1.0;

    const double r0 = l2_residual_laplacian(p, rhs, nx_tot, ny_tot, nz_tot, ng, dx, dy, dz);

    poisson_jacobi_c(rhs.data(), nx_tot, ny_tot, nz_tot, ng, dx, dy, dz, iters, p.data());

    const double r1 = l2_residual_laplacian(p, rhs, nx_tot, ny_tot, nz_tot, ng, dx, dy, dz);

    // Expect a meaningful reduction
    REQUIRE(r1 < 0.7 * r0);
}
