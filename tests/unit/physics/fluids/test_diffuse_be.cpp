
#include <algorithm>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <vector>
#include "kernels_fluids.h"

using Catch::Approx;

static inline size_t idx(int I, int J, int K, int nx_tot, int ny_tot)
{
    return static_cast<size_t>(I + nx_tot * (J + ny_tot * K));
}

static inline bool is_interior(int I, int J, int K, int nx_tot, int ny_tot, int nz_tot, int ng)
{
    return (I >= ng && I < nx_tot - ng && J >= ng && J < ny_tot - ng && K >= ng && K < nz_tot - ng);
}

static void copy_halos(std::vector<double>& dst, const std::vector<double>& src, int nx_tot,
                       int ny_tot, int nz_tot, int ng)
{
    for (int K = 0; K < nz_tot; ++K)
        for (int J = 0; J < ny_tot; ++J)
            for (int I = 0; I < nx_tot; ++I)
                if (!is_interior(I, J, K, nx_tot, ny_tot, nz_tot, ng))
                    dst[idx(I, J, K, nx_tot, ny_tot)] = src[idx(I, J, K, nx_tot, ny_tot)];
}

TEST_CASE("BE Jacobi sweep reduces diffusion residual", "[fluids][diffuse][be]")
{
    const int nx = 10, ny = 7, nz = 5, ng = 1;
    const int nx_tot = nx + 2 * ng, ny_tot = ny + 2 * ng, nz_tot = nz + 2 * ng;
    const double dx = 1.0, dy = 1.0, dz = 1.0, dt = 0.5;

    const size_t N = static_cast<size_t>(nx_tot) * ny_tot * nz_tot;
    std::vector<double> u_rhs(N), v_rhs(N), w_rhs(N), u_iter(N), v_iter(N), w_iter(N),
        u_next(N, 0.0), v_next(N, 0.0), w_next(N, 0.0), nu_eff(N, 1.0e-2);

    // Nontrivial initial iterate / RHS
    for (int K = 0; K < nz_tot; ++K)
        for (int J = 0; J < ny_tot; ++J)
            for (int I = 0; I < nx_tot; ++I)
            {
                const size_t q = idx(I, J, K, nx_tot, ny_tot);
                u_rhs[q] = 0.1 * I + 0.05 * J - 0.02 * K;
                v_rhs[q] = -0.3 + 0.03 * I + 0.04 * J + 0.02 * K;
                w_rhs[q] = 0.2 - 0.01 * I + 0.02 * J - 0.03 * K;
                u_iter[q] = 0.0;
                v_iter[q] = 0.0;
                w_iter[q] = 0.0;
            }

    double res2_before = -1.0, rhs2_before = -1.0;
    diffuse_velocity_be_residual_c(u_rhs.data(), v_rhs.data(), w_rhs.data(), u_iter.data(),
                                   v_iter.data(), w_iter.data(), nu_eff.data(), nx_tot, ny_tot,
                                   nz_tot, ng, dx, dy, dz, dt, &res2_before, &rhs2_before);

    diffuse_velocity_be_sweep_c(u_rhs.data(), v_rhs.data(), w_rhs.data(), u_iter.data(),
                                v_iter.data(), w_iter.data(), nu_eff.data(), nx_tot, ny_tot, nz_tot,
                                ng, dx, dy, dz, dt, u_next.data(), v_next.data(), w_next.data());

    copy_halos(u_next, u_iter, nx_tot, ny_tot, nz_tot, ng);
    copy_halos(v_next, v_iter, nx_tot, ny_tot, nz_tot, ng);
    copy_halos(w_next, w_iter, nx_tot, ny_tot, nz_tot, ng);

    double res2_after = -1.0, rhs2_after = -1.0;
    diffuse_velocity_be_residual_c(u_rhs.data(), v_rhs.data(), w_rhs.data(), u_next.data(),
                                   v_next.data(), w_next.data(), nu_eff.data(), nx_tot, ny_tot,
                                   nz_tot, ng, dx, dy, dz, dt, &res2_after, &rhs2_after);

    REQUIRE(rhs2_after == Approx(rhs2_before));
    REQUIRE(res2_after < res2_before);
}

TEST_CASE("Red-black GS updates only the requested color", "[fluids][diffuse][gs]")
{
    const int nx = 9, ny = 7, nz = 5, ng = 1;
    const int nx_tot = nx + 2 * ng, ny_tot = ny + 2 * ng, nz_tot = nz + 2 * ng;
    const double dx = 1.0, dy = 1.0, dz = 1.0, dt = 0.1;

    const size_t N = static_cast<size_t>(nx_tot) * ny_tot * nz_tot;
    std::vector<double> u(N), v(N), w(N), u_rhs(N), v_rhs(N), w_rhs(N), nu_eff(N, 1.0e-3);

    // Initialize with a strongly varying pattern so some updates are non-trivial
    for (int K = 0; K < nz_tot; ++K)
        for (int J = 0; J < ny_tot; ++J)
            for (int I = 0; I < nx_tot; ++I)
            {
                const size_t q = idx(I, J, K, nx_tot, ny_tot);
                u[q] = 1.0 + I + 10 * J + 100 * K;
                v[q] = -2.0 + 2 * I - 5 * J + 3 * K;
                w[q] = 0.5 - I + 4 * J - 2 * K;
                u_rhs[q] = u[q];
                v_rhs[q] = v[q];
                w_rhs[q] = w[q];
            }

    // Keep a copy for comparison
    auto u_before = u, v_before = v, w_before = w;

    // Update color 0 (red)
    diffuse_velocity_be_gs_color_c(u.data(), v.data(), w.data(), u_rhs.data(), v_rhs.data(),
                                   w_rhs.data(), nu_eff.data(), nx_tot, ny_tot, nz_tot, ng, dx, dy,
                                   dz, dt,
                                   /*color=*/0);

    // Check parity behavior
    bool any_updated = false;
    for (int K = 1; K < nz_tot - 1; ++K)
        for (int J = 1; J < ny_tot - 1; ++J)
            for (int I = 1; I < nx_tot - 1; ++I)
            {
                const int parity = (I - ng) + (J - ng) + (K - ng);
                const size_t q = idx(I, J, K, nx_tot, ny_tot);
                if (parity % 2 == 0)
                {
                    if (std::abs(u[q] - u_before[q]) > 1e-15 ||
                        std::abs(v[q] - v_before[q]) > 1e-15 ||
                        std::abs(w[q] - w_before[q]) > 1e-15)
                    {
                        any_updated = true; // red cells should change
                    }
                }
                else
                {
                    // black cells must be untouched
                    REQUIRE(u[q] == Approx(u_before[q]).margin(0));
                    REQUIRE(v[q] == Approx(v_before[q]).margin(0));
                    REQUIRE(w[q] == Approx(w_before[q]).margin(0));
                }
            }
    REQUIRE(any_updated); // ensure the test actually exercised updates
}
