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

TEST_CASE("Pressure gradient correction updates velocities exactly", "[fluids][corrector]")
{
    const int nx = 7, ny = 5, nz = 4, ng = 1;
    const int nx_tot = nx + 2 * ng, ny_tot = ny + 2 * ng, nz_tot = nz + 2 * ng;
    const double rho = 1.7, dt = 0.3;

    const size_t N = static_cast<size_t>(nx_tot) * ny_tot * nz_tot;

    std::vector<double> u(N), v(N), w(N);
    std::vector<double> dpx(N), dpy(N), dpz(N);

    // Seed fields with simple, deterministic patterns (not constant)
    for (int K = 0; K < nz_tot; ++K)
        for (int J = 0; J < ny_tot; ++J)
            for (int I = 0; I < nx_tot; ++I)
            {
                const size_t q = idx(I, J, K, nx_tot, ny_tot);
                u[q] = 0.1 * I + 0.2 * J + 0.3 * K;
                v[q] = 1.0 + 0.05 * I - 0.04 * J + 0.02 * K;
                w[q] = -0.3 + 0.07 * I + 0.01 * J - 0.06 * K;
                dpx[q] = 0.9 - 0.1 * I + 0.05 * J - 0.03 * K;
                dpy[q] = -0.2 + 0.02 * I + 0.03 * J + 0.01 * K;
                dpz[q] = 0.5 - 0.04 * I + 0.02 * J + 0.01 * K;
            }

    // Save originals for expected update
    auto u0 = u, v0 = v, w0 = w;

    correct_velocity_c(u.data(), v.data(), w.data(), dpx.data(), dpy.data(), dpz.data(), nx_tot,
                       ny_tot, nz_tot, ng, rho, dt);

    const double fac = dt / rho;
    bool bad = false;
    for (int K = ng; K < nz + ng && !bad; ++K)
        for (int J = ng; J < ny + ng && !bad; ++J)
            for (int I = ng; I < nx + ng; ++I)
            {
                size_t q = idx(I, J, K, nx_tot, ny_tot);
                auto exp_u = u0[q] - fac * dpx[q];
                auto exp_v = v0[q] - fac * dpy[q];
                auto exp_w = w0[q] - fac * dpz[q];

                if (std::abs(u[q] - exp_u) > 1e-12 || std::abs(v[q] - exp_v) > 1e-12 ||
                    std::abs(w[q] - exp_w) > 1e-12)
                {
                    INFO("First mismatch at I=" << I << " J=" << J << " K=" << K << " q=" << q);
                    INFO("u diff=" << (u[q] - exp_u) << " v diff=" << (v[q] - exp_v)
                                   << " w diff=" << (w[q] - exp_w));
                    bad = true;
                    break;
                }
            }
    REQUIRE(!bad);
}
