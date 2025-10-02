// test_correct_velocity_varrho.cpp
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <vector>
#include "kernels_fluids.h"

using Catch::Approx;

static inline size_t idx(int I, int J, int K, int nx_tot, int ny_tot)
{
    return (size_t) I + (size_t) nx_tot * ((size_t) J + (size_t) ny_tot * (size_t) K);
}

TEST_CASE("Variable-ρ corrector matches constant-ρ when ρ is uniform",
          "[fluids][corrector][varrho]")
{
    const int nx = 9, ny = 7, nz = 6, ng = 1;
    const int nx_tot = nx + 2 * ng, ny_tot = ny + 2 * ng, nz_tot = nz + 2 * ng;
    const double dt = 0.4, rho_s = 2.0;
    const size_t N = (size_t) nx_tot * ny_tot * nz_tot;

    std::vector<double> u1(N), v1(N), w1(N), u2(N), v2(N), w2(N);
    std::vector<double> dpx(N), dpy(N), dpz(N);
    std::vector<double> rho(N, rho_s);

    for (int K = 0; K < nz_tot; ++K)
        for (int J = 0; J < ny_tot; ++J)
            for (int I = 0; I < nx_tot; ++I)
            {
                const auto q = idx(I, J, K, nx_tot, ny_tot);
                u1[q] = u2[q] = 0.11 * I - 0.07 * J + 0.02 * K;
                v1[q] = v2[q] = -0.3 + 0.05 * I + 0.01 * J - 0.02 * K;
                w1[q] = w2[q] = 0.2 - 0.03 * I + 0.02 * J + 0.04 * K;
                dpx[q] = 0.9 - 0.1 * I + 0.02 * J - 0.03 * K;
                dpy[q] = -0.2 + 0.02 * I + 0.03 * J + 0.01 * K;
                dpz[q] = 0.5 - 0.04 * I + 0.02 * J + 0.01 * K;
            }

    correct_velocity_c(u1.data(), v1.data(), w1.data(), dpx.data(), dpy.data(), dpz.data(), nx_tot,
                       ny_tot, nz_tot, ng, rho_s, dt);
    correct_velocity_varrho_c(u2.data(), v2.data(), w2.data(), dpx.data(), dpy.data(), dpz.data(),
                              nx_tot, ny_tot, nz_tot, ng, rho.data(), dt);

    for (int K = ng; K < nz + ng; ++K)
        for (int J = ng; J < ny + ng; ++J)
            for (int I = ng; I < nx + ng; ++I)
            {
                const auto q = idx(I, J, K, nx_tot, ny_tot);
                REQUIRE(u2[q] == Approx(u1[q]).margin(1e-12));
                REQUIRE(v2[q] == Approx(v1[q]).margin(1e-12));
                REQUIRE(w2[q] == Approx(w1[q]).margin(1e-12));
            }
}
