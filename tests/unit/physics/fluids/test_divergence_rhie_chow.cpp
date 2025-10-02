// test_divergence_rhie_chow.cpp
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

TEST_CASE("Rhie–Chow divergence equals central divergence for uniform p", "[fluids][rc][div]")
{
    const int nx = 12, ny = 8, nz = 6, ng = 1;
    const int nx_tot = nx + 2 * ng, ny_tot = ny + 2 * ng, nz_tot = nz + 2 * ng;
    const double dx = 1.0, dy = 1.0, dz = 1.0, dt = 0.5;

    const size_t N = (size_t) nx_tot * ny_tot * nz_tot;
    std::vector<double> u(N, 0.1), v(N, -0.07), w(N, 0.03), p(N, 2.5), div_c(N, 0.0),
        div_rc(N, 0.0), rho(N, 1.7);

    divergence_c(u.data(), v.data(), w.data(), nx_tot, ny_tot, nz_tot, ng, dx, dy, dz,
                 div_c.data());
    divergence_rhie_chow_c(u.data(), v.data(), w.data(), p.data(), rho.data(), nx_tot, ny_tot,
                           nz_tot, ng, dx, dy, dz, dt, div_rc.data());

    for (int K = ng; K < nz + ng; ++K)
        for (int J = ng; J < ny + ng; ++J)
            for (int I = ng; I < nx + ng; ++I)
            {
                const auto q = idx(I, J, K, nx_tot, ny_tot);
                REQUIRE(div_rc[q] == Approx(div_c[q]).margin(1e-14));
            }
}

TEST_CASE("Rhie–Chow breaks checkerboard pressure null space", "[fluids][rc][div][checker]")
{
    const int nx = 18, ny = 10, nz = 8, ng = 1;
    const int nx_tot = nx + 2 * ng, ny_tot = ny + 2 * ng, nz_tot = nz + 2 * ng;
    const double dx = 1.0, dy = 1.0, dz = 1.0, dt = 0.2;
    const size_t N = (size_t) nx_tot * ny_tot * nz_tot;

    std::vector<double> u(N, 0.0), v(N, 0.0), w(N, 0.0), div_c(N, 0.0), div_rc(N, 0.0), p(N, 0.0);

    auto parity = [&](int I, int J, int K) { return (I + J + K) & 1; };
    for (int K = 0; K < nz_tot; ++K)
        for (int J = 0; J < ny_tot; ++J)
            for (int I = 0; I < nx_tot; ++I)
                p[idx(I, J, K, nx_tot, ny_tot)] = parity(I, J, K) ? +1.0 : -1.0;

    // central divergence of zeros is zero
    divergence_c(u.data(), v.data(), w.data(), nx_tot, ny_tot, nz_tot, ng, dx, dy, dz,
                 div_c.data());

    // Case A: constant-ρ field
    std::vector<double> rhoA(N, 1.0);
    divergence_rhie_chow_c(u.data(), v.data(), w.data(), p.data(), rhoA.data(), nx_tot, ny_tot,
                           nz_tot, ng, dx, dy, dz, dt, div_rc.data());
    double norm_c = 0.0, norm_rc = 0.0;
    for (int K = ng; K < nz + ng; ++K)
        for (int J = ng; J < ny + ng; ++J)
            for (int I = ng; I < nx + ng; ++I)
            {
                const auto q = idx(I, J, K, nx_tot, ny_tot);
                norm_c += div_c[q] * div_c[q];
                norm_rc += div_rc[q] * div_rc[q];
            }
    REQUIRE(norm_c == Approx(0.0).margin(0));
    REQUIRE(norm_rc > 0.0);

    // Case B: spatially varying ρ field
    std::vector<double> rhoB(N, 1.0);
    for (int K = 0; K < nz_tot; ++K)
        for (int J = 0; J < ny_tot; ++J)
            for (int I = 0; I < nx_tot; ++I)
            {
                const auto q = idx(I, J, K, nx_tot, ny_tot);
                rhoB[q] = 1.0 + 0.25 * std::sin(0.1 * I) * std::cos(0.2 * J) * std::cos(0.3 * K);
            }
    std::fill(div_rc.begin(), div_rc.end(), 0.0);
    divergence_rhie_chow_c(u.data(), v.data(), w.data(), p.data(), rhoB.data(), nx_tot, ny_tot,
                           nz_tot, ng, dx, dy, dz, dt, div_rc.data());
    double norm_rc_var = 0.0;
    for (int K = ng; K < nz + ng; ++K)
        for (int J = ng; J < ny + ng; ++J)
            for (int I = ng; I < nx + ng; ++I)
                norm_rc_var +=
                    div_rc[idx(I, J, K, nx_tot, ny_tot)] * div_rc[idx(I, J, K, nx_tot, ny_tot)];
    REQUIRE(norm_rc_var > 0.0);
}
