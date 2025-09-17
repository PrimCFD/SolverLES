#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <vector>
#include "kernels_fluids.h"

using Catch::Approx;

TEST_CASE("div(grad p) matches discrete Laplacian for quadratic p", "[fluids][grad][div]")
{
    const int nx = 10, ny = 8, nz = 6, ng = 2;
    const int nx_tot = nx + 2 * ng, ny_tot = ny + 2 * ng, nz_tot = nz + 2 * ng;
    const double dx = 1.0, dy = 1.0, dz = 1.0;

    std::vector<double> p(nx_tot * ny_tot * nz_tot, 0.0), dpx(p.size(), 0.0), dpy(p.size(), 0.0),
        dpz(p.size(), 0.0), lap(p.size(), 0.0);

    auto idx = [=](int I, int J, int K) { return I + nx_tot * (J + ny_tot * K); };

    // Quadratic: p = A i^2 + B j^2 + C k^2  (indices stand in for coords with dx=dy=dz=1)
    const double A = 1.5, B = 0.5, C = 0.25;
    for (int K = 0; K < nz_tot; ++K)
        for (int J = 0; J < ny_tot; ++J)
            for (int I = 0; I < nx_tot; ++I)
            {
                p[idx(I, J, K)] = A * I * I + B * J * J + C * K * K;
            }

    gradp_c(p.data(), nx_tot, ny_tot, nz_tot, ng, dx, dy, dz, dpx.data(), dpy.data(), dpz.data());
    divergence_c(dpx.data(), dpy.data(), dpz.data(), nx_tot, ny_tot, nz_tot, ng, dx, dy, dz,
                 lap.data());

    // For quadratic with uniform spacing, discrete Laplacian is constant: 2A+2B+2C
    const double expect = 2.0 * (A + B + C);

    for (int K = ng + 1; K < nz + ng - 1; ++K)
        for (int J = ng + 1; J < ny + ng - 1; ++J)
            for (int I = ng + 1; I < nx + ng - 1; ++I)
            {
                REQUIRE(lap[idx(I, J, K)] == Approx(expect).epsilon(1e-12));
            }
}
