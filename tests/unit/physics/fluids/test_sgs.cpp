#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <vector>
#include "kernels_fluids.h"

using Catch::Approx;

TEST_CASE("Smagorinsky SGS on constant-gradient velocity", "[fluids][sgs]")
{
    const int nx = 8, ny = 6, nz = 4, ng = 2;
    const int nx_tot = nx + 2 * ng, ny_tot = ny + 2 * ng, nz_tot = nz + 2 * ng;
    const double dx = 1.0, dy = 1.0, dz = 1.0, Cs = 0.16;

    std::vector<double> u(nx_tot * ny_tot * nz_tot, 0.0), v(nx_tot * ny_tot * nz_tot, 0.0),
        w(nx_tot * ny_tot * nz_tot, 0.0), nu_t(nx_tot * ny_tot * nz_tot, -1.0);

    // Linear field: dudx=a, dvdy=b, dwdz=c; all cross-derivs = 0
    const double a = 1.0, b = 2.0, c = 3.0;
    auto idx = [=](int I, int J, int K) { return I + nx_tot * (J + ny_tot * K); };

    for (int K = 0; K < nz_tot; ++K)
        for (int J = 0; J < ny_tot; ++J)
            for (int I = 0; I < nx_tot; ++I)
            {
                u[idx(I, J, K)] = a * dx * I;
                v[idx(I, J, K)] = b * dy * J;
                w[idx(I, J, K)] = c * dz * K;
            }

    sgs_smagorinsky_c(u.data(), v.data(), w.data(), nx_tot, ny_tot, nz_tot, ng, dx, dy, dz, Cs,
                      nu_t.data());

    // Expected Smag: sqrt(2*(a^2+b^2+c^2)); Δ=(dx*dy*dz)^(1/3)=1
    const double Smag = std::sqrt(2.0 * (a * a + b * b + c * c));
    const double expect = (Cs * Cs) * Smag;

    // Check a few interior points away from ghosts
    for (int K = ng + 1; K < nz + ng - 1; ++K)
        for (int J = ng + 1; J < ny + ng - 1; ++J)
            for (int I = ng + 1; I < nx + ng - 1; ++I)
            {
                REQUIRE(nu_t[idx(I, J, K)] == Approx(expect).epsilon(1e-12));
            }
}
