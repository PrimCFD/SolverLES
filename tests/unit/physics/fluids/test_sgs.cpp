#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <vector>
#include "MacOps.hpp"
using namespace numerics::kernels;

using Catch::Approx;

TEST_CASE("Smagorinsky SGS on constant-gradient velocity (MAC)", "[fluids][sgs]")
{
    const int nx = 8, ny = 6, nz = 4, ng = 2;
    const double dx = 1.0, dy = 1.0, dz = 1.0, Cs = 0.16;

    // centers
    const int nxc_tot = nx + 2 * ng, nyc_tot = ny + 2 * ng, nzc_tot = nz + 2 * ng;

    // face totals (MAC)
    const int nxu_tot = (nx + 1) + 2 * ng, nyu_tot = ny + 2 * ng, nzu_tot = nz + 2 * ng;
    const int nxv_tot = nx + 2 * ng, nyv_tot = (ny + 1) + 2 * ng, nzv_tot = nz + 2 * ng;
    const int nxw_tot = nx + 2 * ng, nyw_tot = ny + 2 * ng, nzw_tot = (nz + 1) + 2 * ng;

    auto idx3 = [](int i, int nx, int j, int ny, int k) { return i + nx * (j + ny * k); };

    std::vector<double> u(nxu_tot * nyu_tot * nzu_tot, 0.0);
    std::vector<double> v(nxv_tot * nyv_tot * nzv_tot, 0.0);
    std::vector<double> w(nxw_tot * nyw_tot * nzw_tot, 0.0);
    std::vector<double> nu_t(nxc_tot * nyc_tot * nzc_tot, -1.0);

    // Target gradients
    const double a = 1.0, b = 2.0, c = 3.0;

    // Fill as linear functions at FACE locations
    for (int K = 0; K < nzu_tot; ++K)
        for (int J = 0; J < nyu_tot; ++J)
            for (int I = 0; I < nxu_tot; ++I)
                u[idx3(I, nxu_tot, J, nyu_tot, K)] = a * dx * ((I - ng) - 0.5); // x-face

    for (int K = 0; K < nzv_tot; ++K)
        for (int J = 0; J < nyv_tot; ++J)
            for (int I = 0; I < nxv_tot; ++I)
                v[idx3(I, nxv_tot, J, nyv_tot, K)] = b * dy * ((J - ng) - 0.5); // y-face

    for (int K = 0; K < nzw_tot; ++K)
        for (int J = 0; J < nyw_tot; ++J)
            for (int I = 0; I < nxw_tot; ++I)
                w[idx3(I, nxw_tot, J, nyw_tot, K)] = c * dz * ((K - ng) - 0.5); // z-face

    // Call new MAC kernel signature
    sgs_smagorinsky(u.data(), v.data(), w.data(), nxu_tot, nyu_tot, nzu_tot, nxv_tot, nyv_tot,
                          nzv_tot, nxw_tot, nyw_tot, nzw_tot, nxc_tot, nyc_tot, nzc_tot, ng, dx, dy,
                          dz, Cs, nu_t.data());

    // Expected Smag remains √[2(a²+b²+c²)], Δ=(dx*dy*dz)^(1/3)=1
    const double Smag = std::sqrt(2.0 * (a * a + b * b + c * c));
    const double expect = (Cs * Cs) * Smag;

    auto cidx = [&](int I, int J, int K) { return idx3(I, nxc_tot, J, nyc_tot, K); };

    // Check interior centers
    for (int K = ng + 1; K < nz + ng - 1; ++K)
        for (int J = ng + 1; J < ny + ng - 1; ++J)
            for (int I = ng + 1; I < nx + ng - 1; ++I)
                REQUIRE(nu_t[cidx(I, J, K)] == Catch::Approx(expect).epsilon(1e-12));
}
