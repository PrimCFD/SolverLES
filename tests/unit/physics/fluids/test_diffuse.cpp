#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <vector>
#include "kernels_fluids.h"

using Catch::Approx;

TEST_CASE("Diffusion step preserves constant velocity field (interior)", "[fluids][diffuse]")
{
    const int nx = 12, ny = 7, nz = 5, ng = 2;
    const int nx_tot = nx + 2 * ng, ny_tot = ny + 2 * ng, nz_tot = nz + 2 * ng;
    const double dx = 1.0, dy = 1.0, dz = 1.0, dt = 0.3;

    std::vector<double> u(nx_tot * ny_tot * nz_tot, 3.14), v(nx_tot * ny_tot * nz_tot, 3.14),
        w(nx_tot * ny_tot * nz_tot, 3.14), nu_eff(nx_tot * ny_tot * nz_tot, 1.0e-3),
        us(u.size(), 0.0), vs(v.size(), 0.0), ws(w.size(), 0.0);

    diffuse_velocity_fe_c(u.data(), v.data(), w.data(), nu_eff.data(), nx_tot, ny_tot, nz_tot, ng,
                          dx, dy, dz, dt, us.data(), vs.data(), ws.data());

    auto idx = [=](int I, int J, int K) { return I + nx_tot * (J + ny_tot * K); };

    // Check interior only; ghosts are owned by BC/exchange elsewhere.
    for (int K = ng; K < nz + ng; ++K)
        for (int J = ng; J < ny + ng; ++J)
            for (int I = ng; I < nx + ng; ++I)
            {
                REQUIRE(us[idx(I, J, K)] == Approx(3.14));
                REQUIRE(vs[idx(I, J, K)] == Approx(3.14));
                REQUIRE(ws[idx(I, J, K)] == Approx(3.14));
            }
}
