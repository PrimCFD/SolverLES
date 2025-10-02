#include <algorithm>
#include <catch2/catch_all.hpp>
#include <vector>
#include "kernels_fluids.h"

using Catch::Approx;

// Helpers for MAC totals
static inline int center_tot(int n, int ng)
{
    return n + 2 * ng;
}
static inline int face_tot(int n, int ng)
{
    return n + 1 + 2 * ng;
}

TEST_CASE(
    "MAC correct_velocity (variable rho at centers) matches constant case when rho is uniform",
    "[mac][correct][varrho]")
{
    const int nx = 5, ny = 4, nz = 3, ng = 1;
    const double dx = 1.0, dy = 1.0, dz = 1.0;
    const double rho0 = 3.0, dt = 0.1;

    const int nxc = center_tot(nx, ng), nyc = center_tot(ny, ng), nzc = center_tot(nz, ng);
    const int nxu = face_tot(nx, ng), nyu = center_tot(ny, ng), nzu = center_tot(nz, ng);
    const int nxv = center_tot(nx, ng), nyv = face_tot(ny, ng), nzv = center_tot(nz, ng);
    const int nxw = center_tot(nx, ng), nyw = center_tot(ny, ng), nzw = face_tot(nz, ng);

    const std::size_t Nu = (std::size_t) nxu * nyu * nzu;
    const std::size_t Nv = (std::size_t) nxv * nyv * nzv;
    const std::size_t Nw = (std::size_t) nxw * nyw * nzw;
    const std::size_t Np = (std::size_t) nxc * nyc * nzc;

    // Simple initial velocities
    std::vector<double> u1(Nu, 0.25), v1(Nv, -0.4), w1(Nw, 0.6);
    std::vector<double> u2 = u1, v2 = v1, w2 = w1;

    // Pressure at centers: linear field -> constant gradients
    std::vector<double> p(Np, 0.0);
    auto idxC = [&](int i, int j, int k)
    {
        return (std::size_t) i +
               (std::size_t) nxc * ((std::size_t) j + (std::size_t) nyc * (std::size_t) k);
    };
    for (int k = 0; k < nzc; ++k)
        for (int j = 0; j < nyc; ++j)
            for (int i = 0; i < nxc; ++i)
                p[idxC(i, j, k)] = 1.5 * (i - ng) - 0.75 * (j - ng) + 0.25 * (k - ng);

    // Face gradients
    std::vector<double> dpx_u(Nu, 0.0), dpy_v(Nv, 0.0), dpz_w(Nw, 0.0);
    gradp_faces_c(p.data(), nxc, nyc, nzc, ng, dx, dy, dz, dpx_u.data(), nxu, nyu, nzu,
                  dpy_v.data(), nxv, nyv, nzv, dpz_w.data(), nxw, nyw, nzw);

    // Uniform rho at centers for this test (exercises varrho code path with trivial averaging)
    std::vector<double> rho_c(Np, rho0);

    // A) constant-ρ path
    correct_velocity_mac_c(u1.data(), v1.data(), w1.data(), dpx_u.data(), dpy_v.data(),
                           dpz_w.data(), nxu, nyu, nzu, nxv, nyv, nzv, nxw, nyw, nzw, nxc, nyc, nzc,
                           ng, rho0, dt);

    // B) variable-ρ path
    correct_velocity_varrho_mac_c(u2.data(), v2.data(), w2.data(), dpx_u.data(), dpy_v.data(),
                                  dpz_w.data(), nxu, nyu, nzu, nxv, nyv, nzv, nxw, nyw, nzw, nxc,
                                  nyc, nzc, ng, rho_c.data(), dt);

    // With uniform ρ, both paths must match bitwise (allow tiny FP slack)
    for (std::size_t i = 0; i < Nu; ++i)
        REQUIRE(u1[i] == Approx(u2[i]));
    for (std::size_t i = 0; i < Nv; ++i)
        REQUIRE(v1[i] == Approx(v2[i]));
    for (std::size_t i = 0; i < Nw; ++i)
        REQUIRE(w1[i] == Approx(w2[i]));
}
