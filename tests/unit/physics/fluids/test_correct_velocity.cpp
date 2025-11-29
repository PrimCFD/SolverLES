#include <algorithm>
#include <catch2/catch_all.hpp>
#include <vector>
#include "MacOps.hpp"
using namespace numerics::kernels;

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

TEST_CASE("MAC correct_velocity (constant rho) matches variable-rho path for uniform rho",
          "[mac][correct]")
{
    const int nx = 7, ny = 6, nz = 5, ng = 1;
    const double dx = 1.0, dy = 1.0, dz = 1.0;
    const double rho0 = 2.5, dt = 0.2;

    const int nxc = center_tot(nx, ng), nyc = center_tot(ny, ng), nzc = center_tot(nz, ng);
    const int nxu = face_tot(nx, ng), nyu = center_tot(ny, ng), nzu = center_tot(nz, ng);
    const int nxv = center_tot(nx, ng), nyv = face_tot(ny, ng), nzv = center_tot(nz, ng);
    const int nxw = center_tot(nx, ng), nyw = center_tot(ny, ng), nzw = face_tot(nz, ng);

    const std::size_t Nu = (std::size_t) nxu * nyu * nzu;
    const std::size_t Nv = (std::size_t) nxv * nyv * nzv;
    const std::size_t Nw = (std::size_t) nxw * nyw * nzw;
    const std::size_t Np = (std::size_t) nxc * nyc * nzc;

    // Random-ish initial velocities (deterministic patterns)
    std::vector<double> u1(Nu), v1(Nv), w1(Nw);
    for (std::size_t i = 0; i < Nu; ++i)
        u1[i] = 0.1 + 0.001 * double(i % 101);
    for (std::size_t i = 0; i < Nv; ++i)
        v1[i] = -0.2 + 0.002 * double(i % 73);
    for (std::size_t i = 0; i < Nw; ++i)
        w1[i] = 0.3 - 0.003 * double(i % 59);

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
                p[idxC(i, j, k)] = 3.0 * (i - ng) + 2.0 * (j - ng) - 1.0 * (k - ng);

    // Face gradients
    std::vector<double> dpx_u(Nu, 0.0), dpy_v(Nv, 0.0), dpz_w(Nw, 0.0);
    grad_p_faces(p.data(), nxc, nyc, nzc, ng, dx, dy, dz, dpx_u.data(), nxu, nyu, nzu,
                  dpy_v.data(), nxv, nyv, nzv, dpz_w.data(), nxw, nyw, nzw);

    // Constant-rho corrector
    correct_velocity_const_rho(u1.data(), v1.data(), w1.data(), dpx_u.data(), dpy_v.data(),
                           dpz_w.data(), nxu, nyu, nzu, nxv, nyv, nzv, nxw, nyw, nzw, nxc, nyc, nzc,
                           ng, rho0, dt);

    // Variable-rho with uniform centers
    std::vector<double> rho_c(Np, rho0);
    correct_velocity_varrho(u2.data(), v2.data(), w2.data(), dpx_u.data(), dpy_v.data(),
                                  dpz_w.data(), nxu, nyu, nzu, nxv, nyv, nzv, nxw, nyw, nzw, nxc,
                                  nyc, nzc, ng, rho_c.data(), dt);

    for (std::size_t i = 0; i < Nu; ++i)
        REQUIRE(u1[i] == Approx(u2[i]).margin(1e-13));
    for (std::size_t i = 0; i < Nv; ++i)
        REQUIRE(v1[i] == Approx(v2[i]).margin(1e-13));
    for (std::size_t i = 0; i < Nw; ++i)
        REQUIRE(w1[i] == Approx(w2[i]).margin(1e-13));
}
