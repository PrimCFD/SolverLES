#include "MacOps.hpp"
#include <algorithm>
#include <catch2/catch_all.hpp>
#include <cmath>
#include <vector>
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

TEST_CASE("MAC gradp_faces + divergence_mac equals the discrete 7-point Laplacian on centers",
          "[mac][grad][div]")
{
    // Slightly larger grid + wider ghosts to keep boundary effects away
    const int nx = 10, ny = 9, nz = 8, ng = 2;
    const double dx = 1.0, dy = 1.0, dz = 1.0;

    // Totals
    const int nxc = center_tot(nx, ng), nyc = center_tot(ny, ng), nzc = center_tot(nz, ng);
    const int nxu = face_tot(nx, ng), nyu = center_tot(ny, ng), nzu = center_tot(nz, ng);
    const int nxv = center_tot(nx, ng), nyv = face_tot(ny, ng), nzv = center_tot(nz, ng);
    const int nxw = center_tot(nx, ng), nyw = center_tot(ny, ng), nzw = face_tot(nz, ng);

    const std::size_t Np = (std::size_t) nxc * nyc * nzc;
    const std::size_t Nu = (std::size_t) nxu * nyu * nzu;
    const std::size_t Nv = (std::size_t) nxv * nyv * nzv;
    const std::size_t Nw = (std::size_t) nxw * nyw * nzw;

    // Smooth polynomial field at centers (includes cross terms)
    std::vector<double> p(Np, 0.0);

    auto idxC = [&](int i, int j, int k)
    {
        return (std::size_t) i +
               (std::size_t) nxc * ((std::size_t) j + (std::size_t) nyc * (std::size_t) k);
    };

    for (int k = 0; k < nzc; ++k)
        for (int j = 0; j < nyc; ++j)
            for (int i = 0; i < nxc; ++i)
            {
                const double x = double(i - ng);
                const double y = double(j - ng);
                const double z = double(k - ng);
                p[idxC(i, j, k)] = 1.7 * x * x - 0.8 * y * y + 0.5 * z * z + 0.1 * x * y -
                                   0.2 * y * z + 0.05 * x * z;
            }

    // Face gradients via kernel
    std::vector<double> dpx_u(Nu, 0.0), dpy_v(Nv, 0.0), dpz_w(Nw, 0.0);
    grad_p_faces(p.data(), nxc, nyc, nzc, ng, dx, dy, dz, dpx_u.data(), nxu, nyu, nzu, dpy_v.data(),
                 nxv, nyv, nzv, dpz_w.data(), nxw, nyw, nzw);

    // div(grad p) via MAC
    std::vector<double> mac_lap(Np, 0.0);
    divergence(dpx_u.data(), dpy_v.data(), dpz_w.data(), nxu, nyu, nzu, nxv, nyv, nzv, nxw, nyw,
               nzw, nxc, nyc, nzc, ng, dx, dy, dz, mac_lap.data());

    // Discrete 7-point Laplacian at centers
    std::vector<double> disc_lap(Np, 0.0);
    const double ax = 1.0 / (dx * dx), ay = 1.0 / (dy * dy), az = 1.0 / (dz * dz);

    for (int k = ng; k < nzc - ng; ++k)
        for (int j = ng; j < nyc - ng; ++j)
            for (int i = ng; i < nxc - ng; ++i)
            {
                const auto c = idxC(i, j, k);
                const auto ip = idxC(i + 1, j, k), im = idxC(i - 1, j, k);
                const auto jp = idxC(i, j + 1, k), jm = idxC(i, j - 1, k);
                const auto kp = idxC(i, j, k + 1), km = idxC(i, j, k - 1);

                disc_lap[c] = ax * (p[ip] - 2.0 * p[c] + p[im]) +
                              ay * (p[jp] - 2.0 * p[c] + p[jm]) + az * (p[kp] - 2.0 * p[c] + p[km]);
            }

    // Compare on a deeper interior slab to fully decouple from BC handling
    const int pad = ng + 1; // skip an extra interior shell
    double linf = 0.0;

    for (int k = pad; k < nzc - pad; ++k)
        for (int j = pad; j < nyc - pad; ++j)
            for (int i = pad; i < nxc - pad; ++i)
            {
                const auto c = idxC(i, j, k);
                linf = std::max(linf, std::abs(mac_lap[c] - disc_lap[c]));
            }

    // Allow a tiny numerical slack for FP roundoff
    REQUIRE(linf == Approx(0.0).margin(1e-12));
}
