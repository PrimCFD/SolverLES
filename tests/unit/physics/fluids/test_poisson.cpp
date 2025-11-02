#include <algorithm>
#include <array>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <string>
#include <vector>

#include "PressurePoisson.hpp" // fluids::make_poisson, IAction
#include "Program.hpp"         // KV, BcTable (parsed inside make_poisson)
#include "master/FieldCatalog.hpp"
#include "master/RunContext.hpp"
#include "master/Views.hpp" // MeshTileView / AnyFieldView box

#include "test_petsc_guard.hpp"
static PetscTestGuard _petsc_guard;

using Catch::Approx;
using core::master::FieldCatalog;
using core::master::MeshTileView;
using core::master::RunContext;
using core::master::plugin::KV;

static inline size_t cidx(int I, int J, int K, int nxc_tot, int nyc_tot)
{
    return (size_t) I + (size_t) nxc_tot * ((size_t) J + (size_t) nyc_tot * (size_t) K);
}
static inline size_t aidx(int I, int J, int K, int nx_tot, int ny_tot)
{
    return (size_t) I + (size_t) nx_tot * ((size_t) J + (size_t) ny_tot * (size_t) K);
}

static inline std::array<std::ptrdiff_t, 3> strides_bytes(int nx, int ny)
{
    const std::ptrdiff_t s0 = (std::ptrdiff_t) sizeof(double);
    const std::ptrdiff_t s1 = (std::ptrdiff_t) nx * s0;
    const std::ptrdiff_t s2 = (std::ptrdiff_t) nx * (std::ptrdiff_t) ny * s0;
    return {s0, s1, s2};
}

// Build MAC face velocities from cell-centered p* so that
// (1/dt) div(u*) = div(∇p*)  (β=1)
static void build_faces_from_pstar(const std::vector<double>& pstar, std::vector<double>& u,
                                   std::vector<double>& v, std::vector<double>& w, int nx, int ny,
                                   int nz, int ng, double dx, double dy, double dz, double dt)
{
    const int nxc_tot = nx + 2 * ng, nyc_tot = ny + 2 * ng, nzc_tot = nz + 2 * ng;
    const int nxu_tot = nxc_tot + 1, nyu_tot = nyc_tot, nzu_tot = nzc_tot;
    const int nxv_tot = nxc_tot, nyv_tot = nyc_tot + 1, nzv_tot = nzc_tot;
    const int nxw_tot = nxc_tot, nyw_tot = nyc_tot, nzw_tot = nzc_tot + 1;

    // u at (i+1/2,j,k)
    for (int K = 0; K < nzu_tot; ++K)
        for (int J = 0; J < nyu_tot; ++J)
            for (int I = 0; I < nxu_tot; ++I)
            {
                if (I == 0 || I == nxu_tot - 1)
                {
                    u[aidx(I, J, K, nxu_tot, nyu_tot)] = 0.0;
                    continue;
                }
                const int ic = I - 1, jc = J, kc = K;
                const double pc = pstar[cidx(ic, jc, kc, nxc_tot, nyc_tot)];
                const double pn = pstar[cidx(ic + 1, jc, kc, nxc_tot, nyc_tot)];
                u[aidx(I, J, K, nxu_tot, nyu_tot)] = dt * (pn - pc) / dx;
            }

    // v at (i,j+1/2,k)
    for (int K = 0; K < nzv_tot; ++K)
        for (int J = 0; J < nyv_tot; ++J)
            for (int I = 0; I < nxv_tot; ++I)
            {
                if (J == 0 || J == nyv_tot - 1)
                {
                    v[aidx(I, J, K, nxv_tot, nyv_tot)] = 0.0;
                    continue;
                }
                const int ic = I, jc = J - 1, kc = K;
                const double pc = pstar[cidx(ic, jc, kc, nxc_tot, nyc_tot)];
                const double pn = pstar[cidx(ic, jc + 1, kc, nxc_tot, nyc_tot)];
                v[aidx(I, J, K, nxv_tot, nyv_tot)] = dt * (pn - pc) / dy;
            }

    // w at (i,j,k+1/2)
    for (int K = 0; K < nzw_tot; ++K)
        for (int J = 0; J < nyw_tot; ++J)
            for (int I = 0; I < nxw_tot; ++I)
            {
                if (K == 0 || K == nzw_tot - 1)
                {
                    w[aidx(I, J, K, nxw_tot, nyw_tot)] = 0.0;
                    continue;
                }
                const int ic = I, jc = J, kc = K - 1;
                const double pc = pstar[cidx(ic, jc, kc, nxc_tot, nyc_tot)];
                const double pn = pstar[cidx(ic, jc, kc + 1, nxc_tot, nyc_tot)];
                w[aidx(I, J, K, nxw_tot, nyw_tot)] = dt * (pn - pc) / dz;
            }
}

static double interior_mean(const std::vector<double>& a, int nx, int ny, int nz, int ng,
                            int nxc_tot, int nyc_tot)
{
    long double s = 0.0L;
    size_t n = 0;
    for (int K = ng; K < nz + ng; ++K)
        for (int J = ng; J < ny + ng; ++J)
            for (int I = ng; I < nx + ng; ++I)
            {
                s += a[cidx(I, J, K, nxc_tot, nyc_tot)];
                ++n;
            }
    return (double) (s / (long double) n);
}

TEST_CASE("MG Poisson reproduces manufactured p* (β=1) using FieldCatalog::register_scalar",
          "[fluids][poisson][mg]")
{
    const int nx = 16, ny = 16, nz = 16, ng = 1;
    const int nxc_tot = nx + 2 * ng, nyc_tot = ny + 2 * ng, nzc_tot = nz + 2 * ng;
    const int nxu_tot = nxc_tot + 1, nyu_tot = nyc_tot, nzu_tot = nzc_tot;
    const int nxv_tot = nxc_tot, nyv_tot = nyc_tot + 1, nzv_tot = nzc_tot;
    const int nxw_tot = nxc_tot, nyw_tot = nyc_tot, nzw_tot = nzc_tot + 1;

    const double dx = 1.0, dy = 1.0, dz = 1.0, dt = 1.0;

    std::vector<double> p(nxc_tot * nyc_tot * nzc_tot, 0.0);
    std::vector<double> pstar(p.size(), 0.0);
    std::vector<double> u(nxu_tot * nyu_tot * nzu_tot, 0.0);
    std::vector<double> v(nxv_tot * nyv_tot * nzv_tot, 0.0);
    std::vector<double> w(nxw_tot * nyw_tot * nzw_tot, 0.0);
    std::vector<double> rho(nxc_tot * nyc_tot * nzc_tot, 1.0); // β=1 → ρ=1

    // Manufactured p*
    for (int K = 0; K < nzc_tot; ++K)
        for (int J = 0; J < nyc_tot; ++J)
            for (int I = 0; I < nxc_tot; ++I)
            {
                pstar[cidx(I, J, K, nxc_tot, nyc_tot)] = std::cos(2 * M_PI * (I - 0.5) / nx) +
                                                         0.3 * std::cos(2 * M_PI * (J - 0.5) / ny) +
                                                         0.2 * std::cos(2 * M_PI * (K - 0.5) / nz);
            }

    build_faces_from_pstar(pstar, u, v, w, nx, ny, nz, ng, dx, dy, dz, dt);

    // One-tile view for [0,nx)×[0,ny)×[0,nz)
    MeshTileView tile{};
    FieldCatalog fields;

    // Register
    fields.register_scalar("p", p.data(), sizeof(double), {nxc_tot, nyc_tot, nzc_tot},
                           strides_bytes(nxc_tot, nyc_tot), core::master::Stagger::Cell);
    fields.register_scalar("u", u.data(), sizeof(double), {nxu_tot, nyu_tot, nzu_tot},
                           strides_bytes(nxu_tot, nyu_tot), core::master::Stagger::IFace);
    fields.register_scalar("v", v.data(), sizeof(double), {nxv_tot, nyv_tot, nzv_tot},
                           strides_bytes(nxv_tot, nyv_tot), core::master::Stagger::JFace);
    fields.register_scalar("w", w.data(), sizeof(double), {nxw_tot, nyw_tot, nzw_tot},
                           strides_bytes(nxw_tot, nyw_tot), core::master::Stagger::KFace);
    fields.register_scalar("rho", rho.data(), sizeof(double), {nxc_tot, nyc_tot, nzc_tot},
                           strides_bytes(nxc_tot, nyc_tot), core::master::Stagger::Cell);

    // Build MG Poisson action
    KV kv{
        {"dx", std::to_string(dx)},
        {"dy", std::to_string(dy)},
        {"dz", std::to_string(dz)},
        {"rho", "1.0"},
        {"iters", "50"},
        {"div_tol", "1e-10"},
    };

    core::master::RunContext rc{};
    auto poisson = fluids::make_poisson(kv, rc);

    poisson->execute(tile, fields, dt);

    // Compare p vs p* up to a constant
    const double mp = interior_mean(p, nx, ny, nz, ng, nxc_tot, nyc_tot);
    const double mps = interior_mean(pstar, nx, ny, nz, ng, nxc_tot, nyc_tot);

    long double num2 = 0.0L, den2 = 0.0L;
    for (int K = ng; K < nz + ng; ++K)
        for (int J = ng; J < ny + ng; ++J)
            for (int I = ng; I < nx + ng; ++I)
            {
                const size_t c = cidx(I, J, K, nxc_tot, nyc_tot);
                const long double d = (p[c] - mp) - (pstar[c] - mps); // Matches sign convention
                num2 += d * d;
                const long double b = (pstar[c] - mps);
                den2 += b * b + 1e-30L;
            }
    const double rel = std::sqrt((double) (num2 / den2));
    REQUIRE(rel < 5e-6);
}
