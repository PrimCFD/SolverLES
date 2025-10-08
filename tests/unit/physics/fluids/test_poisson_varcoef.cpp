#include <algorithm>
#include <array>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <string>
#include <vector>

#include "PressurePoisson.hpp"
#include "Program.hpp"
#include "master/FieldCatalog.hpp"
#include "master/RunContext.hpp"
#include "master/Views.hpp"

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
static inline double harm(double a, double b)
{
    const double aa = std::max(a, 1e-300), bb = std::max(b, 1e-300);
    return 2.0 / (1.0 / aa + 1.0 / bb);
}

static void build_faces_from_pstar_varbeta(const std::vector<double>& pstar,
                                           const std::vector<double>& beta_c,
                                           std::vector<double>& u, std::vector<double>& v,
                                           std::vector<double>& w, int nx, int ny, int nz, int ng,
                                           double dx, double dy, double dz, double dt)
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
                const double betaf = harm(beta_c[cidx(ic, jc, kc, nxc_tot, nyc_tot)],
                                          beta_c[cidx(ic + 1, jc, kc, nxc_tot, nyc_tot)]);
                u[aidx(I, J, K, nxu_tot, nyu_tot)] = dt * betaf * (pn - pc) / dx;
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
                const double betaf = harm(beta_c[cidx(ic, jc, kc, nxc_tot, nyc_tot)],
                                          beta_c[cidx(ic, jc + 1, kc, nxc_tot, nyc_tot)]);
                v[aidx(I, J, K, nxv_tot, nyv_tot)] = dt * betaf * (pn - pc) / dy;
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
                const double betaf = harm(beta_c[cidx(ic, jc, kc, nxc_tot, nyc_tot)],
                                          beta_c[cidx(ic, jc, kc + 1, nxc_tot, nyc_tot)]);
                w[aidx(I, J, K, nxw_tot, nyw_tot)] = dt * betaf * (pn - pc) / dz;
            }
}

// Discrete div(β∇p) at centers (for residual checks)
static void apply_divbgrad(const std::vector<double>& beta_c, const std::vector<double>& p,
                           std::vector<double>& out, int nx, int ny, int nz, int ng, double dx,
                           double dy, double dz)
{
    const int nxc_tot = nx + 2 * ng, nyc_tot = ny + 2 * ng, nzc_tot = nz + 2 * ng;
    // zero all (safe for ghosts too)
    std::fill(out.begin(), out.end(), 0.0);

    for (int K = ng; K < nz + ng; ++K)
        for (int J = ng; J < ny + ng; ++J)
            for (int I = ng; I < nx + ng; ++I)
            {
                const size_t c = cidx(I, J, K, nxc_tot, nyc_tot);
                const size_t W = cidx(I - 1, J, K, nxc_tot, nyc_tot);
                const size_t E = cidx(I + 1, J, K, nxc_tot, nyc_tot);
                const size_t S = cidx(I, J - 1, K, nxc_tot, nyc_tot);
                const size_t N = cidx(I, J + 1, K, nxc_tot, nyc_tot);
                const size_t B = cidx(I, J, K - 1, nxc_tot, nyc_tot);
                const size_t T = cidx(I, J, K + 1, nxc_tot, nyc_tot);

                const double txW = harm(beta_c[W], beta_c[c]) / (dx * dx);
                const double txE = harm(beta_c[E], beta_c[c]) / (dx * dx);
                const double tyS = harm(beta_c[S], beta_c[c]) / (dy * dy);
                const double tyN = harm(beta_c[N], beta_c[c]) / (dy * dy);
                const double tzB = harm(beta_c[B], beta_c[c]) / (dz * dz);
                const double tzT = harm(beta_c[T], beta_c[c]) / (dz * dz);

                out[c] = txE * (p[E] - p[c]) + txW * (p[W] - p[c]) + tyN * (p[N] - p[c]) +
                         tyS * (p[S] - p[c]) + tzT * (p[T] - p[c]) + tzB * (p[B] - p[c]);
            }
}

TEST_CASE("MG variable-β Poisson: A p ≈ rhs and β→1 consistency", "[fluids][poisson][mg][varcoef]")
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
    std::vector<double> beta(nxc_tot * nyc_tot * nzc_tot, 1.0);
    std::vector<double> rho(beta.size(), 1.0); // rho = 1/beta

    // Variable β at centers (smooth)
    for (int K = 0; K < nzc_tot; ++K)
        for (int J = 0; J < nyc_tot; ++J)
            for (int I = 0; I < nxc_tot; ++I)
            {
                const double x = (I - 0.5 - ng) / (double) nx;
                const double y = (J - 0.5 - ng) / (double) ny;
                const double z = (K - 0.5 - ng) / (double) nz;
                const double b = 1.0 + 0.4 * std::sin(2 * M_PI * x) * std::cos(2 * M_PI * y) +
                                 0.2 * std::cos(2 * M_PI * z);
                beta[cidx(I, J, K, nxc_tot, nyc_tot)] = b;
                rho[cidx(I, J, K, nxc_tot, nyc_tot)] = 1.0 / std::max(b, 1e-12);
            }

    // Manufactured p*
    for (int K = 0; K < nzc_tot; ++K)
        for (int J = 0; J < nyc_tot; ++J)
            for (int I = 0; I < nxc_tot; ++I)
            {
                pstar[cidx(I, J, K, nxc_tot, nyc_tot)] = std::sin(2 * M_PI * (I - 0.5) / nx) +
                                                         0.3 * std::cos(2 * M_PI * (J - 0.5) / ny) +
                                                         0.2 * std::cos(2 * M_PI * (K - 0.5) / nz);
            }

    // Build MAC faces so (1/dt)div(u*) = div(β∇p*)
    build_faces_from_pstar_varbeta(pstar, beta, u, v, w, nx, ny, nz, ng, dx, dy, dz, dt);

    MeshTileView tile{};
    tile.box.lo[0] = 0;
    tile.box.lo[1] = 0;
    tile.box.lo[2] = 0;
    tile.box.hi[0] = nx;
    tile.box.hi[1] = ny;
    tile.box.hi[2] = nz;

    FieldCatalog fields;

    // Register fields via FieldCatalog::register_scalar (byte strides)
    fields.register_scalar("p", p.data(), sizeof(double), {nxc_tot, nyc_tot, nzc_tot},
                           strides_bytes(nxc_tot, nyc_tot));
    fields.register_scalar("u", u.data(), sizeof(double), {nxu_tot, nyu_tot, nzu_tot},
                           strides_bytes(nxu_tot, nyu_tot));
    fields.register_scalar("v", v.data(), sizeof(double), {nxv_tot, nyv_tot, nzv_tot},
                           strides_bytes(nxv_tot, nyv_tot));
    fields.register_scalar("w", w.data(), sizeof(double), {nxw_tot, nyw_tot, nzw_tot},
                           strides_bytes(nxw_tot, nyw_tot));
    fields.register_scalar("rho", rho.data(), sizeof(double), {nxc_tot, nyc_tot, nzc_tot},
                           strides_bytes(nxc_tot, nyc_tot));

    // Build MG Poisson (variable density)
    KV kv{
        {"dx", std::to_string(dx)},
        {"dy", std::to_string(dy)},
        {"dz", std::to_string(dz)},
        {"iters", "60"},
        {"div_tol", "5e-11"}
        // leave pressure BCs unset → natural Neumann; solver should handle gauge (zero-mean)
    };
    core::master::RunContext rc{};
    if (PetscTestGuard::petsc_uses_mpi())
    {
        rc.mpi_comm = const_cast<void*>(_petsc_guard.mpi_comm_ptr()); // else leave nullptr
    }
    auto poisson = fluids::make_poisson(kv, rc);

    poisson->execute(tile, fields, dt);

    // Residual check: A p ≈ rhs_exact := div(β∇p*)
    std::vector<double> rhs_exact(p.size(), 0.0), Ap(p.size(), 0.0);
    apply_divbgrad(beta, pstar, rhs_exact, nx, ny, nz, ng, dx, dy, dz);
    apply_divbgrad(beta, p, Ap, nx, ny, nz, ng, dx, dy, dz);

    long double num2 = 0.0L, den2 = 0.0L;
    for (int K = ng; K < nz + ng; ++K)
        for (int J = ng; J < ny + ng; ++J)
            for (int I = ng; I < nx + ng; ++I)
            {
                const size_t c = cidx(I, J, K, nxc_tot, nyc_tot);
                const long double d = (long double) Ap[c] - (long double) rhs_exact[c];
                num2 += d * d;
                den2 += (long double) rhs_exact[c] * (long double) rhs_exact[c] + 1e-30L;
            }
    REQUIRE(std::sqrt((double) (num2 / den2)) < 8e-6);

    // Consistency probe: β→1 reproduces constant-coefficient operator
    std::fill(beta.begin(), beta.end(), 1.0);
    apply_divbgrad(beta, p, Ap, nx, ny, nz, ng, dx, dy, dz);
    apply_divbgrad(beta, pstar, rhs_exact, nx, ny, nz, ng, dx, dy, dz);

    long double diff2 = 0.0L, base2 = 0.0L;
    for (int K = ng; K < nz + ng; ++K)
        for (int J = ng; J < ny + ng; ++J)
            for (int I = ng; I < nx + ng; ++I)
            {
                const size_t c = cidx(I, J, K, nxc_tot, nyc_tot);
                const long double d = (long double) Ap[c] - (long double) rhs_exact[c];
                diff2 += d * d;
                base2 += (long double) rhs_exact[c] * (long double) rhs_exact[c] + 1e-30L;
            }
    REQUIRE(std::sqrt((double) (diff2 / base2)) < 1e-10);
}
