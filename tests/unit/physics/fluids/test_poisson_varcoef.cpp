// test_poisson_varcoef.cpp
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <vector>
#include "kernels_fluids.h"

using Catch::Approx;

static inline size_t idx(int I, int J, int K, int nx_tot, int ny_tot)
{
    return (size_t) I + (size_t) nx_tot * ((size_t) J + (size_t) ny_tot * (size_t) K);
}

TEST_CASE("Variable-coefficient Poisson reduces residual and matches CC when beta=1",
          "[fluids][poisson][varcoef]")
{
    const int nx = 12, ny = 9, nz = 7, ng = 1;
    const int nx_tot = nx + 2 * ng, ny_tot = ny + 2 * ng, nz_tot = nz + 2 * ng;
    const double dx = 1.0, dy = 1.0, dz = 1.0;
    const size_t N = (size_t) nx_tot * ny_tot * nz_tot;

    std::vector<double> beta(N), p(N), rhs(N), tmp(N);

    // beta(x,y,z)
    for (int K = 0; K < nz_tot; ++K)
        for (int J = 0; J < ny_tot; ++J)
            for (int I = 0; I < nx_tot; ++I)
            {
                const auto q = idx(I, J, K, nx_tot, ny_tot);
                beta[q] = 1.0 + 0.3 * std::sin(2 * M_PI * (I - 0.5) / nx) *
                                    std::cos(2 * M_PI * (J - 0.5) / ny) *
                                    std::cos(2 * M_PI * (K - 0.5) / nz);
            }

    auto apply_divbgrad = [&](const std::vector<double>& pin, std::vector<double>& out)
    {
        const double ax = 1.0 / (dx * dx), ay = 1.0 / (dy * dy), az = 1.0 / (dz * dz);
        for (int K = ng; K < nz + ng; ++K)
            for (int J = ng; J < ny + ng; ++J)
                for (int I = ng; I < nx + ng; ++I)
                {
                    const size_t c = idx(I, J, K, nx_tot, ny_tot);
                    const size_t ip = idx(I + 1, J, K, nx_tot, ny_tot),
                                 im = idx(I - 1, J, K, nx_tot, ny_tot);
                    const size_t jp = idx(I, J + 1, K, nx_tot, ny_tot),
                                 jm = idx(I, J - 1, K, nx_tot, ny_tot);
                    const size_t kp = idx(I, J, K + 1, nx_tot, ny_tot),
                                 km = idx(I, J, K - 1, nx_tot, ny_tot);
                    const double be = 0.5 * (beta[c] + beta[ip]), bw = 0.5 * (beta[c] + beta[im]);
                    const double bn = 0.5 * (beta[c] + beta[jp]), bs = 0.5 * (beta[c] + beta[jm]);
                    const double bt = 0.5 * (beta[c] + beta[kp]), bb = 0.5 * (beta[c] + beta[km]);
                    out[c] = ax * (be * pin[ip] - (be + bw) * pin[c] + bw * pin[im]) +
                             ay * (bn * pin[jp] - (bn + bs) * pin[c] + bs * pin[jm]) +
                             az * (bt * pin[kp] - (bt + bb) * pin[c] + bb * pin[km]);
                }
    };

    // Manufactured p*, then build rhs = div(beta grad p*)
    std::vector<double> pstar(N, 0.0), rhs_exact(N, 0.0);
    for (int K = 0; K < nz_tot; ++K)
        for (int J = 0; J < ny_tot; ++J)
            for (int I = 0; I < nx_tot; ++I)
            {
                const auto q = idx(I, J, K, nx_tot, ny_tot);
                pstar[q] = std::sin(2 * M_PI * (I - 0.5) / nx) +
                           0.5 * std::cos(2 * M_PI * (J - 0.5) / ny) +
                           0.25 * std::cos(2 * M_PI * (K - 0.5) / nz);
            }
    apply_divbgrad(pstar, rhs_exact);

    // Solve with var-coef Jacobi
    std::vector<double> p_var(N, 0.0);
    for (int k = 0; k < 150; ++k)
        poisson_jacobi_varcoef_c(rhs_exact.data(), beta.data(), nx_tot, ny_tot, nz_tot, ng, dx, dy,
                                 dz, 1, p_var.data());

    // Check relative residual ||A p - rhs|| / ||rhs||
    std::vector<double> Ap(N, 0.0);
    apply_divbgrad(p_var, Ap);
    long double res2 = 0.0, rhs2 = 0.0;
    for (int K = ng; K < nz + ng; ++K)
        for (int J = ng; J < ny + ng; ++J)
            for (int I = ng; I < nx + ng; ++I)
            {
                const auto c = idx(I, J, K, nx_tot, ny_tot);
                const long double r = (long double) Ap[c] - (long double) rhs_exact[c];
                res2 += r * r;
                rhs2 += (long double) rhs_exact[c] * (long double) rhs_exact[c] + 1e-30L;
            }
    REQUIRE(std::sqrt((double) (res2 / rhs2)) < 5e-2);

    // Equivalence when beta==1: varcoef == constant-coef operator
    std::vector<double> beta1(N, 1.0), rhs1(N, 0.0), p_cc(N, 0.0), p_var1(N, 0.0);
    // rhs1 = Δ p0 for some p0
    std::vector<double> p0(N, 0.0);
    for (int K = 0; K < nz_tot; ++K)
        for (int J = 0; J < ny_tot; ++J)
            for (int I = 0; I < nx_tot; ++I)
                p0[idx(I, J, K, nx_tot, ny_tot)] =
                    std::cos(2 * M_PI * (I - 0.5) / nx) + 0.3 * std::sin(2 * M_PI * (J - 0.5) / ny);

    auto lap = [&](const std::vector<double>& pin, std::vector<double>& out)
    {
        const double ax = 1.0 / (dx * dx), ay = 1.0 / (dy * dy), az = 1.0 / (dz * dz);
        for (int K = ng; K < nz + ng; ++K)
            for (int J = ng; J < ny + ng; ++J)
                for (int I = ng; I < nx + ng; ++I)
                {
                    const size_t c = idx(I, J, K, nx_tot, ny_tot);
                    const size_t ip = idx(I + 1, J, K, nx_tot, ny_tot),
                                 im = idx(I - 1, J, K, nx_tot, ny_tot);
                    const size_t jp = idx(I, J + 1, K, nx_tot, ny_tot),
                                 jm = idx(I, J - 1, K, nx_tot, ny_tot);
                    const size_t kp = idx(I, J, K + 1, nx_tot, ny_tot),
                                 km = idx(I, J, K - 1, nx_tot, ny_tot);
                    out[c] = ax * (pin[ip] - 2 * pin[c] + pin[im]) +
                             ay * (pin[jp] - 2 * pin[c] + pin[jm]) +
                             az * (pin[kp] - 2 * pin[c] + pin[km]);
                }
    };
    lap(p0, rhs1);

    for (int k = 0; k < 120; ++k)
        poisson_jacobi_varcoef_c(rhs1.data(), beta1.data(), nx_tot, ny_tot, nz_tot, ng, dx, dy, dz,
                                 1, p_var1.data());
    for (int k = 0; k < 120; ++k)
        poisson_jacobi_c(rhs1.data(), nx_tot, ny_tot, nz_tot, ng, dx, dy, dz, 1, p_cc.data());

    std::vector<double> Lvar(N, 0.0), Lcc(N, 0.0);
    lap(p_var1, Lvar);
    lap(p_cc, Lcc);
    long double diff2 = 0.0, base2 = 0.0;
    for (int K = ng; K < nz + ng; ++K)
        for (int J = ng; J < ny + ng; ++J)
            for (int I = ng; I < nx + ng; ++I)
            {
                const auto c = idx(I, J, K, nx_tot, ny_tot);
                const long double d = (long double) Lvar[c] - (long double) Lcc[c];
                diff2 += d * d;
                base2 += (long double) Lcc[c] * (long double) Lcc[c] + 1e-30L;
            }
    REQUIRE(std::sqrt((double) (diff2 / base2)) < 1e-10);
}
