#include <algorithm>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <vector>
#include "kernels_fluids.h"

using Catch::Approx;

static inline int center_tot(int n, int ng)
{
    return n + 2 * ng;
}
static inline int face_tot(int n, int ng)
{
    return n + 1 + 2 * ng;
}
static inline std::size_t idx(int I, int J, int K, int nx, int ny)
{
    return std::size_t(I) + std::size_t(nx) * (std::size_t(J) + std::size_t(ny) * std::size_t(K));
}
static inline bool is_int(int I, int J, int K, int nx, int ny, int nz, int ng)
{
    return (I >= ng && I < nx - ng && J >= ng && J < ny - ng && K >= ng && K < nz - ng);
}
static void copy_halos_face(std::vector<double>& dst, const std::vector<double>& src, int nx,
                            int ny, int nz, int ng)
{
    for (int K = 0; K < nz; ++K)
        for (int J = 0; J < ny; ++J)
            for (int I = 0; I < nx; ++I)
                if (!is_int(I, J, K, nx, ny, nz, ng))
                    dst[idx(I, J, K, nx, ny)] = src[idx(I, J, K, nx, ny)];
}

TEST_CASE("BE Jacobi sweep reduces diffusion residual", "[fluids][diffuse][be]")
{
    const int nx = 10, ny = 7, nz = 5, ng = 1;
    const int nxc_tot = center_tot(nx, ng), nyc_tot = center_tot(ny, ng),
              nzc_tot = center_tot(nz, ng);
    const int nxu_tot = face_tot(nx, ng), nyu_tot = nyc_tot, nzu_tot = nzc_tot;
    const int nxv_tot = nxc_tot, nyv_tot = face_tot(ny, ng), nzv_tot = nzc_tot;
    const int nxw_tot = nxc_tot, nyw_tot = nyc_tot, nzw_tot = face_tot(nz, ng);
    const double dx = 1.0, dy = 1.0, dz = 1.0, dt = 0.5;

    std::vector<double> u_rhs((size_t) nxu_tot * nyu_tot * nzu_tot),
        v_rhs((size_t) nxv_tot * nyv_tot * nzv_tot), w_rhs((size_t) nxw_tot * nyw_tot * nzw_tot);
    std::vector<double> u_iter = u_rhs, v_iter = v_rhs, w_iter = w_rhs;
    std::vector<double> u_next(u_rhs.size(), 0.0), v_next(v_rhs.size(), 0.0),
        w_next(w_rhs.size(), 0.0);
    std::vector<double> nu_eff((size_t) nxc_tot * nyc_tot * nzc_tot, 1.0e-2);

    auto init_face = [](std::vector<double>& A, int nx, int ny, int nz)
    {
        for (int K = 0; K < nz; ++K)
            for (int J = 0; J < ny; ++J)
                for (int I = 0; I < nx; ++I)
                {
                    const size_t q = idx(I, J, K, nx, ny);
                    A[q] = 0.1 * I + 0.05 * J - 0.02 * K;
                }
    };
    init_face(u_rhs, nxu_tot, nyu_tot, nzu_tot);
    init_face(v_rhs, nxv_tot, nyv_tot, nzv_tot);
    init_face(w_rhs, nxw_tot, nyw_tot, nzw_tot);
    std::fill(u_iter.begin(), u_iter.end(), 0.0);
    std::fill(v_iter.begin(), v_iter.end(), 0.0);
    std::fill(w_iter.begin(), w_iter.end(), 0.0);

    double res2_before = -1.0, rhs2_before = -1.0;
    diffuse_velocity_be_residual_mac_c(
        u_rhs.data(), v_rhs.data(), w_rhs.data(), u_iter.data(), v_iter.data(), w_iter.data(),
        nu_eff.data(), nxc_tot, nyc_tot, nzc_tot, nxu_tot, nyu_tot, nzu_tot, nxv_tot, nyv_tot,
        nzv_tot, nxw_tot, nyw_tot, nzw_tot, ng, dx, dy, dz, dt, &res2_before, &rhs2_before);

    diffuse_velocity_be_sweep_mac_c(u_rhs.data(), v_rhs.data(), w_rhs.data(), u_iter.data(),
                                    v_iter.data(), w_iter.data(), nu_eff.data(), nxc_tot, nyc_tot,
                                    nzc_tot, nxu_tot, nyu_tot, nzu_tot, nxv_tot, nyv_tot, nzv_tot,
                                    nxw_tot, nyw_tot, nzw_tot, ng, dx, dy, dz, dt, u_next.data(),
                                    v_next.data(), w_next.data());

    copy_halos_face(u_next, u_iter, nxu_tot, nyu_tot, nzu_tot, ng);
    copy_halos_face(v_next, v_iter, nxv_tot, nyv_tot, nzv_tot, ng);
    copy_halos_face(w_next, w_iter, nxw_tot, nyw_tot, nzw_tot, ng);

    double res2_after = -1.0, rhs2_after = -1.0;
    diffuse_velocity_be_residual_mac_c(
        u_rhs.data(), v_rhs.data(), w_rhs.data(), u_next.data(), v_next.data(), w_next.data(),
        nu_eff.data(), nxc_tot, nyc_tot, nzc_tot, nxu_tot, nyu_tot, nzu_tot, nxv_tot, nyv_tot,
        nzv_tot, nxw_tot, nyw_tot, nzw_tot, ng, dx, dy, dz, dt, &res2_after, &rhs2_after);

    REQUIRE(rhs2_after == Approx(rhs2_before));
    REQUIRE(res2_after < res2_before);
}

TEST_CASE("Red-black GS updates only the requested color", "[fluids][diffuse][gs]")
{
    const int nx = 9, ny = 7, nz = 5, ng = 1;
    const int nxc_tot = center_tot(nx, ng), nyc_tot = center_tot(ny, ng),
              nzc_tot = center_tot(nz, ng);
    const int nxu_tot = face_tot(nx, ng), nyu_tot = nyc_tot, nzu_tot = nzc_tot;
    const int nxv_tot = nxc_tot, nyv_tot = face_tot(ny, ng), nzv_tot = nzc_tot;
    const int nxw_tot = nxc_tot, nyw_tot = nyc_tot, nzw_tot = face_tot(nz, ng);
    const double dx = 1.0, dy = 1.0, dz = 1.0, dt = 0.1;

    std::vector<double> u((size_t) nxu_tot * nyu_tot * nzu_tot),
        v((size_t) nxv_tot * nyv_tot * nzv_tot), w((size_t) nxw_tot * nyw_tot * nzw_tot);
    auto u_rhs = u, v_rhs = v, w_rhs = w;
    std::vector<double> nu_eff((size_t) nxc_tot * nyc_tot * nzc_tot, 1.0e-3);

    // Initialize each face field on its own MAC lattice to avoid any out-of-bounds writes.
    // u (x-faces)
    for (int K = 0; K < nzu_tot; ++K)
        for (int J = 0; J < nyu_tot; ++J)
            for (int I = 0; I < nxu_tot; ++I)
            {
                const size_t q = idx(I, J, K, nxu_tot, nyu_tot);
                u[q] = 1.0 + I + 10 * J + 100 * K;
                u_rhs[q] = u[q];
            }
    // v (y-faces)
    for (int K = 0; K < nzv_tot; ++K)
        for (int J = 0; J < nyv_tot; ++J)
            for (int I = 0; I < nxv_tot; ++I)
            {
                const size_t q = idx(I, J, K, nxv_tot, nyv_tot);
                v[q] = -2.0 + 2 * I - 5 * J + 3 * K;
                v_rhs[q] = v[q];
            }
    // w (z-faces)
    for (int K = 0; K < nzw_tot; ++K)
        for (int J = 0; J < nyw_tot; ++J)
            for (int I = 0; I < nxw_tot; ++I)
            {
                const size_t q = idx(I, J, K, nxw_tot, nyw_tot);
                w[q] = 0.5 - I + 4 * J - 2 * K;
                w_rhs[q] = w[q];
            }

    // Keep a copy for comparison
    auto u_before = u, v_before = v, w_before = w;

    // Update color 0 (red) on MAC grid
    diffuse_velocity_be_gs_color_mac_c(u.data(), v.data(), w.data(), u_rhs.data(), v_rhs.data(),
                                       w_rhs.data(), nu_eff.data(), nxc_tot, nyc_tot, nzc_tot,
                                       nxu_tot, nyu_tot, nzu_tot, nxv_tot, nyv_tot, nzv_tot,
                                       nxw_tot, nyw_tot, nzw_tot, ng, dx, dy, dz, dt, /*color=*/0);

    // Check parity behavior (per-component to respect each face-grid extent)
    bool any_updated = false;

    // u (x-faces): interior
    for (int K = 1; K < nzu_tot - 1; ++K)
        for (int J = 1; J < nyu_tot - 1; ++J)
            for (int I = 1; I < nxu_tot - 1; ++I)
            {
                const int parity = (I - ng) + (J - ng) + (K - ng);
                const size_t q = idx(I, J, K, nxu_tot, nyu_tot);
                if ((parity & 1) == 0)
                { // red
                    if (std::abs(u[q] - u_before[q]) > 1e-15)
                        any_updated = true;
                }
                else
                {
                    REQUIRE(u[q] == Approx(u_before[q]).margin(0));
                }
            }

    // v (y-faces): interior
    for (int K = 1; K < nzv_tot - 1; ++K)
        for (int J = 1; J < nyv_tot - 1; ++J)
            for (int I = 1; I < nxv_tot - 1; ++I)
            {
                const int parity = (I - ng) + (J - ng) + (K - ng);
                const size_t q = idx(I, J, K, nxv_tot, nyv_tot);
                if ((parity & 1) == 0)
                { // red
                    if (std::abs(v[q] - v_before[q]) > 1e-15)
                        any_updated = true;
                }
                else
                {
                    REQUIRE(v[q] == Approx(v_before[q]).margin(0));
                }
            }

    // w (z-faces): interior
    for (int K = 1; K < nzw_tot - 1; ++K)
        for (int J = 1; J < nyw_tot - 1; ++J)
            for (int I = 1; I < nxw_tot - 1; ++I)
            {
                const int parity = (I - ng) + (J - ng) + (K - ng);
                const size_t q = idx(I, J, K, nxw_tot, nyw_tot);
                if ((parity & 1) == 0)
                { // red
                    if (std::abs(w[q] - w_before[q]) > 1e-15)
                        any_updated = true;
                }
                else
                {
                    REQUIRE(w[q] == Approx(w_before[q]).margin(0));
                }
            }

    REQUIRE(any_updated); // ensure the test actually exercised updates
}
