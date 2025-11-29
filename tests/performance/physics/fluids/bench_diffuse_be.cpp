
#include "simple_bench.hpp"
#include <random>
#include <vector>
#include "MacOps.hpp"
using namespace numerics::kernels;

int main()
{
    const int nx = 64, ny = 64, nz = 64, ng = 2;
    const int nxc_tot = nx + 2 * ng, nyc_tot = ny + 2 * ng, nzc_tot = nz + 2 * ng;
    const int nxu_tot = (nx + 1) + 2 * ng, nyu_tot = nyc_tot, nzu_tot = nzc_tot;
    const int nxv_tot = nxc_tot, nyv_tot = (ny + 1) + 2 * ng, nzv_tot = nzc_tot;
    const int nxw_tot = nxc_tot, nyw_tot = nyc_tot, nzw_tot = (nz + 1) + 2 * ng;
    const double dx = 1.0, dy = 1.0, dz = 1.0, dt = 1e-3;

    std::vector<double> u_rhs((size_t) nxu_tot * nyu_tot * nzu_tot),
        v_rhs((size_t) nxv_tot * nyv_tot * nzv_tot), w_rhs((size_t) nxw_tot * nyw_tot * nzw_tot);
    std::vector<double> u0 = u_rhs, v0 = v_rhs, w0 = w_rhs;
    std::vector<double> u1 = u_rhs, v1 = v_rhs, w1 = w_rhs;
    std::vector<double> nu_eff((size_t) nxc_tot * nyc_tot * nzc_tot, 1.0e-3);

    std::mt19937 rng(2025);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    for (auto& x : u_rhs)
        x = dist(rng);
    u0 = u_rhs;
    for (auto& x : v_rhs)
        x = dist(rng);
    v0 = v_rhs;
    for (auto& x : w_rhs)
        x = dist(rng);
    w0 = w_rhs;

    // --- Jacobi sweep ping-pong ---
    auto [mean_jac, std_jac] = bench::run(
        [&]
        {
            diffuse_be_jacobi_sweep(u_rhs.data(), v_rhs.data(), w_rhs.data(), u0.data(),
                                            v0.data(), w0.data(), nu_eff.data(), nxc_tot, nyc_tot,
                                            nzc_tot, nxu_tot, nyu_tot, nzu_tot, nxv_tot, nyv_tot,
                                            nzv_tot, nxw_tot, nyw_tot, nzw_tot, ng, dx, dy, dz, dt,
                                            u1.data(), v1.data(), w1.data());
            // swap iterates
            std::swap(u0, u1);
            std::swap(v0, v1);
            std::swap(w0, w1);
        });

    // Rough streaming estimate per sweep
    const double Nfaces = double(nxu_tot) * nyu_tot * nzu_tot +
                          double(nxv_tot) * nyv_tot * nzv_tot + double(nxw_tot) * nyw_tot * nzw_tot;
    double bytes_jac = (7.0 * Nfaces / 3.0 + 3.0 * Nfaces / 3.0) * sizeof(double);

    bench::report("fluids_be_jacobi_sweep_mac_64^3", mean_jac, std_jac, bytes_jac);

    // --- Red/black GS (one "iteration" = red + black) ---
    // Start from u0/v0/w0 field
    auto [mean_gs, std_gs] = bench::run(
        [&]
        {
            diffuse_be_rbgs_color(
                u0.data(), v0.data(), w0.data(), u_rhs.data(), v_rhs.data(), w_rhs.data(),
                nu_eff.data(), nxc_tot, nyc_tot, nzc_tot, nxu_tot, nyu_tot, nzu_tot, nxv_tot,
                nyv_tot, nzv_tot, nxw_tot, nyw_tot, nzw_tot, ng, dx, dy, dz, dt, /*color=*/0);
            diffuse_be_rbgs_color(
                u0.data(), v0.data(), w0.data(), u_rhs.data(), v_rhs.data(), w_rhs.data(),
                nu_eff.data(), nxc_tot, nyc_tot, nzc_tot, nxu_tot, nyu_tot, nzu_tot, nxv_tot,
                nyv_tot, nzv_tot, nxw_tot, nyw_tot, nzw_tot, ng, dx, dy, dz, dt, /*color=*/1);
        });

    // Each color updates ~half interior faces; rough estimate per color ~(7N/2 reads + 3N/2 writes)
    double bytes_gs = (7.0 + 3.0) * 0.5 * Nfaces * sizeof(double) * 2.0; // red+black

    bench::report("fluids_be_gs_rb_mac_64^3", mean_gs, std_gs, bytes_gs);

    return 0;
}
