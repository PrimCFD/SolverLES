#include "MacOps.hpp"
#include "simple_bench.hpp"
#include <random>
#include <vector>
using namespace numerics::kernels;

static inline size_t idx(int I, int J, int K, int nx_tot, int ny_tot)
{
    return static_cast<size_t>(I + nx_tot * (J + ny_tot * K));
}

int main()
{
    // Interior sizes and ghosts
    const int nx = 128, ny = 64, nz = 64, ng = 1;

    // Cell-centered totals (centers)
    const int nxc_tot = nx + 2 * ng;
    const int nyc_tot = ny + 2 * ng;
    const int nzc_tot = nz + 2 * ng;

    // Face-centered totals (MAC)
    const int nxu_tot = nx + 1 + 2 * ng; // u has +1 along x
    const int nyu_tot = ny + 0 + 2 * ng;
    const int nzu_tot = nz + 0 + 2 * ng;

    const int nxv_tot = nx + 0 + 2 * ng;
    const int nyv_tot = ny + 1 + 2 * ng; // v has +1 along y
    const int nzv_tot = nz + 0 + 2 * ng;

    const int nxw_tot = nx + 0 + 2 * ng;
    const int nyw_tot = ny + 0 + 2 * ng;
    const int nzw_tot = nz + 1 + 2 * ng; // w has +1 along z

    const double dx = 1.0, dy = 1.0, dz = 1.0;
    const double rho = 1.0, dt = 1e-3;

    // Element counts
    const size_t Nc = static_cast<size_t>(nxc_tot) * nyc_tot * nzc_tot;
    const size_t Nu = static_cast<size_t>(nxu_tot) * nyu_tot * nzu_tot;
    const size_t Nv = static_cast<size_t>(nxv_tot) * nyv_tot * nzv_tot;
    const size_t Nw = static_cast<size_t>(nxw_tot) * nyw_tot * nzw_tot;

    // Arrays: face velocities, cell pressure, face-gradients, cell divergence
    std::vector<double> u(Nu), v(Nv), w(Nw);
    std::vector<double> p(Nc, 0.0), dpx_u(Nu), dpy_v(Nv), dpz_w(Nw);
    std::vector<double> div(Nc);

    // Random init
    std::mt19937 rng(1312);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    for (auto& x : u)
        x = dist(rng);
    for (auto& x : v)
        x = dist(rng);
    for (auto& x : w)
        x = dist(rng);
    for (auto& x : p)
        x = dist(rng);

    // grad(p) -> face gradients
    auto [mean_grad, std_grad] = bench::run(
        [&]
        {
            grad_p_faces(p.data(), nxc_tot, nyc_tot, nzc_tot, ng, dx, dy, dz, dpx_u.data(), nxu_tot,
                         nyu_tot, nzu_tot, dpy_v.data(), nxv_tot, nyv_tot, nzv_tot, dpz_w.data(),
                         nxw_tot, nyw_tot, nzw_tot);
        });
    // bytes ~ read p (Nc) + write three face grads (Nu+Nv+Nw)
    double bytes_grad = (static_cast<double>(Nc + Nu + Nv + Nw)) * sizeof(double);
    bench::report("fluids_gradp_faces_128x64x64", mean_grad, std_grad, bytes_grad);

    // div(u) at centers from face velocities
    auto [mean_div, std_div] = bench::run(
        [&]
        {
            divergence(u.data(), v.data(), w.data(), nxu_tot, nyu_tot, nzu_tot, nxv_tot, nyv_tot,
                       nzv_tot, nxw_tot, nyw_tot, nzw_tot, nxc_tot, nyc_tot, nzc_tot, ng, dx, dy,
                       dz, div.data());
        });
    // bytes ~ read u/v/w faces (Nu+Nv+Nw) + write div centers (Nc)
    double bytes_div = (static_cast<double>(Nu + Nv + Nw + Nc)) * sizeof(double);
    bench::report("fluids_divergence_mac_128x64x64", mean_div, std_div, bytes_div);

    // corrector: u -= (dt/rho) * gradp on faces (constant rho path)
    auto [mean_corr, std_corr] = bench::run(
        [&]
        {
            correct_velocity_const_rho(u.data(), v.data(), w.data(), dpx_u.data(), dpy_v.data(),
                                       dpz_w.data(), nxu_tot, nyu_tot, nzu_tot, nxv_tot, nyv_tot,
                                       nzv_tot, nxw_tot, nyw_tot, nzw_tot, nxc_tot, nyc_tot,
                                       nzc_tot, ng, rho, dt);
        });
    // bytes ~ inout u/v/w faces (Nu+Nv+Nw) + read grads (Nu+Nv+Nw)
    double bytes_corr = (static_cast<double>(2 * (Nu + Nv + Nw))) * sizeof(double);
    bench::report("fluids_correct_velocity_mac_128x64x64", mean_corr, std_corr, bytes_corr);

    return 0;
}
