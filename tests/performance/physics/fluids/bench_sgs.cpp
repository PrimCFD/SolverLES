#include "simple_bench.hpp"
#include <random>
#include <vector>
#include "MacOps.hpp"
using namespace numerics::kernels;

int main()
{
    // Interior size and ghosts (need ng>=2 for wide stencils in practice)
    const int nx = 64, ny = 64, nz = 64, ng = 2;
    const double dx = 1.0, dy = 1.0, dz = 1.0, Cs = 0.16;

    // Center totals (cell-centered)
    const int nxc_tot = nx + 2 * ng;
    const int nyc_tot = ny + 2 * ng;
    const int nzc_tot = nz + 2 * ng;

    // MAC face totals
    const int nxu_tot = nx + 1 + 2 * ng, nyu_tot = ny + 0 + 2 * ng, nzu_tot = nz + 0 + 2 * ng;
    const int nxv_tot = nx + 0 + 2 * ng, nyv_tot = ny + 1 + 2 * ng, nzv_tot = nz + 0 + 2 * ng;
    const int nxw_tot = nx + 0 + 2 * ng, nyw_tot = ny + 0 + 2 * ng, nzw_tot = nz + 1 + 2 * ng;

    // Allocate MAC velocities and center nu_t
    std::vector<double> u((size_t) nxu_tot * nyu_tot * nzu_tot);
    std::vector<double> v((size_t) nxv_tot * nyv_tot * nzv_tot);
    std::vector<double> w((size_t) nxw_tot * nyw_tot * nzw_tot);
    std::vector<double> nu_t((size_t) nxc_tot * nyc_tot * nzc_tot, 0.0);

    std::mt19937 rng(12345);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    for (auto& x : u)
        x = dist(rng);
    for (auto& x : v)
        x = dist(rng);
    for (auto& x : w)
        x = dist(rng);

    auto [mean, stddev] = bench::run(
        [&]
        {
            // MAC-face kernel signature:
            // sgs_smagorinsky(u, v, w,
            //                       nxu, nyu, nzu, nxv, nyv, nzv, nxw, nyw, nzw,
            //                       nxc, nyc, nzc, ng, dx, dy, dz, Cs, nu_t_centers)
            sgs_smagorinsky(u.data(), v.data(), w.data(), nxu_tot, nyu_tot, nzu_tot, nxv_tot,
                                  nyv_tot, nzv_tot, nxw_tot, nyw_tot, nzw_tot, nxc_tot, nyc_tot,
                                  nzc_tot, ng, dx, dy, dz, Cs, nu_t.data());
        });

    // Rough I/O traffic: read u/v/w faces + write center nu_t
    const double bytes =
        ((double) nxu_tot * nyu_tot * nzu_tot + (double) nxv_tot * nyv_tot * nzv_tot +
         (double) nxw_tot * nyw_tot * nzw_tot + (double) nxc_tot * nyc_tot * nzc_tot) *
        sizeof(double);
    bench::report("fluids_sgs_mac_64^3", mean, stddev, bytes);
    return 0;
}