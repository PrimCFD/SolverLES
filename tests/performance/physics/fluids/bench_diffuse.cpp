#include "simple_bench.hpp"
#include <random>
#include <vector>
#include "kernels_fluids.h"

int main()
{
    const int nx = 64, ny = 64, nz = 64, ng = 2;
    const int nxc_tot = nx + 2 * ng, nyc_tot = ny + 2 * ng, nzc_tot = nz + 2 * ng;
    const int nxu_tot = (nx + 1) + 2 * ng, nyu_tot = nyc_tot, nzu_tot = nzc_tot;
    const int nxv_tot = nxc_tot, nyv_tot = (ny + 1) + 2 * ng, nzv_tot = nzc_tot;
    const int nxw_tot = nxc_tot, nyw_tot = nyc_tot, nzw_tot = (nz + 1) + 2 * ng;
    const double dx = 1.0, dy = 1.0, dz = 1.0, dt = 1e-3;

    std::vector<double> u((size_t) nxu_tot * nyu_tot * nzu_tot),
        v((size_t) nxv_tot * nyv_tot * nzv_tot), w((size_t) nxw_tot * nyw_tot * nzw_tot),
        nu_eff((size_t) nxc_tot * nyc_tot * nzc_tot, 1.0e-3), us(u.size()), vs(v.size()),
        ws(w.size());

    std::mt19937 rng(42);
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
            diffuse_velocity_fe_mac_c(u.data(), nxu_tot, nyu_tot, nzu_tot, v.data(), nxv_tot,
                                      nyv_tot, nzv_tot, w.data(), nxw_tot, nyw_tot, nzw_tot,
                                      nu_eff.data(), nxc_tot, nyc_tot, nzc_tot, ng, dx, dy, dz, dt,
                                      us.data(), vs.data(), ws.data());
        });

    const size_t Nfaces = (size_t) nxu_tot * nyu_tot * nzu_tot +
                          (size_t) nxv_tot * nyv_tot * nzv_tot +
                          (size_t) nxw_tot * nyw_tot * nzw_tot;
    bench::report("fluids_diffuse_mac_64^3", mean, stddev,
                  /*bytes (rough):*/ (Nfaces * 3 + (size_t) nxc_tot * nyc_tot * nzc_tot) *
                      sizeof(double));
    return 0;
}
