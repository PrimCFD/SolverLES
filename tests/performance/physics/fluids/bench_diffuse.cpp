#include "simple_bench.hpp"
#include <random>
#include <vector>
#include "kernels_fluids.h"

int main()
{
    const int nx = 64, ny = 64, nz = 64, ng = 2;
    const int nx_tot = nx + 2 * ng, ny_tot = ny + 2 * ng, nz_tot = nz + 2 * ng;
    const double dx = 1.0, dy = 1.0, dz = 1.0, dt = 1e-3;

    std::vector<double> u(nx_tot * ny_tot * nz_tot), v(nx_tot * ny_tot * nz_tot),
        w(nx_tot * ny_tot * nz_tot), nu_eff(nx_tot * ny_tot * nz_tot, 1.0e-3), us(u.size()),
        vs(v.size()), ws(w.size());

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
            diffuse_velocity_fe_c(u.data(), v.data(), w.data(), nu_eff.data(), nx_tot, ny_tot,
                                  nz_tot, ng, dx, dy, dz, dt, us.data(), vs.data(), ws.data());
        });

    bench::report("fluids_diffuse_64^3", mean, stddev,
                  /*bytes=*/(nx_tot * ny_tot * nz_tot * 6 * sizeof(double)));
    return 0;
}
