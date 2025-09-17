#include "simple_bench.hpp"
#include <random>
#include <vector>
#include "kernels_fluids.h"

int main()
{
    const int nx = 64, ny = 64, nz = 64, ng = 2;
    const int nx_tot = nx + 2 * ng, ny_tot = ny + 2 * ng, nz_tot = nz + 2 * ng;
    const double dx = 1.0, dy = 1.0, dz = 1.0, Cs = 0.16;

    std::vector<double> u(nx_tot * ny_tot * nz_tot), v(nx_tot * ny_tot * nz_tot),
        w(nx_tot * ny_tot * nz_tot), nu_t(nx_tot * ny_tot * nz_tot, 0.0);

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
            sgs_smagorinsky_c(u.data(), v.data(), w.data(), nx_tot, ny_tot, nz_tot, ng, dx, dy, dz,
                              Cs, nu_t.data());
        });

    bench::report("fluids_sgs_64^3", mean, stddev,
                  /*bytes=*/(nx_tot * ny_tot * nz_tot * 4 * sizeof(double)));
    return 0;
}
