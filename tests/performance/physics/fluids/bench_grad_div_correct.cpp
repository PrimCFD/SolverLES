
#include "simple_bench.hpp"
#include <random>
#include <vector>
#include "kernels_fluids.h"

static inline size_t idx(int I, int J, int K, int nx_tot, int ny_tot)
{
    return static_cast<size_t>(I + nx_tot * (J + ny_tot * K));
}

int main()
{
    const int nx = 128, ny = 64, nz = 64, ng = 1;
    const int nx_tot = nx + 2 * ng, ny_tot = ny + 2 * ng, nz_tot = nz + 2 * ng;
    const double dx = 1.0, dy = 1.0, dz = 1.0;
    const double rho = 1.0, dt = 1e-3;
    const size_t N = static_cast<size_t>(nx_tot) * ny_tot * nz_tot;

    std::vector<double> u(N), v(N), w(N);
    std::vector<double> p(N, 0.0), dpx(N), dpy(N), dpz(N);
    std::vector<double> div(N);

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

    // grad(p)
    auto [mean_grad, std_grad] = bench::run(
        [&] {
            gradp_c(p.data(), nx_tot, ny_tot, nz_tot, ng, dx, dy, dz, dpx.data(), dpy.data(),
                    dpz.data());
        });
    double bytes_grad = (1.0 + 3.0) * N * sizeof(double);
    bench::report("fluids_gradp_128x64x64", mean_grad, std_grad, bytes_grad);

    // div(u)
    auto [mean_div, std_div] = bench::run(
        [&]
        {
            divergence_c(u.data(), v.data(), w.data(), nx_tot, ny_tot, nz_tot, ng, dx, dy, dz,
                         div.data());
        });
    double bytes_div = (3.0 + 1.0) * N * sizeof(double);
    bench::report("fluids_divergence_128x64x64", mean_div, std_div, bytes_div);

    // corrector
    auto [mean_corr, std_corr] = bench::run(
        [&]
        {
            correct_velocity_c(u.data(), v.data(), w.data(), dpx.data(), dpy.data(), dpz.data(),
                               nx_tot, ny_tot, nz_tot, ng, rho, dt);
        });
    double bytes_corr = (3.0 + 3.0) * N * sizeof(double);
    bench::report("fluids_correct_velocity_128x64x64", mean_corr, std_corr, bytes_corr);

    return 0;
}
