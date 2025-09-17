
#include "simple_bench.hpp"
#include <random>
#include <vector>
#include "kernels_fluids.h"

int main()
{
    const int nx = 64, ny = 64, nz = 64, ng = 1;
    const int nx_tot = nx + 2 * ng, ny_tot = ny + 2 * ng, nz_tot = nz + 2 * ng;
    const double dx = 1.0, dy = 1.0, dz = 1.0;
    const int iters = 10;

    const size_t N = static_cast<size_t>(nx_tot) * ny_tot * nz_tot;
    std::vector<double> rhs(N, 0.0), p(N, 0.0);

    // random sparse RHS
    std::mt19937 rng(777);
    for (int k = ng; k < nz_tot - ng; ++k)
        for (int j = ng; j < ny_tot - ng; ++j)
            for (int i = ng; i < nx_tot - ng; ++i)
            {
                size_t c = static_cast<size_t>(i + nx_tot * (j + ny_tot * k));
                if (((i + j + k) % 23) == 0)
                    rhs[c] = 1.0;
            }

    auto [mean, stddev] = bench::run(
        [&]
        { poisson_jacobi_c(rhs.data(), nx_tot, ny_tot, nz_tot, ng, dx, dy, dz, iters, p.data()); });

    // Estimate: each iter reads 6 neighbors + rhs and writes p; rough streaming bytes per iter: ~8N
    double bytes = 8.0 * N * sizeof(double) * iters;
    bench::report("fluids_poisson_jacobi_64^3x10", mean, stddev, bytes);
    return 0;
}
