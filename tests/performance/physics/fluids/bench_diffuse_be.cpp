
#include "simple_bench.hpp"
#include <random>
#include <vector>
#include "kernels_fluids.h"

int main()
{
    const int nx = 64, ny = 64, nz = 64, ng = 2;
    const int nx_tot = nx + 2 * ng, ny_tot = ny + 2 * ng, nz_tot = nz + 2 * ng;
    const double dx = 1.0, dy = 1.0, dz = 1.0, dt = 1e-3;

    const size_t N = static_cast<size_t>(nx_tot) * ny_tot * nz_tot;
    std::vector<double> u_rhs(N), v_rhs(N), w_rhs(N);
    std::vector<double> u0(N), v0(N), w0(N);
    std::vector<double> u1(N), v1(N), w1(N);
    std::vector<double> nu_eff(N, 1.0e-3);

    std::mt19937 rng(2025);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    for (size_t i = 0; i < N; ++i)
    {
        u_rhs[i] = u0[i] = dist(rng);
        v_rhs[i] = v0[i] = dist(rng);
        w_rhs[i] = w0[i] = dist(rng);
    }

    // --- Jacobi sweep ping-pong ---
    auto [mean_jac, std_jac] = bench::run(
        [&]
        {
            diffuse_velocity_be_sweep_c(u_rhs.data(), v_rhs.data(), w_rhs.data(), u0.data(),
                                        v0.data(), w0.data(), nu_eff.data(), nx_tot, ny_tot, nz_tot,
                                        ng, dx, dy, dz, dt, u1.data(), v1.data(), w1.data());
            // swap iterates
            std::swap(u0, u1);
            std::swap(v0, v1);
            std::swap(w0, w1);
        });

    // Estimate memory traffic per sweep (reads u_rhs,v_rhs,w_rhs,nu + 6 neighbors of
    // u_iter/v_iter/w_iter; write u_next/v_next/w_next). We keep a simple streaming estimate:
    double bytes_jac = (3.0 + 1.0 + 3.0) * N * sizeof(double); // ~7N reads + 3N writes

    bench::report("fluids_be_jacobi_sweep_64^3", mean_jac, std_jac, bytes_jac);

    // --- Red/black GS (one "iteration" = red + black) ---
    // Start from u0/v0/w0 field
    auto [mean_gs, std_gs] = bench::run(
        [&]
        {
            // red then black, in-place
            diffuse_velocity_be_gs_color_c(u0.data(), v0.data(), w0.data(), u_rhs.data(),
                                           v_rhs.data(), w_rhs.data(), nu_eff.data(), nx_tot,
                                           ny_tot, nz_tot, ng, dx, dy, dz, dt, /*color=*/0);
            diffuse_velocity_be_gs_color_c(u0.data(), v0.data(), w0.data(), u_rhs.data(),
                                           v_rhs.data(), w_rhs.data(), nu_eff.data(), nx_tot,
                                           ny_tot, nz_tot, ng, dx, dy, dz, dt, /*color=*/1);
        });

    // Each color updates ~half interior; rough streaming estimate per color ~ (7N/2 reads + 3N/2
    // writes)
    double bytes_gs = (7.0 + 3.0) * 0.5 * N * sizeof(double) * 2.0; // red+black

    bench::report("fluids_be_gs_rb_64^3", mean_gs, std_gs, bytes_gs);

    return 0;
}
