#include "lin.hpp"
#include <algorithm>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <limits>
#include <mpi.h>
#include <vector>
#include "kernels_fluids.h" // Fortran C-bindings

using Catch::Approx;

static inline void allreduce_max(double& v)
{
    double g = 0.0;
    MPI_Allreduce(&v, &g, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    v = g;
}
static inline void allreduce_min(double& v)
{
    double g = 0.0;
    MPI_Allreduce(&v, &g, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    v = g;
}

TEST_CASE("SGS Smagorinsky is globally constant for uniform velocity", "[mpi][fluids][sgs]")
{
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Local domain + halos
    const int ng = 1;
    const int nx = 10, ny = 8, nz = 6; // interior per-rank
    const int nx_tot = nx + 2 * ng, ny_tot = ny + 2 * ng, nz_tot = nz + 2 * ng;

    const double dx = 1.0, dy = 1.0, dz = 1.0;
    const double Cs = 0.17;

    const std::size_t Ntot = static_cast<std::size_t>(nx_tot) * ny_tot * nz_tot;
    std::vector<double> u(Ntot), v(Ntot), w(Ntot), nu_t(Ntot, -999.0);

    // Fill entire array (including halos) uniformly so ghosts are valid without exchanges
    std::fill(u.begin(), u.end(), 3.14);
    std::fill(v.begin(), v.end(), -2.0);
    std::fill(w.begin(), w.end(), 0.5);

    // Compute model viscosity on INTERIOR only
    sgs_smagorinsky_c(u.data(), v.data(), w.data(), nx_tot, ny_tot, nz_tot, ng, dx, dy, dz, Cs,
                      nu_t.data());

    bool any_written = false;
    for (int K = 0; K < nz; ++K)
        for (int J = 0; J < ny; ++J)
            for (int I = 0; I < nx; ++I)
                any_written |= (nu_t[lin(I + ng, J + ng, K + ng, nx_tot, ny_tot)] != -999.0);
    REQUIRE(any_written);

    // Local stats over interior
    double lmin = std::numeric_limits<double>::infinity();
    double lmax = -std::numeric_limits<double>::infinity();
    double lerr = 0.0; // max |nu_t|

    for (int K = 0; K < nz; ++K)
        for (int J = 0; J < ny; ++J)
            for (int I = 0; I < nx; ++I)
            {
                const std::size_t c = lin(I + ng, J + ng, K + ng, nx_tot, ny_tot);
                lmin = std::min(lmin, nu_t[c]);
                lmax = std::max(lmax, nu_t[c]);
                lerr = std::max(lerr, std::abs(nu_t[c]));
            }

    // Global aggregates
    double gmin = lmin, gmax = lmax, gerr = lerr;
    allreduce_min(gmin);
    allreduce_max(gmax);
    allreduce_max(gerr);

    if (rank == 0)
    {
        INFO("gmin=" << gmin << ", gmax=" << gmax << ", gerr=" << gerr);
    }

    REQUIRE(gerr <= 1e-12);
    REQUIRE(std::abs(gmax - gmin) <= 1e-12);
}