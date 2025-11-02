#include "memory/MpiBox.hpp"
#include "simple_bench.hpp"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <mpi.h>
#include <string>
#include <vector>

// PETSc/MPI lifecycle
#include "test_petsc_guard.hpp"
#include <petsc.h>

#include "PressurePoisson.hpp" // fluids::make_poisson
#include "Program.hpp"         // KV, BcTable consumed inside make_poisson
#include "master/FieldCatalog.hpp"
#include "master/RunContext.hpp"
#include "master/Views.hpp"

// Byte-stride helper
static inline std::array<std::ptrdiff_t, 3> strides_bytes(int nx, int ny)
{
    const std::ptrdiff_t s0 = (std::ptrdiff_t) sizeof(double);
    const std::ptrdiff_t s1 = (std::ptrdiff_t) nx * s0;
    const std::ptrdiff_t s2 = (std::ptrdiff_t) nx * (std::ptrdiff_t) ny * s0;
    return {s0, s1, s2};
}

// Build MAC faces so (1/dt) div(u*) = div(∇p*) with β=1
static void build_faces_from_pstar(const std::vector<double>& pstar, std::vector<double>& u,
                                   std::vector<double>& v, std::vector<double>& w, int nx, int ny,
                                   int nz, int ng, double dx, double dy, double dz, double dt)
{
    const int nxc_tot = nx + 2 * ng, nyc_tot = ny + 2 * ng, nzc_tot = nz + 2 * ng;
    const int nxu_tot = nxc_tot + 1, nyu_tot = nyc_tot, nzu_tot = nzc_tot;
    const int nxv_tot = nxc_tot, nyv_tot = nyc_tot + 1, nzv_tot = nzc_tot;
    const int nxw_tot = nxc_tot, nyw_tot = nyc_tot, nzw_tot = nzc_tot + 1;
    auto cidx = [&](int I, int J, int K)
    { return (size_t) I + (size_t) nxc_tot * ((size_t) J + (size_t) nyc_tot * (size_t) K); };
    auto aidx = [&](int I, int J, int K, int nxT, int nyT)
    { return (size_t) I + (size_t) nxT * ((size_t) J + (size_t) nyT * (size_t) K); };
    // u(i+1/2,j,k)
    for (int K = 0; K < nzu_tot; ++K)
        for (int J = 0; J < nyu_tot; ++J)
            for (int I = 0; I < nxu_tot; ++I)
            {
                if (I == 0 || I == nxu_tot - 1)
                {
                    u[aidx(I, J, K, nxu_tot, nyu_tot)] = 0.0;
                    continue;
                }
                const int ic = I - 1, jc = J, kc = K;
                const double pc = pstar[cidx(ic, jc, kc)];
                const double pe = pstar[cidx(ic + 1, jc, kc)];
                u[aidx(I, J, K, nxu_tot, nyu_tot)] = -(dt / dx) * (pe - pc);
            }
    // v(i,j+1/2,k)
    for (int K = 0; K < nzv_tot; ++K)
        for (int J = 0; J < nyv_tot; ++J)
            for (int I = 0; I < nxv_tot; ++I)
            {
                if (J == 0 || J == nyv_tot - 1)
                {
                    v[aidx(I, J, K, nxv_tot, nyv_tot)] = 0.0;
                    continue;
                }
                const int ic = I, jc = J - 1, kc = K;
                const double ps = pstar[cidx(ic, jc, kc)];
                const double pn = pstar[cidx(ic, jc + 1, kc)];
                v[aidx(I, J, K, nxv_tot, nyv_tot)] = -(dt / dy) * (pn - ps);
            }
    // w(i,j,k+1/2)
    for (int K = 0; K < nzw_tot; ++K)
        for (int J = 0; J < nyw_tot; ++J)
            for (int I = 0; I < nxw_tot; ++I)
            {
                if (K == 0 || K == nzw_tot - 1)
                {
                    w[aidx(I, J, K, nxw_tot, nyw_tot)] = 0.0;
                    continue;
                }
                const int ic = I, jc = J, kc = K - 1;
                const double pb = pstar[cidx(ic, jc, kc)];
                const double pt = pstar[cidx(ic, jc, kc + 1)];
                w[aidx(I, J, K, nxw_tot, nyw_tot)] = -(dt / dz) * (pt - pb);
            }
}

int main(int argc, char** argv)
{

    PetscTestGuard petsc_guard(argc, argv); // Honours runtime flags

    // Problem size & geometry (overridable via PETSc options -nx/-ny/-nz)
    int nx = 64, ny = 64, nz = 64;
    const int ng = 1;

    // Allow override via PETSc options: -nx/-ny/-nz
    {
        PetscBool setx = PETSC_FALSE, sety = PETSC_FALSE, setz = PETSC_FALSE;
        (void) PetscOptionsGetInt(NULL, NULL, "-nx", &nx, &setx);
        (void) PetscOptionsGetInt(NULL, NULL, "-ny", &ny, &sety);
        (void) PetscOptionsGetInt(NULL, NULL, "-nz", &nz, &setz);
        (void) PetscPrintf(PETSC_COMM_WORLD, "[bench] grid chosen: nx=%d ny=%d nz=%d\n", nx, ny,
                           nz);
    }

    const int nxc_tot = nx + 2 * ng, nyc_tot = ny + 2 * ng, nzc_tot = nz + 2 * ng;
    const int nxu_tot = nxc_tot + 1, nyu_tot = nyc_tot, nzu_tot = nzc_tot;
    const int nxv_tot = nxc_tot, nyv_tot = nyc_tot + 1, nzv_tot = nzc_tot;
    const int nxw_tot = nxc_tot, nyw_tot = nyc_tot, nzw_tot = nzc_tot + 1;
    const double dx = 1.0, dy = 1.0, dz = 1.0, dt = 1.0;

    // Allocate fields like the tests
    std::vector<double> p((size_t) nxc_tot * nyc_tot * nzc_tot, 0.0);
    std::vector<double> rho((size_t) nxc_tot * nyc_tot * nzc_tot, 1.0);
    std::vector<double> u((size_t) nxu_tot * nyu_tot * nzu_tot, 0.0);
    std::vector<double> v((size_t) nxv_tot * nyv_tot * nzv_tot, 0.0);
    std::vector<double> w((size_t) nxw_tot * nyw_tot * nzw_tot, 0.0);
    std::vector<double> pstar = p; // manufactured p*

    // Manufactured p*
    for (int K = 0; K < nzc_tot; ++K)
        for (int J = 0; J < nyc_tot; ++J)
            for (int I = 0; I < nxc_tot; ++I)
                pstar[(size_t) I +
                      (size_t) nxc_tot * ((size_t) J + (size_t) nyc_tot * (size_t) K)] =
                    std::cos(2 * M_PI * (I - 0.5) / nx) +
                    0.3 * std::cos(2 * M_PI * (J - 0.5) / ny) +
                    0.2 * std::cos(2 * M_PI * (K - 0.5) / nz);

    // Build MAC faces so (1/dt)div(u*) = div(∇p*) with β=1
    build_faces_from_pstar(pstar, u, v, w, nx, ny, nz, ng, dx, dy, dz, dt);

    // One-tile view like the tests
    core::master::MeshTileView tile{};

    core::master::FieldCatalog fields;
    fields.register_scalar("p", p.data(), sizeof(double), {nxc_tot, nyc_tot, nzc_tot},
                           strides_bytes(nxc_tot, nyc_tot), core::master::Stagger::Cell);
    fields.register_scalar("u", u.data(), sizeof(double), {nxu_tot, nyu_tot, nzu_tot},
                           strides_bytes(nxu_tot, nyu_tot), core::master::Stagger::IFace);
    fields.register_scalar("v", v.data(), sizeof(double), {nxv_tot, nyv_tot, nzv_tot},
                           strides_bytes(nxv_tot, nyv_tot), core::master::Stagger::JFace);
    fields.register_scalar("w", w.data(), sizeof(double), {nxw_tot, nyw_tot, nzw_tot},
                           strides_bytes(nxw_tot, nyw_tot), core::master::Stagger::KFace);
    fields.register_scalar("rho", rho.data(), sizeof(double), {nxc_tot, nyc_tot, nzc_tot},
                           strides_bytes(nxc_tot, nyc_tot), core::master::Stagger::Cell);

    // Build the MG Poisson action
    core::master::plugin::KV kv{{"dx", std::to_string(dx)},
                                {"dy", std::to_string(dy)},
                                {"dz", std::to_string(dz)},
                                {"rho", "1.0"},
                                {"iters", "50"},
                                {"div_tol", "1e-10"}};

    // --- Build an MPI Cartesian communicator just like the app does ---
    core::master::RunContext rc{};
    {
        int world_size = 1, world_rank = 0;
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

        // Heuristic: try to form a near-cube 3D grid of ranks
        int dims[3] = {0, 0, 0};
        MPI_Dims_create(world_size, 3, dims); // fills dims in-place
        int periods[3] = {0, 0, 0};           // non-periodic for this bench

        MPI_Comm cart_comm = MPI_COMM_NULL;
        MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, /*reorder*/ 1, &cart_comm);
        if (cart_comm == MPI_COMM_NULL)
            cart_comm = MPI_COMM_WORLD; // fallback

        // Stash the MPI topology/comm into the runtime context so the solver uses it
        rc.mpi_comm = mpi_box(cart_comm);
        rc.world_size = world_size;
        rc.world_rank = world_rank;
    }

    // Time a full action execute (this runs MG inside PressurePoisson)
    // >>> Build the solver ONCE <<<
    auto poisson = fluids::make_poisson(kv, rc);

    // Warm-up once to pay any one-time setup (coarse Pmats via Galerkin, etc.)
    poisson->execute(tile, fields, dt);

    // Allow quick override of iterations: BENCH_ITERS (default 5 for heavy tests)
    int iters = 5;
    if (const char* s = std::getenv("BENCH_ITERS"))
    {
        try
        {
            iters = std::max(1, std::stoi(std::string{s}));
        }
        catch (...)
        {
        }
    }

    // MPI-correct timing: global max across ranks per iteration, then stats
    auto [mean_us, stddev_us] = bench::run_mpi_max(
        MPI_COMM_WORLD,
        [&]
        {
            poisson->execute(tile, fields, dt); // no reconstruction
        },
        iters);
    bench::report_root(MPI_COMM_WORLD, "fluids_poisson_petsc_mg_64^3", mean_us, stddev_us, 0.0);
    return 0;
}