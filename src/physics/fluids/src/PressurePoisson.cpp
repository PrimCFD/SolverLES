#include "PressurePoisson.hpp"
#include "master/FieldCatalog.hpp"
#include "master/HaloOps.hpp"
#include <cstdlib>
#include <stdexcept>
#include "kernels_fluids.h"
#ifdef HAVE_MPI
#include <mpi.h>
#endif
#include <algorithm>

using namespace core::master;

namespace fluids
{

static inline double to_d(const std::string& s, double dflt)
{
    char* e = nullptr;
    double v = std::strtod(s.c_str(), &e);
    return (e && *e == 0) ? v : dflt;
}
static inline int to_i(const std::string& s, int dflt)
{
    char* e = nullptr;
    long v = std::strtol(s.c_str(), &e, 10);
    return (e && *e == 0) ? (int) v : dflt;
}

PressurePoisson::PressurePoisson(double rho, double dx, double dy, double dz, int iters,
                                 bool use_rhie_chow, void* mpi_comm)
    : rho_(rho), dx_(dx), dy_(dy), dz_(dz), iters_(iters), use_rhie_chow_(use_rhie_chow),
      mpi_comm_(mpi_comm)
{
    info_.name = "pressure_poisson";
    info_.phases = plugin::Phase::PostExchange;
}

std::shared_ptr<plugin::IAction> make_poisson(const plugin::KV& kv,
                                              const core::master::RunContext& rc)
{
    auto lower = [](std::string s)
    {
        std::transform(s.begin(), s.end(), s.begin(),
                       [](unsigned char c) { return std::tolower(c); });
        return s;
    };

    auto get = [&](const char* k, const char* dflt) -> std::string
    {
        if (auto it = kv.find(k); it != kv.end())
            return it->second;
        return dflt;
    };

    auto to_b = [&](const std::string& s, bool dflt) -> bool
    {
        auto v = lower(s);
        if (v == "1" || v == "true" || v == "on" || v == "yes")
            return true;
        if (v == "0" || v == "false" || v == "off" || v == "no")
            return false;
        return dflt;
    };

    return std::make_shared<PressurePoisson>(
        to_d(get("rho", "1.0"), 1.0), to_d(get("dx", "1.0"), 1.0), to_d(get("dy", "1.0"), 1.0),
        to_d(get("dz", "1.0"), 1.0), to_i(get("p_iters", "80"), 80),
        to_b(get("rhie_chow", "false"), false), rc.mpi_comm);
}

void PressurePoisson::execute(const MeshTileView& tile, FieldCatalog& fields, double dt)
{
    (void) dt;
    if (!fields.contains("u") || !fields.contains("v") || !fields.contains("w") ||
        !fields.contains("p"))
        throw std::runtime_error("[fluids.poisson] u/v/w/p required.");

    auto vu = fields.view("u");
    auto vv = fields.view("v");
    auto vw = fields.view("w");
    auto vp = fields.view("p");

    const int nx_tot = vu.extents[0], ny_tot = vu.extents[1], nz_tot = vu.extents[2];
    const int nx = tile.box.hi[0] - tile.box.lo[0];
    const int ng = (nx_tot - nx) / 2;

    if ((int) div_.size() != nx_tot * ny_tot * nz_tot)
        div_.assign((std::size_t) nx_tot * ny_tot * nz_tot, 0.0);

    core::mesh::Mesh mesh_like = *tile.mesh;
    // Ensure local & ng match the current field
    mesh_like.local = {nx_tot - 2 * ng, ny_tot - 2 * ng, nz_tot - 2 * ng};
    mesh_like.ng = ng;

    core::master::exchange_named_fields(fields, mesh_like, mpi_comm_, {"u", "v", "w"});

    // Build RHS using either plain divergence or Rhie–Chow face velocities.
    // If rho is a field => variable-coefficient Poisson with beta=1/rho and RHS=(1/dt)*div(u*).
    // Else constant-coefficient Poisson with RHS=(rho/dt)*div(u*).
    const std::size_t N = (std::size_t) nx_tot * ny_tot * nz_tot;

    const bool have_rho_field = fields.contains("rho");
    const double eps = 1.0e-300;

    const double* u = static_cast<const double*>(vu.host_ptr);
    const double* v = static_cast<const double*>(vv.host_ptr);
    const double* w = static_cast<const double*>(vw.host_ptr);
    double* p_ptr = static_cast<double*>(vp.host_ptr);

    std::vector<double> rho_buf; // used either as beta source or as const-ρ array for RC
    const double* rho_ptr = nullptr;
    if (have_rho_field)
    {
        auto vr = fields.view("rho");
        rho_ptr = static_cast<const double*>(vr.host_ptr);
    }
    else
    {
        // synthesize constant-ρ buffer if RC is requested
        if (use_rhie_chow_)
        {
            rho_buf.assign(N, rho_);
            rho_ptr = rho_buf.data();
        }
    }

    if (use_rhie_chow_ && rho_ptr)
    {
        // RC divergence (uses p and ρ)
        divergence_rhie_chow_c(u, v, w, p_ptr, rho_ptr, nx_tot, ny_tot, nz_tot, ng, dx_, dy_, dz_,
                               dt, div_.data());
    }
    else
    {
        // Plain central divergence
        divergence_c(u, v, w, nx_tot, ny_tot, nz_tot, ng, dx_, dy_, dz_, div_.data());
    }

    if (have_rho_field)
    {
        // variable-coefficient Poisson: ∇·(β∇p)= (1/dt) div(u*), with β=1/ρ
        for (double& x : div_)
            x *= (1.0 / dt);
    }
    else
    {
        // constant-coefficient: Δp = (ρ/dt) div(u*)
        for (double& x : div_)
            x *= (rho_ / dt);
    }

    // Solve Laplacian(p) = RHS with Jacobi, exchanging halos each sweep
    if (!tile.mesh)
        throw std::runtime_error("[fluids.poisson] tile.mesh is null; Scheduler must set it.");

    // If p_iters == 0, skip solve but ensure halos are consistent for corrector
    if (iters_ == 0)
    {
        std::fill_n(static_cast<double*>(vp.host_ptr),
                    static_cast<std::size_t>(nx_tot) * ny_tot * nz_tot, 0.0);
        core::master::exchange_named_fields(fields, mesh_like, mpi_comm_, {"p"});
        return;
    }

    const int kmax = std::max(1, iters_);
    if (have_rho_field)
    {
        // variable-coefficient: build beta = 1 / rho
        beta_.resize(N);
        for (std::size_t i = 0; i < N; ++i)
            beta_[i] = 1.0 / std::max(rho_ptr[i], eps);
        for (int k = 0; k < kmax; ++k)
        {
            poisson_jacobi_varcoef_c(div_.data(), beta_.data(), nx_tot, ny_tot, nz_tot, ng, dx_,
                                     dy_, dz_, /*iters=*/1, p_ptr);
            core::master::exchange_named_fields(fields, mesh_like, mpi_comm_, {"p"});
        }
    }
    else
    {
        for (int k = 0; k < kmax; ++k)
        {
            poisson_jacobi_c(div_.data(), nx_tot, ny_tot, nz_tot, ng, dx_, dy_, dz_, /*iters=*/1,
                             p_ptr);
            core::master::exchange_named_fields(fields, mesh_like, mpi_comm_, {"p"});
        }
    }

    // Remove global mean (fixes Neumann nullspace), then exchange once more
    {
        const int nx_i = nx_tot - 2 * ng, ny_i = ny_tot - 2 * ng, nz_i = nz_tot - 2 * ng;
        double local_sum = 0.0;
        std::size_t local_cnt = static_cast<std::size_t>(nx_i) * ny_i * nz_i;
        for (int kk = 0; kk < nz_i; ++kk)
            for (int jj = 0; jj < ny_i; ++jj)
                for (int ii = 0; ii < nx_i; ++ii)
                {
                    const std::size_t c =
                        1 + static_cast<std::size_t>(ii + ng) +
                        static_cast<std::size_t>(nx_tot) *
                            (static_cast<std::size_t>(jj + ng) +
                             static_cast<std::size_t>(ny_tot) * static_cast<std::size_t>(kk + ng));
                    local_sum += p_ptr[c];
                }
        double global_sum = local_sum;
        unsigned long long global_cnt = static_cast<unsigned long long>(local_cnt);
#ifdef HAVE_MPI
        if (mpi_comm_)
        {
            MPI_Comm comm = static_cast<MPI_Comm>(mpi_comm_);
            MPI_Allreduce(MPI_IN_PLACE, &global_sum, 1, MPI_DOUBLE, MPI_SUM, comm);
            MPI_Allreduce(MPI_IN_PLACE, &global_cnt, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, comm);
        }
#endif
        const double mean = (global_cnt > 0) ? (global_sum / static_cast<double>(global_cnt)) : 0.0;
        for (int kk = 0; kk < nz_i; ++kk)
            for (int jj = 0; jj < ny_i; ++jj)
                for (int ii = 0; ii < nx_i; ++ii)
                {
                    const std::size_t c =
                        1 + static_cast<std::size_t>(ii + ng) +
                        static_cast<std::size_t>(nx_tot) *
                            (static_cast<std::size_t>(jj + ng) +
                             static_cast<std::size_t>(ny_tot) * static_cast<std::size_t>(kk + ng));
                    p_ptr[c] -= mean;
                }
        core::master::exchange_named_fields(fields, mesh_like, mpi_comm_, {"p"});
    }
}

} // namespace fluids
