#include "PressurePoisson.hpp"
#include "master/FieldCatalog.hpp"
#include "master/HaloOps.hpp"
#include <algorithm>
#include <cstdlib>
#include <stdexcept>
#include "kernels_fluids.h"
#ifdef HAVE_MPI
#include <mpi.h>
#endif

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
                                 void* mpi_comm)
    : rho_(rho), dx_(dx), dy_(dy), dz_(dz), iters_(iters), mpi_comm_(mpi_comm)
{
    info_.name = "pressure_poisson";
    info_.phases = plugin::Phase::PostExchange; // unchanged
}

std::shared_ptr<plugin::IAction> make_poisson(const plugin::KV& kv,
                                              const core::master::RunContext& rc)
{
    auto get = [&](const char* k, const char* dflt) -> std::string
    {
        if (auto it = kv.find(k); it != kv.end())
            return it->second;
        return dflt;
    };

    return std::make_shared<PressurePoisson>(
        to_d(get("rho", "1.0"), 1.0), to_d(get("dx", "1.0"), 1.0), to_d(get("dy", "1.0"), 1.0),
        to_d(get("dz", "1.0"), 1.0), to_i(get("p_iters", "80"), 80), rc.mpi_comm);
}

void PressurePoisson::execute(const MeshTileView& tile, FieldCatalog& fields, double dt)
{
    (void) dt;
    if (!fields.contains("u") || !fields.contains("v") || !fields.contains("w") ||
        !fields.contains("p"))
        throw std::runtime_error("[fluids.poisson] u/v/w/p required.");

    // Views
    auto vu = fields.view("u");
    auto vv = fields.view("v");
    auto vw = fields.view("w");
    auto vp = fields.view("p");

    // Use per-field totals; ng for centers (p) derives from tile extents
    const int nx_c_tot = vp.extents[0], ny_c_tot = vp.extents[1], nz_c_tot = vp.extents[2];
    const int nx_i = tile.box.hi[0] - tile.box.lo[0];
    const int ng = (nx_c_tot - nx_i) / 2;

    // Allocate divergence (center-sized)
    const std::size_t Np = (std::size_t) nx_c_tot * ny_c_tot * nz_c_tot;
    if (div_.size() != Np)
        div_.assign(Np, 0.0);

    // Helper: exchange one field using its own extents (MAC-safe, minimal change)
    auto exchange_one = [&](const char* name, const AnyFieldView& v_like)
    {
        core::mesh::Mesh m = *tile.mesh;
        const int nx_tot = v_like.extents[0], ny_tot = v_like.extents[1],
                  nz_tot = v_like.extents[2];
        const int ng_f = (nx_tot - nx_i) / 2;
        m.local = {nx_tot - 2 * ng_f, ny_tot - 2 * ng_f, nz_tot - 2 * ng_f};
        m.ng = ng_f;
        core::master::exchange_named_fields(fields, m, mpi_comm_, {name});
    };

    // Make sure u/v/w halos are valid before divergence (MAC requires face ghosts)
    exchange_one("u", vu);
    exchange_one("v", vv);
    exchange_one("w", vw);

    // Divergence RHS (MAC faces -> centers)
    const double* u = static_cast<const double*>(vu.host_ptr);
    const double* v = static_cast<const double*>(vv.host_ptr);
    const double* w = static_cast<const double*>(vw.host_ptr);
    double* p = static_cast<double*>(vp.host_ptr);

    divergence_mac_c(u, v, w, vu.extents[0], vu.extents[1], vu.extents[2], vv.extents[0],
                     vv.extents[1], vv.extents[2], vw.extents[0], vw.extents[1], vw.extents[2],
                     nx_c_tot, ny_c_tot, nz_c_tot, ng, dx_, dy_, dz_, div_.data());

    // Scale RHS for Poisson type:
    //   variable-coefficient:   ∇·(β ∇p) = (1/dt) div(u*)
    //   constant-coefficient:         Δp = (ρ/dt) div(u*)
    const bool have_rho = fields.contains("rho");
    if (have_rho)
    {
        for (double& x : div_)
            x *= (1.0 / dt);
    }
    else
    {
        for (double& x : div_)
            x *= (rho_ / dt);
    }

    // Jacobi solve at centers with per-iteration halo exchange for p
    if (!tile.mesh)
        throw std::runtime_error("[fluids.poisson] tile.mesh is null.");
    core::mesh::Mesh mesh_p = *tile.mesh;
    mesh_p.local = {nx_c_tot - 2 * ng, ny_c_tot - 2 * ng, nz_c_tot - 2 * ng};
    mesh_p.ng = ng;

    const int kmax = std::max(1, iters_);
    const double eps = 1.0e-300;
    if (iters_ == 0)
    {
        std::fill_n(p, Np, 0.0);
        core::master::exchange_named_fields(fields, mesh_p, mpi_comm_, {"p"});
        return;
    }

    if (have_rho)
    {
        auto vr = fields.view("rho");
        const double* rho_c = static_cast<const double*>(vr.host_ptr);
        beta_.resize(Np);
        for (std::size_t i = 0; i < Np; ++i)
            beta_[i] = 1.0 / std::max(rho_c[i], eps);

        for (int k = 0; k < kmax; ++k)
        {
            poisson_jacobi_varcoef_c(div_.data(), beta_.data(), nx_c_tot, ny_c_tot, nz_c_tot, ng,
                                     dx_, dy_, dz_, /*iters=*/1, p);
            core::master::exchange_named_fields(fields, mesh_p, mpi_comm_, {"p"});
        }
    }
    else
    {
        for (int k = 0; k < kmax; ++k)
        {
            poisson_jacobi_c(div_.data(), nx_c_tot, ny_c_tot, nz_c_tot, ng, dx_, dy_, dz_,
                             /*iters=*/1, p);
            core::master::exchange_named_fields(fields, mesh_p, mpi_comm_, {"p"});
        }
    }

    // Remove global mean on centers (fix Neumann nullspace), then exchange p once more
    {
        const int nx_i_c = nx_c_tot - 2 * ng, ny_i_c = ny_c_tot - 2 * ng,
                  nz_i_c = nz_c_tot - 2 * ng;
        double local_sum = 0.0;
        std::size_t local_cnt = (std::size_t) nx_i_c * ny_i_c * nz_i_c;

        for (int kk = 0; kk < nz_i_c; ++kk)
            for (int jj = 0; jj < ny_i_c; ++jj)
                for (int ii = 0; ii < nx_i_c; ++ii)
                {
                    const std::size_t c =
                        1 + (std::size_t)(ii + ng) +
                        (std::size_t) nx_c_tot * ((std::size_t)(jj + ng) +
                                                  (std::size_t) ny_c_tot * (std::size_t)(kk + ng));
                    local_sum += p[c];
                }

        double global_sum = local_sum;
        unsigned long long global_cnt = (unsigned long long) local_cnt;
#ifdef HAVE_MPI
        if (mpi_comm_)
        {
            MPI_Comm comm = static_cast<MPI_Comm>(mpi_comm_);
            MPI_Allreduce(MPI_IN_PLACE, &global_sum, 1, MPI_DOUBLE, MPI_SUM, comm);
            MPI_Allreduce(MPI_IN_PLACE, &global_cnt, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, comm);
        }
#endif
        const double mean = (global_cnt > 0) ? (global_sum / (double) global_cnt) : 0.0;
        for (int kk = 0; kk < nz_i_c; ++kk)
            for (int jj = 0; jj < ny_i_c; ++jj)
                for (int ii = 0; ii < nx_i_c; ++ii)
                {
                    const std::size_t c =
                        1 + (std::size_t)(ii + ng) +
                        (std::size_t) nx_c_tot * ((std::size_t)(jj + ng) +
                                                  (std::size_t) ny_c_tot * (std::size_t)(kk + ng));
                    p[c] -= mean;
                }
        core::master::exchange_named_fields(fields, mesh_p, mpi_comm_, {"p"});
    }
}

} // namespace fluids
