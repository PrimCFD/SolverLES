#include "MomentumPredictor.hpp"
#include "master/FieldCatalog.hpp"
#include "master/HaloOps.hpp"
#include <algorithm> // std::transform, std::copy_n
#include <cctype>    // std::tolower
#include <cmath>     // std::sqrt
#include <cstdlib>
#include <cstring> // std::memcpy
#include <stdexcept>
#include "kernels_fluids.h"
#ifdef HAVE_MPI
#include <mpi.h>
#endif

static std::string to_lower(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

// MVP: diffusion-only predictor; convection can be added later

using namespace core::master;

namespace fluids
{

static inline double to_d(const std::string& s, double dflt)
{
    char* e = nullptr;
    double v = std::strtod(s.c_str(), &e);
    return (e && *e == 0) ? v : dflt;
}

Predictor::Predictor(double rho, double nu, double dx, double dy, double dz, Mode mode)
    : rho_(rho), nu_(nu), dx_(dx), dy_(dy), dz_(dz), mode_(mode)
{
    info_.name = "momentum_predictor";
    info_.phases = plugin::Phase::Interior;
}

std::shared_ptr<plugin::IAction> make_predictor(const plugin::KV& kv,
                                                const core::master::RunContext& rc)
{
    auto get = [&](const char* k, const char* dflt) -> std::string
    {
        if (auto it = kv.find(k); it != kv.end())
            return it->second;
        return dflt;
    };

    const std::string scheme = to_lower(get("time_scheme", "fe"));
    const Predictor::Mode mode =
        (scheme == "be" || scheme == "backward_euler") ? Predictor::Mode::BE : Predictor::Mode::FE;

    auto P = std::make_shared<Predictor>(
        to_d(get("rho", "1.0"), 1.0), to_d(get("nu", "1e-3"), 1e-3), to_d(get("dx", "1.0"), 1.0),
        to_d(get("dy", "1.0"), 1.0), to_d(get("dz", "1.0"), 1.0), mode);

    // BE inner solver knobs (default: 50 iters, 1e-8 rtol)
    const int be_iters = (int) std::strtol(get("be_max_iters", "50").c_str(), nullptr, 10);
    const double be_rtol = to_d(get("be_rtol", "1e-8"), 1e-8);
    const std::string be_solver_str = to_lower(get("be_solver", "rbgs")); // "jacobi"|"rbgs"
    Predictor::BESolver be_solver =
        be_solver_str == "jacobi" ? Predictor::BESolver::Jacobi : Predictor::BESolver::RBGS;
    P->set_be_controls(be_iters, be_rtol, rc.mpi_comm,
                       be_solver); // Keep MPI communicator (opaque) for halo_ops
    return P;
}

void Predictor::execute(const MeshTileView& tile, FieldCatalog& fields, double dt)
{
    if (!fields.contains("u") || !fields.contains("v") || !fields.contains("w"))
        throw std::runtime_error("[fluids.predictor] u/v/w required.");

    auto vu = fields.view("u");
    auto vv = fields.view("v");
    auto vw = fields.view("w");

    const int nx_tot = vu.extents[0], ny_tot = vu.extents[1], nz_tot = vu.extents[2];
    const int nx = tile.box.hi[0] - tile.box.lo[0];
    const int ng = (nx_tot - nx) / 2;

    const std::size_t N = static_cast<std::size_t>(nx_tot) * ny_tot * nz_tot;

    // Freeze RHS = u^n, v^n, w^n for this timestep
    std::vector<double> urhs(N), vrhs(N), wrhs(N);
    std::copy_n(static_cast<const double*>(vu.host_ptr), N, urhs.data());
    std::copy_n(static_cast<const double*>(vv.host_ptr), N, vrhs.data());
    std::copy_n(static_cast<const double*>(vw.host_ptr), N, wrhs.data());

    // Effective viscosity array: nu + nu_t (if present)
    std::vector<double> nu_eff(N, nu_);
    if (fields.contains("nu_t"))
    {
        auto vt = fields.view("nu_t");
        const double* t = static_cast<const double*>(vt.host_ptr);
        for (std::size_t q = 0; q < N; ++q)
            nu_eff[q] = nu_ + t[q];
    }

    // Scratch storage for FE or BE next iterate
    if (us_.size() != N)
    {
        us_.assign(N, 0.0);
        vs_.assign(N, 0.0);
        ws_.assign(N, 0.0);
    }

    if (!tile.mesh)
        throw std::runtime_error("[fluids.predictor] tile.mesh is null; Scheduler must set it.");
    core::mesh::Mesh mesh_like = *tile.mesh;
    // Ensure local & ng match the current field
    mesh_like.local = {nx_tot - 2 * ng, ny_tot - 2 * ng, nz_tot - 2 * ng};
    mesh_like.ng = ng;

    double* u_ptr = static_cast<double*>(vu.host_ptr);
    double* v_ptr = static_cast<double*>(vv.host_ptr);
    double* w_ptr = static_cast<double*>(vw.host_ptr);

    if (mode_ == Mode::FE)
    {
        // ---- FE: single explicit update ----
        diffuse_velocity_fe_c(u_ptr, v_ptr, w_ptr, nu_eff.data(), nx_tot, ny_tot, nz_tot, ng, dx_,
                              dy_, dz_, dt, us_.data(), vs_.data(), ws_.data());

        std::memcpy(u_ptr, us_.data(), N * sizeof(double));
        std::memcpy(v_ptr, vs_.data(), N * sizeof(double));
        std::memcpy(w_ptr, ws_.data(), N * sizeof(double));
        core::master::exchange_named_fields(fields, mesh_like, mpi_comm_, {"u", "v", "w"});
        return;
    }

    // ---- BE: host-controlled solve (Jacobi or RBGS), centralized halos ----

    // Initial guess u^(0) := u^n; ensure halos are consistent before first sweep
    core::master::exchange_named_fields(fields, mesh_like, mpi_comm_, {"u", "v", "w"});

    const int kmax = std::max(1, be_max_iters_);
    const double rtol = (be_rtol_ > 0.0 ? be_rtol_ : 1e-8);

    for (int k = 0; k < kmax; ++k)
    {
        if (be_solver_ == BESolver::Jacobi)
        {
            // ---------- Jacobi (as before) ----------
            diffuse_velocity_be_sweep_c(urhs.data(), vrhs.data(), wrhs.data(), u_ptr, v_ptr, w_ptr,
                                        nu_eff.data(), nx_tot, ny_tot, nz_tot, ng, dx_, dy_, dz_,
                                        dt, us_.data(), vs_.data(), ws_.data());
            std::copy_n(us_.data(), N, u_ptr);
            std::copy_n(vs_.data(), N, v_ptr);
            std::copy_n(ws_.data(), N, w_ptr);
            core::master::exchange_named_fields(fields, mesh_like, mpi_comm_, {"u", "v", "w"});
        }
        else
        {
            // ---------- RBGS (in-place color sweeps) ----------
            // Red
            diffuse_velocity_be_gs_color_c(u_ptr, v_ptr, w_ptr, urhs.data(), vrhs.data(),
                                           wrhs.data(), nu_eff.data(), nx_tot, ny_tot, nz_tot, ng,
                                           dx_, dy_, dz_, dt, /*color=*/0);
            core::master::exchange_named_fields(fields, mesh_like, mpi_comm_, {"u", "v", "w"});
            // Black
            diffuse_velocity_be_gs_color_c(u_ptr, v_ptr, w_ptr, urhs.data(), vrhs.data(),
                                           wrhs.data(), nu_eff.data(), nx_tot, ny_tot, nz_tot, ng,
                                           dx_, dy_, dz_, dt, /*color=*/1);
            core::master::exchange_named_fields(fields, mesh_like, mpi_comm_, {"u", "v", "w"});
        }

        // Residual (needs valid halos; we just exchanged)
        double res2 = 0.0, rhs2 = 0.0;
        diffuse_velocity_be_residual_c(urhs.data(), vrhs.data(), wrhs.data(), u_ptr, v_ptr, w_ptr,
                                       nu_eff.data(), nx_tot, ny_tot, nz_tot, ng, dx_, dy_, dz_, dt,
                                       &res2, &rhs2);
#ifdef HAVE_MPI
        if (mpi_comm_)
        {
            MPI_Comm comm = static_cast<MPI_Comm>(mpi_comm_);
            double pair[2] = {res2, rhs2};
            MPI_Allreduce(MPI_IN_PLACE, pair, 2, MPI_DOUBLE, MPI_SUM, comm);
            res2 = pair[0];
            rhs2 = pair[1];
        }
#endif

        if (rhs2 > 0.0 && std::sqrt(res2 / rhs2) < rtol)
            break;
    }
}

} // namespace fluids
