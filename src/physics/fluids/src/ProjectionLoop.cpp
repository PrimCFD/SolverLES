#include "ProjectionLoop.hpp"
#include "MacOps.hpp"
#include "master/FieldCatalog.hpp"
#include "master/HaloOps.hpp"
#include "master/Views.hpp"
#include <algorithm>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <stdexcept>
using namespace numerics::kernels;

#include "memory/MpiBox.hpp"

using namespace core::master;
using core::master::plugin::KV;

namespace fluids
{

static inline double to_d(const std::string& s, double dflt)
{
    char* e = nullptr;
    const double v = std::strtod(s.c_str(), &e);
    return (e && *e == 0) ? v : dflt;
}
static inline int to_i(const std::string& s, int dflt)
{
    char* e = nullptr;
    const long v = std::strtol(s.c_str(), &e, 10);
    return (e && *e == 0) ? int(v) : dflt;
}
static inline std::string lower(std::string x)
{
    std::transform(x.begin(), x.end(), x.begin(), [](unsigned char c) { return std::tolower(c); });
    return x;
}

ProjectionLoop::ProjectionLoop(Options opt, std::shared_ptr<plugin::IAction> sgs,
                               std::shared_ptr<plugin::IAction> bc,
                               std::shared_ptr<plugin::IAction> predictor,
                               std::shared_ptr<plugin::IAction> poisson,
                               std::shared_ptr<plugin::IAction> corrector, void* mpi_comm)
    : ActionBase({"projection_loop", plugin::Phase::PostExchange,
                  /*Access=*/{/*reads*/ {}, /*writes*/ {}, /*halos*/ {}}}),
      opt_(opt), sgs_(std::move(sgs)), bc_(std::move(bc)), pred_(std::move(predictor)),
      psolve_(std::move(poisson)), corr_(std::move(corrector)), mpi_comm_(mpi_comm)
{
}

double ProjectionLoop::compute_div_linf(const core::master::MeshTileView& tile,
                                        core::master::FieldCatalog& fields) const
{
    // Views (u,v,w faces; p centers to infer ng)
    auto vu = fields.view("u");
    auto vv = fields.view("v");
    auto vw = fields.view("w");
    auto vp = fields.view("p");

    // --- MAC shape invariants with single Source of Truth (mesh) ---
    if (!tile.mesh)
        throw std::runtime_error("[fluids.predictor] tile.mesh is null; Scheduler must set it.");
    const auto& mesh = *tile.mesh;
    const int ng = mesh.ng;
    const int nxc = vp.extents[0] - 2 * ng;
    const int nyc = vp.extents[1] - 2 * ng;
    const int nzc = vp.extents[2] - 2 * ng;
    auto check = [](const char* nm, int got, int exp)
    {
        if (got != exp)
        {
            throw std::runtime_error(std::string("[MAC layout] Field '") + nm +
                                     "' interior count is " + std::to_string(got) +
                                     " but expected " + std::to_string(exp) +
                                     " (face-centered must be +1 along its normal).");
        }
    };
    const int nxu = vu.extents[0] - 2 * ng, nyu = vu.extents[1] - 2 * ng,
              nzu = vu.extents[2] - 2 * ng;
    const int nxv = vv.extents[0] - 2 * ng, nyv = vv.extents[1] - 2 * ng,
              nzv = vv.extents[2] - 2 * ng;
    const int nxw = vw.extents[0] - 2 * ng, nyw = vw.extents[1] - 2 * ng,
              nzw = vw.extents[2] - 2 * ng;
    check("u(x-face).nx", nxu, nxc + 1);
    check("u.ny", nyu, nyc);
    check("u.nz", nzu, nzc);
    check("v.nx", nxv, nxc);
    check("v(y-face).ny", nyv, nyc + 1);
    check("v.nz", nzv, nzc);
    check("w.nx", nxw, nxc);
    check("w.ny", nyw, nyc);
    check("w(z-face).nz", nzw, nzc + 1);

    // Center totals (divergence is center-located)
    const int nxc_tot = vp.extents[0];
    const int nyc_tot = vp.extents[1];
    const int nzc_tot = vp.extents[2];

    // Build divergence (center-sized) and compute via MAC operator
    std::vector<double> div((std::size_t) nxc_tot * nyc_tot * nzc_tot, 0.0);

    numerics::kernels::divergence(
        static_cast<const double*>(vu.host_ptr), static_cast<const double*>(vv.host_ptr),
        static_cast<const double*>(vw.host_ptr), vu.extents[0], vu.extents[1], vu.extents[2],
        vv.extents[0], vv.extents[1], vv.extents[2], vw.extents[0], vw.extents[1], vw.extents[2],
        nxc_tot, nyc_tot, nzc_tot, ng, opt_.dx, opt_.dy, opt_.dz, div.data());

    // Linf on interior cells only (centers)
    double linf = 0.0;
    for (int k = 0; k < nzc_tot - 2 * ng; ++k)
    {
        const int kk = k + ng;
        for (int j = 0; j < nyc_tot - 2 * ng; ++j)
        {
            const int jj = j + ng;
            std::size_t c = 1u + (std::size_t) ng + // i = 0+ng
                            (std::size_t) nxc_tot *
                                ((std::size_t) jj + (std::size_t) nyc_tot * (std::size_t) kk);

            // walk i interior contiguously
            for (int i = 0; i < nxc_tot - 2 * ng; ++i, ++c)
                linf = std::max(linf, std::abs(div[c]));
        }
    }

    if (mpi_comm_)
    {
        MPI_Comm comm = mpi_unbox(mpi_comm_);
        double linf_global = 0.0;
        MPI_Allreduce(&linf, &linf_global, 1, MPI_DOUBLE, MPI_MAX, comm);
        linf = linf_global;
    }

    return linf;
}

void ProjectionLoop::execute(const MeshTileView& tile, FieldCatalog& fields, double dt)
{

    // Single source of truth for halos/topology
    if (!tile.mesh)
        throw std::runtime_error("[fluids.predictor] tile.mesh is null; Scheduler must set it.");
    const auto& mesh = *tile.mesh;
    const int ng = mesh.ng;

    // Helper: one pressure-correction sweep (no new predictor)
    auto pressure_correction = [&]()
    {
        psolve_->execute(tile, fields, dt);
        core::master::exchange_named_fields(fields, mesh, mpi_comm_, {"p"});
        corr_->execute(tile, fields, dt);
        core::master::exchange_named_fields(fields, mesh, mpi_comm_, {"u", "v", "w"});

        if (bc_)
            bc_->execute(tile, fields, 0.0);
        core::master::exchange_named_fields(fields, mesh, mpi_comm_, {"p", "u", "v", "w"});
    };

    // --------- IPISO mode (predictor once, then fixed # pressure corrections) ---------
    if (opt_.mode == Mode::IPISO)
    {
        if (bc_)
            bc_->execute(tile, fields, 0.0);
        if (sgs_)
            sgs_->execute(tile, fields, 0.0);
        if (fields.contains("nu_t"))
            core::master::exchange_named_fields(fields, mesh, mpi_comm_, {"nu_t"});
        pred_->execute(tile, fields, dt); // momentum predictor once per step
        if (bc_)
            bc_->execute(tile, fields, 0.0);

        core::master::exchange_named_fields(fields, mesh, mpi_comm_, {"u", "v", "w"});

        // Baseline divergence for relative drop criterion
        double r0 = compute_div_linf(tile, fields);
        if (r0 == 0.0)
            r0 = 1.0;

        const int nCorr = std::max(0, opt_.num_corrections);
        for (int m = 0; m < nCorr; ++m)
        {
            pressure_correction();
            const double r = compute_div_linf(tile, fields);
            if (r / r0 <= opt_.div_tol)
                break; // early exit if sufficiently divergence-free
        }
        return;
    }

    // Legacy SIMPLE loops

    if (opt_.mode == Mode::FE)
    {
        // SIMPLE-style FE: predictor ONCE, then fe_iters pressure corrections.
        if (bc_)
            bc_->execute(tile, fields, 0.0);
        if (sgs_)
            sgs_->execute(tile, fields, 0.0);
        if (fields.contains("nu_t"))
            core::master::exchange_named_fields(fields, mesh, mpi_comm_, {"nu_t"});
        pred_->execute(tile, fields, dt); // diffusion/advection once per time step
        if (bc_)
            bc_->execute(tile, fields, 0.0);

        core::master::exchange_named_fields(fields, mesh, mpi_comm_, {"u", "v", "w"});

        const int n = std::max(1, opt_.fe_iters);
        for (int it = 0; it < n; ++it)
            pressure_correction();
        return;
    }

    // ---------- BE mode ----------
    // Do predictor ONCE (true BE solve), then pressure corrections only.
    if (bc_)
        bc_->execute(tile, fields, 0.0);
    if (sgs_)
        sgs_->execute(tile, fields, 0.0);
    if (fields.contains("nu_t"))
        core::master::exchange_named_fields(fields, mesh, mpi_comm_, {"nu_t"});
    pred_->execute(tile, fields, dt);
    if (bc_)
        bc_->execute(tile, fields, 0.0);

    core::master::exchange_named_fields(fields, mesh, mpi_comm_, {"u", "v", "w"});

    double r0 = compute_div_linf(tile, fields);
    if (r0 == 0.0)
        r0 = 1.0;
    int it = 0;
    while (it < std::max(1, opt_.max_iters))
    {
        pressure_correction();
        ++it;
        const double r = compute_div_linf(tile, fields);
        if (r / r0 <= opt_.rtol)
            break;
    }

    return;
}

std::shared_ptr<plugin::IAction> make_projection_loop(
    const KV& kv, const core::master::RunContext& rc, std::shared_ptr<plugin::IAction> sgs,
    std::shared_ptr<plugin::IAction> bc, std::shared_ptr<plugin::IAction> predictor,
    std::shared_ptr<plugin::IAction> poisson, std::shared_ptr<plugin::IAction> corrector)
{
    auto get = [&](const char* k, const char* dflt) -> std::string
    {
        if (auto it = kv.find(k); it != kv.end())
            return it->second;
        return dflt;
    };
    ProjectionLoop::Options o;
    const std::string scheme = lower(get("scheme", "fe")); // "fe", "be", or "ipiso"
    if (scheme == "be" || scheme == "backward_euler")
        o.mode = ProjectionLoop::Mode::BE;
    else if (scheme == "ipiso" || scheme == "piso")
        o.mode = ProjectionLoop::Mode::IPISO;
    else
        o.mode = ProjectionLoop::Mode::FE;
    o.fe_iters = to_i(get("fe_iters", "1"), 1);
    o.rtol = to_d(get("nonlinear_rtol", "1e-3"), 1e-3);
    o.max_iters = to_i(get("nonlinear_max_iters", "50"), 50);
    o.dx = to_d(get("dx", "1.0"), 1.0);
    o.dy = to_d(get("dy", "1.0"), 1.0);
    o.dz = to_d(get("dz", "1.0"), 1.0);
    // IPISO extras
    o.num_corrections = to_i(get("num_corrections", "2"), 2);
    o.div_tol = to_d(get("div_tol", "1e-7"), 1e-7);
    return std::make_shared<ProjectionLoop>(o, sgs, bc, predictor, poisson, corrector, rc.mpi_comm);
}

} // namespace fluids
