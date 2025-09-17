#include "ProjectionLoop.hpp"
#include "master/FieldCatalog.hpp"
#include "master/HaloOps.hpp"
#include "master/Views.hpp"
#include <algorithm>
#include <cstdlib>
#include <stdexcept>
#include "kernels_fluids.h"
#ifdef HAVE_MPI
#include <mpi.h>
#endif

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

double ProjectionLoop::compute_div_linf(const MeshTileView& tile, FieldCatalog& fields) const
{
    auto vu = fields.view("u");
    auto vv = fields.view("v");
    auto vw = fields.view("w");
    const int nx_tot = vu.extents[0], ny_tot = vu.extents[1], nz_tot = vu.extents[2];
    const int nx = tile.box.hi[0] - tile.box.lo[0];
    const int ng = (nx_tot - nx) / 2;
    // Ensure halos up-to-date for divergence check
    if (!tile.mesh)
        throw std::runtime_error("[fluids.projection] tile.mesh is null; Scheduler must set it.");
    core::mesh::Mesh mesh_like = *tile.mesh;
    mesh_like.local = {nx_tot - 2 * ng, ny_tot - 2 * ng, nz_tot - 2 * ng};
    mesh_like.ng = ng;
    core::master::exchange_named_fields(fields, mesh_like, mpi_comm_, {"u", "v", "w"});
    std::vector<double> div((std::size_t) nx_tot * ny_tot * nz_tot, 0.0);
    divergence_c(static_cast<const double*>(vu.host_ptr), static_cast<const double*>(vv.host_ptr),
                 static_cast<const double*>(vw.host_ptr), nx_tot, ny_tot, nz_tot, ng, opt_.dx,
                 opt_.dy, opt_.dz, div.data());
    double linf = 0.0;
    for (double x : div)
        linf = std::max(linf, std::abs(x));

#ifdef HAVE_MPI
    if (mpi_comm_)
    {
        MPI_Comm comm = static_cast<MPI_Comm>(mpi_comm_);
        double linf_global = 0.0;
        MPI_Allreduce(&linf, &linf_global, 1, MPI_DOUBLE, MPI_MAX, comm);
        linf = linf_global;
    }
#endif

    return linf;
}

void ProjectionLoop::execute(const MeshTileView& tile, FieldCatalog& fields, double dt)
{

    auto vu = fields.view("u");
    const int nx_tot = vu.extents[0], ny_tot = vu.extents[1], nz_tot = vu.extents[2];
    const int nx = tile.box.hi[0] - tile.box.lo[0];
    const int ng = (nx_tot - nx) / 2;

    core::mesh::Mesh mesh_like = *tile.mesh;
    mesh_like.local = {nx_tot - 2 * ng, ny_tot - 2 * ng, nz_tot - 2 * ng};
    mesh_like.ng = ng;

    // Helper: one pressure-correction sweep (no new predictor)
    auto pressure_correction = [&]()
    {
        psolve_->execute(tile, fields, dt);
        corr_->execute(tile, fields, dt);
        core::master::exchange_named_fields(fields, mesh_like, mpi_comm_, {"u", "v", "w"});

        if (bc_)
            bc_->execute(tile, fields, 0.0);
    };

    if (opt_.mode == Mode::FE)
    {
        // SIMPLE-style FE: predictor ONCE, then fe_iters pressure corrections.
        if (bc_)
            bc_->execute(tile, fields, 0.0);
        if (sgs_)
            sgs_->execute(tile, fields, 0.0);
        if (fields.contains("nu_t"))
            core::master::exchange_named_fields(fields, mesh_like, mpi_comm_, {"nu_t"});
        pred_->execute(tile, fields, dt); // diffusion/advection once per time step
        if (bc_)
            bc_->execute(tile, fields, 0.0);

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
        core::master::exchange_named_fields(fields, mesh_like, mpi_comm_, {"nu_t"});
    pred_->execute(tile, fields, dt);
    if (bc_)
        bc_->execute(tile, fields, 0.0);

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
    const std::string scheme = lower(get("time_scheme", "fe")); // "fe" or "be"
    o.mode = (scheme == "be" || scheme == "backward_euler") ? ProjectionLoop::Mode::BE
                                                            : ProjectionLoop::Mode::FE;
    o.fe_iters = to_i(get("fe_iters", "1"), 1);
    o.rtol = to_d(get("nonlinear_rtol", "1e-3"), 1e-3);
    o.max_iters = to_i(get("nonlinear_max_iters", "50"), 50);
    o.dx = to_d(get("dx", "1.0"), 1.0);
    o.dy = to_d(get("dy", "1.0"), 1.0);
    o.dz = to_d(get("dz", "1.0"), 1.0);
    return std::make_shared<ProjectionLoop>(o, sgs, bc, predictor, poisson, corrector, rc.mpi_comm);
}

} // namespace fluids
