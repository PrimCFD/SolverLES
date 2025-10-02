#pragma once
#include "master/RunContext.hpp"
#include "master/plugin/Action.hpp"
#include <memory>
#include <vector>

namespace fluids
{

struct Predictor final : core::master::plugin::IAction
{
    enum class Mode
    {
        FE = 0,
        BE = 1
    };
    enum class BESolver
    {
        Jacobi,
        RBGS
    };
    Predictor(double rho, double nu, double dx, double dy, double dz, Mode mode);

    // --- BE inner-solver controls
    void set_be_controls(int iters, double rtol, void* mpi_comm, BESolver be_solver)
    {
        be_max_iters_ = iters;
        be_rtol_ = rtol;
        mpi_comm_ = mpi_comm;
        be_solver_ = be_solver;
    }

    const core::master::plugin::ActionInfo& info() const noexcept override { return info_; }
    void execute(const core::master::MeshTileView& tile, core::master::FieldCatalog& fields,
                 double dt) override;

    void set_adv_blend(double v) { adv_blend_ = std::max(0.0, std::min(1.0, v)); }
    void set_advect_enabled(bool b) { advect_enabled_ = b; }

  private:
    core::master::plugin::ActionInfo info_;
    double rho_, nu_, dx_, dy_, dz_;
    // scratch arrays for u*, v*, w*
    std::vector<double> us_, vs_, ws_;
    Mode mode_{Mode::FE};
    int be_max_iters_ = 50;
    double be_rtol_ = 1e-8;
    void* mpi_comm_ = nullptr;            // opaque; valid only in MPI builds
    BESolver be_solver_ = BESolver::RBGS; // default: RBGS (faster than Jacobi)
    double adv_blend_ = 0.0;
    bool advect_enabled_ = true;
};

std::shared_ptr<core::master::plugin::IAction> make_predictor(const core::master::plugin::KV& kv,
                                                              const core::master::RunContext& rc);

} // namespace fluids
