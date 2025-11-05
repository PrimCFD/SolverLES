#pragma once
#include "master/RunContext.hpp"
#include "master/plugin/Action.hpp"
#include <deque>
#include <memory>
#include <vector>

namespace fluids
{

struct Predictor final : core::master::plugin::IAction
{
    enum class Mode
    {
        FE = 0,
        BE = 1,
        ABM3 = 2 // 3rd-order Adams–Bashforth on advection with automatic AB1/AB2 startup
    };
    enum class IMPSolver
    {
        Jacobi,
        RBGS
    };

    Predictor(double rho, double nu, double dx, double dy, double dz, Mode mode);

    // --- BE inner-solver controls
    void set_pred_imp_controls(int iters, double rtol, void* mpi_comm, IMPSolver pred_imp_solver)
    {
        pred_imp_max_iters_ = iters;
        pred_imp_rtol_ = rtol;
        mpi_comm_ = mpi_comm;
        pred_imp_solver_ = pred_imp_solver;
    }

    const core::master::plugin::ActionInfo& info() const noexcept override { return info_; }
    void execute(const core::master::MeshTileView& tile, core::master::FieldCatalog& fields,
                 double dt) override;

    // --- KK advection controls ---
    void set_advect_enabled(bool b) { advect_enabled_ = b; }

    void set_imp_order(int k) { imp_order_ = k; }

    // Full reset (both adv & diff histories) — use sparingly
    void reset_abm_history()
    {
        Nu_hist_.clear();
        Nv_hist_.clear();
        Nw_hist_.clear();
        Du_hist_.clear();
        Dv_hist_.clear();
        Dw_hist_.clear();
        adv_hist_len_ = 0;
        diff_hist_len_ = 0;
    }
    // Advection-only reset (keep diffusion history intact)
    void reset_adv_history()
    {
        Nu_hist_.clear();
        Nv_hist_.clear();
        Nw_hist_.clear();
        adv_hist_len_ = 0;
    }
    // Diffusion-only reset (keep advection history intact)
    void reset_diff_history()
    {
        Du_hist_.clear();
        Dv_hist_.clear();
        Dw_hist_.clear();
        diff_hist_len_ = 0;
    }

  private:
    core::master::plugin::ActionInfo info_;
    double rho_, nu_, dx_, dy_, dz_;
    // scratch arrays for u*, v*, w*
    std::vector<double> us_, vs_, ws_;
    Mode mode_{Mode::FE};
    int pred_imp_max_iters_ = 50;
    double pred_imp_rtol_ = 1e-8;
    void* mpi_comm_ = nullptr;                    // opaque; valid only in MPI builds
    IMPSolver pred_imp_solver_ = IMPSolver::RBGS; // default: RBGS (faster than Jacobi)
    int imp_order_ = 3;                           // 1=BE, 2=CN/AM2, 3=AM3

    // KK scheme parameters
    bool advect_enabled_ = true;

    // -------- ABM3 advection/diffusion history (face-centered, MAC layout) --------

    int adv_hist_len_ = 0, diff_hist_len_ = 0; // 0,1,2,3 available entries

    // We store the last up-to-3 tendencies N(u) at n, n-1, n-2 per component.
    // On first two steps, we automatically use AB1 then AB2.
    // Each entry is a full-sized array matching the corresponding face field.
    std::deque<std::vector<double>> Nu_hist_; // newest at front: [N^n, N^{n-1}, N^{n-2}]
    std::deque<std::vector<double>> Nv_hist_;
    std::deque<std::vector<double>> Nw_hist_;

    // ---- history for diffusion operator application (AM3 RHS terms) ----
    // newest at front: [D(u^n), D(u^{n-1}), ...] where D(u)=L u at faces
    std::deque<std::vector<double>> Du_hist_;
    std::deque<std::vector<double>> Dv_hist_;
    std::deque<std::vector<double>> Dw_hist_;

    // Helper to push a new tendency snapshot, trimming to length 3
    void push_adv_hist(const std::vector<double>& Nu, const std::vector<double>& Nv,
                       const std::vector<double>& Nw);

    void push_diff_hist(const std::vector<double>& Du, const std::vector<double>& Dv,
                        const std::vector<double>& Dw);
};

std::shared_ptr<core::master::plugin::IAction> make_predictor(const core::master::plugin::KV& kv,
                                                              const core::master::RunContext& rc);

} // namespace fluids
