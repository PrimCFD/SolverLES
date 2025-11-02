#pragma once
#include <mpi.h>
#include "Program.hpp"
#include "master/RunContext.hpp"
#include "master/plugin/Action.hpp"
#include <cstddef>
#include <memory>
#include <vector>

namespace fluids
{

struct PPImpl;

struct PressurePoisson final : core::master::plugin::IAction
{
    PressurePoisson(double rho, double dx, double dy, double dz, int iters, void* mpi_comm,
                    const BcTable& bcs);

    const core::master::plugin::ActionInfo& info() const noexcept override { return info_; }
    void execute(const core::master::MeshTileView& tile, core::master::FieldCatalog& fields,
                 double dt) override;

    ~PressurePoisson() override;

    // Set target on ||div(u^{n+1})||_2 used to derive KSP tolerances
    void set_div_tol(double v) { user_div_tol_ = v; }

  private:
    core::master::plugin::ActionInfo info_;
    double rho_, dx_, dy_, dz_;
    int iters_;
    std::vector<double> div_; // scratch divergence
    std::vector<double> beta_;
    void* mpi_comm_ = nullptr;
    std::unique_ptr<PPImpl> impl_;
    BcTable bcs_; // KV-driven boundary conditions (shared with ApplyBCs semantics)
    // User divergence target (defaults to 1e-7 if not provided in YAML)
    double user_div_tol_ = 1e-7;
};

std::shared_ptr<core::master::plugin::IAction> make_poisson(const core::master::plugin::KV& kv,
                                                            const core::master::RunContext& rc);

} // namespace fluids
