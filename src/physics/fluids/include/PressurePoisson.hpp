#pragma once
#include "master/RunContext.hpp"
#include "master/plugin/Action.hpp"
#include <cstddef>
#include <memory>
#include <vector>

namespace fluids
{

struct PressurePoisson final : core::master::plugin::IAction
{
    PressurePoisson(double rho, double dx, double dy, double dz, int iters, void* mpi_comm);

    const core::master::plugin::ActionInfo& info() const noexcept override { return info_; }
    void execute(const core::master::MeshTileView& tile, core::master::FieldCatalog& fields,
                 double dt) override;

  private:
    core::master::plugin::ActionInfo info_;
    double rho_, dx_, dy_, dz_;
    int iters_;
    std::vector<double> div_; // scratch divergence
    std::vector<double> beta_;
    void* mpi_comm_ = nullptr;
};

std::shared_ptr<core::master::plugin::IAction> make_poisson(const core::master::plugin::KV& kv,
                                                            const core::master::RunContext& rc);

} // namespace fluids
