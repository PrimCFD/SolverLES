#pragma once
#include "master/plugin/Action.hpp"
#include <memory>

namespace fluids
{

struct Corrector final : core::master::plugin::IAction
{
    Corrector(double rho, double dx, double dy, double dz);

    const core::master::plugin::ActionInfo& info() const noexcept override { return info_; }
    void execute(const core::master::MeshTileView& tile, core::master::FieldCatalog& fields,
                 double dt) override;

  private:
    core::master::plugin::ActionInfo info_;
    double rho_, dx_, dy_, dz_;
};

std::shared_ptr<core::master::plugin::IAction> make_corrector(const core::master::plugin::KV& kv);

} // namespace fluids
