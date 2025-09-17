#pragma once
#include "master/plugin/Action.hpp"
#include <memory>

namespace fluids
{

struct SGS final : core::master::plugin::IAction
{
    SGS(double dx, double dy, double dz, double Cs);

    const core::master::plugin::ActionInfo& info() const noexcept override { return info_; }
    void execute(const core::master::MeshTileView& tile, core::master::FieldCatalog& fields,
                 double /*dt*/) override;

  private:
    core::master::plugin::ActionInfo info_;
    double dx_, dy_, dz_, Cs_;
    // scratch only if nu_t not provided by app
    std::vector<double> scratch_;
    bool wrote_to_field_{false};
};

std::shared_ptr<core::master::plugin::IAction> make_sgs(const core::master::plugin::KV& kv);

} // namespace fluids
