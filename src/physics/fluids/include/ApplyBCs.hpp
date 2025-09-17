#pragma once
#include "Program.hpp" // Params, BcSpec, BcTable, parse_params
#include "master/plugin/Action.hpp"
#include <memory>

namespace fluids
{

struct ApplyBCs final : core::master::plugin::IAction
{
    explicit ApplyBCs(const Params&);

    const core::master::plugin::ActionInfo& info() const noexcept override { return info_; }
    void execute(const core::master::MeshTileView& tile, core::master::FieldCatalog& fields,
                 double /*dt*/) override;

  private:
    core::master::plugin::ActionInfo info_;
    BcTable bcs_;
};

std::shared_ptr<core::master::plugin::IAction> make_apply_bcs(const core::master::plugin::KV& kv);

} // namespace fluids
