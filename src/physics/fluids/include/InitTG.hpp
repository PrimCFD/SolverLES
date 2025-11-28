#pragma once
#include "master/RunContext.hpp"
#include "master/plugin/Action.hpp"
#include <memory>
#include <string>

namespace fluids
{

struct InitTG final : core::master::plugin::IAction
{
    InitTG(double Lx, double Ly, double Lz, double U0, void* mpi_comm);

    const core::master::plugin::ActionInfo& info() const noexcept override { return info_; }
    void execute(const core::master::MeshTileView& tile, core::master::FieldCatalog& fields,
                 double /*dt*/) override;

  private:
    core::master::plugin::ActionInfo info_;
    double Lx_, Ly_, Lz_, U0_;
    void* mpi_comm_{nullptr};
};

std::shared_ptr<core::master::plugin::IAction> make_init_tg(const core::master::plugin::KV& kv,
                                                            const core::master::RunContext& rc);

} // namespace fluids
