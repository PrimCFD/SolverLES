#pragma once
#include "master/RunContext.hpp"
#include "master/plugin/Action.hpp"
#include <memory>
#include <vector>

namespace core
{
namespace master
{
    struct RunContext;
    struct MeshTileView;
    class FieldCatalog;
} // namespace master
} // namespace core

namespace fluids
{

class ProjectionLoop : public core::master::plugin::ActionBase
{
  public:
    enum class Mode
    {
        FE,
        BE
    };

    struct Options
    {
        Mode mode = Mode::FE;
        int fe_iters = 1;              // for FE
        double rtol = 1e-3;            // for BE
        int max_iters = 50;            // for BE
        double dx = 1, dy = 1, dz = 1; // for residual (div) kernel
    };

    ProjectionLoop(Options opt, std::shared_ptr<core::master::plugin::IAction> sgs,
                   std::shared_ptr<core::master::plugin::IAction> bc,
                   std::shared_ptr<core::master::plugin::IAction> predictor,
                   std::shared_ptr<core::master::plugin::IAction> poisson,
                   std::shared_ptr<core::master::plugin::IAction> corrector, void* mpi_comm);

    void execute(const core::master::MeshTileView& tile, core::master::FieldCatalog& fields,
                 double dt) override;

  private:
    Options opt_;
    std::shared_ptr<core::master::plugin::IAction> sgs_, bc_, pred_, psolve_, corr_;
    void* mpi_comm_ = nullptr;
    // scratch for residuals
    std::vector<double> div_;
    double compute_div_linf(const core::master::MeshTileView&, core::master::FieldCatalog&) const;
};

std::shared_ptr<core::master::plugin::IAction>
make_projection_loop(const core::master::plugin::KV& kv, const core::master::RunContext& rc,
                     std::shared_ptr<core::master::plugin::IAction> sgs,
                     std::shared_ptr<core::master::plugin::IAction> bc,
                     std::shared_ptr<core::master::plugin::IAction> predictor,
                     std::shared_ptr<core::master::plugin::IAction> poisson,
                     std::shared_ptr<core::master::plugin::IAction> corrector);

} // namespace fluids
