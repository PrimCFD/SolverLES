#include "Program.hpp"
#include "ProjectionLoop.hpp"
#include "master/RunContext.hpp"
#include "master/plugin/Program.hpp"
#include "master/plugin/Registry.hpp"

using core::master::RunContext;
using core::master::plugin::KV;
using core::master::plugin::Registry;

extern "C" bool physics_register_v1(Registry* R)
{
    if (!R)
        return false;

    // Program: Pressure scheme + Smagorinsky SGS
    R->add_program("les", [](const KV& kv, const RunContext& rc)
                   { return std::make_unique<fluids::LESProgram>(kv, rc); });

    // (Optional) Also expose individual actions for experimentation/debug.
    // Users typically won’t instantiate them directly.
    R->add_action("fluids.init_tg",
                  [](const KV& kv, const RunContext&) { return fluids::make_init_tg(kv); });
    R->add_action("fluids.sgs",
                  [](const KV& kv, const RunContext&) { return fluids::make_sgs(kv); });
    R->add_action("fluids.momentum_predictor", [](const KV& kv, const RunContext& rc)
                  { return fluids::make_predictor(kv, rc); });
    R->add_action("fluids.pressure_poisson",
                  [](const KV& kv, const RunContext& rc) { return fluids::make_poisson(kv, rc); });
    R->add_action("fluids.velocity_corrector",
                  [](const KV& kv, const RunContext&) { return fluids::make_corrector(kv); });

    R->add_action("fluids.projection_loop",
                  [](const KV& kv, const RunContext& rc)
                  {
                      // build a loop with fresh sub-actions; handy for A/B testing
                      return fluids::make_projection_loop(
                          kv, rc, fluids::make_sgs(kv), fluids::make_apply_bcs(kv),
                          fluids::make_predictor(kv, rc), fluids::make_poisson(kv, rc),
                          fluids::make_corrector(kv));
                  });

    return true;
}
