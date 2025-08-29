#pragma once
#include "master/plugin/Action.hpp"
#include <memory>
#include <vector>

/**
 * @file Program.hpp
 * @brief Step planning interface that assembles actions into a plan.
 *
 * @details
 * A program returns a :cpp:struct:`core::master::plugin::StepPlan` for the
 * current step (and may vary by sub-step if needed). The scheduler executes
 * that plan while enforcing halo/BC ordering.
 *
 * @rst
 * .. code-block:: cpp
 *
 *    struct RK3 final : IProgram {
 *      StepPlan plan_step(double dt) override {
 *        StepPlan P{};
 *        // push_back shared_ptr<IAction> for each sub-stage
 *        return P;
 *      }
 *    };
 * @endrst
 */

namespace core::master::plugin
{

// A per-step plan produced by the program
struct StepPlan
{
    std::vector<std::shared_ptr<IAction>> tiled;   // tile-wise actions
    std::vector<std::shared_ptr<IGlobal>> globals; // global actions
};

// A program assembles actions (e.g., “explicit RK3” or “operator-split”).
struct IProgram
{
    virtual ~IProgram() = default;
    virtual StepPlan plan_step(double dt) = 0; // may change across steps if needed
};

} // namespace core::master::plugin
