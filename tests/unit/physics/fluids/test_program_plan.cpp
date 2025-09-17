#include "Program.hpp"
#include "master/RunContext.hpp"
#include "master/plugin/Program.hpp"
#include <catch2/catch_test_macros.hpp>

using namespace core::master::plugin;

TEST_CASE("LESProgram creates expected action sequence", "[fluids][program]")
{
    core::master::plugin::KV kv = {{"time_scheme", "fe"}, {"fe_iters", "3"}, {"dx", "1"},
                                   {"dy", "1"},           {"dz", "1"},       {"rho", "1"},
                                   {"nu", "1e-3"},        {"Cs", "0.16"}};
    core::master::RunContext rc{};
    fluids::LESProgram prog(kv, rc);

    auto plan = prog.plan_step(1.0);
    // LESProgram returns a single loop action (plus optional one-time init if configured).
    REQUIRE(plan.tiled.size() == 1);
}