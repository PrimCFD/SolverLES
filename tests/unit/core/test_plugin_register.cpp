#include <catch2/catch_test_macros.hpp>
#include "master/PluginHost.hpp"
#include "master/RunContext.hpp"
#include "master/plugin/Program.hpp"  // IProgram, StepPlan
#include "master/plugin/Action.hpp"   // KV

using namespace core::master;

TEST_CASE("PluginHost provides builtin noop program", "[plugin][noop]") {
  RunContext rc{};
  PluginHost host;

  plugin::KV cfg; // empty config
  auto prog = host.make_program("noop", cfg, rc);
  REQUIRE(prog != nullptr);

  auto plan = prog->plan_step(1.0);
  REQUIRE(plan.tiled.empty());
  REQUIRE(plan.globals.empty());
}
