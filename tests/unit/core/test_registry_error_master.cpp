#include "master/RunContext.hpp"
#include "master/plugin/Action.hpp" // KV
#include "master/plugin/Registry.hpp"
#include <catch2/catch_test_macros.hpp>
#include <stdexcept>

TEST_CASE("Registry throws on unknown program/action/global keys", "[registry]")
{
    core::master::plugin::Registry reg;
    core::master::RunContext rc{};
    core::master::plugin::KV kv; // empty

    REQUIRE_THROWS_AS(reg.make_program("does_not_exist", kv, rc), std::runtime_error);
    REQUIRE_THROWS_AS(reg.make_action("does_not_exist", kv, rc), std::runtime_error);
    REQUIRE_THROWS_AS(reg.make_global("does_not_exist", kv, rc), std::runtime_error);
}
