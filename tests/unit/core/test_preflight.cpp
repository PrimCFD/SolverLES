#include "master/Views.hpp"
#include "master/io/Preflight.hpp"
#include "master/io/WritePlan.hpp"
#include <catch2/catch_test_macros.hpp>

using namespace core::master;
using namespace core::master::io;

TEST_CASE("Preflight rejects when RAM/Disk insufficient", "[io][preflight]")
{
    // Build a trivial plan: one 32^3 double field
    AnyFieldView v{};
    v.name = "rho";
    v.host_ptr = reinterpret_cast<void*>(0x1);
    v.elem_size = sizeof(double);
    v.extents = {32, 32, 32};
    v.strides = {(std::ptrdiff_t) sizeof(double), (std::ptrdiff_t)(sizeof(double) * 32),
                 (std::ptrdiff_t)(sizeof(double) * 32 * 32)};
    auto plan = build_write_plan(std::span<const AnyFieldView>(&v, 1), /*override=*/0);

    WriterConfig cfg;
    cfg.backend = WriterConfig::Backend::XDMF; // irrelevant here

    // Deliberately tiny limits
    auto [ok1, msg1] = run_preflight(cfg, plan, /*world=*/1,
                                     /*ram*/ 1 << 16, /*disk*/ 1 << 16);
    REQUIRE_FALSE(ok1);

    // Generous limits
    auto [ok2, msg2] = run_preflight(cfg, plan, /*world=*/1,
                                     /*ram*/ 1ull << 30, /*disk*/ 1ull << 30);
    REQUIRE(ok2);
    (void) msg1;
    (void) msg2;
}
