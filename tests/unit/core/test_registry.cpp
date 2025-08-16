#include "MemoryManager.hpp"
#include <catch2/catch_all.hpp>

using core::MemoryManager;

TEST_CASE("MemoryManager registry tracks blocks", "[memory][registry]")
{
    auto& mm = MemoryManager::instance();
    auto* p1 = mm.allocate<double>(100);
    auto* p2 = mm.allocate<double>(50);

    // internal API only visible in tests
    REQUIRE(mm.debug_count() == 2);

    mm.release(p1);
    mm.release(p2);
    REQUIRE(mm.debug_count() == 0);
}