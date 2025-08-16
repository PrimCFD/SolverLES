#include "AlignedAlloc.hpp"
#include "MemoryManager.hpp"
#include <catch2/catch_all.hpp>
#include <cstdint>

using core::MemoryManager;

TEMPLATE_TEST_CASE("Aligned allocations respect HW boundary", "[memory][alignment]", float, double,
                   int64_t)
{
    auto& mm = MemoryManager::instance();
    const std::size_t N = 257; // an odd count
    TestType* ptr = mm.allocate<TestType>(N);
    REQUIRE(reinterpret_cast<uintptr_t>(ptr) % core::detail::HW_ALIGN == 0);

    mm.release(ptr);
}