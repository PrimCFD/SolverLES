#include "memory/MemoryManager.hpp"
#include <catch2/catch_all.hpp>

TEST_CASE("device_ptr/host_ptr invariants", "[memory][mirrors]")
{
    auto& mm = core::memory::MemoryManager::instance();

#ifndef HAVE_CUDA
    SECTION("CPU-only: device_ptr is null")
    {
        double* h = mm.allocate<double>(8);
        REQUIRE(mm.device_ptr(h) == nullptr);
        REQUIRE(mm.host_ptr(h) == h);
        mm.release(h);
    }
#else
    const bool um = mm.using_unified_memory();
    double* h = mm.allocate<double>(16);

    SECTION("Forward mirror lookup")
    {
        void* d = mm.device_ptr(h);
        if (um)
        {
            REQUIRE(d == h); // UM: same pointer
        }
        else
        {
            REQUIRE(d != nullptr);
            REQUIRE(d != h); // explicit mirrors
        }
    }
    SECTION("Reverse mirror lookup")
    {
        void* d = mm.device_ptr(h);
        REQUIRE(mm.host_ptr(h) == h);
        REQUIRE(mm.host_ptr(d) == h); // map device back to host
    }

    mm.release(h);
#endif
}
