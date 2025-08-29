#include "memory/MemoryManager.hpp"
#include <catch2/catch_all.hpp>
#include <cstring>
#include <vector>

#ifdef HAVE_CUDA
#include <cuda_runtime.h>
#endif

TEST_CASE("H2D/D2H path preserves bytes", "[memory][copy]")
{
    auto& mm = core::MemoryManager::instance();
    constexpr std::size_t N = 1 << 20;
    constexpr std::size_t BYTES = N * sizeof(double);

    // Fill with deterministic pattern
    std::vector<unsigned char> ref(BYTES);
    for (std::size_t i = 0; i < BYTES; ++i)
        ref[i] = static_cast<unsigned char>(i * 131u + 7u);

    double* h = mm.allocate<double>(N);
    std::memcpy(h, ref.data(), BYTES);

    // Call the API regardless of CUDA availability:
    mm.to_device(h, BYTES);
#ifdef HAVE_CUDA
    cudaDeviceSynchronize();
#endif

    // Clobber host, bring back
    std::memset(h, 0xCD, BYTES);
    mm.to_host(h, BYTES);
#ifdef HAVE_CUDA
    cudaDeviceSynchronize();
#endif

    // CPU-only: to_device/to_host are no-ops → buffer stays clobbered?
    // No: we copied ref in before clobber; without CUDA, to_host is a no-op too,
    // so h remains 0xCD and would fail. Instead, verify contract:
    // - If CUDA: round-trip reproduces ref.
    // - If no CUDA: to_device/to_host are no-ops → the data after "to_host"
    //   must equal what we last wrote to *host* (which is 0xCD).
#ifndef HAVE_CUDA
    std::vector<unsigned char> cd(BYTES, 0xCD);
    REQUIRE(std::memcmp(h, cd.data(), BYTES) == 0);
#else
    REQUIRE(std::memcmp(h, ref.data(), BYTES) == 0);
#endif

    mm.release(h);
}
