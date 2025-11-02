#include "master/Log.hpp"
#include "memory/MemoryManager.hpp"
#include "simple_bench.hpp"
#include <cstdlib>
#include <iostream>

#ifdef HAVE_CUDA // macro set when CUDAToolkit_FOUND
#include <cuda_runtime.h>
#endif

int main(int argc, char** argv)
{
    // Logger: rank0-only INFO by default (overridable via SOLVER_LOG=...)
    core::master::logx::init({core::master::logx::Level::Info, /*color*/true, /*rank0_only*/true});

    const std::size_t BYTES = (argc > 1)
                                  ? std::strtoull(argv[1], nullptr, 10) // user-supplied bytes
                                  : (1u << 24);                         // default 16 MiB

    const std::size_t N = BYTES / sizeof(double);
    auto& mm = core::memory::MemoryManager::instance();
    double* host = mm.allocate<double>(N);

#ifndef HAVE_CUDA
    std::cout << "bench_copy skipped: CUDA not enabled\n";
    LOGI("bench_copy: CUDA not enabled; test skipped.\n");
    mm.release(host);
    return 0; // test passes (skipped)
#else
    // Warm-up once (avoids first-touch penalty on UM)
    mm.to_device(host, BYTES);
    cudaDeviceSynchronize();
    mm.to_host(host, BYTES);
    cudaDeviceSynchronize();

    auto [mean, stddev] = bench::run(
        [&]
        {
            mm.to_device(host, BYTES);
            cudaDeviceSynchronize(); // ensure transfer finished
            mm.to_host(host, BYTES);
            cudaDeviceSynchronize();
        });

    // 2×BYTES moved per iteration (H→D + D→H)
    bench::report("copy_roundtrip_" + std::to_string(BYTES >> 20) + "MiB", mean, stddev,
                  2.0 * BYTES);

    mm.release(host);
    return 0;
#endif
}
