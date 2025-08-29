#include "memory/MemoryManager.hpp"
#include "simple_bench.hpp"

int main()
{
    constexpr std::size_t N = 1 << 20; // 1 Mi doubles
    constexpr std::size_t BYTES = N * sizeof(double);

    auto& mm = core::memory::MemoryManager::instance();
    auto [mean, stddev] = bench::run(
        [&]
        {
            double* p = mm.allocate<double>(N);
            mm.release(p);
        });

    bench::report("allocate_release_1Mi_doubles", mean, stddev, BYTES);
    return 0;
}