#include "memory/MemoryManager.hpp"
#include "mesh/Field.hpp"
#include "mesh/Mesh.hpp"
#include "simple_bench.hpp"

int main()
{
    const int n = 64; // 64³ interior cells
    core::mesh::Mesh M{.local = {n, n, n}, .ng = 2};
    auto& mm = core::memory::MemoryManager::instance();
    double* raw = mm.allocate<double>(M.volume_with_ghosts());
    core::mesh::Field<double> F{raw, M.extents(), M.ng};

    auto [mean, stddev] = bench::run(
        [&]
        {
            double acc = 0.0;
            for (int k = 0; k < n; ++k)
                for (int j = 0; j < n; ++j)
                    for (int i = 0; i < n; ++i)
                        acc += F(i, j, k);
            asm volatile("" ::"g"(acc)); // inhibit optimisation
        });

    bench::report("indexer_64cubed", mean, stddev);
    mm.release(raw);
    return 0;
}