#include "memory/MemoryManager.hpp"
#include "mesh/Field.hpp"
#include "mesh/Mesh.hpp"
#include <catch2/catch_all.hpp>

using namespace core;

TEST_CASE("Field indexer matches manual formula", "[field]")
{
    Mesh M{.local = {4, 5, 6}, .ng = 2};
    auto e = M.extents(); // 8,9,10
    auto& mm = MemoryManager::instance();
    double* raw = mm.allocate<double>(M.volume_with_ghosts());
    Field<double> f{raw, e, M.ng};

    REQUIRE(&f(0, 0, 0) == raw + (2 * e[1] + 2) * e[0] + 2); // front-lower-left
    REQUIRE(&f(3, 4, 5) == raw + ((5 + 2) * e[1] + (4 + 2)) * e[0] + (3 + 2));

    mm.release(raw);
}