#pragma once
#include <array>
#include <cstddef>

namespace core
{

struct Mesh
{
    std::array<int, 3> local{}; // interior sizes
    int ng = 0;

    std::array<int, 3> extents() const
    {
        return {local[0] + 2 * ng, local[1] + 2 * ng, local[2] + 2 * ng};
    }
    std::size_t volume_with_ghosts() const
    {
        const auto e = extents();
        return static_cast<std::size_t>(e[0]) * static_cast<std::size_t>(e[1]) *
               static_cast<std::size_t>(e[2]);
    }
};

} // namespace core
