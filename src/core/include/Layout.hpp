#pragma once
#include <cstddef>

namespace layout
{

// Linearization with totals INCLUDING ghosts
struct Indexer3D
{
    int nx{}, ny{}, nz{};
    inline std::size_t operator()(int i, int j, int k) const noexcept
    {
        return static_cast<std::size_t>((k * ny + j) * nx + i);
    }
};

} // namespace layout
