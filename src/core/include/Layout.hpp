#pragma once
#include <array>
#include <cstddef>

namespace layout
{

// Flat 3‑D indexer generator – header‑only, constexpr‑friendly.
struct Indexer3D
{
    std::size_t nx, ny, nz; // with ghosts
    constexpr std::size_t operator()(int i, int j, int k) const noexcept
    {
        return (static_cast<std::size_t>(k) * ny + j) * nx + i;
    }
};

constexpr inline Indexer3D make_indexer(const std::array<int, 3>& ext)
{
    return Indexer3D{static_cast<std::size_t>(ext[0]), static_cast<std::size_t>(ext[1]),
                     static_cast<std::size_t>(ext[2])};
}

} // namespace layout