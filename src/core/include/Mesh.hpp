#pragma once
#include <array>
#include <cstddef>

/**
* @file Mesh.hpp
* @ingroup memory
* @brief Local domain geometry and ghost width.
*
* Holds the interior cell counts and the ghost width \c ng. Utility methods
* return ghost‑inclusive extents and total volume used by memory allocation.
* The mesh also carries the Cartesian communicator used by halo exchange.
*/

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
