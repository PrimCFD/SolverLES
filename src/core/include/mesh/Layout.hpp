#pragma once
#include <cstddef>

/**
 * @file Layout.hpp
 * @ingroup memory
 * @brief Indexing utilities and layout strategies for structured 3‑D grids.
 *
 * Provides \c layout::Indexer3D which maps ( \p i, \p j, \p k ) to linear addresses inside
 * ghost‑padded arrays using SoA storage. The indexer is trivially inlinable and
 * shared between C++ and Fortran kernels to guarantee identical addressing.
 *
 * @rst
 *.. code-block:: cpp
 *
 *   // 0‑based, ghosts included
 *   inline std::size_t idx(int i,int j,int k) const noexcept {
 *   return (std::size_t(k)*Ny_tot + j)*Nx_tot + i;
 *   }
 * @endrst
 *
 * Alternative strategies (e.g., AoSoA) can be provided behind the same API.
 */

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
