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

namespace core::mesh
{

struct Mesh
{
    std::array<int, 3> local{};     // interior sizes (this tile)
    std::array<int, 3> global{};    // interior sizes (whole domain)
    std::array<int, 3> global_lo{}; // this tile’s global lower-left-back interior index
    int ng = 0;
    std::array<bool, 3> periodic{false, false, false}; // x,y,z periodicity

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

    // Ghost-inclusive totals for cell-centered arrays
    inline std::array<int, 3> totals_cell() const noexcept { return extents(); }

    // Ghost-inclusive totals for a face-aligned array along axis a=0(x),1(y),2(z).
    // MAC per-tile rule: interior count is (N + 1) along the normal axis.
    inline std::array<int, 3> totals_face_axis(int a) const noexcept
    {
        std::array<int, 3> t = {local[0], local[1], local[2]};
        t[a] += 1; // +1 interior along face normal
        t[0] += 2 * ng;
        t[1] += 2 * ng;
        t[2] += 2 * ng;
        return t;
    }

    // Desired MPI process grid (m,n,p) for DMDA.
    // If any entry is 0, the solver will call MPI_Dims_create() to fill it.
    // If all are >0, they are used verbatim.
    std::array<int,3> proc_grid{0,0,0};
};

} // namespace core::mesh
