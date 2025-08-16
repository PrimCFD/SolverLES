#pragma once
#include <array>
#ifdef HAVE_MPI
#include <mpi.h>
#endif

namespace core
{

struct Mesh
{
    std::array<int, 3> local{}; // interior cell counts
    int ng = 0;                 // ghost width
#ifdef HAVE_MPI
    MPI_Comm cart = MPI_COMM_NULL;
#endif
    constexpr std::array<int, 3> extents() const
    {
        return {local[0] + 2 * ng, local[1] + 2 * ng, local[2] + 2 * ng};
    }
    constexpr std::size_t volume_with_ghosts() const
    {
        auto e = extents();
        return static_cast<std::size_t>(e[0]) * e[1] * e[2];
    }
};

} // namespace core