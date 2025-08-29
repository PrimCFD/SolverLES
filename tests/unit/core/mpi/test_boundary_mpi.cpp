#include <catch2/catch_test_macros.hpp>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include "mesh/Boundary.hpp"

TEST_CASE("is_physical_face sentinel under single-rank MPI", "[boundary][mpi]")
{
#ifdef HAVE_MPI
    int init = 0;
    MPI_Initialized(&init);
    REQUIRE(init == 1);

    int size = 0, rank = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    REQUIRE(size >= 1);

    // A missing neighbor is represented by MPI_PROC_NULL (negative).
    REQUIRE(core::mesh::is_physical_face(MPI_PROC_NULL));

    // This is *not* a physical face sentinel; a neighbor rank (even self)
    // is a real neighbor in decomposition terms.
    REQUIRE_FALSE(core::mesh::is_physical_face(rank));
#else
    // Without MPI we still expect negative to indicate "physical"
    REQUIRE(core::mesh::is_physical_face(-1));
    REQUIRE_FALSE(core::mesh::is_physical_face(0));
#endif
}
