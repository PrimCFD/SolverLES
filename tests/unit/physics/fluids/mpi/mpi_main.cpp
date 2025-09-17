#define CATCH_CONFIG_RUNNER
#include <catch2/catch_session.hpp>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

int main(int argc, char* argv[])
{
#ifdef HAVE_MPI
    int init = 0;
    MPI_Initialized(&init);
    if (!init)
    {
        int provided = 0;
        MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    }
#endif

    const int rc = Catch::Session().run(argc, argv);

#ifdef HAVE_MPI
    int fin = 0;
    MPI_Finalized(&fin);
    if (!fin)
    {
        MPI_Barrier(MPI_COMM_WORLD); // keep output tidy
        MPI_Finalize();
    }
#endif
    return rc;
}