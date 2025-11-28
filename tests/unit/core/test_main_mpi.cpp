#define CATCH_CONFIG_RUNNER
#include <catch2/catch_all.hpp>
#include <mpi.h>

struct MpiSession
{
    MpiSession(int& argc, char**& argv)
    {
        int inited = 0;
        MPI_Initialized(&inited);
        if (!inited)
            MPI_Init(&argc, &argv);
    }
    ~MpiSession()
    {
        int finalized = 0;
        MPI_Finalized(&finalized);
        if (!finalized)
            MPI_Finalize();
    }
};

int main(int argc, char** argv)
{
    MpiSession mpi(argc, argv);
    return Catch::Session().run(argc, argv);
}