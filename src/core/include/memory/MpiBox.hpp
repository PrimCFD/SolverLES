#pragma once

#ifdef HAVE_MPI
#include <mpi.h>

inline void* mpi_box(MPI_Comm in)
{
    MPI_Comm* p = new MPI_Comm(MPI_COMM_NULL);
    // Create a process-group-equivalent communicator with our MPI
    if (in != MPI_COMM_NULL)
    {
        MPI_Comm_dup(in, p);
    }
    return p;
}
inline MPI_Comm mpi_unbox(const void* p)
{
    return p ? *reinterpret_cast<const MPI_Comm*>(p) : MPI_COMM_NULL;
}
inline void mpi_box_free(void* p)
{
    if (!p)
        return;
    MPI_Comm* pc = reinterpret_cast<MPI_Comm*>(p);
    if (*pc != MPI_COMM_NULL)
    {
        MPI_Comm tmp = *pc;
        MPI_Comm_free(&tmp);
    }
    delete pc;
}
#else

inline void* mpi_box(int)
{
    return nullptr;
}
inline int mpi_unbox(const void*)
{
    return 0;
}
inline void mpi_box_free(void*) {}
#endif