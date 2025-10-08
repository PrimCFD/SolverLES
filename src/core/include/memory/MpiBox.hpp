#pragma once
#ifdef HAVE_MPI
#include <mpi.h>
inline void* mpi_box(MPI_Comm c)
{
    return static_cast<void*>(new MPI_Comm(c));
}
inline MPI_Comm mpi_unbox(const void* p)
{
    return p ? *reinterpret_cast<const MPI_Comm*>(p) : MPI_COMM_NULL;
}
inline void mpi_box_free(void* p)
{
    delete reinterpret_cast<const MPI_Comm*>(p);
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