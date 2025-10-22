#pragma once
#include <petscsys.h>

// Only include/operate on MPI if (a) the project was built with MPI and
// (b) PETSc itself was compiled with MPI support.
#if defined(HAVE_MPI) && defined(PETSC_HAVE_MPI)
#include <mpi.h>
#define PETSC_TESTS_USE_MPI 1
#else
#define PETSC_TESTS_USE_MPI 0
#endif

struct PetscTestGuard
{
    // Ctor that forwards argc/argv so PETSc can parse CLI options.
    PetscTestGuard(int& argc, char**& argv)
    {
#if PETSC_TESTS_USE_MPI
        int mpi_inited = 0;
        MPI_Initialized(&mpi_inited);
        if (!mpi_inited)
        {
            int prov = 0;
            MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &prov);
            owns_mpi_ = true;
        }
#endif

        PetscBool pinit = PETSC_FALSE;
        PetscInitialized(&pinit);
        if (!pinit)
        {
            // Forward argc/argv so PETSc reads flags like -ksp_view, -nx/-ny/-nz, etc.
            PetscInitialize(&argc, &argv, nullptr, nullptr);
            owns_petsc_ = true;
        }
    }
    
    PetscTestGuard()
    {
#if PETSC_TESTS_USE_MPI
        int mpi_inited = 0;
        MPI_Initialized(&mpi_inited);
        if (!mpi_inited)
        {
            int prov = 0;
            MPI_Init_thread(nullptr, nullptr, MPI_THREAD_FUNNELED, &prov);
            owns_mpi_ = true;
        }
#endif

        PetscBool pinit = PETSC_FALSE;
        PetscInitialized(&pinit);
        if (!pinit)
        {
            // No options file: tests don’t rely on it; pass nullptrs.
            PetscInitialize(nullptr, nullptr, nullptr, nullptr);
            owns_petsc_ = true;
        }
    }

    ~PetscTestGuard()
    {
        // Finalize PETSc first
        PetscBool pfin = PETSC_FALSE;
        PetscFinalized(&pfin);
        if (!pfin && owns_petsc_)
        {
            PetscFinalize();
        }

#if PETSC_TESTS_USE_MPI
        int mpi_fin = 0;
        MPI_Finalized(&mpi_fin);
        if (!mpi_fin && owns_mpi_)
        {
            MPI_Finalize();
        }
#endif
    }

    // --- Helpers for tests ---
    // True iff PETSc was compiled with MPI and we can/should provide an MPI_Comm.
    static constexpr bool petsc_uses_mpi() noexcept
    {
#if PETSC_TESTS_USE_MPI
        return true;
#else
        return false;
#endif
    }

#if PETSC_TESTS_USE_MPI
    // Returns &MPI_COMM_WORLD if MPI is initialized; otherwise nullptr.
    const void* mpi_comm_ptr() const noexcept
    {
        int init = 0;
        MPI_Initialized(&init);
        return init ? reinterpret_cast<const void*>(&MPI_COMM_WORLD) : nullptr;
    }
#else
    const void* mpi_comm_ptr() const noexcept { return nullptr; }
#endif

  private:
    bool owns_petsc_ = false;
#if PETSC_TESTS_USE_MPI
    bool owns_mpi_ = false;
#endif
};
