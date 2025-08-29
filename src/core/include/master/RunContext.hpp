#pragma once
#include <cstddef>
#include <memory>

/**
 * @file RunContext.hpp
 * @brief Runtime handles shared across the orchestration layer.
 *
 * @details
 * Carries low-level resources owned by the application: MPI communicator,
 * device stream/queue, and a pointer to the MemoryManager. These are passed
 * by reference throughout the core so plugins can use the same stream and
 * residency policies as the scheduler.
 *
 * @rst
 * .. tip::
 *    Treat :cpp:struct:`core::master::RunContext` as immutable during a run.
 *    If you need to replace the stream or communicator, rebuild the
 *    :cpp:class:`core::master::Scheduler`.
 * @endrst
 */

namespace core
{
namespace memory
{
    class MemoryManager;
}

namespace master
{

    /// Runtime handles shared across the run (MPI, device stream, memory mgr).
    struct RunContext
    {
        // Opaque handles so headers stay light; define in the .cpp if needed.
        void* mpi_comm{nullptr};      // e.g., MPI_Comm*
        void* device_stream{nullptr}; // e.g., cudaStream_t
        core::memory::MemoryManager* mem{nullptr};

        bool cuda_aware_mpi{false};
        int world_rank{0};
        int world_size{1};
    };

} // namespace master
} // namespace core
