#pragma once
#include "memory/MemoryManager.hpp"
#include <cstddef>
#include <memory>

/**
 * @file RunContext.hpp
 * @brief Runtime handles shared across the orchestration layer.
 *
 * @details
 * `RunContext` carries low-level resources owned by the application: the MPI communicator,
 * a device stream/queue, and a pointer to the MemoryManager. These are passed by reference
 * throughout the core so plugins and writers use the **same stream** and residency policies
 * as the scheduler.
 *
 * The fields are intentionally **opaque** (stored as `void*`) to keep public headers light;
 * feature toggles (`HAVE_CUDA`) guard the actual includes and casts inside `.cpp`.
 *
 * @rst
 * .. code-block:: cpp
 *
 *   core::master::RunContext rc{};
 *   rc.mpi_comm     = reinterpret_cast<void*>(&comm);
 *   rc.device_stream= reinterpret_cast<void*>(cudaStream); // guarded by HAVE_CUDA
 *   rc.mem          = &core::memory::MemoryManager::instance();
 * @endrst
 *
 * @note Treat RunContext as **immutable during a run**. If you need to replace the stream or
 * communicator, rebuild the :cpp:class:`core::master::Scheduler` with a fresh context.
 */

namespace core::master
{

/// Runtime handles shared across the run (MPI, device stream, memory mgr).
struct RunContext
{
    // Opaque handles so headers stay light
    void* mpi_comm{nullptr};      // e.g., MPI_Comm*
    void* device_stream{nullptr}; // e.g., cudaStream_t
    core::memory::MemoryManager* mem{nullptr};

    bool cuda_aware_mpi{false};
    int world_rank{0};
    int world_size{1};
};

} // namespace core::master
