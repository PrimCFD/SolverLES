#pragma once
#include "AlignedAlloc.hpp"
#include <cstddef>
#include <cstdint>
#include <mutex>
#include <unordered_map>

/**
 * @defgroup memory Memory & Data‑Layout Subsystem
 * @brief Ownership, transfers, layout and halo exchange for structured grids.
 *
 * The memory layer provides:
 * - **Structure‑of‑Arrays (SoA)** global layout for stencil‑friendly access.
 * - **Ghost/halo padding** around local domains for branch‑free kernels.
 * - **Explicit host↔device mirrors** or **CUDA Unified Memory** (UM) paths.
 * - **64‑byte host alignment** (cache‑/SIMD‑friendly) and optional pinned host buffers.
 * - **Async transfers & prefetch** with streams; overlap with MPI halo exchange.
 *
 * Compile‑time feature toggles:
 * - \c HAVE_CUDA – enable GPU code paths.
 * - \c USE_CUDA_UM – use cudaMallocManaged + prefetch instead of explicit mirrors.
 * - \c USE_PINNED_HOST – allocate page‑locked host staging buffers for faster H2D/D2H.
 *
 * @note All raw pointers returned by this layer are **owned** by the MemoryManager.
 * Non‑owning *views* (e.g., Field<T>) expose typed access without taking ownership.
 * @see MemoryManager, Field, layout::Indexer3D, Mesh, HaloExchange
 */

/**
 * @page memory_overview Memory & Data‑Layout Overview
 * @ingroup memory
 *
 * ## Goals
 * - Keep hot loops pointer‑fast and vectorizable on CPU and coalesced on GPU.
 * - Provide one ownership point (RAII) and leak‑free lifetime.
 * - Make halo exchange contiguous and overlap it with compute.
 *
 * ## Layout
 * **SoA** per primitive: \c rho[], \c ux[], \c uy[], \c uz[], \c p[], … allocated with ghosts.
 * Linear index (0‑based, ghosts included):
 * `idx(i,j,k) = (k * Ny_tot + j) * Nx_tot + i`.
 *
 * ## Alignment
 * Host allocations are aligned (≥64 B). On GPU, UM or explicit mirrors are used.
 * Optional **pinned host** buffers improve PCIe/NVLink throughput.
 *
 * ## CPU/GPU
 * - **Unified Memory**: one pointer, on‑demand migration + optional prefetch.
 * - **Explicit mirrors**: \c to_device() / \c to_host() async copies into device buffers.
 *
 * ## Halos
 * Each local field is \c (Nx+2*ng)×(Ny+2*ng)×(Nz+2*ng). Six face slabs are exchanged with
 * neighbors. Ghost padding turns boundary updates into unit‑stride copies and simplifies kernels.
 */

/**
 * @file MemoryManager.hpp
 * @ingroup memory
 * @brief Singleton RAII owner of all solver buffers (host and device).
 *
 * The MemoryManager owns every allocated block and tracks it in an internal
 * registry. It supports two back‑ends:
 * 1. **Unified Memory (UM)** – single pointer via \c cudaMallocManaged with
 * optional \c cudaMemPrefetchAsync.
 * 2. **Explicit mirrors** – separate host/device buffers with async H2D/D2H.
 *
 * ### Thread‑safety
 * All public methods are thread‑safe. The registry is protected by a mutex.
 *
 * ### Typical use
 * @rst
 *.. code-block:: cpp
 *
 *   auto& mm = MemoryManager::instance();
 *   // Allocate SoA arrays with ghosts
 *   double* rho = mm.allocate<double>(nx_tot*ny_tot*nz_tot);
 *   // Transfer before a GPU kernel when not using UM
 *   mm.to_device(rho, bytes, stream);
 *   // On shutdown or scope end, release once
 *   mm.release(rho);
 * @endrst
 *
 * @note Non‑owning views (Field<T>) should wrap the returned pointer for safer access.
 */

/**
 * @class MemoryManager
 * @ingroup memory
 * @brief Centralized allocator + transfer engine.
 *
 * @par Features
 * - \c allocate\<T\>(n) – aligned host allocation; may also create a device mirror.
 * - \c to_device(ptr, bytes, stream) / \c to_host(ptr, bytes, stream) – async transfer or prefetch.
 * - \c device_ptr(host_ptr) – returns device mirror when using explicit mirrors.
 * - \c host_ptr(maybe_device) – returns owning host pointer for a device mirror.
 * - \c release(ptr) – frees host/device memory and unregisters the block.
 * - \c using_unified_memory() – compile‑time/runtime policy check.
 *
 * @warning Pointers are invalid after \c release. Views must not outlive the block.
 */

#ifdef HAVE_CUDA
#include <cuda_runtime.h>
#endif

namespace core
{

enum class AllocKind : uint8_t
{
    HostAligned,
    HostPinned,
    Device,
    Unified
};

struct Block
{
    void* host = nullptr;   // host-visible pointer
    void* device = nullptr; // device mirror when !UM (UM: same as host)
    std::size_t bytes = 0;
    AllocKind kind = AllocKind::HostAligned;
};

class MemoryManager
{
  public:
    static MemoryManager& instance();

    template <class T> T* allocate(std::size_t n);

    void release(void* host_ptr) noexcept;

    // Explicit H<->D copies. No-ops (or prefetch) in UM / CPU-only builds.
    void to_device(const void* host_ptr, std::size_t bytes, void* stream = nullptr);
    void to_host(const void* host_ptr, std::size_t bytes, void* stream = nullptr);

    // Query mirrors
    void* device_ptr(const void* host_ptr) const noexcept;
    void* host_ptr(const void* maybe_device) const noexcept;

    bool using_unified_memory() const noexcept;

    // Introspection for tests
    std::size_t debug_count() const noexcept;

    ~MemoryManager();

  private:
    MemoryManager() = default;
    MemoryManager(const MemoryManager&) = delete;
    MemoryManager& operator=(const MemoryManager&) = delete;

    Block* find_block_unlocked(const void* host_ptr) noexcept;
    const Block* find_block_unlocked(const void* host_ptr) const noexcept;

    mutable std::mutex mtx_;
    std::unordered_map<const void*, Block> registry_; // keyed by host pointer
};

// ---- template implementation ----

template <class T> T* MemoryManager::allocate(std::size_t n)
{
    const std::size_t bytes = n * sizeof(T);
    Block blk{};

#ifdef HAVE_CUDA
    if (using_unified_memory())
    {
        blk.kind = AllocKind::Unified;
        blk.host = core::detail::aligned_malloc(bytes); // UM via cudaMallocManaged
        blk.device = blk.host;
    }
    else
    {
#ifdef USE_PINNED_HOST
        // Prefer pinned host memory for faster H2D/D2H and MPI staging
        void* h = nullptr;
        if (cudaHostAlloc(&h, bytes, cudaHostAllocDefault) != cudaSuccess)
            throw std::bad_alloc{};
        blk.kind = AllocKind::HostPinned;
        blk.host = h;
#else
        blk.kind = AllocKind::HostAligned;
        blk.host = core::detail::aligned_malloc(bytes);
#endif
        void* d = nullptr;
        if (cudaMalloc(&d, bytes) != cudaSuccess)
        {
#ifdef USE_PINNED_HOST
            if (blk.host)
                cudaFreeHost(blk.host);
#else
            core::detail::aligned_free(blk.host);
#endif
            throw std::bad_alloc{};
        }
        blk.device = d;
    }
#else
    // CPU-only build
    blk.kind = AllocKind::HostAligned;
    blk.host = core::detail::aligned_malloc(bytes);
    blk.device = nullptr;
#endif

    blk.bytes = bytes;
    {
        std::lock_guard<std::mutex> lk(mtx_);
        registry_.emplace(blk.host, blk);
    }
    return static_cast<T*>(blk.host);
}

} // namespace core
