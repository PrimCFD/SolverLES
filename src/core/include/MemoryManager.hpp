#pragma once
#include "AlignedAlloc.hpp"
#include <cstddef>
#include <cstdint>
#include <mutex>
#include <unordered_map>

// Optional toggles expected from the build system:
// - HAVE_CUDA
// - USE_CUDA_UM
// - USE_PINNED_HOST  (optional: prefer pinned host mem when not using UM & CUDA is available)

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
