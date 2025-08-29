#include "memory/MemoryManager.hpp"
#include <algorithm>
#include <new>

namespace core::memory
{

MemoryManager& MemoryManager::instance()
{
    static MemoryManager mm;
    return mm;
}

bool MemoryManager::using_unified_memory() const noexcept
{
#ifdef HAVE_CUDA
#ifdef USE_CUDA_UM
    return true;
#else
    return false;
#endif
#else
    return false;
#endif
}

std::size_t MemoryManager::debug_count() const noexcept
{
    std::lock_guard<std::mutex> lk(mtx_);
    return registry_.size();
}

Block* MemoryManager::find_block_unlocked(const void* host_ptr) noexcept
{
    auto it = registry_.find(host_ptr);
    if (it == registry_.end())
        return nullptr;
    return &it->second;
}

const Block* MemoryManager::find_block_unlocked(const void* host_ptr) const noexcept
{
    auto it = registry_.find(host_ptr);
    if (it == registry_.end())
        return nullptr;
    return &it->second;
}

void* MemoryManager::device_ptr(const void* host_ptr) const noexcept
{
    std::lock_guard<std::mutex> lk(mtx_);
    auto b = find_block_unlocked(host_ptr);
    return b ? b->device : nullptr;
}

void* MemoryManager::host_ptr(const void* maybe_device) const noexcept
{
    std::lock_guard<std::mutex> lk(mtx_);
    // First try key match (host pointer)
    if (auto it = registry_.find(maybe_device); it != registry_.end())
        return it->second.host;
    // Slow path: search mirrors
    for (auto& kv : registry_)
    {
        if (kv.second.device == maybe_device)
            return kv.second.host;
    }
    return nullptr;
}

void MemoryManager::to_device(const void* host_ptr, std::size_t bytes, void* stream)
{
#ifdef HAVE_CUDA
    if (using_unified_memory())
    {
        // Optional prefetch to current device for performance
        int dev = 0;
        cudaGetDevice(&dev);
        cudaMemPrefetchAsync(const_cast<void*>(host_ptr), bytes, dev,
                             static_cast<cudaStream_t>(stream));
        return;
    }

    std::lock_guard<std::mutex> lk(mtx_);
    auto b = find_block_unlocked(host_ptr);
    if (!b || !b->device)
        return;
    cudaMemcpyAsync(b->device, b->host, bytes, cudaMemcpyHostToDevice,
                    static_cast<cudaStream_t>(stream));
#else
    (void) host_ptr;
    (void) bytes;
    (void) stream;
#endif
}

void MemoryManager::to_host(const void* host_ptr, std::size_t bytes, void* stream)
{
#ifdef HAVE_CUDA
    if (using_unified_memory())
    {
        cudaMemPrefetchAsync(const_cast<void*>(host_ptr), bytes, cudaCpuDeviceId,
                             static_cast<cudaStream_t>(stream));
        return;
    }

    std::lock_guard<std::mutex> lk(mtx_);
    auto b = find_block_unlocked(host_ptr);
    if (!b || !b->device)
        return;
    cudaMemcpyAsync(b->host, b->device, bytes, cudaMemcpyDeviceToHost,
                    static_cast<cudaStream_t>(stream));
#else
    (void) host_ptr;
    (void) bytes;
    (void) stream;
#endif
}

void MemoryManager::release(void* host_ptr) noexcept
{
    if (!host_ptr)
        return;
    std::lock_guard<std::mutex> lk(mtx_);
    auto it = registry_.find(host_ptr);
    if (it == registry_.end())
        return;

    Block blk = it->second;
#ifdef HAVE_CUDA
    if (using_unified_memory())
    {
        if (blk.host)
            cudaFree(blk.host);
    }
    else
    {
        if (blk.device)
            cudaFree(blk.device);
#ifdef USE_PINNED_HOST
        if (blk.host)
            cudaFreeHost(blk.host);
#else
        core::memory::aligned_free(blk.host);
#endif
    }
#else
    core::memory::aligned_free(blk.host);
#endif
    registry_.erase(it);
}

MemoryManager::~MemoryManager()
{
    std::lock_guard<std::mutex> lk(mtx_);
    for (auto& kv : registry_)
    {
        Block& blk = kv.second;
#ifdef HAVE_CUDA
        if (using_unified_memory())
        {
            if (blk.host)
                cudaFree(blk.host);
        }
        else
        {
            if (blk.device)
                cudaFree(blk.device);
#ifdef USE_PINNED_HOST
            if (blk.host)
                cudaFreeHost(blk.host);
#else
            core::memory::aligned_free(blk.host);
#endif
        }
#else
        core::memory::aligned_free(blk.host);
#endif
    }
    registry_.clear();
}

} // namespace core::memory
