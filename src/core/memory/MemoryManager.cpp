#include "MemoryManager.hpp"
#include "AlignedAlloc.hpp"
#ifdef HAVE_CUDA
#include <cuda_runtime.h>
#endif
#include <algorithm>

namespace core
{

MemoryManager& MemoryManager::instance()
{
    static MemoryManager mgr;
    return mgr;
}

std::size_t MemoryManager::debug_count() const noexcept
{
    std::lock_guard<std::mutex> lock(mtx_);
    return registry_.size();
}

MemoryManager::~MemoryManager()
{
    for (auto& [ptr, blk] : registry_)
    {
        detail::aligned_free(ptr);
    }
}

void MemoryManager::release(void* p) noexcept
{
    if (!p)
        return;
    std::lock_guard<std::mutex> lock(mtx_);
    auto it = registry_.find(p);
    if (it != registry_.end())
    {
        detail::aligned_free(p);
        registry_.erase(it);
    }
}

void MemoryManager::to_device(void* p, std::size_t bytes)
{
#ifdef HAVE_CUDA
#ifndef USE_CUDA_UM
    cudaMemcpy(p, p, bytes, cudaMemcpyHostToDevice); // placeholder – implement real copy
#else
    (void) p;
    (void) bytes; // no‑op under UM
#endif
#else
    (void) p;
    (void) bytes;
#endif
}

void MemoryManager::to_host(void* p, std::size_t bytes)
{
#ifdef HAVE_CUDA
#ifndef USE_CUDA_UM
    cudaMemcpy(p, p, bytes, cudaMemcpyDeviceToHost); // placeholder
#else
    (void) p;
    (void) bytes;
#endif
#else
    (void) p;
    (void) bytes;
#endif
}

bool MemoryManager::on_device(const void* /*p*/) const noexcept
{
#ifdef USE_CUDA_UM
    return true; // UM pages migrate; true is safe default
#else
    return false; // host only by default
#endif
}

} // namespace core