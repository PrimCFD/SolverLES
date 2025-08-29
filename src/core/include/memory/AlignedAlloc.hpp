#pragma once
#include <cstddef>
#include <cstdlib>
#include <new>

/**
 * @file AlignedAlloc.hpp
 * @ingroup memory
 * @brief Low‑level aligned allocation helpers for host memory.
 *
 * Provides \c aligned_malloc(bytes, alignment)  and \c aligned_free(ptr) used by the
 * memory subsystem to satisfy cache‑line/SIMD alignment (≥64 B by default).
 * When \c USE_PINNED_HOST is defined and CUDA is available, host buffers may be
 * allocated with \c cudaHostAlloc instead (see MemoryManager for policy).
 */

#ifdef HAVE_CUDA
#include <cuda_runtime.h>
#endif

namespace core::memory
{

#if defined(USE_CUDA_UM)
// Managed allocations are at least 256B aligned on most GPUs; 128 is safe+wide for host SIMD.
inline constexpr std::size_t HW_ALIGN = 128;
#else
inline constexpr std::size_t HW_ALIGN = 64;
#endif

inline void* aligned_malloc(std::size_t bytes, std::size_t alignment = HW_ALIGN)
{
#if defined(USE_CUDA_UM) && defined(HAVE_CUDA)
    // Unified Memory path: one pointer valid on host and device
    void* p = nullptr;
    if (cudaMallocManaged(&p, bytes, cudaMemAttachGlobal) != cudaSuccess)
    {
        throw std::bad_alloc{};
    }
    return p;

#else
// Host-only path (or UM requested but CUDA not available)
#if defined(_MSC_VER)
    void* p = _aligned_malloc(bytes, alignment);
    if (!p)
        throw std::bad_alloc{};
    return p;
#else
    // std::aligned_alloc requires size multiple of alignment
    if (alignment == 0 || (alignment & (alignment - 1)) != 0)
        alignment = HW_ALIGN;
    std::size_t padded = ((bytes + alignment - 1) / alignment) * alignment;
    void* p = std::aligned_alloc(alignment, padded);
    if (!p)
        throw std::bad_alloc{};
    return p;
#endif
#endif
}

inline void aligned_free(void* p) noexcept
{
#if defined(USE_CUDA_UM) && defined(HAVE_CUDA)
    if (p)
        cudaFree(p);
#else
#if defined(_MSC_VER)
    _aligned_free(p);
#else
    std::free(p);
#endif
#endif
}

} // namespace core::memory
