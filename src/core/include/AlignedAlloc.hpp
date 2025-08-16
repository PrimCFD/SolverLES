#pragma once
#include <cstddef>
#include <cstdlib>
#include <stdexcept>
#ifdef HAVE_CUDA
#include <cuda_runtime.h>
#endif

namespace core::detail
{

#ifdef USE_CUDA_UM
constexpr std::size_t HW_ALIGN = 128; // covers CPU 64 B as well
#else
constexpr std::size_t HW_ALIGN = 64;
#endif

inline void* aligned_malloc(std::size_t bytes)
{
#ifdef USE_CUDA_UM
    void* p = nullptr;
#ifdef HAVE_CUDA
    if (cudaMallocManaged(&p, bytes) != cudaSuccess)
        p = nullptr;
#else
    (void) bytes;
    p = nullptr; // UM requested but CUDA absent
#endif
    if (!p)
        throw std::bad_alloc{};
    return p;
#else
    void* p = std::aligned_alloc(HW_ALIGN, ((bytes + HW_ALIGN - 1) & ~(HW_ALIGN - 1)));
    if (!p)
        throw std::bad_alloc{};
    return p;
#endif
}

inline void aligned_free(void* p) noexcept
{
#ifdef USE_CUDA_UM
#ifdef HAVE_CUDA
    cudaFree(p);
#else
    std::free(p);
#endif
#else
    std::free(p);
#endif
}

} // namespace core::detail