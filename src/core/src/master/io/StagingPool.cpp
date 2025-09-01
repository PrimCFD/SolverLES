#include "master/io/StagingPool.hpp"
#include <cstdlib>
#include <cstring>
#include <new>

namespace core::master::io
{

static void* page_aligned_alloc(std::size_t bytes)
{
    void* p = nullptr;
#if defined(_MSC_VER)
    p = _aligned_malloc(bytes, 4096);
    if (!p)
        throw std::bad_alloc{};
#else
    if (posix_memalign(&p, 4096, bytes))
        throw std::bad_alloc{};
#endif
    return p;
}

void StagingPool::reserve(std::size_t count, std::size_t bytes)
{
    std::lock_guard<std::mutex> lock(mtx_);
    if (buffers_.size() < count)
        buffers_.resize(count);
    for (auto& b : buffers_)
    {
        if (b.cap >= bytes)
            continue;
        // allocate new and swap
        void* p = page_aligned_alloc(bytes);
        b.ptr.reset(p);
        b.cap = bytes;
    }
}

/// \cond DOXYGEN_EXCLUDE
void StagingPool::Block::Deleter::operator()(void* p) const noexcept
{
#if defined(_MSC_VER)
    _aligned_free(p);
#else
    free(p);
#endif
}
/// \encond

StagingPool::~StagingPool() = default;

} // namespace core::master::io