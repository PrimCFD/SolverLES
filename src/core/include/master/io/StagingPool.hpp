#pragma once
#include <cstddef>
#include <memory>
#include <mutex>
#include <vector>

/**
 * @file StagingPool.hpp
 * @brief Reusable page-aligned host buffers sized to the largest field.
 *
 * @details
 * Writers request N buffers, each at least `bytes` capacity. Implemented with page-aligned
 * allocations; integrate with `MemoryManager` for pinned allocations if needed.
 *
 * Thread-safe for distinct buffer indices.
 */

namespace core::master::io
{

// A reusable pool of host buffers sized to the largest field we need to stage.
// Uses page-aligned allocations; integrate with your MemoryManager for pinned allocations.
class StagingPool
{
  public:
    explicit StagingPool() = default;
    ~StagingPool();

    // Ensure we have N buffers, each at least `bytes` capacity.
    void reserve(std::size_t count, std::size_t bytes);

    // Borrow buffer i (0<=i<count). Thread-safe for distinct i.
    void* data(std::size_t i) noexcept { return buffers_[i].ptr.get(); }
    std::size_t capacity(std::size_t i) const noexcept { return buffers_[i].cap; }

  private:
    struct Block
    {

        struct Deleter
        {
            void operator()(void* p) const noexcept;
        };

        std::unique_ptr<void, Deleter> ptr{nullptr};
        std::size_t cap = 0;
        
    };

    std::vector<Block> buffers_;
    std::mutex mtx_;
};

} // namespace core::master::io