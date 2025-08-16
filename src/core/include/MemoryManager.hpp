#pragma once
#include "AlignedAlloc.hpp"
#include <cstddef>
#include <mutex>
#include <unordered_map>
#ifdef HAVE_MPI
#include <mpi.h>
#endif

namespace core
{

class MemoryManager
{
  public:
    // Meyers singleton – thread‑safe in C++11 and above
    static MemoryManager& instance();

    // For tests: number of tracked allocations
    std::size_t debug_count() const noexcept;

    template <class T> T* allocate(std::size_t n);

    void release(void* p) noexcept;

    // Host/device helpers become no‑ops if Unified Memory is enabled or
    // CUDA not present.
    void to_device(void* p, std::size_t bytes);
    void to_host(void* p, std::size_t bytes);
    bool on_device(const void* p) const noexcept;

  private:
    MemoryManager() = default;
    ~MemoryManager();

    struct Block
    {
        void* ptr;
        std::size_t bytes;
        bool device;
    };

    std::unordered_map<void*, Block> registry_;
    mutable std::mutex mtx_;
};

template <class T> inline T* MemoryManager::allocate(std::size_t n)
{
    std::lock_guard<std::mutex> lock(mtx_);
    const std::size_t bytes = n * sizeof(T);
    void* raw = core::detail::aligned_malloc(bytes);
    registry_[raw] = {raw, bytes, /*device*/ false};
    return static_cast<T*>(raw);
}

} // namespace core