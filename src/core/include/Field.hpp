#pragma once
#include "Layout.hpp"
#include <array>
#include <cstddef>
#include <span>

namespace core
{

template <class T> class Field
{
    T* data_ = nullptr;        // host (or UM) pointer
    std::array<int, 3> ext_{}; // totals including ghosts
    int ng_ = 0;
    layout::Indexer3D idx_{};

  public:
    Field() = default;
    Field(T* p, std::array<int, 3> e, int ng)
        : data_(p), ext_(e), ng_(ng), idx_{layout::Indexer3D{e[0], e[1], e[2]}}
    {
    }

    inline T& operator()(int i, int j, int k) noexcept
    {
        return data_[idx_(i + ng_, j + ng_, k + ng_)];
    }
    inline const T& operator()(int i, int j, int k) const noexcept
    {
        return data_[idx_(i + ng_, j + ng_, k + ng_)];
    }

    std::span<T> span() noexcept
    {
        return {data_, static_cast<std::size_t>(ext_[0]) * static_cast<std::size_t>(ext_[1]) *
                           static_cast<std::size_t>(ext_[2])};
    }

    // Raw typed pointer accessors
    T* raw() noexcept { return data_; }
    const T* raw() const noexcept { return data_; }

    const std::array<int, 3>& extents() const noexcept { return ext_; }
    int ng() const noexcept { return ng_; }
};

} // namespace core
