#pragma once
#include "Layout.hpp"
#include <array>
#include <span>

namespace core
{

template <class T> class Field
{
  public:
    Field() = default;

    Field(T* raw, std::array<int, 3> ext_with_ghosts, int ng)
        : data_{raw}, ext_{ext_with_ghosts}, ng_{ng}, idx_{layout::make_indexer(ext_with_ghosts)}
    {
    }

    constexpr T& operator()(int i, int j, int k) noexcept
    {
        return data_[idx_(i + ng_, j + ng_, k + ng_)];
    }
    constexpr const T& operator()(int i, int j, int k) const noexcept
    {
        return data_[idx_(i + ng_, j + ng_, k + ng_)];
    }

    [[nodiscard]] void* c_ptr() noexcept { return static_cast<void*>(data_); }

    std::span<T> span() noexcept
    {
        return {data_, static_cast<std::size_t>(ext_[0] * ext_[1] * ext_[2])};
    }

  private:
    T* data_ = nullptr; // non‑owning
    std::array<int, 3> ext_{};
    int ng_ = 0;
    layout::Indexer3D idx_{};
};

} // namespace core