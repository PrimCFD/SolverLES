#pragma once
#include "Layout.hpp"
#include <array>
#include <cstddef>
#include <span>

/**
 * @file Field.hpp
 * @ingroup memory
 * @brief Non‑owning, typed view into SoA data with ghost padding.
 *
 * \c Field\<T\> wraps a raw pointer and the ghost‑inclusive extents. It does not
 * allocate or free memory; ownership remains with MemoryManager. Designed for
 * zero‑overhead access in inner loops and ISO_C_BINDING interop.
 *
 * @tparam T arithmetic element type (e.g., double, float)
 *
 *@rst
 *.. code-block:: cpp
 *
 *   Field<double> rho(rho_ptr, {nx_tot, ny_tot, nz_tot}, ng);
 *   double val = rho(i,j,k); // adds ghost offset, uses Indexer3D
 *   void* cptr = rho.c_ptr(); // pass to Fortran bind(C) kernel
 * @endrst
 *
 * @note Bounds are not checked in release builds. Prefer unit tests to validate extents.
 */

namespace core::mesh
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

} // namespace core::mesh
