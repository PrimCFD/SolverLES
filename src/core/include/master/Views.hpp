#pragma once
#include "mesh/Mesh.hpp"
#include <array>
#include <cstddef>
#include <cstdint>
#include <string>

/**
 * @file Views.hpp
 * @brief ABI-safe POD views shared with plugins and writers.
 *
 * @details
 * The views in this header are **trivial** (no STL containers except `std::string`) so they can
 * cross shared-library boundaries without ABI pitfalls. Writers and plugins read from these views
 * but never free the underlying storage, which is owned by the application/core.
 *
 * - :cpp:struct:`core::master::AnyFieldView` describes a 3-D field with **byte strides**.
 * - :cpp:struct:`core::master::MeshTileView` identifies a local sub-box and carries the stream.
 *
 * @rst
 * .. code-block:: cpp
 *
 *   AnyFieldView v{ "rho", rho_ptr, sizeof(double), {nx,ny,nz},
 *                   { 8, nx*8, nx*ny*8 } };
 *   MeshTileView tile{ Box3i{{0,0,0},{nx,ny,nz}}, rc.device_stream };
 * @endrst
 */

namespace core::master
{

struct Box3i
{
    std::array<int, 3> lo{0, 0, 0}, hi{0, 0, 0};
};

struct AnyFieldView
{
    std::string name;
    void* host_ptr{nullptr};
    std::size_t elem_size{0};
    std::array<int, 3> extents{0, 0, 0};
    std::array<std::ptrdiff_t, 3> strides{0, 0, 0};
};

template <class T> struct FieldView
{
    std::string name;
    T* data{nullptr};
    std::array<int, 3> extents{0, 0, 0};
    std::array<std::ptrdiff_t, 3> strides{0, 0, 0};
};

struct MeshTileView
{
    Box3i box; // interior window (no halos)
    void* stream{nullptr};
    const core::mesh::Mesh* mesh = nullptr;
};

inline AnyFieldView make_interior_copy(const AnyFieldView& v, int ng)
{
    AnyFieldView out = v; // copy metadata
    const int nx_i = v.extents[0] - 2 * ng;
    const int ny_i = v.extents[1] - 2 * ng;
    const int nz_i = v.extents[2] - 2 * ng;

    // Basic sanity (optional: turn into asserts)
    if (nx_i <= 0 || ny_i <= 0 || nz_i <= 0)
        return out; // leave unchanged

    out.extents = {nx_i, ny_i, nz_i};

    // byte-wise pointer bump: p + ng*sx + ng*sy + ng*sz
    auto* b = static_cast<unsigned char*>(v.host_ptr);
    b += std::size_t(ng) * std::size_t(v.strides[0]);
    b += std::size_t(ng) * std::size_t(v.strides[1]);
    b += std::size_t(ng) * std::size_t(v.strides[2]);
    out.host_ptr = b;

    return out;
}

} // namespace core::master
