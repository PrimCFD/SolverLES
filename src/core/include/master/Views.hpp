#pragma once
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
};

} // namespace core::master
