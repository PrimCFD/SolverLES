#pragma once
#include <array>
#include <cstddef>
#include <cstdint>
#include <string>

/**
 * @file Views.hpp
 * @brief Plain Old Data (POD) views used across plugin and writer ABIs.
 *
 * @details
 * These structs are purposely trivial (no STL/Eigen) so they can cross
 * shared-library boundaries without ABI pitfalls. Plugins may retrieve
 * device mirrors of host pointers via the MemoryManager using the host
 * address stored in these views.
 *
 * @rst
 * .. code-block:: cpp
 *
 *    // Access rules:
 *    // - "host_ptr" belongs to the core; plugins must not free it.
 *    // - "stream" comes from RunContext; use it for kernels and copies.
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
