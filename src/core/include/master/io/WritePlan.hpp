#pragma once
#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <unordered_map>
#include <vector>

/**
 * @file WritePlan.hpp
 * @brief Per-field layout plan synthesized from catalog views.
 *
 * @details
 * A :cpp:struct:`WritePlan` captures logical extents, **byte** strides, element size after
 * precision policy, contiguity flags, and total byte counts per field. Writers use the plan to
 * decide between **direct** writes and **staging**, and to size the :cpp:class:`StagingPool`.
 *
 * @rst
 * .. code-block:: cpp
 *
 *   WritePlan P = build_write_plan(fc.selected_for_output(), 4); // override precision
 *   for (const FieldPlan& f : P.fields) { // allocate staging, choose path
 *     // ...
 *   }
 * @endrst
 */

namespace core::master
{
struct AnyFieldView;
struct MeshTileView;
} // namespace core::master

namespace core::master::io
{

struct FieldShape
{
    std::string name;
    std::size_t elem_size = 0;             // bytes per element *after* precision policy
    int nx = 0, ny = 0, nz = 0;            // logical extents to write
    std::ptrdiff_t sx = 0, sy = 0, sz = 0; // byte strides in source memory
    bool contiguous = false;
};

struct FieldPlan
{
    FieldShape shape;
    std::size_t bytes = 0; // nx*ny*nz*elem_size
};

struct WritePlan
{
    std::vector<FieldPlan> fields;
    std::unordered_map<std::string, std::size_t> index_of;
    const FieldPlan* find(std::string_view name) const
    {
        auto it = index_of.find(std::string(name));
        return it == index_of.end() ? nullptr : &fields[it->second];
    }
};

// Build directly from the catalog's selected views
WritePlan build_write_plan(std::span<const AnyFieldView> selected,
                           std::size_t elem_size_override /*0=keep*/);

} // namespace core::master::io