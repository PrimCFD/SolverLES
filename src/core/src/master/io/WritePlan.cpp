#include "master/io/WritePlan.hpp"
#include "master/Views.hpp" // AnyFieldView
#include <algorithm>

namespace core::master::io
{

static inline bool is_contiguous(const AnyFieldView& v, std::size_t s)
{
    const auto nx = v.extents[0], ny = v.extents[1], nz = v.extents[2];
    const auto sx = v.strides[0], sy = v.strides[1], sz = v.strides[2];
    return sx == (std::ptrdiff_t) s && sy == (std::ptrdiff_t)(s * nx) &&
           sz == (std::ptrdiff_t)(s * nx * ny);
}

WritePlan build_write_plan(std::span<const AnyFieldView> selected, std::size_t elem_size_override)
{
    WritePlan plan;
    plan.fields.reserve(selected.size());
    for (const AnyFieldView& v : selected)
    {
        FieldPlan fp{};
        fp.shape.name = v.name;
        fp.shape.elem_size = elem_size_override ? elem_size_override : v.elem_size;
        fp.shape.nx = v.extents[0];
        fp.shape.ny = v.extents[1];
        fp.shape.nz = v.extents[2];
        fp.shape.sx = v.strides[0];
        fp.shape.sy = v.strides[1];
        fp.shape.sz = v.strides[2];
        fp.shape.contiguous = is_contiguous(v, fp.shape.elem_size);
        fp.bytes = std::size_t(fp.shape.nx) * std::size_t(fp.shape.ny) * std::size_t(fp.shape.nz) *
                   fp.shape.elem_size;

        plan.index_of.emplace(fp.shape.name, plan.fields.size());
        plan.fields.push_back(fp);
    }
    return plan;
}

} // namespace core::master::io