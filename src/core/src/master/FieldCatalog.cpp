#include "master/FieldCatalog.hpp"
#include <stdexcept>

using namespace core::master;

void FieldCatalog::register_scalar(std::string name, void* host_ptr, std::size_t elem_size,
                                   std::array<int, 3> extents,
                                   std::array<std::ptrdiff_t, 3> strides)
{
    if (idx_.count(name))
        throw std::runtime_error("Field already registered: " + name);
    AnyFieldView v;
    v.name = name;
    v.host_ptr = host_ptr;
    v.elem_size = elem_size;
    v.extents = extents;
    v.strides = strides;
    idx_[v.name] = views_.size();
    views_.push_back(v);
}

bool FieldCatalog::contains(std::string_view name) const noexcept
{
    return idx_.find(std::string(name)) != idx_.end();
}

AnyFieldView FieldCatalog::view(std::string_view name) const
{
    auto it = idx_.find(std::string(name));
    if (it == idx_.end())
        throw std::runtime_error("Unknown field: " + std::string(name));
    return views_[it->second];
}

void FieldCatalog::select_for_output(std::string_view name)
{
    auto v = view(name);
    selected_views_.push_back(v);
}
