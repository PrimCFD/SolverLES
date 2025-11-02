#include "master/FieldCatalog.hpp"
#include "memory/MemoryManager.hpp"
#include "mesh/Mesh.hpp"
#include <algorithm>
#include <stdexcept>

using namespace core::master;

void FieldCatalog::register_scalar(std::string name, void* host_ptr, std::size_t elem_size,
                                   std::array<int, 3> extents,
                                   std::array<std::ptrdiff_t, 3> strides, Stagger stag)
{
    if (idx_.count(name))
        throw std::runtime_error("Field already registered: " + name);
    AnyFieldView v;
    v.name = name;
    v.host_ptr = host_ptr;
    v.elem_size = elem_size;
    v.extents = extents;
    v.strides = strides;
    v.stagger = stag;
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

static inline std::array<std::ptrdiff_t, 3> row_major_strides_bytes(const std::array<int, 3>& e,
                                                                    std::size_t elem)
{
    return {(std::ptrdiff_t) elem, (std::ptrdiff_t) elem * e[0],
            (std::ptrdiff_t) elem * e[0] * e[1]};
}

// ---------------------- creation API -----------------
void FieldCatalog::create_center_scalar(std::string name, const core::mesh::Mesh& mesh)
{
    if (idx_.count(name))
        throw std::runtime_error("Field already registered: " + name);
    const auto e = mesh.totals_cell();
    const std::size_t N = (std::size_t) e[0] * e[1] * e[2];
    auto& mm = core::memory::MemoryManager::instance();
    double* p = mm.allocate<double>(N);
    std::fill_n(p, N, 0.0);
    owned_[name] = OwnedBlock{p, N * sizeof(double)};
    auto sb = row_major_strides_bytes(e, sizeof(double));
    register_scalar(std::move(name), p, sizeof(double), e, sb, Stagger::Cell);
}

void FieldCatalog::create_face_scalar(std::string name, int axis, const core::mesh::Mesh& mesh)
{
    if (axis < 0 || axis > 2)
        throw std::runtime_error("create_face_scalar: bad axis");
    if (idx_.count(name))
        throw std::runtime_error("Field already registered: " + name);
    const auto e = mesh.totals_face_axis(axis);
    const std::size_t N = (std::size_t) e[0] * e[1] * e[2];
    auto& mm = core::memory::MemoryManager::instance();
    double* p = mm.allocate<double>(N);
    std::fill_n(p, N, 0.0);
    owned_[name] = OwnedBlock{p, N * sizeof(double)};
    auto sb = row_major_strides_bytes(e, sizeof(double));
    const Stagger stag = (axis == 0)   ? Stagger::IFace
                         : (axis == 1) ? Stagger::JFace
                                       : Stagger::KFace;
    register_scalar(std::move(name), p, sizeof(double), e, sb, stag);
}
