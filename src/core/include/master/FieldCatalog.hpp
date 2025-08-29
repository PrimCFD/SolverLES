#pragma once
#include "master/Views.hpp"
#include <span>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

/**
 * @file FieldCatalog.hpp
 * @brief Registry of fields and their POD views for plugins and I/O.
 *
 * @details
 * The catalog maps logical names to views and maintains a selection set
 * for output. Ownership of storage remains in the application/core.
 *
 * @rst
 * .. code-block:: cpp
 *
 *    FieldCatalog fc;
 *    fc.register_scalar("rho", rho.data(), sizeof(double), {nx,ny,nz}, {1,sx,sy});
 *    fc.select_for_output("rho");
 *
 * @endrst
 */

namespace core::master
{

class FieldCatalog
{
  public:
    void register_scalar(std::string name, void* host_ptr, std::size_t elem_size,
                         std::array<int, 3> extents, std::array<std::ptrdiff_t, 3> strides);

    bool contains(std::string_view name) const noexcept;
    AnyFieldView view(std::string_view name) const;

    void select_for_output(std::string_view name);
    std::span<const AnyFieldView> all_views() const noexcept { return views_; }
    std::span<const AnyFieldView> selected_for_output() const noexcept { return selected_views_; }

  private:
    std::vector<AnyFieldView> views_;
    std::unordered_map<std::string, std::size_t> idx_;
    std::vector<AnyFieldView> selected_views_;
};

} // namespace core::master
