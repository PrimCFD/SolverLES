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
 * The catalog maps **logical names** to :cpp:struct:`AnyFieldView` and maintains a **selection**
 * used by writers. Ownership of storage remains with the application/core.
 *
 * @rst
 * .. code-block:: cpp
 *
 *   FieldCatalog fc;
 *   fc.register_scalar("rho", rho_ptr, sizeof(double), {nx,ny,nz},
 *                      { 8, nx*8, nx*ny*8 });
 *   fc.select_for_output("rho");
 *
 *   for (const auto& v : fc.selected_for_output()) { // pass to writer
 *       // writer.write_view(v);
 *   }
 *
 * @endrst
 */

namespace core::master
{

class FieldCatalog
{
  public:
    void register_scalar(std::string name, void* host_ptr, std::size_t elem_size,
                         std::array<int, 3> extents, std::array<std::ptrdiff_t, 3> strides,
                         Stagger stag);

    bool contains(std::string_view name) const noexcept;
    AnyFieldView view(std::string_view name) const;

    void select_for_output(std::string_view name);

    // These allocate storage with MemoryManager and register the field with correct extents and
    // staggering.
    // * Cell-centered scalar (double)
    void create_center_scalar(std::string name, const core::mesh::Mesh& mesh);
    // * Face-centered scalar (double) along axis a=0(x),1(y),2(z)  => interior size is N+1 on 'a'
    void create_face_scalar(std::string name, int axis, const core::mesh::Mesh& mesh);

    std::span<const AnyFieldView> all_views() const noexcept { return views_; }
    std::span<const AnyFieldView> selected_for_output() const noexcept { return selected_views_; }

  private:
    std::vector<AnyFieldView> views_;
    std::unordered_map<std::string, std::size_t> idx_;
    std::vector<AnyFieldView> selected_views_;
    struct OwnedBlock
    {
        void* ptr{};
        std::size_t bytes{};
    };
    std::unordered_map<std::string, OwnedBlock> owned_; // core-owned allocations
};

} // namespace core::master
