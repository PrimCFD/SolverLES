#include "master/FieldCatalog.hpp"
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <vector>

using namespace core::master;

TEST_CASE("FieldCatalog view properties and bad lookups", "[fieldcatalog]")
{
    FieldCatalog fc;
    std::vector<float> a(16, 0.f);

    fc.register_scalar("T", a.data(), sizeof(float), std::array<int, 3>{2, 2, 4},
                       std::array<std::ptrdiff_t, 3>{1, 2, 8});
    REQUIRE(fc.contains("T"));

    auto v = fc.view("T");
    REQUIRE(v.elem_size == sizeof(float));
    REQUIRE(v.extents == std::array<int, 3>{2, 2, 4});
    REQUIRE(v.strides == std::array<std::ptrdiff_t, 3>{1, 2, 8});

    // Selection appends (no dedup in FieldCatalog)
    const auto before = fc.selected_for_output().size();
    fc.select_for_output("T"); // duplicates allowed
    REQUIRE(fc.selected_for_output().size() == before + 1);

    // Bad name fails
    REQUIRE_FALSE(fc.contains("bad"));
    REQUIRE_THROWS(fc.view("bad"));
}
