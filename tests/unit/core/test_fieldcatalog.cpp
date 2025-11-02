#include "master/FieldCatalog.hpp"
#include <catch2/catch_test_macros.hpp>
#include <vector>

using namespace core::master;

TEST_CASE("FieldCatalog basic registration + selection", "[fieldcatalog]")
{
    FieldCatalog fc;
    std::vector<double> a(8, 0.0);

    fc.register_scalar("a", a.data(), sizeof(double), {2, 2, 2}, {1, 2, 4},
                       core::master::Stagger::Cell);
    REQUIRE(fc.contains("a"));

    auto v = fc.view("a");
    REQUIRE(v.elem_size == sizeof(double));
    fc.select_for_output("a");
    REQUIRE_FALSE(fc.selected_for_output().empty());
}
