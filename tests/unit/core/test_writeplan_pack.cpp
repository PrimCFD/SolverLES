#include "master/Views.hpp" // AnyFieldView
#include "master/io/WritePlan.hpp"
#include "memory/MemoryManager.hpp"
#include <catch2/catch_all.hpp>

using core::master::AnyFieldView;
using core::master::io::build_write_plan;
using core::master::io::WritePlan;
using core::memory::MemoryManager;

static AnyFieldView make_view(const char* name, void* base, size_t elem, int nx, int ny, int nz,
                              ptrdiff_t sx, ptrdiff_t sy, ptrdiff_t sz)
{
    AnyFieldView v{};
    v.name = name;
    v.host_ptr = base;
    v.elem_size = elem;
    v.extents[0] = nx;
    v.extents[1] = ny;
    v.extents[2] = nz;
    v.strides[0] = sx;
    v.strides[1] = sy;
    v.strides[2] = sz;
    return v;
}

TEST_CASE("WritePlan contiguous mapping", "[io][plan]")
{
    auto& mm = MemoryManager::instance();
    const int nx = 8, ny = 4, nz = 3;
    const size_t elem = sizeof(double);
    const size_t N = size_t(nx) * ny * nz;
    double* buf = mm.allocate<double>(N);

    AnyFieldView v = make_view("rho", buf, elem, nx, ny, nz, elem, elem * nx, elem * nx * ny);
    std::vector<AnyFieldView> sel{v};
    auto plan = build_write_plan(std::span<const AnyFieldView>(sel.data(), sel.size()), 0);

    REQUIRE(plan.fields.size() == 1);
    REQUIRE(plan.fields[0].shape.contiguous);
    REQUIRE(plan.fields[0].bytes == N * elem);

    mm.release(buf);
}

TEST_CASE("WritePlan strided mapping", "[io][plan]")
{
    auto& mm = MemoryManager::instance();
    const int nx = 7, ny = 5, nz = 2;
    const size_t elem = sizeof(float);
    // add padding in Y and Z
    const ptrdiff_t sx = elem, sy = elem * (nx + 3), sz = sy * (ny + 2);
    const size_t span_bytes = size_t((nx - 1) * sx + (ny - 1) * sy + (nz - 1) * sz) + elem;

    unsigned char* raw =
        reinterpret_cast<unsigned char*>(mm.allocate<std::byte>((span_bytes + elem - 1) / elem));
    AnyFieldView v = make_view("T", raw, elem, nx, ny, nz, sx, sy, sz);

    std::vector<AnyFieldView> sel{v};
    auto plan = build_write_plan(std::span<const AnyFieldView>(sel.data(), sel.size()), 0);

    REQUIRE_FALSE(plan.fields[0].shape.contiguous);
    REQUIRE(plan.fields[0].bytes == size_t(nx) * ny * nz * elem);

    mm.release(reinterpret_cast<std::byte*>(raw));
}
