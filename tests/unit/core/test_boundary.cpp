#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <numeric>
#include <vector>

#include "mesh/Boundary.hpp"

using core::mesh::Array3DView;
using core::mesh::Axis;
using core::mesh::BCOp;
using core::mesh::MirrorMask;

using Catch::Approx;

namespace
{
template <class T>
Array3DView<T> make_view(std::vector<T>& buf, int nx, int ny, int nz, int hx, int hy, int hz)
{
    const std::size_t total = static_cast<std::size_t>(nx + 2 * hx) *
                              static_cast<std::size_t>(ny + 2 * hy) *
                              static_cast<std::size_t>(nz + 2 * hz);
    REQUIRE(buf.size() >= total);
    return Array3DView<T>{buf.data(), nx, ny, nz, hx, hy, hz};
}

inline double pattern(int i, int j, int k)
{
    return i + 10.0 * j + 100.0 * k;
}
} // namespace

TEST_CASE("Dirichlet face fill", "[boundary][dirichlet]")
{
    const int nx = 4, ny = 3, nz = 2, hx = 2, hy = 1, hz = 1;
    std::vector<double> a((nx + 2 * hx) * (ny + 2 * hy) * (nz + 2 * hz), -99.0);
    auto A = make_view(a, nx, ny, nz, hx, hy, hz);

    // initialize interior
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
                A.at(i, j, k) = pattern(i, j, k);

    const double c = 7.25;
    core::mesh::apply_scalar_bc(A, Axis::I, -1, BCOp::Dirichlet, c); // minus-I face

    // Check ghost layers on minus-I side (i = -1, -2)
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
        {
            REQUIRE(A.at(-1, j, k) == Approx(c));
            REQUIRE(A.at(-2, j, k) == Approx(c));
        }
}

TEST_CASE("Neumann zero (copy interior)", "[boundary][neumann0]")
{
    const int nx = 5, ny = 4, nz = 3, hx = 1, hy = 2, hz = 1;
    std::vector<double> a((nx + 2 * hx) * (ny + 2 * hy) * (nz + 2 * hz), -99.0);
    auto A = make_view(a, nx, ny, nz, hx, hy, hz);

    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
                A.at(i, j, k) = pattern(i, j, k);

    core::mesh::apply_scalar_bc(A, Axis::J, +1, BCOp::NeumannZero); // plus-J face

    for (int k = 0; k < nz; ++k)
        for (int i = 0; i < nx; ++i)
        {
            REQUIRE(A.at(i, ny, k) == Approx(A.at(i, ny - 1, k)));
            REQUIRE(A.at(i, ny + 1, k) == Approx(A.at(i, ny - 1, k)));
        }
}

TEST_CASE("Extrapolate1 (linear)", "[boundary][extrapolate1]")
{
    const int nx = 3, ny = 3, nz = 5, hx = 1, hy = 1, hz = 2;
    std::vector<double> a((nx + 2 * hx) * (ny + 2 * hy) * (nz + 2 * hz), -99.0);
    auto A = make_view(a, nx, ny, nz, hx, hy, hz);

    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
                A.at(i, j, k) = pattern(i, j, k);

    // minus-K face extrapolation
    core::mesh::apply_scalar_bc(A, Axis::K, -1, BCOp::Extrapolate1);

    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i)
        {
            const double f0 = A.at(i, j, 0);
            const double f1 = A.at(i, j, 1);
            const double d = f0 - f1;
            REQUIRE(A.at(i, j, -1) == Approx(f0 + 1 * d));
            REQUIRE(A.at(i, j, -2) == Approx(f0 + 2 * d));
        }
}

TEST_CASE("Vector mirror with sign mask", "[boundary][mirror-vector]")
{
    const int nx = 5, ny = 2, nz = 2, hx = 2, hy = 1, hz = 1;
    std::vector<double> ux((nx + 2 * hx) * (ny + 2 * hy) * (nz + 2 * hz), -9.0);
    std::vector<double> uy(ux.size(), -9.0);
    std::vector<double> uz(ux.size(), -9.0);

    auto Ux = make_view(ux, nx, ny, nz, hx, hy, hz);
    auto Uy = make_view(uy, nx, ny, nz, hx, hy, hz);
    auto Uz = make_view(uz, nx, ny, nz, hx, hy, hz);

    // interior: simple pattern per component
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
            {
                Ux.at(i, j, k) = 1.0 + pattern(i, j, k); // x-comp
                Uy.at(i, j, k) = 2.0 + pattern(i, j, k); // y-comp
                Uz.at(i, j, k) = 3.0 + pattern(i, j, k); // z-comp
            }

    // Mirror across minus-I face with no-slip mask: (-1,+1,+1)
    core::mesh::apply_mirror_vector(Ux, Uy, Uz, Axis::I, -1, MirrorMask{{-1, 1, 1}});

    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
        {
            // layer l=1
            REQUIRE(Ux.at(-1, j, k) == Approx(-Ux.at(0, j, k)));
            REQUIRE(Uy.at(-1, j, k) == Approx(+Uy.at(0, j, k)));
            REQUIRE(Uz.at(-1, j, k) == Approx(+Uz.at(0, j, k)));
            // layer l=2 (mirrors interior l-1 -> 1)
            REQUIRE(Ux.at(-2, j, k) == Approx(-Ux.at(1, j, k)));
            REQUIRE(Uy.at(-2, j, k) == Approx(+Uy.at(1, j, k)));
            REQUIRE(Uz.at(-2, j, k) == Approx(+Uz.at(1, j, k)));
        }
}
