#pragma once
#include "Field.hpp" // <-- new: delegate all indexing to Field<T>
#include <array>
#include <cstddef>
#include <cstdint>
#include <type_traits>
#include <vector>

/**
 * @file Boundary.hpp
 * @ingroup memory
 * @brief Face-wise ghost–cell boundary operators for structured 3-D fields.
 *
 * Provides scalar and vector boundary condition (BC) kernels that fill ghost layers
 * on the 6 axis-aligned faces of a cartesian block. Indexing is delegated to
 * `Field<T>` (via the lightweight `Array3DView<T>` adapter) so there is a single
 * source of truth for memory layout (`layout::Indexer3D` is only accessed by `Field`).
 *
 * @section boundary-assumptions Numerical assumptions & conventions
 *
 * - **Ghost layers:** uniform width :code:`ng >= 1` on all axes (I,J,K). All functions
 *   fill **every** ghost layer :code:`l = 1..ng`. Anisotropic ghosts are not required
 *   at the moment; can be introduced by extending `Field` (BC code stays unchanged).
 * - **Grid spacing:** unit spacing is assumed along the normal when computing
 *   one-sided differences; i.e. :math:`\Delta n = 1`. If using physical
 *   spacing :math:`\Delta n \neq 1`, need to be sclaed outside these helpers.
 * - **Faces vs. edges/corners:** functions write **faces only**. If applying BCs on
 *   multiple axes, edges/corners are the result of the last face write ("last write wins").
 *   If need for a special edge/corner policy for wide stencils, do an explicit pass.
 * - **Domains with neighbors:** periodicity and inter-rank neighbors are handled by the
 *   halo exchange. Apply these BCs only on **physical faces** (no neighbor).
 * - **Preconditions:** interior extent along the operated axis must satisfy
 *   :code:`n_axis >= 1` (>=2 for `Extrapolate1`). The backing storage must be contiguous.
 *
 * @section boundary-ops Boundary operators (semantics & numerical order)
 *
 * - **Dirichlet** (`BCOp::Dirichlet`) — set a constant boundary value:
 *   :math:`f_{\text{ghost}} = c`. No derivative constraints are implied.
 *
 * - **NeumannZero** (`BCOp::NeumannZero`) — zero normal derivative,
 *   implemented as **constant extrapolation** using the nearest interior value:
 *   :math:`f(-l) = f(0)` on a minus face (and analogously on plus faces).
 *   This enforces :math:`\partial_n f = 0` with **first-order** accuracy at the face.
 *
 * - **Extrapolate1** (`BCOp::Extrapolate1`) — **linear** (first-order one-sided)
 *   extrapolation using the outward gradient estimated from the first two interior cells:
 *   on a minus face, :math:`d = f(0) - f(1)`, then
 *   :math:`f(-l) = f(0) + l\,d,\; l=1..ng`. Requires at least 2 interior cells along
 *   the operated axis. Accuracy is set by the gradient estimate (**first order**).
 *
 * - **Mirror (vector)** (`BCOp::Mirror`) — reflect interior values across the face with
 *   per-component parity. For a face with outward normal :math:`\hat{n}`, the mirrored
 *   value is :math:`\mathbf{u}_\text{ghost}(l) = \mathbf{S}\,\mathbf{u}_\text{int}(l-1)`
 *   where :math:`\mathbf{S}=\mathrm{diag}(s_x,s_y,s_z)` with :math:`s_i\in\{-1,+1\}`.
 *   Typical slip wall in x: :code:`MirrorMask{ {-1,+1,+1} }` (flip normal, keep tangentials).
 *   For scalars, `Mirror` is ignored in the dispatcher.
 *
 * - **Periodic** (`BCOp::Periodic`) — **no-op** here; handled by the halo exchange.
 *
 * @section boundary-complexity Complexity & threading
 *
 * Work is :math:`\mathcal{O}(\text{face\_area} \times ng)` per face. Loops are written
 * to stride contiguously in the fastest (I) dimension inside the innermost loop.
 * Threading over the two in-face indices is safe provided the same cells are not written
 * concurrently by multiple faces.
 *
 * @rst
 *.. code-block:: cpp
 *
 *   using core::mesh::Axis;
 *   using core::mesh::BCOp;
 *   using core::mesh::Array3DView;
 *   using core::mesh::MirrorMask;
 *
 *   // Build views from existing Field<T> (uniform ng inferred)
 *   Array3DView<double> rho(rho_field);
 *   Array3DView<double> ux (ux_field), uy(uy_field), uz(uz_field);
 *
 *   // Scalar BCs
 *   apply_scalar_bc(rho, Axis::I, -1, BCOp::NeumannZero);   // x- face: ∂n f = 0 (1st order)
 *   apply_scalar_bc(rho, Axis::K, +1, BCOp::Extrapolate1);  // z+ face: linear extrapolation
 *   apply_scalar_bc(rho, Axis::J, -1, BCOp::Dirichlet, 0.0);// y- face: f = 0
 *
 *   // Vector mirror (slip wall on x faces)
 *   MirrorMask slip_x{{-1, +1, +1}};
 *   apply_mirror_vector(ux, uy, uz, Axis::I, -1, slip_x);   // x- face
 *   apply_mirror_vector(ux, uy, uz, Axis::I, +1, slip_x);   // x+ face
 * @endrst
 *
 * @note These routines operate on **faces only**.
 */

namespace core::mesh
{

// ---------------------------------------------
// Lightweight 3D view that delegates to Field<T>
// (row-major: i fastest, then j, then k)
// ---------------------------------------------
template <class T> struct Array3DView
{
    static_assert(!std::is_const_v<T>, "Use a non-const T for Array3DView.");

    // Delegation target
    Field<T>* f = nullptr;

    // Interior extents and halos kept for geometry loops
    int nx = 0, ny = 0, nz = 0;
    int hx = 0, hy = 0, hz = 0;

    Array3DView() = default;

    // Build from a Field<T> (uniform ghosts inferred from Field)
    explicit Array3DView(Field<T>& field) : f(&field)
    {
        const int ng = field.ng();
        const auto e = field.extents(); // totals including ghosts
        nx = e[0] - 2 * ng;
        ny = e[1] - 2 * ng;
        nz = e[2] - 2 * ng;
        hx = hy = hz = ng;
    }

    // Build from a Field<T> with explicit interior extents (still uniform ghosts)
    Array3DView(Field<T>& field, int nx_in, int ny_in, int nz_in, int ng)
        : f(&field), nx(nx_in), ny(ny_in), nz(nz_in), hx(ng), hy(ng), hz(ng)
    {
    }

    // number of cells including halos per axis
    int sx() const noexcept { return nx + 2 * hx; }
    int sy() const noexcept { return ny + 2 * hy; }
    int sz() const noexcept { return nz + 2 * hz; }

    // linear index for (i,j,k) computed via Field (no Layout math here)
    inline std::size_t index(int i, int j, int k) const noexcept
    {
        // Pointer arithmetic relative to Field base; relies on Field::operator()
        return static_cast<std::size_t>(&((*f)(i, j, k)) - f->raw());
    }

    inline T& at(int i, int j, int k) noexcept { return (*f)(i, j, k); }
    inline const T& at(int i, int j, int k) const noexcept { return (*f)(i, j, k); }
};

// ---------------------------------------------
// Geometry helpers
// ---------------------------------------------
enum class Axis : uint8_t
{
    I = 0,
    J = 1,
    K = 2
};

// face_sign: -1 for "minus" face, +1 for "plus" face
constexpr inline bool is_minus(int face_sign) noexcept
{
    return face_sign < 0;
}
constexpr inline bool is_plus(int face_sign) noexcept
{
    return face_sign > 0;
}

// For MPI integration: in typical decompositions, a physical face means
// "no neighbor" — i.e., neighbor rank is MPI_PROC_NULL (usually negative).
// This helper avoids pulling in <mpi.h> from the header.
constexpr inline bool is_physical_face(int neighbor_rank_like) noexcept
{
    return neighbor_rank_like < 0; // true for MPI_PROC_NULL or any negative sentinel
}

// ---------------------------------------------
// Boundary condition vocabulary
// ---------------------------------------------
enum class BCOp : uint8_t
{
    Dirichlet,    // set ghost = value
    NeumannZero,  // copy nearest interior along normal
    Extrapolate1, // linear extrapolation using first interior gradient
    Mirror,       // reflect interior across face (vector: with sign mask)
    Periodic      // handled by halo exchange; not applied here
};

// Per-component sign for Mirror on vectors (1 or -1)
struct MirrorMask
{
    std::array<int8_t, 3> sign{1, 1, 1};
};

// ---------------------------------------------
// Core per-face operators (scalar)
// ---------------------------------------------
template <class T> inline void apply_dirichlet(Array3DView<T> a, Axis ax, int face_sign, T value)
{
    const int H = (ax == Axis::I ? a.hx : (ax == Axis::J ? a.hy : a.hz));
    if (H == 0)
        return;

    switch (ax)
    {
    case Axis::I:
    {
        if (is_minus(face_sign))
        {
            for (int k = 0; k < a.nz; ++k)
                for (int j = 0; j < a.ny; ++j)
                    for (int l = 1; l <= a.hx; ++l)
                        a.at(-l, j, k) = value;
        }
        else
        {
            for (int k = 0; k < a.nz; ++k)
                for (int j = 0; j < a.ny; ++j)
                    for (int l = 1; l <= a.hx; ++l)
                        a.at(a.nx + (l - 1), j, k) = value;
        }
        break;
    }
    case Axis::J:
    {
        if (is_minus(face_sign))
        {
            for (int k = 0; k < a.nz; ++k)
                for (int i = 0; i < a.nx; ++i)
                    for (int l = 1; l <= a.hy; ++l)
                        a.at(i, -l, k) = value;
        }
        else
        {
            for (int k = 0; k < a.nz; ++k)
                for (int i = 0; i < a.nx; ++i)
                    for (int l = 1; l <= a.hy; ++l)
                        a.at(i, a.ny + (l - 1), k) = value;
        }
        break;
    }
    case Axis::K:
    {
        if (is_minus(face_sign))
        {
            for (int j = 0; j < a.ny; ++j)
                for (int i = 0; i < a.nx; ++i)
                    for (int l = 1; l <= a.hz; ++l)
                        a.at(i, j, -l) = value;
        }
        else
        {
            for (int j = 0; j < a.ny; ++j)
                for (int i = 0; i < a.nx; ++i)
                    for (int l = 1; l <= a.hz; ++l)
                        a.at(i, j, a.nz + (l - 1)) = value;
        }
        break;
    }
    }
}

template <class T> inline void apply_neumann_zero(Array3DView<T> a, Axis ax, int face_sign)
{
    const int H = (ax == Axis::I ? a.hx : (ax == Axis::J ? a.hy : a.hz));
    if (H == 0)
        return;

    switch (ax)
    {
    case Axis::I:
    {
        if (is_minus(face_sign))
        {
            for (int k = 0; k < a.nz; ++k)
                for (int j = 0; j < a.ny; ++j)
                    for (int l = 1; l <= a.hx; ++l)
                        a.at(-l, j, k) = a.at(0, j, k);
        }
        else
        {
            for (int k = 0; k < a.nz; ++k)
                for (int j = 0; j < a.ny; ++j)
                    for (int l = 1; l <= a.hx; ++l)
                        a.at(a.nx + (l - 1), j, k) = a.at(a.nx - 1, j, k);
        }
        break;
    }
    case Axis::J:
    {
        if (is_minus(face_sign))
        {
            for (int k = 0; k < a.nz; ++k)
                for (int i = 0; i < a.nx; ++i)
                    for (int l = 1; l <= a.hy; ++l)
                        a.at(i, -l, k) = a.at(i, 0, k);
        }
        else
        {
            for (int k = 0; k < a.nz; ++k)
                for (int i = 0; i < a.nx; ++i)
                    for (int l = 1; l <= a.hy; ++l)
                        a.at(i, a.ny + (l - 1), k) = a.at(i, a.ny - 1, k);
        }
        break;
    }
    case Axis::K:
    {
        if (is_minus(face_sign))
        {
            for (int j = 0; j < a.ny; ++j)
                for (int i = 0; i < a.nx; ++i)
                    for (int l = 1; l <= a.hz; ++l)
                        a.at(i, j, -l) = a.at(i, j, 0);
        }
        else
        {
            for (int j = 0; j < a.ny; ++j)
                for (int i = 0; i < a.nx; ++i)
                    for (int l = 1; l <= a.hz; ++l)
                        a.at(i, j, a.nz + (l - 1)) = a.at(i, j, a.nz - 1);
        }
        break;
    }
    }
}

template <class T> inline void apply_extrapolate1(Array3DView<T> a, Axis ax, int face_sign)
{
    const int H = (ax == Axis::I ? a.hx : (ax == Axis::J ? a.hy : a.hz));
    if (H == 0)
        return;

    switch (ax)
    {
    case Axis::I:
    {
        if (is_minus(face_sign))
        {
            for (int k = 0; k < a.nz; ++k)
                for (int j = 0; j < a.ny; ++j)
                {
                    const T f0 = a.at(0, j, k);
                    const T f1 = a.at(1, j, k);
                    const T d = f0 - f1; // outward gradient
                    for (int l = 1; l <= a.hx; ++l)
                        a.at(-l, j, k) = f0 + static_cast<T>(l) * d;
                }
        }
        else
        {
            for (int k = 0; k < a.nz; ++k)
                for (int j = 0; j < a.ny; ++j)
                {
                    const T f0 = a.at(a.nx - 1, j, k);
                    const T f1 = a.at(a.nx - 2, j, k);
                    const T d = f0 - f1;
                    for (int l = 1; l <= a.hx; ++l)
                        a.at(a.nx + (l - 1), j, k) = f0 + static_cast<T>(l) * d;
                }
        }
        break;
    }
    case Axis::J:
    {
        if (is_minus(face_sign))
        {
            for (int k = 0; k < a.nz; ++k)
                for (int i = 0; i < a.nx; ++i)
                {
                    const T f0 = a.at(i, 0, k);
                    const T f1 = a.at(i, 1, k);
                    const T d = f0 - f1;
                    for (int l = 1; l <= a.hy; ++l)
                        a.at(i, -l, k) = f0 + static_cast<T>(l) * d;
                }
        }
        else
        {
            for (int k = 0; k < a.nz; ++k)
                for (int i = 0; i < a.nx; ++i)
                {
                    const T f0 = a.at(i, a.ny - 1, k);
                    const T f1 = a.at(i, a.ny - 2, k);
                    const T d = f0 - f1;
                    for (int l = 1; l <= a.hy; ++l)
                        a.at(i, a.ny + (l - 1), k) = f0 + static_cast<T>(l) * d;
                }
        }
        break;
    }
    case Axis::K:
    {
        if (is_minus(face_sign))
        {
            for (int j = 0; j < a.ny; ++j)
                for (int i = 0; i < a.nx; ++i)
                {
                    const T f0 = a.at(i, j, 0);
                    const T f1 = a.at(i, j, 1);
                    const T d = f0 - f1;
                    for (int l = 1; l <= a.hz; ++l)
                        a.at(i, j, -l) = f0 + static_cast<T>(l) * d;
                }
        }
        else
        {
            for (int j = 0; j < a.ny; ++j)
                for (int i = 0; i < a.nx; ++i)
                {
                    const T f0 = a.at(i, j, a.nz - 1);
                    const T f1 = a.at(i, j, a.nz - 2);
                    const T d = f0 - f1;
                    for (int l = 1; l <= a.hz; ++l)
                        a.at(i, j, a.nz + (l - 1)) = f0 + static_cast<T>(l) * d;
                }
        }
        break;
    }
    }
}

// ---------------------------------------------
// Mirror for vector fields U=(ux,uy,uz)
// Mirror means reflect across the face; mask controls component parity.
// For minus face: ghost(-l) = sign * interior(l-1)
// For plus  face: ghost(n + (l-1)) = sign * interior(n - l)
// ---------------------------------------------
template <class T>
inline void apply_mirror_vector(Array3DView<T> ux, Array3DView<T> uy, Array3DView<T> uz, Axis ax,
                                int face_sign, MirrorMask mask)
{
    auto mirror_scalar = [&](Array3DView<T> a, int sgn_normal)
    {
        switch (ax)
        {
        case Axis::I:
        {
            if (is_minus(face_sign))
            {
                for (int k = 0; k < a.nz; ++k)
                    for (int j = 0; j < a.ny; ++j)
                        for (int l = 1; l <= a.hx; ++l)
                            a.at(-l, j, k) = static_cast<T>(sgn_normal) * a.at(l - 1, j, k);
            }
            else
            {
                for (int k = 0; k < a.nz; ++k)
                    for (int j = 0; j < a.ny; ++j)
                        for (int l = 1; l <= a.hx; ++l)
                            a.at(a.nx + (l - 1), j, k) =
                                static_cast<T>(sgn_normal) * a.at(a.nx - l, j, k);
            }
            break;
        }
        case Axis::J:
        {
            if (is_minus(face_sign))
            {
                for (int k = 0; k < a.nz; ++k)
                    for (int i = 0; i < a.nx; ++i)
                        for (int l = 1; l <= a.hy; ++l)
                            a.at(i, -l, k) = static_cast<T>(sgn_normal) * a.at(i, l - 1, k);
            }
            else
            {
                for (int k = 0; k < a.nz; ++k)
                    for (int i = 0; i < a.nx; ++i)
                        for (int l = 1; l <= a.hy; ++l)
                            a.at(i, a.ny + (l - 1), k) =
                                static_cast<T>(sgn_normal) * a.at(i, a.ny - l, k);
            }
            break;
        }
        case Axis::K:
        {
            if (is_minus(face_sign))
            {
                for (int j = 0; j < a.ny; ++j)
                    for (int i = 0; i < a.nx; ++i)
                        for (int l = 1; l <= a.hz; ++l)
                            a.at(i, j, -l) = static_cast<T>(sgn_normal) * a.at(i, j, l - 1);
            }
            else
            {
                for (int j = 0; j < a.ny; ++j)
                    for (int i = 0; i < a.nx; ++i)
                        for (int l = 1; l <= a.hz; ++l)
                            a.at(i, j, a.nz + (l - 1)) =
                                static_cast<T>(sgn_normal) * a.at(i, j, a.nz - l);
            }
            break;
        }
        }
    };

    // Interpret which component is "normal" to the face and apply signs
    if (ax == Axis::I)
    {
        mirror_scalar(ux, mask.sign[0]); // normal component
        mirror_scalar(uy, mask.sign[1]); // tangential
        mirror_scalar(uz, mask.sign[2]); // tangential
    }
    else if (ax == Axis::J)
    {
        mirror_scalar(ux, mask.sign[0]);
        mirror_scalar(uy, mask.sign[1]);
        mirror_scalar(uz, mask.sign[2]);
    }
    else
    { // K
        mirror_scalar(ux, mask.sign[0]);
        mirror_scalar(uy, mask.sign[1]);
        mirror_scalar(uz, mask.sign[2]);
    }
}

// ---------------------------------------------
// Convenience dispatcher for scalar BCs
// ---------------------------------------------
template <class T>
inline void apply_scalar_bc(Array3DView<T> a, Axis ax, int face_sign, BCOp op, T value = T{})
{
    switch (op)
    {
    case BCOp::Dirichlet:
        apply_dirichlet(a, ax, face_sign, value);
        break;
    case BCOp::NeumannZero:
        apply_neumann_zero(a, ax, face_sign);
        break;
    case BCOp::Extrapolate1:
        apply_extrapolate1(a, ax, face_sign);
        break;
    case BCOp::Mirror: /* not meaningful for scalar; ignore */
        break;
    case BCOp::Periodic: /* handled by halo exchange */
        break;
    }
}

} // namespace core::mesh
