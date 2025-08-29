#pragma once
#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>
#include <type_traits>

namespace mesh {

// ---------------------------------------------
// Lightweight 3D view for contiguous arrays with halos
// (row-major: i fastest, then j, then k)
// ---------------------------------------------
template <class T>
struct Array3DView {
  static_assert(!std::is_const_v<T>, "Use a non-const T for Array3DView.");
  T* data = nullptr;          // base pointer
  int nx = 0, ny = 0, nz = 0; // interior extents
  int hx = 0, hy = 0, hz = 0; // halo width per axis

  // number of cells including halos per axis
  int sx() const noexcept { return nx + 2 * hx; }
  int sy() const noexcept { return ny + 2 * hy; }
  int sz() const noexcept { return nz + 2 * hz; }

  // linear index for (i,j,k) where i∈[-hx, nx+hx-1], etc.
  inline std::size_t index(int i, int j, int k) const noexcept {
    const int ii = i + hx;
    const int jj = j + hy;
    const int kk = k + hz;
    return static_cast<std::size_t>(ii) +
           static_cast<std::size_t>(sx()) *
               (static_cast<std::size_t>(jj) +
                static_cast<std::size_t>(sy()) * static_cast<std::size_t>(kk));
  }

  inline T& at(int i, int j, int k) noexcept { return data[index(i, j, k)]; }
  inline const T& at(int i, int j, int k) const noexcept { return data[index(i, j, k)]; }
};

// ---------------------------------------------
// Geometry helpers
// ---------------------------------------------
enum class Axis : uint8_t { I = 0, J = 1, K = 2 };

// face_sign: -1 for "minus" face, +1 for "plus" face
constexpr inline bool is_minus(int face_sign) noexcept { return face_sign < 0; }
constexpr inline bool is_plus (int face_sign) noexcept { return face_sign > 0; }

// For MPI integration: in typical decompositions, a physical face means
// "no neighbor" — i.e., neighbor rank is MPI_PROC_NULL (usually negative).
// This helper avoids pulling in <mpi.h> from the header.
constexpr inline bool is_physical_face(int neighbor_rank_like) noexcept {
  return neighbor_rank_like < 0; // true for MPI_PROC_NULL or any negative sentinel
}

// ---------------------------------------------
// Boundary condition vocabulary
// ---------------------------------------------
enum class BCOp : uint8_t {
  Dirichlet,     // set ghost = value
  NeumannZero,   // copy nearest interior along normal
  Extrapolate1,  // linear extrapolation using first interior gradient
  Mirror,        // reflect interior across face (vector: with sign mask)
  Periodic       // handled by halo exchange; not applied here
};

// Per-component sign for Mirror on vectors (1 or -1)
struct MirrorMask {
  std::array<int8_t, 3> sign {1,1,1};
};

// ---------------------------------------------
// Core per-face operators (scalar)
// ---------------------------------------------
template <class T>
inline void apply_dirichlet(Array3DView<T> a, Axis ax, int face_sign, T value)
{
  const int H = (ax == Axis::I ? a.hx : (ax == Axis::J ? a.hy : a.hz));
  if (H == 0) return;

  switch (ax) {
    case Axis::I: {
      if (is_minus(face_sign)) {
        for (int k = 0; k < a.nz; ++k)
          for (int j = 0; j < a.ny; ++j)
            for (int l = 1; l <= a.hx; ++l)
              a.at(-l, j, k) = value;
      } else {
        for (int k = 0; k < a.nz; ++k)
          for (int j = 0; j < a.ny; ++j)
            for (int l = 1; l <= a.hx; ++l)
              a.at(a.nx + (l - 1), j, k) = value;
      }
      break;
    }
    case Axis::J: {
      if (is_minus(face_sign)) {
        for (int k = 0; k < a.nz; ++k)
          for (int i = 0; i < a.nx; ++i)
            for (int l = 1; l <= a.hy; ++l)
              a.at(i, -l, k) = value;
      } else {
        for (int k = 0; k < a.nz; ++k)
          for (int i = 0; i < a.nx; ++i)
            for (int l = 1; l <= a.hy; ++l)
              a.at(i, a.ny + (l - 1), k) = value;
      }
      break;
    }
    case Axis::K: {
      if (is_minus(face_sign)) {
        for (int j = 0; j < a.ny; ++j)
          for (int i = 0; i < a.nx; ++i)
            for (int l = 1; l <= a.hz; ++l)
              a.at(i, j, -l) = value;
      } else {
        for (int j = 0; j < a.ny; ++j)
          for (int i = 0; i < a.nx; ++i)
            for (int l = 1; l <= a.hz; ++l)
              a.at(i, j, a.nz + (l - 1)) = value;
      }
      break;
    }
  }
}

template <class T>
inline void apply_neumann_zero(Array3DView<T> a, Axis ax, int face_sign)
{
  const int H = (ax == Axis::I ? a.hx : (ax == Axis::J ? a.hy : a.hz));
  if (H == 0) return;

  switch (ax) {
    case Axis::I: {
      if (is_minus(face_sign)) {
        for (int k = 0; k < a.nz; ++k)
          for (int j = 0; j < a.ny; ++j)
            for (int l = 1; l <= a.hx; ++l)
              a.at(-l, j, k) = a.at(0, j, k);
      } else {
        for (int k = 0; k < a.nz; ++k)
          for (int j = 0; j < a.ny; ++j)
            for (int l = 1; l <= a.hx; ++l)
              a.at(a.nx + (l - 1), j, k) = a.at(a.nx - 1, j, k);
      }
      break;
    }
    case Axis::J: {
      if (is_minus(face_sign)) {
        for (int k = 0; k < a.nz; ++k)
          for (int i = 0; i < a.nx; ++i)
            for (int l = 1; l <= a.hy; ++l)
              a.at(i, -l, k) = a.at(i, 0, k);
      } else {
        for (int k = 0; k < a.nz; ++k)
          for (int i = 0; i < a.nx; ++i)
            for (int l = 1; l <= a.hy; ++l)
              a.at(i, a.ny + (l - 1), k) = a.at(i, a.ny - 1, k);
      }
      break;
    }
    case Axis::K: {
      if (is_minus(face_sign)) {
        for (int j = 0; j < a.ny; ++j)
          for (int i = 0; i < a.nx; ++i)
            for (int l = 1; l <= a.hz; ++l)
              a.at(i, j, -l) = a.at(i, j, 0);
      } else {
        for (int j = 0; j < a.ny; ++j)
          for (int i = 0; i < a.nx; ++i)
            for (int l = 1; l <= a.hz; ++l)
              a.at(i, j, a.nz + (l - 1)) = a.at(i, j, a.nz - 1);
      }
      break;
    }
  }
}

template <class T>
inline void apply_extrapolate1(Array3DView<T> a, Axis ax, int face_sign)
{
  const int H = (ax == Axis::I ? a.hx : (ax == Axis::J ? a.hy : a.hz));
  if (H == 0) return;

  switch (ax) {
    case Axis::I: {
      if (is_minus(face_sign)) {
        for (int k = 0; k < a.nz; ++k)
          for (int j = 0; j < a.ny; ++j) {
            const T f0 = a.at(0, j, k);
            const T f1 = a.at(1, j, k);
            const T d  = f0 - f1; // outward gradient
            for (int l = 1; l <= a.hx; ++l)
              a.at(-l, j, k) = f0 + static_cast<T>(l) * d;
          }
      } else {
        for (int k = 0; k < a.nz; ++k)
          for (int j = 0; j < a.ny; ++j) {
            const T f0 = a.at(a.nx - 1, j, k);
            const T f1 = a.at(a.nx - 2, j, k);
            const T d  = f0 - f1;
            for (int l = 1; l <= a.hx; ++l)
              a.at(a.nx + (l - 1), j, k) = f0 + static_cast<T>(l) * d;
          }
      }
      break;
    }
    case Axis::J: {
      if (is_minus(face_sign)) {
        for (int k = 0; k < a.nz; ++k)
          for (int i = 0; i < a.nx; ++i) {
            const T f0 = a.at(i, 0, k);
            const T f1 = a.at(i, 1, k);
            const T d  = f0 - f1;
            for (int l = 1; l <= a.hy; ++l)
              a.at(i, -l, k) = f0 + static_cast<T>(l) * d;
          }
      } else {
        for (int k = 0; k < a.nz; ++k)
          for (int i = 0; i < a.nx; ++i) {
            const T f0 = a.at(i, a.ny - 1, k);
            const T f1 = a.at(i, a.ny - 2, k);
            const T d  = f0 - f1;
            for (int l = 1; l <= a.hy; ++l)
              a.at(i, a.ny + (l - 1), k) = f0 + static_cast<T>(l) * d;
          }
      }
      break;
    }
    case Axis::K: {
      if (is_minus(face_sign)) {
        for (int j = 0; j < a.ny; ++j)
          for (int i = 0; i < a.nx; ++i) {
            const T f0 = a.at(i, j, 0);
            const T f1 = a.at(i, j, 1);
            const T d  = f0 - f1;
            for (int l = 1; l <= a.hz; ++l)
              a.at(i, j, -l) = f0 + static_cast<T>(l) * d;
          }
      } else {
        for (int j = 0; j < a.ny; ++j)
          for (int i = 0; i < a.nx; ++i) {
            const T f0 = a.at(i, j, a.nz - 1);
            const T f1 = a.at(i, j, a.nz - 2);
            const T d  = f0 - f1;
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
inline void apply_mirror_vector(Array3DView<T> ux,
                                Array3DView<T> uy,
                                Array3DView<T> uz,
                                Axis ax, int face_sign, MirrorMask mask)
{
  auto mirror_scalar = [&](Array3DView<T> a, int sgn_normal) {
    switch (ax) {
      case Axis::I: {
        if (is_minus(face_sign)) {
          for (int k = 0; k < a.nz; ++k)
            for (int j = 0; j < a.ny; ++j)
              for (int l = 1; l <= a.hx; ++l)
                a.at(-l, j, k) = static_cast<T>(sgn_normal) * a.at(l - 1, j, k);
        } else {
          for (int k = 0; k < a.nz; ++k)
            for (int j = 0; j < a.ny; ++j)
              for (int l = 1; l <= a.hx; ++l)
                a.at(a.nx + (l - 1), j, k) =
                    static_cast<T>(sgn_normal) * a.at(a.nx - l, j, k);
        }
        break;
      }
      case Axis::J: {
        if (is_minus(face_sign)) {
          for (int k = 0; k < a.nz; ++k)
            for (int i = 0; i < a.nx; ++i)
              for (int l = 1; l <= a.hy; ++l)
                a.at(i, -l, k) = static_cast<T>(sgn_normal) * a.at(i, l - 1, k);
        } else {
          for (int k = 0; k < a.nz; ++k)
            for (int i = 0; i < a.nx; ++i)
              for (int l = 1; l <= a.hy; ++l)
                a.at(i, a.ny + (l - 1), k) =
                    static_cast<T>(sgn_normal) * a.at(i, a.ny - l, k);
        }
        break;
      }
      case Axis::K: {
        if (is_minus(face_sign)) {
          for (int j = 0; j < a.ny; ++j)
            for (int i = 0; i < a.nx; ++i)
              for (int l = 1; l <= a.hz; ++l)
                a.at(i, j, -l) = static_cast<T>(sgn_normal) * a.at(i, j, l - 1);
        } else {
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
  if (ax == Axis::I) {
    mirror_scalar(ux, mask.sign[0]); // normal component
    mirror_scalar(uy, mask.sign[1]); // tangential
    mirror_scalar(uz, mask.sign[2]); // tangential
  } else if (ax == Axis::J) {
    mirror_scalar(ux, mask.sign[0]);
    mirror_scalar(uy, mask.sign[1]);
    mirror_scalar(uz, mask.sign[2]);
  } else { // K
    mirror_scalar(ux, mask.sign[0]);
    mirror_scalar(uy, mask.sign[1]);
    mirror_scalar(uz, mask.sign[2]);
  }
}

// ---------------------------------------------
// Convenience dispatcher for scalar BCs
// ---------------------------------------------
template <class T>
inline void apply_scalar_bc(Array3DView<T> a, Axis ax, int face_sign,
                            BCOp op, T value = T{})
{
  switch (op) {
    case BCOp::Dirichlet:    apply_dirichlet(a, ax, face_sign, value); break;
    case BCOp::NeumannZero:  apply_neumann_zero(a, ax, face_sign);     break;
    case BCOp::Extrapolate1: apply_extrapolate1(a, ax, face_sign);     break;
    case BCOp::Mirror:       /* not meaningful for scalar; ignore */   break;
    case BCOp::Periodic:     /* handled by halo exchange */            break;
  }
}

} // namespace mesh
