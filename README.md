# FVM–LES–Plasma Solver

> **Status:** _Initial draft._

---

## 1&nbsp;· Project Vision

The objective is to build an **open, modular finite‑volume LES solver for thermal–plasma flows** that lets researchers drop in new sub‑grid‑scale (SGS) models with minimal friction while still scaling to modern CPU & GPU clusters. The project is intended to research new SGS models for thermal plasma jets simulations accounting for steep density and temperature gradients.

---

## 2&nbsp;· Repository Layout (high‑level)

```text
fvm-les-plasma/
├─ CMakeLists.txt            # Root super‑build
├─ .clang-format             # C++ style
├─ .fprettifyrc              # Fortran style
├─ .gitignorev
├─ LICENSE
├─ README.md                 # Newcomer entry point
│
├─ cmake/                    # CMake helper modules
│   ├─ FetchKokkos.cmake
│   ├─ FetchCGNS.cmake
│   └─ CompileOptions.cmake
│
├─ extern/                   # Third‑party sources (mirrored, no edits)
│   ├─ petsc/
│   └─ cgns/
│
├─ src/                      # Solver source code
│   ├─ core/                 # C++ runtime (orchestration only)
│   ├─ gui/                  # Qt front‑end
│   ├─ plugins/              # Hot‑swappable physics modules
│   │   ├─ flux/             # ‑ IFluxScheme implementations
│   │   ├─ sgs/              # ‑ ISGSModel implementations
│   │   └─ time/             # ‑ ITimeStepper implementations
│   ├─ kernels/              # Shared Fortran math kernels
│   └─ bindings/             # C/Fortran interop helpers
│
├─ tests/                    # Unit, regression & perf tests
├─ examples/                 # Tiny run‑ready cases (< 60 s)
├─ benchmarks/               # Large scaling cases (LFS)
│
├─ docs/                     # Sphinx + architecture notes
├─ scripts/                  # Dev‑ops & helper scripts
└─ .github/                  # CI definitions
```

---

## 3&nbsp;· Directory Details

| Path | Role | Notes |
|------|------|-------|
| **src/core** | Owns mesh, fields, time loop, I/O; _no physics_. | Written in modern C++ (C++20). |
| **src/plugins/flux** | Spatial discretisation schemes. | Each sub‑dir → one shared lib (e.g. `libweno.so`). |
| **src/plugins/sgs** | LES SGS models. | Research mostly done here. |
| **src/plugins/time** | Time‑integration schemes. | Explicit RK, implicit, etc. |
| **src/kernels** | Optimised Fortran math; reused by plug‑ins. | Vectorised / GPU‑offloaded via >OpenMP/CUDA Fortran. |
| **extern** | Vendored: CGNS, PETSc. | Pulled via `FetchContent`; do **not** modify in‑tree. |
| **tests** | Unit, regression CGNS, etc | CI runs these on every PR. |
| **examples** | Minimal decks: Taylor–Green, Sod, etc. | Should finish < 1 minute serial. |
| **benchmarks** | 256³+ meshes for scaling. | Stored via Git‑LFS or external URL. |

---

## 4&nbsp;· Build & Install Quick Start

```bash
# Clone and configure (CPU only, Debug)
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build -j8

# Run the hello‑mesh example
./build/bin/solver examples/hello_mesh.yaml
```

*See `/docs/developer_guide/build.md` for GPU flags and cluster hints.*

---

## 5&nbsp;· Continuous Integration

* **linux.yml** – GCC 13 / Clang 18 matrix; runs unit & regression tests.
* **style.yml** – clang‑format, fprettify, cmake‑lint.
* **docs.yml** – builds Sphinx docs, pushes to `gh-pages`.

---

## 6&nbsp;· License

Distributed under the **MIT License** as outlined in `LICENSE`.

