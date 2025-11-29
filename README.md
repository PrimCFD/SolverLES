<p align="center">
  <img src="assets/logo/KolmoPlasLogo512.png" alt="KolmoPlas logo" width="160">
</p>

# KolmoPlas

> **Status:** _Early stage (hydrodynamics)_

[![Docs](https://img.shields.io/badge/docs-online-blue)](https://primcfd.github.io/SolverLES/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[![Linux CI](https://github.com/PrimCFD/SolverLES/actions/workflows/linux.yml/badge.svg?branch=main)](https://github.com/PrimCFD/SolverLES/actions/workflows/linux.yml)
[![Docs](https://github.com/PrimCFD/SolverLES/actions/workflows/docs.yml/badge.svg?branch=main)](https://github.com/PrimCFD/SolverLES/actions/workflows/docs.yml)
[![Style](https://github.com/PrimCFD/SolverLES/actions/workflows/style.yml/badge.svg?branch=main)](https://github.com/PrimCFD/SolverLES/actions/workflows/style.yml)


---

## 1&nbsp;· Overview

The project is an **open, modular finite‑volume LES solver for thermal–plasma flows** that lets researchers drop in new sub‑grid‑scale (SGS) models with minimal friction while still scaling to modern CPU & GPU clusters. The project is intended to research new SGS models for thermal plasma jets simulations accounting for steep density and temperature gradients.

---

## 2&nbsp;· Repository Layout (high‑level)

```text
SolverLES/

├─ .github/workflows/        # CI definitions
├─ assets/                   # Doc and visual assets
├─ cmake/                    # CMake helper modules
├─ docs/                     # Sphinx/Doxygen furo doc
├─ examples/                 # Tiny example run‑ready cases
├─ extern/                   # Third‑party sources (tar balls mirrored, no edits)
├─ scripts/                  # Dev‑ops & helper scripts
├─ src/                      # Solver source code
│   ├─ core/                 # C++ runtime static lib (orchestration, memory management)
│   ├─ gui/                  # Qt/VTK front‑end
│   ├─ physics/              # Hot‑swappable physics modules dynamic lib
│   ├─ kernels/              # Shared Fortran math kernels
│   ├─ bindings/             # C/Fortran interop helpers
│   ├─ ipc/                  # Inter-process communication GUI/Solver
│   └─ apps/                 # Executable entry points
├─ tests/                    # Unit, regression & perf tests
├─ .clang-format             # C++ style
├─ .cmake-format.yml         # CMake style
├─ .fprettify.yml            # Fortran style
├─ .gitignore
├─ CMakeLists.txt            # Root super‑build
├─ LICENSE
└─ README.md                 # This file
```

---

## 3&nbsp;· Directory Details

| Path | Role | Notes |
|------|------|-------|
| **src/core** | Owns mesh, fields, time loop, I/O; _no physics_. | Written in modern C++ (C++20). |
| **src/physics** | Physics shared library plugins | Each physics → one shared lib (e.g. `libfluids.so`). |
| **src/physics/kernels** | Optimised Fortran math; reused by plug‑ins. | Vectorised / GPU‑offloaded via OpenMP/CUDA Fortran. |
| **extern** | Vendored: CGNS, PETSc, etc. | Pulled via `FetchContent`/`ExternalProject`; do **not** modify in‑tree. |
| **tests** | Unit, regression CGNS, etc | CI runs these on every PR. |
| **examples** | Minimal decks: coarse Taylor–Green, etc. | Should finish < 1 minute serial. |

---

## 4&nbsp;· Quick Start (using `scripts/` helpers)

This project ships convenience scripts in `scripts/` for reliable, repeatable developer workflows (builds, docs, MPI, offline vendor cache, cleaning, formatting). See the cheatsheet in that folder for details.  

### 4.1 Prerequisites
- **CMake ≥ 3.24**, a C/C++/Fortran toolchain; **Ninja** is recommended.
- **Parrallel stack** (OpenMP) technically overridable but no gains.
- Optional/when needed:
  - **MPI stack** (Open MPI) for MPI builds and tests.
  - **Doxygen** and **Sphinx** (`sphinx-build`) for docs.
  - Formatters: `clang-format`, `fprettify`, `cmake-format`. (if using a Python virtual environment, run `source venv/bin/activate`)


> Tip: CI uses recent gcc on Linux, might be tweaking required for clang; matching that locally avoids surprises (see §5 CI). System gcc + OpenMP typically is the HPC setup.

### 4.2 Fast path (CPU, Release)
```bash
# Configure + build (Release by default) into ./build
./scripts/build.sh

# Run the hello-mesh example
./build/bin/solver examples/hello_mesh.yaml
```

### 4.3 Common switches
```bash
# Debug build into a custom directory
BUILD_DIR=build-debug CMAKE_BUILD_TYPE=Debug ./scripts/build.sh

# Extra CMake options (examples)
EXTRA_CMAKE_ARGS="-DBUILD_EXAMPLES=ON -DBUILD_GUI=OFF" ./scripts/build.sh
```

### 4.4 Offline / reproducible third-party cache
Pre-download third-party sources into `extern/` and build fully offline later:
```bash
# Populate extern/ with reproducible archives + MANIFEST.prefetch + SHA256SUMS
./scripts/prefetch_third_party.sh

# Then force an offline build (uses the cached archives)
OFFLINE=1 ./scripts/build.sh
```

### 4.5 MPI quickstart (laptop-friendly)
```bash
# Prepare a consistent MPI launcher env (auto-detects vendor/launcher)
source scripts/mpi_env.sh auto

# Launch with N ranks (works across srun/mpirun/mpiexec)
mpi_exec 4 ./build/bin/solver examples/hello_mesh.yaml

# Or run MPI enqbled test cases end-to-end
./scripts/run_regression_tests.sh
```

### 4.6 Docs build & local preview
```bash
# Build Doxygen XML + Sphinx HTML
./scripts/build_docs.sh

# Serve locally and open a browser tab
./scripts/build_docs.sh --serve --open
```

### 4.7 Cleaning up build trees & vendor cache
```bash
# Clean CMake/Ninja artifacts (keeps build-docs/)
./scripts/clean_build.sh --all -y

# Optionally wipe vendored archives/sources under extern/
./scripts/clean_extern.sh -y
```

### 4.8 Format sources (C/C++, Fortran, CMake)
```bash
# Best-effort (ignores missing tools)
./scripts/format_all.sh

# Fail if any formatter is missing or errors
STRICT=1 ./scripts/format_all.sh
```

---

## 5&nbsp;· Continuous Integration
[![Linux CI](https://github.com/PrimCFD/SolverLES/actions/workflows/linux.yml/badge.svg?branch=main)](https://github.com/PrimCFD/SolverLES/actions/workflows/linux.yml) [![Docs](https://github.com/PrimCFD/SolverLES/actions/workflows/docs.yml/badge.svg?branch=main)](https://github.com/PrimCFD/SolverLES/actions/workflows/docs.yml) [![Style](https://github.com/PrimCFD/SolverLES/actions/workflows/style.yml/badge.svg?branch=main)](https://github.com/PrimCFD/SolverLES/actions/workflows/style.yml)

All CI is always built from clean slate (containerized) to check the whole pipeline on Hardware with missing tools

* **linux.yml** – GCC 13 / Clang 18 matrix; runs unit, performance & regression tests (GPU tests need hosted runner).
* **style.yml** – clang‑format, fprettify, cmake‑lint.
* **docs.yml** – builds Sphinx docs, pushes to `gh-pages`. (run on GitHub manually)

---

## 6&nbsp;· License

Distributed under the **MIT License** as outlined in `LICENSE`.

