# `scripts/` Developer Workflow Cheatsheet

## TL;DR table

| Script | Purpose | Typical use |
|---|---|---|
| `build.sh` | Configure & build the project via CMake (env‑driven). | `./scripts/build.sh` |
| `build_docs.sh` | Build clean Sphinx HTML (+ Doxygen XML) and optionally serve locally. | `./scripts/build_docs.sh --serve --open` |
| `clean_build.sh` | Remove CMake/Ninja build trees (`build-*`, keeps docs). | `./scripts/clean_build.sh --all -y` |
| `clean_extern.sh` | Wipe vendored third‑party sources under `extern/` (keeps `extern/README.md`). | `./scripts/clean_extern.sh -y` |
| `format_all.sh` | Run code formatters for C/C++, Fortran, and CMake (skips vendor/build). | `./scripts/format_all.sh` |
| `mpi_env.sh` | Sourceable helpers to set MPI launcher/env; provides `mpi_exec`. | `source scripts/mpi_env.sh auto` |
| `prefetch_third_party.sh` | Pre‑download third‑party sources and create reproducible archives in `extern/`. | `./scripts/prefetch_third_party.sh` |
| `push_ci.sh` | Push the current branch to GitHub and runs the CI on a locally hosted runner (tmux terminal). | `./scripts/push_ci.sh` |
| `run_unit_tests.sh` | Build (optional) and run **unit** tests via CTest; write JUnit XML. | `./scripts/run_unit_tests.sh` |
| `run_perf_tests.sh` | Build (optional) and run **perf** tests via CTest; write JUnit XML. | `./scripts/run_perf_tests.sh` |
| `run_regression_tests.sh` | Configure with MPI and run **MPI** integration tests, laptop/cluster friendly. | `./scripts/run_regression_tests.sh` |

> Unless stated otherwise, run scripts from the repo root. Many respect `BUILD_DIR` if you want custom build folders.

---

## Common prerequisites

- **CMake ≥ 3.24** and a C/C++/Fortran toolchain. Ninja recommended.
- **Python 3** (only for `build_docs.sh --serve`).
- Optional formatters: `clang-format`, `fprettify`, `cmake-format` (`STRICT=1` makes formatting failures fatal).
- For MPI workflows: an MPI stack (Open MPI/MPICH) and, on clusters, a launcher such as `srun`.
- For docs: `doxygen` and `sphinx-build` (the script will report if missing).

---

## `build.sh` — one‑stop CMake build

**What it does**  
- Picks a generator (prefers Ninja), configures CMake in `BUILD_DIR` and builds with parallel jobs.  
- Respects and forwards FetchContent & third‑party hints; can work **offline** if archives are present in `extern/`.  
- Filters noisy third‑party subbuild logs while still surfacing real errors.

**Key env vars (all optional)**  
- `BUILD_DIR` — build directory (default: `build`).  
- `CMAKE_BUILD_TYPE` — `Debug|Release|RelWithDebInfo|MinSizeRel` (default: `Release`).  
- `BUILD_TESTS` — `ON|OFF` (default: `ON`).  
- `MPIEXEC_PREFLAGS` — forwarded to CMake/CTest (e.g. `--oversubscribe`).  
- `CMAKE_GENERATOR` — override generator (e.g. `Ninja`).  
- `CMAKE_TOOLCHAIN_FILE` — forwarded as‑is.  
- `EXTRA_CMAKE_ARGS` — extra space‑separated CMake args.  
- `OFFLINE=1` — force disconnected `FetchContent` if you’ve already mirrored deps in `extern/`.
- `NPROCS` — override auto‑detected parallel job count.

**Examples**  
```bash
# Default Release build with Ninja
./scripts/build.sh

# Debug build into a custom dir
BUILD_DIR=build-debug CMAKE_BUILD_TYPE=Debug ./scripts/build.sh

# Pass OpenMPI oversubscribe
MPIEXEC_PREFLAGS=--oversubscribe ./scripts/build.sh

# Extra CMake toggles
EXTRA_CMAKE_ARGS="-DBUILD_EXAMPLES=ON -DBUILD_GUI=OFF" ./scripts/build.sh
```

---

## `build_docs.sh` — docs builder/previewer

**What it does**  
- Builds only the documentation: Doxygen XML → Sphinx HTML (no full CMake build).  
- Uses fixed paths under `build-docs/` to match `docs/conf.py`.  
- Can start a local HTTP server for quick preview.

**Flags**  
- `--clean` — remove previous docs build.  
- `--force` — clean + rebuild (CI‑friendly).  
- `--no-doxygen` — skip Doxygen step, reuse existing XML.  
- `--serve [--port N]` — serve HTML locally (default port 8000).  
- `--open` — try to open a browser tab after `--serve`.  
- `--quiet` — less verbose output.

**Examples**  
```bash
# Build docs once
./scripts/build_docs.sh

# Force clean + rebuild (useful in CI)
./scripts/build_docs.sh --force

# Rebuild and preview locally
./scripts/build_docs.sh --serve --open --port 8001
```

---

## `clean_build.sh` — clean build trees

**What it does**  
- Deletes CMake/Ninja artifacts in selected build directories.  
- Skips `build-docs/`

**Flags**  
- `-a, --all` — clean all `build-*` directories found.  
- `-k, --keep-deps` — preserve `_deps/` (downloaded third‑party sources) inside each build dir.  
- `-n, --dry-run` — show what would be removed.  
- `-y, --yes` — do not prompt.  
- `-v, --verbose` — more logging.  
- `-h, --help` — usage.

**Examples**  
```bash
# Clean a specific dir
./scripts/clean_build.sh build-debug

# Preview what would be deleted
./scripts/clean_build.sh --all -n

# Clean multiple builds but keep vendor downloads
./scripts/clean_build.sh --keep-deps build-regression build-perf -y
```

---

## `clean_extern.sh` — wipe vendored sources

**What it does**  
- Removes everything under `extern/` **except** `extern/README.md`.  
- Supports dry‑run and prompt‑less modes.

**Flags**  
- `-y, --yes` — do not prompt.  
- `-n, --dry-run` — list what would be removed.

**Examples**  
```bash
# See what would be deleted
./scripts/clean_extern.sh -n

# Nuke vendored archives/sources (keeps README)
./scripts/clean_extern.sh -y
```

---

## `format_all.sh` — code formatting

**What it does**  
- Runs the available formatters over project sources, skipping vendor/build trees:  
  - C/C++ via `clang-format`  
  - Fortran via `fprettify`  
  - CMake via `cmake-format`  
- Ignores missing tools unless `STRICT=1`.

**Env**  
- `STRICT=1` — fail the script if any formatter is missing or errors.

**Examples**  
```bash
# Best‑effort formatting
./scripts/format_all.sh

# CI‑style strict mode
STRICT=1 ./scripts/format_all.sh
```

---

## `mpi_env.sh` — sourceable MPI helpers

**What it does**  
- Detects MPI vendor/launcher and exports CMake/CTest‑friendly variables:  
  `MPIEXEC_EXECUTABLE`, `MPIEXEC_NUMPROC_FLAG`, `MPIEXEC_PREFLAGS`, `MPIEXEC_POSTFLAGS`.  
- Provides a convenience function `mpi_exec <np> <cmd …>` that calls the right launcher (`srun`, `mpirun`, `mpiexec`).  
- Modes:  
  - `auto` — detect cluster vs laptop, choose safe defaults.  
  - `emulate` — force laptop‑friendly settings (loopback/TCP/oversubscribe for Open MPI).  
  - `cluster` — disable emulation and favor scheduler/HCAs.

**Usage**  
```bash
# Source it (MUST be sourced)
source scripts/mpi_env.sh auto

# Then launch binaries consistently
mpi_exec 4 ./build/bin/your_mpi_program
```

---

## `prefetch_third_party.sh` — offline vendor cache

**What it does**  
- Configures CMake in a temporary dir with `-DPREFETCH_THIRD_PARTY=ON`, letting the project's CMake download declared third‑party sources.  
- Packages each fetched source as a reproducible archive into `extern/`, and writes `MANIFEST.prefetch` and `SHA256SUMS`.  
- Speeds up CI and enables fully offline builds when paired with `OFFLINE=1` in `build.sh`.

**Example**  
```bash
./scripts/prefetch_third_party.sh
# Archives now in extern/, plus MANIFEST.prefetch and SHA256SUMS
```

---

## Test runners

### `push_ci.sh`

**What it does**
- Pushes changes to GitHub  
- Configures tmux terminal with the actions user and actions-runner/run.sh  
- Waits on and runs CI config /.workflow/linux.yml  

**Example**  
```bash
./scripts/push_ci.sh
```

### `run_unit_tests.sh`
- **What it does:** Optionally builds (via `build.sh`) then runs `ctest -L unit`. Creates JUnit XML under `$BUILD_DIR/test-reports/unit/` (fallback created if CTest/Catch2 doesn’t emit one).  
- **Env:**  
  - `BUILD_DIR` (default: `build-unit`)  
  - `SKIP_BUILD` (`0|1`, default `0`)  
  - `CTEST_PARALLEL_LEVEL` (defaults to CPU count)  
  - `CTEST_TIMEOUT` (default: `900`)  
  - `REPORT_DIR` (default: `$BUILD_DIR/test-reports/unit`)  
- **Example:**  
  ```bash
  ./scripts/run_unit_tests.sh
  SKIP_BUILD=1 ./scripts/run_unit_tests.sh
  ```

### `run_perf_tests.sh`
- **What it does:** Same shape as unit, but runs `ctest -L perf` and writes to `$BUILD_DIR/test-reports/perf/`.  
- **Env:** `BUILD_DIR` (default `build-perf`), `SKIP_BUILD`, `CTEST_PARALLEL_LEVEL`, `CTEST_TIMEOUT`, `REPORT_DIR`.  
- **Example:**  
  ```bash
  CMAKE_BUILD_TYPE=Release ./scripts/run_perf_tests.sh
  ```

### `run_regression_tests.sh`
- **What it does:** Sources `mpi_env.sh`, configures with MPI and runs cases with the solver. Emits JUnit to `$BUILD_DIR/test-reports/regression/` (with a safe fallback if no launcher is found).  
- **Examples:**  
  ```bash
  ./scripts/run_regression_tests.sh
  BUILD_DIR=build-regression ./scripts/run_regression_tests.sh
  ```

---

## Tips

- Use `BUILD_DIR=…` consistently to keep separate trees for `unit`, `perf`, and `regression`.  
- `clean_build.sh --keep-deps` is handy when you don’t want to re‑download third‑party code.  
- Pair `prefetch_third_party.sh` with `OFFLINE=1 ./scripts/build.sh` for air‑gapped or flaky‑network environments.
