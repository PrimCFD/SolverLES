# `tests/` — Workflow, conventions, and how to add tests

This README explains how our tests are organized, how to run them locally and in CI, and exactly how to add **unit**, **MPI**, and **performance** tests. It reflects the current layout and CMake wiring so you can be productive immediately.

---

## Directory layout (current)

```text
tests/
├─ unit/
│  └─ core/
│     ├─ CMakeLists.txt            # builds a single Catch2 binary (unit_core)
│     ├─ test_*.cpp                # regular unit tests (Catch2)
│     ├─ gpu/                      # optional CUDA-aware unit tests
│     │  └─ test_*.cpp
│     └─ mpi/
│        ├─ mpi_main.cpp           # Catch2 runner that init/finalizes MPI
│        └─ test_*.cpp             # MPI-aware test sources
└─ performance/
   └─ core/
      ├─ CMakeLists.txt            # builds small bench binaries
      ├─ bench_*.cpp               # perf micro-benchmarks
      └─ simple_bench.hpp          # helpers (timer, mean/stddev, report)
```

> More suites will be added as other subsystems appear (e.g. `tests/unit/gui`, `tests/regression`, etc).

---

## TL;DR

* **Run everything (CPU):**

  ```bash
  ./scripts/run_unit_tests.sh
  CMAKE_BUILD_TYPE=Release ./scripts/run_perf_tests.sh
  ./scripts/run_mpi_tests.sh --mode emulate --np 4
  ```
* **Labels:** CTest labels are `unit`, `perf`, and `mpi`.
* **JUnit reports:**

  * Unit → `build-unit/test-reports/unit/`
  * Perf → `build-perf/test-reports/perf/`
  * MPI  → `build-mpi/test-reports/mpi/`
* **CUDA/GPU:** Tests that require CUDA compile with `ENABLE_CUDA=ON` and will auto‑skip at runtime when unavailable.

---

## How things are wired (CMake + Catch2 + CTest)

### Unit tests (CPU + optional CUDA)

* `tests/unit/core/CMakeLists.txt` compiles a single binary **`unit_core`** from all `test_*.cpp` sources.
* We use a helper macro to auto‑register every `TEST_CASE` in that binary as a CTest case:

  ```cmake
  register_catch_tests(unit_core
    LABEL unit
    TEST_PREFIX "core::"
    REPORTER junit
    OUTPUT_DIR "${CMAKE_BINARY_DIR}/test-reports/unit"
    OUTPUT_PREFIX "catch-"
    OUTPUT_SUFFIX ".xml")
  ```
* Optional CUDA unit tests live under `tests/unit/core/gpu/`. When built with `-DENABLE_CUDA=ON`, we compile a separate `unit_gpu` and mark it with the same `unit` label; individual cases typically guard with `#ifdef HAVE_CUDA` and `SUCCEED("skipped")` when CUDA is off.

### MPI tests

* MPI tests use a dedicated Catch2 runner (`mpi/mpi_main.cpp`) which initializes/finalizes MPI once per process.
* Each MPI group is built as its **own binary** (e.g. `test_halo`, `test_boundary_mpi`, `test_io`). We then add concrete CTest invocations with an explicit number of ranks via an `add_mpi_test(name exe np)` helper.
* All MPI tests carry the `mpi` label and write JUnit XML into `test-reports/mpi/`.

### Performance tests

* Perf tests are tiny standalone **bench\_**\* executables (e.g. `bench_allocate`, `bench_copy`).
* Each bench is registered as a CTest test with the `perf` label via `add_perf_test(name exe [args…])`.
* Reports go to `test-reports/perf/`. These tests are designed to be short and repeatable; they print one line of stats (name, mean, stddev, payload size).

---

## Local workflows

### 1) Unit tests (CPU)

```bash
./scripts/run_unit_tests.sh
```

* Builds into `build-unit/` by default (override with `BUILD_DIR`).
* Parallel by default (`CTEST_PARALLEL_LEVEL` to override).

### 2) Performance tests (CPU)

```bash
CMAKE_BUILD_TYPE=Release ./scripts/run_perf_tests.sh
```

* Prefer `Release` for meaningful numbers.

### 3) MPI tests (laptop‑friendly)

```bash
./scripts/run_mpi_tests.sh --mode emulate --np 4
```

* `--mode emulate` picks safe Open MPI flags for localhost; on clusters use `--mode cluster` and your scheduler launcher.

### 4) GPU (optional)

* Build with `-DENABLE_CUDA=ON` and a CUDA toolchain in PATH.
* Perf: `CMAKE_BUILD_TYPE=Release ./scripts/run_perf_tests.sh`
* Some unit tests auto‑skip when `HAVE_CUDA` isn’t defined at compile time.

> All runners support `SKIP_BUILD=1` to reuse an existing build, and respect `OFFLINE=1` if you have `extern/` cached.

---

## How to add a **unit** test (CPU)

1. Create a new file under `tests/unit/core/`, named `test_<topic>.cpp`.

2. Write a Catch2 test:

```cpp
#include <catch2/catch_test_macros.hpp>

TEST_CASE("my feature works", "[category][core]") {
  // Arrange/Act/Assert
}
```

3. Add your file to the `unit_core` sources in `tests/unit/core/CMakeLists.txt`:

```cmake
add_executable(unit_core
  test_alignment.cpp
  # … existing tests …
  test_my_feature.cpp)  # ← add this line

target_link_libraries(unit_core PRIVATE core Catch2::Catch2WithMain)
```

4. Build & run:

```bash
./scripts/run_unit_tests.sh
# or: ctest -L unit --test-dir build-unit
```

**Tips**

* Use **tags** in the second argument to `TEST_CASE` (`[memory]`, `[io]`, `[mesh]`, etc.).
* Prefer deterministic inputs; avoid sleeping/timeouts inside unit tests.

---

## How to add an **MPI** test

1. Put your test source under `tests/unit/core/mpi/` (e.g. `test_exchange.cpp`). Use standard Catch2 macros in the file (no explicit `main`).

2. Reference the common MPI runner (`mpi_main.cpp`) and create a new test binary in `tests/unit/core/CMakeLists.txt`:

```cmake
add_executable(test_exchange mpi/mpi_main.cpp mpi/test_exchange.cpp)
target_link_libraries(test_exchange PRIVATE core Catch2::Catch2 MPI::MPI_CXX)
```

3. Register one or more invocations (with ranks) via the helper:

```cmake
add_mpi_test(core_exchange_2rank test_exchange 2)
add_mpi_test(core_exchange_4rank test_exchange 4)
```

4. Build & run locally:

```bash
./scripts/run_mpi_tests.sh --mode emulate --np 4
# or: ctest -L mpi --test-dir build-mpi
```

**Guidelines**

* Keep per‑rank output minimal; the runner does an `MPI_Barrier` before finalize to tidy logs.
* Assert expectations that hold for any rank topology you register (e.g. 2 and 4 ranks).

---

## How to add a **performance** test

1. Create a `bench_<topic>.cpp` in `tests/performance/core/` that uses `simple_bench.hpp`:

```cpp
#include "simple_bench.hpp"
int main(){
  auto [mean, stddev] = bench::run([&]{ /* code under test */ });
  bench::report("my_bench", mean, stddev, /*bytes=*/0);
  return 0;
}
```

2. Register it in `tests/performance/core/CMakeLists.txt`:

```cmake
add_executable(bench_my_stuff bench_my_stuff.cpp simple_bench.hpp)
target_link_libraries(bench_my_stuff PRIVATE core)
add_perf_test(perf_my_stuff bench_my_stuff)
```

3. Run:

```bash
CMAKE_BUILD_TYPE=Release ./scripts/run_perf_tests.sh
# or: ctest -L perf --test-dir build-perf
```

**GPU perf**

* If your bench needs CUDA, guard it with `#ifdef HAVE_CUDA`, link `CUDA::cudart`, and add `target_compile_definitions(bench_x PRIVATE HAVE_CUDA)` when `ENABLE_CUDA` is on. The test should print a clear "skipped" message when CUDA is absent.

---

## CI workflow (what runs on PRs)

* **Prefetch:** Third‑party sources are cached to `extern/` in a warm‑up job, then downloaded by test jobs.
* **CPU matrix:**

  * *Unit (CPU):* runs `./scripts/run_unit_tests.sh`.
  * *Perf (CPU):* runs `./scripts/run_perf_tests.sh` with `CMAKE_BUILD_TYPE=Release`.
  * *MPI (emulated):* installs Open MPI and runs `./scripts/run_mpi_tests.sh --mode emulate --np 4 --keep-going`.
* **GPU perf:** Runs perf suite inside `nvidia/cuda:12.4.1-devel-ubuntu22.04` on a self‑hosted runner with GPUs.
* All jobs upload JUnit XML from the `test-reports/*` folders as artifacts.

**Reproducing CI locally**

```bash
# Optional, for offline builds
./scripts/prefetch_third_party.sh && OFFLINE=1 ./scripts/build.sh

# Then run the same scripts as CI
./scripts/run_unit_tests.sh
CMAKE_BUILD_TYPE=Release ./scripts/run_perf_tests.sh
./scripts/run_mpi_tests.sh --mode emulate --np 4
```

---

## Conventions

* **Filenames:** `test_<topic>.cpp`, `bench_<topic>.cpp`
* **Namespaces/headers:** Prefer including public headers (`master/...`, `memory/...`, `mesh/...`).
* **Determinism:** Avoid timing‑based assertions in perf tests; keep iteration counts small.
* **Output:** Perf tests should print a single line with a stable key (used by dashboards).
* **Skips:** Use `#ifdef HAVE_CUDA` and `SUCCEED("CUDA not enabled; test skipped")` for GPU tests.
* **Labels:** Unit/MPI/Perf are separated via CTest labels; keep unit tests fast (<1s) and MPI tests modest in size.

---

## Troubleshooting

* **No tests found:** Re‑run CMake if you added new sources to a target but didn’t re‑configure.
* **MPI launcher errors:** Use `--mode emulate` locally; on clusters, source `scripts/mpi_env.sh cluster` and rerun.
* **Missing JUnit XML:** The registration macros generate them; ensure you didn’t bypass the helper in CMake.
* **GPU not detected:** Verify `-DENABLE_CUDA=ON`, `nvcc --version`, and that your test/bench defines `HAVE_CUDA`.

---

## Examples in the tree

* **Alignment (unit):** checks allocator alignment on various types.
* **Boundary MPI (mpi):** verifies `is_physical_face` behavior with/without MPI.
* **IO under CUDA (unit/gpu):** minimal XDMF/HDF5 write path with device residency.
* **Alloc perf (perf):** measure allocate/release variance and report mean/stddev.

These are good starting points when creating
