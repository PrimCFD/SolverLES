# `extern/` Cached tarballs for offline & reproducible builds

This directory holds **prefetched source archives** for third‑party dependencies. When an archive is present here, the build consumes it; when absent, CMake falls back to upstream URLs.

**Relevant pieces in this repo**

* `scripts/prefetch_third_party.sh` — populates and packs third‑party sources into `extern/`.
* `scripts/clean_extern.sh` — cleans `extern/` except this file.
* `scripts/clean_build.sh` — cleans `build-*/` can keep `deps_/` in targeted or all `build-*/` directories.
* `scripts/build.sh` — configures and builds, preferring archives in `extern/` and supporting a fully offline mode.
* `cmake/PrefetchDependencies.cmake` + `cmake/Fetch*.cmake` — dependency logic (HDF5, CGNS, OpenBLAS, PETSc, Catch2, yaml-cpp).

---

## Quick start (from repo root)

**1) Prefetch once (online) → create/update cache**

```bash
scripts/prefetch_third_party.sh
# ⇒ writes extern/*-src.tgz and updates extern/MANIFEST.prefetch + extern/SHA256SUMS
```

**2) Build using the cache (works offline)**

```bash
OFFLINE=1 scripts/build.sh
```

**3) Clean third‑party state**

```bash
scripts/clean_build.sh  # removes staged sources/installs; next build reuses extern/*
```

---

## How it works (dependency build logic)

1. **Prefetch**

   * The script configures CMake with `-DPREFETCH_THIRD_PARTY=ON` (via `cmake/PrefetchDependencies.cmake`).
   * Each `cmake/Fetch<Name>.cmake` runs a **populate‑only** step and registers a `<name>_src` target.
   * Populated source trees under `.prefetch/_deps/<name>-src/` are packed as `extern/<name>-src.tgz`.
   * `extern/MANIFEST.prefetch` lists what was captured; `extern/SHA256SUMS` records integrity hashes.

2. **Build (cache‑first)**

   * `scripts/build.sh` scans `extern/` for archives whose filenames contain one of:
     **Catch2**, **CGNS**, **HDF5**, **OpenBLAS**, **PETSc**, **yaml-cpp**; formats supported: `.tar.gz`, `.tgz`, `.tar.xz`, `.tar.bz2`, `.zip`.
   * For each hit, the archive is unpacked to `build/_deps/<name>-src/` and the corresponding
     `FETCHCONTENT_SOURCE_DIR_<NAME>` is set so CMake uses the **local sources** instead of downloading.
   * If **all** required deps are satisfied from `extern/` (or export `OFFLINE=1`), the build passes
     `-DFETCHCONTENT_FULLY_DISCONNECTED=ON` (and `...UPDATES_DISCONNECTED=ON`) to **forbid network access**.

3. **Staged installs for config‑packages**

   * Some deps are **built & installed into a local prefix** under `build/_deps/*-install/` during configure
     so dependents can `find_package(...)` without touching the system:

     * `HDF5::hdf5` with `HDF5_DIR` set to `build/_deps/hdf5-install/cmake`.
     * `CGNS::cgns` with `CGNS_DIR` set to `build/_deps/cgns-install/lib/cmake/CGNS`.
     * `OpenBLAS::OpenBLAS` and convenience `BLAS::BLAS`, `LAPACK::LAPACK` interface targets.
     * `PETSc` staged under `build/_deps/petsc-install` with `PETSC_DIR` pointing there.
   * These prefixes are also appended to `CMAKE_PREFIX_PATH` so regular `find_package(Config)` works.

---

## Typical workflows

**Fresh clone → one online prefetch → offline build forever**

```bash
scripts/prefetch_third_party.sh
OFFLINE=1 scripts/build.sh
```

**Air‑gapped or CI without internet**

```bash
# On a connected host
scripts/prefetch_third_party.sh
( cd extern && tar -cf ../extern-cache.tar . )

# Transfer extern-cache.tar alongside the repo on the offline host, then:
mkdir -p extern && tar -xf extern-cache.tar -C extern
OFFLINE=1 scripts/build.sh
```

**Bump a dependency**

```bash
# Adjust the version/pin in the appropriate cmake/Fetch*.cmake
scripts/prefetch_third_party.sh
scripts/clean_build.sh --build
OFFLINE=1 scripts/build.sh
```

---

## Adding or replacing cached archives manually

* Drop an archive in `extern/` whose **filename contains the package token** (`hdf5`, `cgns`, `openblas`, `petsc`, `catch2`, `yaml-cpp`).
* Supported formats: `.tar.gz`, `.tgz`, `.tar.xz`, `.tar.bz2`, `.zip`.
* Prefer using `scripts/prefetch_third_party.sh` so `MANIFEST.prefetch` and `SHA256SUMS` stay in sync.

**Example layout**

```
extern/
  hdf5-src.tgz
  cgns-src.tgz
  openblas-src.tgz
  petsc-src.tgz
  catch2-src.tgz
  yaml-cpp-src.tgz
  MANIFEST.prefetch
  SHA256SUMS
```

---

## Troubleshooting

* **CMake still tries to download**

  1. Verify the expected archives exist in `extern/` and match the package token.
  2. Clear stale state: `scripts/clean_build.sh`.
  3. Re‑run with `OFFLINE=1` to hard‑fail on any attempted network access.

* **HDF5/CGNS/OpenBLAS/PETSc/yaml-cpp not found at link time**

  * Remove the affected staged directories under `build/_deps/*-install/` and rebuild so package configs are regenerated.

* **Multiple archives for the same package**

  * Keep **one** archive per package in `extern/` to avoid ambiguity during the scan.

---

## Directory map

```
extern/
  <name>-src.tgz            # cached source snapshots created by the prefetch script
  MANIFEST.prefetch         # human‑readable list of cached items
  SHA256SUMS                # integrity checksums
build/_deps/
  <name>-src/               # unpacked sources used for the build
  hdf5-install/             # local install prefixes for config‑packages (examples)
  cgns-install/
  openblas-install/
  petsc-install/
```

---

## Related environment knobs (used by `scripts/build.sh`)

* `OFFLINE=1` — force `FetchContent` into fully disconnected mode.
* `BUILD_DIR` — custom build directory (default: `build`).
* `CMAKE_BUILD_TYPE` — `Debug`/`Release`/`RelWithDebInfo`/`MinSizeRel` (default: `Release`).
* `EXTRA_CMAKE_ARGS` — extra flags forwarded to CMake (space‑separated).
* MPI & CUDA toggles (`ENABLE_MPI`, `MPIEXEC_*`, `ENABLE_CUDA`, `USE_CUDA_UM`) are honored if relevant to the build.
