# extern/ — cached tarballs for offline & reproducible builds

This directory holds **prefetched source archives** for third‑party deps. When present,
the build uses them; when absent, it falls back to upstream URLs.

**Scripts this README matches:** `scripts/prefetch_third_party.sh`, `scripts/build.sh`,
`scripts/clean_extern.sh`, `scripts/clean_build.sh`.

---

## Quick commands (from repo root)

**Prefetch once (online) → create/update cache**
```bash
scripts/prefetch_third_party.sh
# ⇒ writes extern/*-src.tgz and updates extern/MANIFEST.prefetch + extern/SHA256SUMS
```

**Build fully offline using the cache**
```bash
OFFLINE=1 scripts/build.sh
```

**Run * tests (builds if needed)**
```bash
scripts/run_*_tests.sh
```

**Clean build dir but keep downloaded sources (_deps)**
```bash
scripts/clean_build.sh --keep-deps
```

**Wipe extern/ cache (keep this README)**
```bash
scripts/clean_extern.sh
```

---

## What happens under the hood

- **Prefetch:** `scripts/prefetch_third_party.sh` configures CMake with
  `-DPREFETCH_THIRD_PARTY=ON`, populates `_deps/<name>-src`, then packs each source tree
  into `extern/<name>-src.tgz`. It also writes `MANIFEST.prefetch` and `SHA256SUMS`.
- **Build:** `scripts/build.sh` scans `extern/` for archives matching `{tar.gz,tgz,tar.xz,tar.bz2,zip}`
  that contain **Catch2, CGNS, HDF5, PETSc** in the filename. For each match, it extracts
  into `build/_deps/<name>-src` and sets the corresponding
  `FETCHCONTENT_SOURCE_DIR_<NAME>` so CMake uses **local sources**.
- **Offline toggle:** If **all** deps are satisfied from `extern/` (or if `OFFLINE=1`), the
  build sets `-DFETCHCONTENT_FULLY_DISCONNECTED=ON` to prevent any network access.
- **Config‑package deps:** Some deps (e.g. **HDF5**) are locally **installed into a small
  staging prefix** inside `build/_deps/…-install/` during configure so dependents like
  **CGNS** can `find_package` them without touching the system.

---

## Tips & troubleshooting

- **Add/replace archives:** Any filename containing the package name works (e.g. `hdf5-foo.tgz`);
  supported formats: `.tar.gz`, `.tgz`, `.tar.xz`, `.tar.bz2`, `.zip`.
- **Force offline:** `OFFLINE=1 scripts/build.sh`.
- **CMake tries to download:** ensure the expected archives exist in `extern/` **and**
  remove stale build cache: `scripts/clean_build.sh --keep-deps` (or `rm -rf build/_deps`),
  then rebuild offline.
- **HDF5/CGNS not found:** clean `_deps` and reconfigure so the staged HDF5 package dir
  is recreated and passed to CGNS.

---

## Directory map

```
extern/
  <name>-src.tgz            # cached source snapshots (created by prefetch script)
  MANIFEST.prefetch         # human‑readable list
  SHA256SUMS                # checksums
build/_deps/
  <name>-src/               # unpacked sources used for the build
  hdf5-install/             # local install prefix for config packages (example)
```