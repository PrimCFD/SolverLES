#!/usr/bin/env bash
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$here/.."
prefetch_dir="${repo_root}/.prefetch"
extern_dir="${repo_root}/extern"
manifest="${extern_dir}/MANIFEST.prefetch"
checksums="${extern_dir}/SHA256SUMS"

cleanup() { rm -rf "${prefetch_dir}"; }
trap cleanup EXIT

echo "▶️  Prefetching third-party sources (one-time, online)…"

cmake -S "${repo_root}" -B "${prefetch_dir}" \
      -DPREFETCH_THIRD_PARTY=ON \
      -DCMAKE_BUILD_TYPE=Release

# Populate _deps/*-src without building full libs
cmake --build "${prefetch_dir}" --target prefetch-archives --parallel

mkdir -p "${extern_dir}"
rm -f "${manifest}" "${checksums}"
echo "# Third-party source snapshots created on $(date -u +%Y-%m-%dT%H:%M:%SZ)" > "${manifest}"

echo "📦  Packing to extern/*.tgz with checksums…"
for d in "${prefetch_dir}"/_deps/*-src; do
  [[ -d "$d" ]] || continue
  base="$(basename "$d")"
  tgz="${extern_dir}/${base}.tgz"
  # Reproducible tar: stable sort & numeric owner
  tar --sort=name --owner=0 --group=0 --numeric-owner -czf "${tgz}" \
    -C "$(dirname "$d")" "$(basename "$d")"
  echo "${base}.tgz  ← ${d}" >> "${manifest}"
  if command -v sha256sum >/dev/null; then
    ( cd "${extern_dir}" && sha256sum "${base}.tgz" ) >> "${checksums}"
  else
    openssl dgst -sha256 "${tgz}" | awk '{print $2"  '"${base}.tgz"'"}' >> "${checksums}"
  fi
done

echo "🧾  Manifest:  ${manifest}"
echo "🔐  Checksums: ${checksums}"
echo "✅  Archives ready in ${extern_dir}/"

# Fetch PETSc third-party tarballs for OFFLINE configure
# We ask PETSc's configure to print required URLs, then we download into a cache dir.
petsc_src="${prefetch_dir}/_deps/petsc-src"
pkg_cache="${extern_dir}/.petsc-downloads"
urls_file="${pkg_cache}/urls.txt"
mkdir -p "${pkg_cache}"

if [[ -d "${petsc_src}" ]]; then
  echo "🌐  Collecting PETSc package URLs for offline builds…"
  pushd "${petsc_src}" >/dev/null
  # Ask configure to print only the URLs needed for our requested downloads.
  # (PETSc docs: --with-packages-download-dir makes configure list the URLs.) 
  # We request just MPICH + reference BLAS/LAPACK to guarantee an MPI build offline.
  set +e
  CFG_OUT="$(python3 ./configure \
      --with-packages-download-dir="${pkg_cache}" \
      --download-mpich \
      --with-debugging=0 \
      --prefix=/tmp/prefetch-petsc 2>&1)"
  set -e
  # Extract, sanitize, and keep only MPICH (we use our own OpenBLAS for BLAS/LAPACK).

  echo "${CFG_OUT}" \
    | grep -Eo 'https?://[^[:space:]]+' \
    | sed -E "s/[\"',)]+$//" \
    | grep -Ei '/mpich[^/]*\.(tar\.(gz|bz2|xz)|tgz)$' \
    | sort -u > "${urls_file}"

  echo "⬇️   Downloading PETSc package tarballs into ${pkg_cache}/ …"
  while read -r url; do
    [[ -z "${url}" ]] && continue
    # Prefer PETSc mirror if available (safe, stable host)
    case "${url}" in
      https://github.com/pmodels/mpich/releases/*)
        mirror="https://web.cels.anl.gov/projects/petsc/download/externalpackages/$(basename "${url}")"
        # If mirror exists, switch to it
        if curl -fsI "${mirror}" >/dev/null; then url="${mirror}"; fi
        ;;
    esac
    fname="${pkg_cache}/$(basename "${url}")"
    if [[ ! -f "${fname}" ]]; then
      echo "  - ${url}"
      curl -L --retry 3 --fail -o "${fname}.part" "${url}"
      mv "${fname}.part" "${fname}"
    fi
  done < "${urls_file}"
  popd >/dev/null

  echo "📚  Cached PETSc downloads:"
  ls -1 "${pkg_cache}" | sed 's/^/  • /'
  echo "✅  PETSc offline cache ready at ${pkg_cache}/"
else
  echo "⚠️  Skipping PETSc offline cache: sources not found at ${petsc_src}"
fi
