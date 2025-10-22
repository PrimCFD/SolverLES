#!/usr/bin/env bash
set -euo pipefail

# Optional: vendor-agnostic prefetch configuration
# PREFETCH_MPI controls whether to cache an MPI vendor for PETSc offline builds.
# Values: "mpich" to cache MPICH tarballs; empty/anything else = no MPI cached.
# You can also pass --mpi=mpich on the command line.
PREFETCH_MPI="${PREFETCH_MPI:-}"
# Enable PETSc/MUMPS download URL discovery (parallel LU stack)
PREFETCH_MUMPS="${PREFETCH_MUMPS:-1}"

for arg in "$@"; do
  case "$arg" in
    --mpi=mpich) PREFETCH_MPI="mpich" ;;
    --mpi=none|--mpi=system|--mpi="") PREFETCH_MPI="" ;;
    --mumps) PREFETCH_MUMPS="1" ;;
    -h|--help)
      cat <<EOF
Usage: $(basename "$0") [--mpi=mpich|none] [--mumps]
  --mpi=mpich  Cache MPICH tarballs for PETSc offline configure
  --mpi=none   Do not cache an MPI (default); PETSc will use system MPI at build time
  --mumps      Also cache MUMPS + ScaLAPACK/BLACS + ParMETIS/METIS (+PT-Scotch/Scotch)
Environment:
  PREFETCH_MPI=mpich   Same as --mpi=mpich
  PREFETCH_MUMPS=1     Same as --mumps
EOF
      exit 0;;
  esac
done

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
  # Optionally cache MPI vendor tarballs for offline PETSc configure.
  # When PREFETCH_MPI=mpich, we request MPICH; otherwise we let PETSc use system MPI later.
  if [[ "${PREFETCH_MPI}" == "mpich" || -n "${PREFETCH_MUMPS}" ]]; then
    echo "⚙️   PETSc configure (URL discovery)"
    CFG_ARGS=(--with-packages-download-dir="${pkg_cache}" --with-debugging=0 --prefix=/tmp/prefetch-petsc)
    [[ "${PREFETCH_MPI}" == "mpich" ]] && CFG_ARGS+=(--download-mpich)
    if [[ -n "${PREFETCH_MUMPS}" ]]; then
      # Ask PETSc to include parallel LU stack URLs:
      # MUMPS + ScaLAPACK/BLACS + ParMETIS/METIS (+ optional PT-Scotch/Scotch)
      # (PETSc knows and publishes the exact tarball URLs during configure.)  <!-- docs: install & MUMPS -->
      CFG_ARGS+=(--download-mumps --download-scalapack --download-blacs \
                 --download-parmetis --download-metis)
    fi
    echo "     configure args: ${CFG_ARGS[*]}"
    # Save full configure output for debugging and grep from file (more reliable than a shell var)
    CFG_LOG="${pkg_cache}/configure_url_discovery.log"
    python3 ./configure "${CFG_ARGS[@]}" > "${CFG_LOG}" 2>&1 || true
  else
    echo "⚙️   PETSc configure (URL discovery) without extra downloads"
    CFG_LOG="${pkg_cache}/configure_url_discovery.log"
    python3 ./configure --with-packages-download-dir="${pkg_cache}" --with-debugging=0 --prefix=/tmp/prefetch-petsc > "${CFG_LOG}" 2>&1 || true
  fi

  set -e

  # Extract candidate URLs from PETSc configure output:
  #  - Allow Bitbucket (pkg-metis / pkg-parmetis), GitHub (Reference-ScaLAPACK),
  #    INRIA GitLab (scotch), MUMPS site, and PETSc mirrors.
  #  - Do NOT require the package name to appear in the filename (ScaLAPACK uses commit-hash tarballs).
  : > "${urls_file}"
  grep -Eo 'https?://[^[:space:]'"'"'"]+\.(tar\.(gz|bz2|xz)|tgz)(\?[^[:space:]'"'"'"]*)?' "${CFG_LOG}" \
  | sed -E 's/[)",]+$//' \
  | grep -Ei '(bitbucket\.org/petsc/pkg-(metis|parmetis)|github\.com/Reference-ScaLAPACK|mumps-solver\.org|/projects/petsc/(mirror|download)/externalpackages/)' \
  | sort -u > "${urls_file}" || true

  echo "📝  Configure log saved: ${CFG_LOG}"

  set -e

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
