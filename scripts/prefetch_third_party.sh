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
