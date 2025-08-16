#!/usr/bin/env bash
set -euo pipefail

# Remove everything in extern/ except extern/README.md
# Options: -y/--yes (no prompt), -n/--dry-run (show what would be deleted)

YES=0
DRY=0

usage() {
  cat <<'EOF'
Usage: scripts/clean-extern.sh [-y|--yes] [-n|--dry-run]

Deletes all files and directories under extern/, except extern/README.md.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -y|--yes)    YES=1 ;;
    -n|--dry-run) DRY=1 ;;
    -h|--help)   usage; exit 0 ;;
    *) echo "Unknown option: $1"; usage; exit 2 ;;
  esac
  shift
done

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"
# Prefer git to find repo root; otherwise assume scripts/ is one level below root
repo_root="$(git -C "$script_dir" rev-parse --show-toplevel 2>/dev/null || echo "$script_dir/..")"
extern_dir="${repo_root}/extern"

# Safety checks
if [[ ! -d "$extern_dir" ]]; then
  echo "❌ extern directory not found at: $extern_dir"
  exit 1
fi

# Build deletion list (depth-first so directories get removed after their content)
# Keep ONLY the top-level extern/README.md
mapfile -d '' TO_DELETE < <(
  find "$extern_dir" -mindepth 1 -depth ! -path "$extern_dir/README.md" -print0
)

if (( ${#TO_DELETE[@]} == 0 )); then
  echo "✅ Nothing to delete. extern/ is already clean (README.md preserved)."
  exit 0
fi

echo "Will remove ${#TO_DELETE[@]} item(s) under: $extern_dir"
if (( DRY )); then
  printf '  %s\n' "${TO_DELETE[@]}"
  echo "🔎 Dry run only — nothing deleted."
  exit 0
fi

if (( ! YES )); then
  read -r -p "Proceed with deletion (preserving extern/README.md)? [y/N] " ans
  case "${ans,,}" in
    y|yes) ;;
    *) echo "Aborted."; exit 0 ;;
  esac
fi

# Delete
printf '%s\0' "${TO_DELETE[@]}" | xargs -0 rm -rf --
echo "🧹 Done. Preserved: $extern_dir/README.md"
