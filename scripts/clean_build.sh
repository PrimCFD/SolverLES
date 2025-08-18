#!/usr/bin/env bash
# Thorough CMake/Ninja build cleaner for this repo
# Works on build-*/ dirs except build-docs/ (docs-only)
# Usage:
#   ./clean.sh [-a|--all] [-k|--keep-deps] [-n|--dry-run] [-y|--yes] [-v|--verbose] [BUILD_DIR ...]
#
# If no BUILD_DIR is provided and --all is not set, falls back to $BUILD_DIR if present.
# Examples:
#   ./clean.sh build-debug
#   ./clean.sh --keep-deps build-mpi build-perf
#   ./clean.sh --all -n     # preview what would be removed

set -Eeuo pipefail

usage() {
  sed -n '2,30p' "$0" | sed 's/^# \{0,1\}//'
}

# ---- options ----
DRY_RUN=0
KEEP_DEPS=0
ALL=0
YES=0
VERBOSE=0

log() { [[ $VERBOSE -eq 1 ]] && printf '%s\n' "$*"; }
say() { printf '%s\n' "$*"; }
warn() { printf '⚠️  %s\n' "$*" >&2; }
die() { printf '❌ %s\n' "$*" >&2; exit 1; }

rmrf() {
  if [[ $DRY_RUN -eq 1 ]]; then
    printf '[dry-run] rm -rf %q\n' "$@"
  else
    rm -rf -- "$@"
  fi
}

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"
if git rev-parse --show-toplevel &>/dev/null; then
  repo_root="$(git rev-parse --show-toplevel)"
else
  repo_root="$script_dir"
fi

# Consider a dir "safe" to clean if it looks like a CMake/Ninja build tree
is_safe_build_dir() {
  local d="$1"
  [[ -d "$d" ]] || return 1
  [[ -f "$d/CMakeCache.txt" ]] && return 0
  [[ -f "$d/build.ninja" ]] && return 0
  [[ -f "$d/Makefile" ]] && return 0
  [[ -d "$d/CMakeFiles" ]] && return 0
  return 1
}

discover_build_dirs() {
  # Search a couple of common places/patterns (depth-limited)
  find "$repo_root" \
    -maxdepth 2 -type d \( \
        -name 'build' -o \
        -name 'build-*' -o \
        -name 'cmake-build-*' \
      \) 2>/dev/null | while read -r d; do
        if is_safe_build_dir "$d"; then printf '%s\n' "$d"; fi
      done
}

clean_one_dir() {
  local dir="$1"
  [[ -d "$dir" ]] || { warn "Skip: $dir (not a directory)"; return; }

  # Safety guardrails
  local abs
  abs="$(cd "$dir" && pwd)"
  [[ -n "$abs" ]] || die "Cannot resolve path for $dir"
  [[ "$abs" != "/" ]] || die "Refusing to clean '/'"
  is_safe_build_dir "$abs" || { warn "Skip: $abs (doesn't look like a build dir)"; return; }

  say "🧹 Cleaning: $abs"
  if [[ $KEEP_DEPS -eq 1 ]]; then
    # Keep build/_deps (downloaded third-party sources), remove everything else.
    # Also prune common heavy subtrees inside the build while preserving _deps.
    mapfile -t top_level < <(find "$abs" -mindepth 1 -maxdepth 1 -printf '%f\n')
    for entry in "${top_level[@]}"; do
      if [[ "$entry" == "_deps" ]]; then
        log "keep: $abs/_deps/"
        continue
      fi
      rmrf "$abs/$entry"
    done
    # Remove typical CMake/Ninja files that might be hidden or recreated next configure
    rmrf "$abs"/{CMakeCache.txt,CTestTestfile.cmake,cmake_install.cmake,.ninja_log,.ninja_deps} 2>/dev/null || true
  else
    # Nuke the entire build directory
    rmrf "$abs"
    return
  fi

  # Extra pass: deep-clean well-known subtrees if they exist but slipped through
  for p in CMakeFiles Testing test-reports bin lib lib64 include share \
           obj .objs .cache .cmake-api ; do
    [[ -e "$abs/$p" ]] && rmrf "$abs/$p"
  done

  # Yank common build system files regardless of depth (keeps _deps intact)
  find "$abs" -path "$abs/_deps" -prune -o -type f \( \
      -name 'CMakeCache.txt' -o -name 'cmake_install.cmake' -o -name 'CTestTestfile.cmake' -o \
      -name 'build.ninja' -o -name '.ninja_log' -o -name '.ninja_deps' -o -name 'Makefile' -o \
      -name 'compile_commands.json' \
    \) -print0 2>/dev/null | xargs -0 -r ${DRY_RUN:+echo} rm -f

  say "✅ Done: $abs"
}

# ---- parse args ----
TARGETS=()
if [[ $# -eq 0 ]]; then
  if [[ -n "${BUILD_DIR:-}" ]]; then
    TARGETS+=("$BUILD_DIR")
  else
    ALL=1
  fi
fi

while [[ $# -gt 0 ]]; do
  case "$1" in
    -a|--all) ALL=1 ;;
    -k|--keep-deps) KEEP_DEPS=1 ;;
    -n|--dry-run) DRY_RUN=1 ;;
    -y|--yes) YES=1 ;;
    -v|--verbose) VERBOSE=1 ;;
    -h|--help) usage; exit 0 ;;
    *) TARGETS+=("$1") ;;
  esac
  shift
done

if [[ $ALL -eq 1 ]]; then
  mapfile -t found < <(discover_build_dirs)
  TARGETS+=("${found[@]}")
fi

# De-dupe targets
if [[ ${#TARGETS[@]} -eq 0 ]]; then
  warn "No build directories found."
  exit 0
fi
mapfile -t TARGETS < <(printf '%s\n' "${TARGETS[@]}" | awk 'NF' | sort -u)

say "Targets:"
for t in "${TARGETS[@]}"; do say "  - $t"; done

if [[ $YES -ne 1 && $DRY_RUN -ne 1 ]]; then
  suffix=""
  if [[ $KEEP_DEPS -eq 1 ]]; then suffix=" (keeping _deps)"; fi
  read -r -p "Proceed with deletion${suffix}? [y/N] " reply
  [[ "${reply,,}" == "y" || "${reply,,}" == "yes" ]] || { say "Aborted."; exit 1; }
fi

for t in "${TARGETS[@]}"; do clean_one_dir "$t"; done
