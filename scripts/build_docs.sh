#!/usr/bin/env bash
set -euo pipefail

# -------------------------------------------------
# build-docs.sh — lightweight docs builder/previewer
# -------------------------------------------------
# Features
#  • Build *only* the documentation (Doxygen XML + Sphinx HTML)
#  • Avoid configuring or building the full CMake project
#  • Usable locally and in CI/GitHub Actions
#  • Minimal rebuilds where possible
#
# Usage
#  ./scripts/build-docs.sh [--clean] [--force] [--no-doxygen] [--serve [--port N]] [--open]
#  CI example: ./scripts/build-docs.sh --force

here_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Find repo root by walking up
repo_root="$here_dir"
while [[ ! -f "$repo_root/CMakeLists.txt" && "$repo_root" != "/" ]]; do
  repo_root="$(dirname "$repo_root")"
done
[[ -f "$repo_root/CMakeLists.txt" ]] || { echo "❌  CMakeLists.txt not found above $here_dir"; exit 1; }

# Fixed build locations to match conf.py
build_root="$repo_root/build-docs"
docs_src="$repo_root/docs"
docs_build_dir="$build_root/docs"
html_dir="$docs_build_dir/html"
doctrees_dir="$docs_build_dir/doctrees"
doxygen_out_dir="$docs_build_dir/_build/doxygen"
doxygen_xml_dir="$doxygen_out_dir/xml"

# Inputs
force=false
clean=false
serve=false
open_browser=false
port=8000
run_doxygen=true
quiet=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    --force)        force=true; shift;;
    --clean)        clean=true; shift;;
    --serve)        serve=true; shift;;
    --port)         port="$2"; shift 2;;
    --open)         open_browser=true; shift;;
    --no-doxygen)   run_doxygen=false; shift;;
    -q|--quiet)     quiet=true; shift;;
    -h|--help)
      sed -n '1,80p' "$0" | sed -n '1,40p'
      exit 0;;
    *) echo "Unknown option: $1"; exit 2;;
  esac
done

log() { $quiet || echo -e "$*"; }

# Basic checks
[[ -d "$docs_src" ]] || { echo "❌  docs/ not found at $docs_src"; exit 1; }
command -v sphinx-build >/dev/null 2>&1 || { echo "❌  sphinx-build not found. Try: pip install sphinx breathe exhale"; exit 1; }
if $run_doxygen; then
  command -v doxygen >/dev/null 2>&1 || { echo "❌  doxygen not found. Install it (e.g., apt-get install doxygen)"; exit 1; }
fi

# Clean
if $clean; then
  log "🧹  Cleaning $docs_build_dir"
  rm -rf "$docs_build_dir"
fi

mkdir -p "$docs_build_dir" "$html_dir" "$doctrees_dir" "$doxygen_out_dir"

# Locate Doxyfile template (prefer docs/Doxyfile.in, fallback to repo root)
if [[ -f "$docs_src/Doxyfile.in" ]]; then
  doxy_in="$docs_src/Doxyfile.in"
elif [[ -f "$repo_root/Doxyfile.in" ]]; then
  doxy_in="$repo_root/Doxyfile.in"
else
  echo "❌  Doxyfile.in not found (looked in docs/ and repo root)"; exit 1
fi

doxyfile="$docs_build_dir/Doxyfile"

# Decide if Doxygen needs to run
need_doxy=false
stamp="$doxygen_out_dir/.stamp"
if $force; then
  need_doxy=true
elif ! $run_doxygen; then
  need_doxy=false
elif [[ ! -d "$doxygen_xml_dir" ]]; then
  need_doxy=true
elif [[ ! -f "$stamp" ]]; then
  need_doxy=true
else
  # Rebuild if any source/docs files are newer than the last run
  if find "$repo_root/src" -type f \( -name '*.hpp' -o -name '*.h' -o -name '*.cpp' -o -name '*.f90' \) -newer "$stamp" -print -quit | grep -q .; then
    need_doxy=true
  elif find "$doxy_in" -newer "$stamp" -print -quit | grep -q .; then
    need_doxy=true
  else
    need_doxy=false
  fi
fi

# Configure Doxyfile without invoking CMake (keep it in docs_build_dir)
if $run_doxygen && $need_doxy; then
  log "📝  Generating Doxyfile → $doxyfile"
  cp "$doxy_in" "$doxyfile"
  # Normalize paths (no trailing slashes)
  src_dir="$repo_root/src"
  sed -E -i \
    -e "s#^INPUT[[:space:]]*=.*#INPUT = ${src_dir//#/\\#}#" \
    -e "s#^STRIP_FROM_PATH[[:space:]]*=.*#STRIP_FROM_PATH = ${src_dir//#/\\#}#" \
    -e "s#^OUTPUT_DIRECTORY[[:space:]]*=.*#OUTPUT_DIRECTORY = ${doxygen_out_dir//#/\\#}#" \
    "$doxyfile"

  log "🔧  Running Doxygen (this scans headers only; no compile)"
  (cd "$docs_build_dir" && doxygen "$doxyfile")
  date +%s > "$stamp"
else
  log "⏭️  Skipping Doxygen (no changes detected)"
fi

# Sphinx build
log "📚  Building Sphinx HTML → $html_dir"
strict_flags=()
[[ "${SPHINX_STRICT:-0}" == "1" ]] && strict_flags=(-W)
# Parallel: let Sphinx decide (auto)
if $quiet; then
  sphinx-build -q -j auto -b html -d "$doctrees_dir" "$docs_src" "$html_dir" "${strict_flags[@]}"
else
  sphinx-build -j auto -b html -d "$doctrees_dir" "$docs_src" "$html_dir" "${strict_flags[@]}"
fi

log "✅  Docs ready: $html_dir"

# Preview server (optional)
if $serve; then
  command -v python3 >/dev/null 2>&1 || { echo "❌  python3 required for --serve"; exit 1; }
  log "🌐  Serving docs at http://localhost:$port (Ctrl+C to stop)"
  if $open_browser; then
    # Try to open a browser without failing the script if not available
    (xdg-open "http://localhost:$port" 2>/dev/null || open "http://localhost:$port" 2>/dev/null || true) &
  fi
  cd "$html_dir"
  exec python3 -m http.server "$port"
fi
