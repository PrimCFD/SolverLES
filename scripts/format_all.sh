#!/usr/bin/env bash
# Format C/C++, Fortran, and CMake files in THIS repo (skip vendor/build dirs)
set -euo pipefail

cd "$(dirname "$0")/.."

STRICT=${STRICT:-0}

# Patterns to prune (match subtree with the trailing *)
PRUNE_GLOBS=(
  './.git*' './.venv*' './extern*' './.prefetch*'
  './build*' './_deps*' './third_party*' './vendor*'
)

# Build a safe prune expression only if we have patterns
PRUNE_EXPR=()
if ((${#PRUNE_GLOBS[@]} > 0)); then
  PRUNE_EXPR=( '(' )
  first=1
  for pat in "${PRUNE_GLOBS[@]}"; do
    if (( first )); then
      PRUNE_EXPR+=( -path "$pat" )
      first=0
    else
      PRUNE_EXPR+=( -o -path "$pat" )
    fi
  done
  PRUNE_EXPR+=( ')' -prune -o )
fi

# -------- C/C++ --------
if command -v clang-format >/dev/null; then
  echo "🔧 Formatting C/C++ with clang-format…"
  find . "${PRUNE_EXPR[@]}" -type f \
       \( -name '*.c' -o -name '*.cc' -o -name '*.cpp' -o -name '*.h' -o -name '*.hpp' \) \
       -exec clang-format -i {} +
else
  echo "⚠️  clang-format not found"; [[ "$STRICT" == "1" ]] && exit 1
fi

# -------- Fortran --------
if command -v fprettify >/dev/null; then
  echo "🧮 Formatting Fortran with fprettify…"
  find . "${PRUNE_EXPR[@]}" -type f -name '*.f90' \
       -exec fprettify -c .fprettify.yml {} \;
else
  echo "⚠️  fprettify not found"; [[ "$STRICT" == "1" ]] && exit 1
fi

# -------- CMake --------
if command -v cmake-format >/dev/null; then
  echo "🛠️  Formatting CMake with cmake-format…"
  # Only project CMake (skip vendor/build trees). Format file-by-file.
  while IFS= read -r -d '' f; do
    if ! cmake-format -i "$f"; then
      echo "⚠️  cmake-format failed on $f — skipping"
      [[ "$STRICT" == "1" ]] && exit 1
    fi
  done < <(find . "${PRUNE_EXPR[@]}" -type f \( -name 'CMakeLists.txt' -o -name '*.cmake' \) -print0)
else
  echo "⚠️  cmake-format not found"; [[ "$STRICT" == "1" ]] && exit 1
fi

echo "✅ Formatting done"
