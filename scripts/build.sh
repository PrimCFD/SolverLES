#!/usr/bin/env bash
set -euo pipefail

# -------------------------------------------------
# build.sh — centralized, env-driven CMake build
# -------------------------------------------------
# Env (all optional):
#   BUILD_DIR                build dir (relative to repo root unless absolute) [default: build]
#   CMAKE_BUILD_TYPE         Debug/Release/RelWithDebInfo/MinSizeRel [default: Release]
#   BUILD_TESTS              ON/OFF (configure tests) [default: ON]
#   ENABLE_MPI               ON/OFF (toggle MPI paths) [default: OFF]
#   MPIEXEC_PREFLAGS         forwarded to CMake (e.g., --oversubscribe for OpenMPI)
#   CMAKE_GENERATOR          override generator (e.g., "Ninja")
#   CMAKE_TOOLCHAIN_FILE     forwarded if set
#   EXTRA_CMAKE_ARGS         extra cmake args (space-separated)
#   OFFLINE                  if "1", force disconnected FetchContent
#
# Examples:
#   CMAKE_BUILD_TYPE=Debug ./build.sh
#   BUILD_DIR=build-mpi ENABLE_MPI=ON ./build.sh
#   EXTRA_CMAKE_ARGS="-DUSE_SYSTEM_HDF5=OFF" ./build.sh

# --- Locate repo root by walking up to nearest CMakeLists.txt
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$here"
while [[ ! -f "${repo_root}/CMakeLists.txt" && "$repo_root" != "/" ]]; do
  repo_root="$(dirname "$repo_root")"
done
[[ -f "${repo_root}/CMakeLists.txt" ]] || { echo "❌  CMakeLists.txt not found above $here"; exit 1; }

# --- Resolve build dir (abs/rel)
if [[ -n "${BUILD_DIR:-}" ]]; then
  if [[ "$BUILD_DIR" = /* ]]; then
    build_dir="$BUILD_DIR"
  else
    build_dir="${repo_root}/${BUILD_DIR}"
  fi
else
  build_dir="${repo_root}/build"
fi
mkdir -p "${build_dir}"

deps_dir="${build_dir}/_deps"
extern_dir="${repo_root}/extern"
mkdir -p "${deps_dir}" "${extern_dir}"

# --- CPU count (portable)
NPROCS=${NPROCS:-$(
  command -v nproc >/dev/null && nproc || \
  getconf _NPROCESSORS_ONLN 2>/dev/null || \
  sysctl -n hw.logicalcpu 2>/dev/null || \
  echo 2
)}

# --- Defaults (overridable)
CMAKE_BUILD_TYPE="${CMAKE_BUILD_TYPE:-Release}"
BUILD_TESTS="${BUILD_TESTS:-ON}"
ENABLE_MPI="${ENABLE_MPI:-OFF}"
ENABLE_CUDA="${ENABLE_CUDA:-OFF}"
USE_CUDA_UM="${USE_CUDA_UM:-OFF}" 

# --- Prefer Ninja unless overridden
if [[ -n "${CMAKE_GENERATOR:-}" ]]; then
  generator="${CMAKE_GENERATOR}"
elif command -v ninja >/dev/null 2>&1 || command -v ninja-build >/dev/null 2>&1; then
  generator="Ninja"
else
  generator="Unix Makefiles"
fi

# --- Core CMake args (keep original project knobs)
cmake_args=(
  -S "${repo_root}"
  -B "${build_dir}"
  -G "${generator}"
  -DCMAKE_BUILD_TYPE="${CMAKE_BUILD_TYPE}"
  -DFETCHCONTENT_BASE_DIR="${deps_dir}"
  -DUSE_SYSTEM_CGNS=OFF
  -DUSE_SYSTEM_HDF5=OFF
  -DUSE_SYSTEM_PETSC=OFF
  -DBUILD_TESTS="${BUILD_TESTS}"
  -DENABLE_MPI="${ENABLE_MPI}"
  -DENABLE_CUDA="${ENABLE_CUDA}"
  -DUSE_CUDA_UM="${USE_CUDA_UM}"
)
[[ -n "${MPIEXEC_NUMPROC_FLAG:-}" ]] && cmake_args+=( -DMPIEXEC_NUMPROC_FLAG="${MPIEXEC_NUMPROC_FLAG}" )
[[ -n "${MPIEXEC_PREFLAGS:-}"    ]] && cmake_args+=( -DMPIEXEC_PREFLAGS="${MPIEXEC_PREFLAGS}" )
[[ -n "${MPIEXEC_POSTFLAGS:-}"   ]] && cmake_args+=( -DMPIEXEC_POSTFLAGS="${MPIEXEC_POSTFLAGS}" )
[[ -n "${MPIEXEC_EXECUTABLE:-}"  ]] && cmake_args+=( -DMPIEXEC_EXECUTABLE="${MPIEXEC_EXECUTABLE}" )
[[ -n "${CMAKE_TOOLCHAIN_FILE:-}" ]] && cmake_args+=( -DCMAKE_TOOLCHAIN_FILE="${CMAKE_TOOLCHAIN_FILE}" )

# --- Scan extern/ for pre-fetched archives; wire up FETCHCONTENT_SOURCE_DIR_* like before
declare -A fetch_var=(
  [Catch2]=FETCHCONTENT_SOURCE_DIR_CATCH2
  [CGNS]=FETCHCONTENT_SOURCE_DIR_CGNS
  [HDF5]=FETCHCONTENT_SOURCE_DIR_HDF5
  [PETSc]=FETCHCONTENT_SOURCE_DIR_PETSC
)
packages=(Catch2 CGNS HDF5 PETSc)

echo "🔍  Scanning extern/ for pre-fetched archives…"
all_offline=true
for pkg in "${packages[@]}"; do
  lc_pkg=${pkg,,}
  shopt -s nullglob nocaseglob
  tarfiles=( "${extern_dir}"/*${pkg}*.{tar.gz,tar.xz,tar.bz2,tgz,zip} )
  shopt -u nullglob nocaseglob

  src_dir="${deps_dir}/${lc_pkg}-src"
  if ((${#tarfiles[@]})); then
    archive="${tarfiles[0]}"
    echo "   • $pkg ← $(basename "$archive")"
    if [[ ! -d "${src_dir}" || "$archive" -nt "${src_dir}" ]]; then
      rm -rf "${src_dir}"; mkdir -p "${src_dir}"
      case "$archive" in
        *.tar.gz|*.tgz) tar -xzf "$archive" --strip-components=1 -C "${src_dir}";;
        *.tar.xz)       tar -xJf "$archive" --strip-components=1 -C "${src_dir}";;
        *.tar.bz2)      tar -xjf "$archive" --strip-components=1 -C "${src_dir}";;
        *.zip)
          command -v unzip >/dev/null || { echo "❌ unzip required for $archive"; exit 1; }
          unzip -q "$archive" -d "${src_dir}"
          first=$(find "${src_dir}" -mindepth 1 -maxdepth 1 -type d | head -n1)
          [[ -d $first ]] && mv "$first"/* "${src_dir}/" && rmdir "$first"
          ;;
        *) echo "❌ Unknown archive format: $archive"; exit 1;;
      esac
    fi
    cmake_args+=( -D"${fetch_var[$pkg]}"="${src_dir}" )
  else
    echo "   • $pkg – no local archive → will fetch online"
    all_offline=false
  fi
done

# --- Disconnected mode
if [[ "${OFFLINE:-0}" == "1" || "$all_offline" == true ]]; then
  cmake_args+=( -DFETCHCONTENT_FULLY_DISCONNECTED=ON -DFETCHCONTENT_UPDATES_DISCONNECTED=ON )
fi

# --- Extra ad-hoc flags (space-separated)
extra_args=()
if [[ -n "${EXTRA_CMAKE_ARGS:-}" ]]; then
  # shellcheck disable=SC2206  # intentional splitting
  extra_args=(${EXTRA_CMAKE_ARGS})
fi

# --- Log summary
echo "📁 Repo root:        ${repo_root}"
echo "🏗️  Build dir:        ${build_dir}"
echo "🧰 Generator:         ${generator}"
echo "🔧 Build type:        ${CMAKE_BUILD_TYPE}"
echo "🧪 Build tests:       ${BUILD_TESTS}"
echo "🧵 Parallel jobs:     ${NPROCS}"
[[ -n "${MPIEXEC_PREFLAGS:-}" ]] && echo "🚀 MPI pref.:          ${MPIEXEC_PREFLAGS}"
[[ -n "${EXTRA_CMAKE_ARGS:-}" ]] && echo "➕ Extra CMake args:   ${EXTRA_CMAKE_ARGS}"
[[ "${OFFLINE:-0}" == "1" || "$all_offline" == true ]] && echo "📦 FetchContent:       disconnected"

# --- Configure & build
echo -e "\n⚙️  Configuring CMake…"
cmake "${cmake_args[@]}" "${extra_args[@]}"

echo -e "\n🛠️  Building…"
# Preserve non-zero exit code through the pipe
set -o pipefail

if [[ "${FILTER_DEPS_LOG:-1}" == "1" ]]; then
  cmake --build "${build_dir}" -j"${NPROCS}" 2>&1 | awk '
    BEGIN { IGNORECASE = 1 }
    # Hide lines from third-party subbuilds unless they contain errors
    /\/_deps\// && $0 !~ /(error:|fatal error|undefined reference|ld:)/ { next }
    { print; fflush() }
  '
else
  cmake --build "${build_dir}" -j"${NPROCS}"
fi

echo -e "\n✅  Build finished at ${build_dir}"
