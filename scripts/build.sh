#!/usr/bin/env bash
set -euo pipefail

# -------------------------------------------------
# build.sh — centralized, env-driven CMake build
# -------------------------------------------------
# Env (all optional):
#   BUILD_DIR                build dir (relative to repo root unless absolute) [default: build]
#   CMAKE_BUILD_TYPE         Debug/Release/RelWithDebInfo/MinSizeRel [default: Release]
#   BUILD_TESTS              ON/OFF (configure tests) [default: ON]
#   MPIEXEC_PREFLAGS         forwarded to CMake (e.g., --oversubscribe for OpenMPI)
#   CMAKE_GENERATOR          override generator (e.g., "Ninja")
#   CMAKE_TOOLCHAIN_FILE     forwarded if set
#   EXTRA_CMAKE_ARGS         extra cmake args (space-separated)
#   OFFLINE                  if "1", force disconnected FetchContent
#
# Examples:
#   CMAKE_BUILD_TYPE=Debug ./build.sh
#   BUILD_DIR=build-regression ./build.sh
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

# --- OpenMP runtime defaults (honor user overrides if already set)
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-$NPROCS}"
export OMP_PLACES="${OMP_PLACES:-cores}"
export OMP_PROC_BIND="${OMP_PROC_BIND:-close}"
export OMP_DYNAMIC="${OMP_DYNAMIC:-FALSE}"
# Avoid nested BLAS threads stepping on OpenMP
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"

# --- Defaults (overridable)
CMAKE_BUILD_TYPE="${CMAKE_BUILD_TYPE:-Release}"
BUILD_TESTS="${BUILD_TESTS:-OFF}"
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
  -DUSE_SYSTEM_YAML_CPP=OFF
  -DUSE_SYSTEM_BLAS=OFF
  -DUSE_SYSTEM_CATCH2=OFF
  -DBUILD_TESTS="${BUILD_TESTS}"
  -DENABLE_CUDA="${ENABLE_CUDA}"
  -DUSE_CUDA_UM="${USE_CUDA_UM}"
  -DCMAKE_MESSAGE_LOG_LEVEL=WARNING -Wno-dev
)

# If the user points to a system PETSc, prefer it
if [[ -n "${PETSC_DIR:-}" ]]; then
  cmake_args+=( -DUSE_SYSTEM_PETSC=ON )
  [[ -n "${PETSC_ARCH:-}" ]] && cmake_args+=( -DPETSC_ARCH="${PETSC_ARCH}" )
fi

 # Strong but *overridable* defaults
 : "${OPT_LEVEL:=3}"          # allow OPT_LEVEL=2 etc.
 : "${ARCH_FLAGS:=-march=native}"  # allow ARCH_FLAGS="" for portability builds
 rel_flags="-O${OPT_LEVEL} -DNDEBUG ${ARCH_FLAGS}"
 if [[ -z "${CMAKE_C_FLAGS_RELEASE:-}" ]];      then cmake_args+=( -DCMAKE_C_FLAGS_RELEASE="${rel_flags}" ); fi
 if [[ -z "${CMAKE_CXX_FLAGS_RELEASE:-}" ]];    then cmake_args+=( -DCMAKE_CXX_FLAGS_RELEASE="${rel_flags}" ); fi
 if [[ -z "${CMAKE_Fortran_FLAGS_RELEASE:-}" ]]; then cmake_args+=( -DCMAKE_Fortran_FLAGS_RELEASE="${rel_flags}" ); fi

[[ -n "${MPIEXEC_NUMPROC_FLAG:-}" ]] && cmake_args+=( -DMPIEXEC_NUMPROC_FLAG="${MPIEXEC_NUMPROC_FLAG}" )
[[ -n "${MPIEXEC_PREFLAGS:-}"    ]] && cmake_args+=( -DMPIEXEC_PREFLAGS="${MPIEXEC_PREFLAGS}" )
[[ -n "${MPIEXEC_POSTFLAGS:-}"   ]] && cmake_args+=( -DMPIEXEC_POSTFLAGS="${MPIEXEC_POSTFLAGS}" )
[[ -n "${MPIEXEC_EXECUTABLE:-}"  ]] && cmake_args+=( -DMPIEXEC_EXECUTABLE="${MPIEXEC_EXECUTABLE}" )
[[ -n "${CMAKE_TOOLCHAIN_FILE:-}" ]] && cmake_args+=( -DCMAKE_TOOLCHAIN_FILE="${CMAKE_TOOLCHAIN_FILE}" )

# Honor explicit MPI wrappers from environment
if [[ -n "${MPI_C_COMPILER:-}" ]]; then
  cmake_args+=( -DMPI_C_COMPILER="${MPI_C_COMPILER}" )
fi
if [[ -n "${MPI_CXX_COMPILER:-}" ]]; then
  cmake_args+=( -DMPI_CXX_COMPILER="${MPI_CXX_COMPILER}" )
fi
if [[ -n "${MPI_Fortran_COMPILER:-}" ]]; then
  cmake_args+=( -DMPI_Fortran_COMPILER="${MPI_Fortran_COMPILER}" )
fi

# --- Scan extern/ for pre-fetched archives; wire up FETCHCONTENT_SOURCE_DIR_* 
declare -A fetch_var=(
  [Catch2]=FETCHCONTENT_SOURCE_DIR_CATCH2
  [CGNS]=FETCHCONTENT_SOURCE_DIR_CGNS
  [HDF5]=FETCHCONTENT_SOURCE_DIR_HDF5
  [OpenBLAS]=FETCHCONTENT_SOURCE_DIR_OPENBLAS
  [PETSc]=FETCHCONTENT_SOURCE_DIR_PETSC
  [yaml-cpp]=FETCHCONTENT_SOURCE_DIR_YAML_CPP 
)
packages=(Catch2 CGNS HDF5 OpenBLAS PETSc yaml-cpp)

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

# If PETSC_DIR is set by the user (system PETSc), make it win
if [[ -n "${PETSC_DIR:-}" ]]; then
  export PKG_CONFIG_PATH="${PETSC_DIR}/lib/pkgconfig:${PETSC_DIR}/${PETSC_ARCH:-}/lib/pkgconfig:${PKG_CONFIG_PATH:-}"
  export CMAKE_PREFIX_PATH="${PETSC_DIR}:${CMAKE_PREFIX_PATH:-}"
fi

# If using vendored PETSc, expose its pkg-config after it’s built:
# (on re-configures, it may already exist; adding path early is harmless)
if [[ -z "${PETSC_DIR:-}" ]]; then
  export PKG_CONFIG_PATH="${BUILD_DIR:-build}/_deps/petsc-install/lib/pkgconfig:${PKG_CONFIG_PATH:-}"
  export CMAKE_PREFIX_PATH="${BUILD_DIR:-build}/_deps/petsc-install:${CMAKE_PREFIX_PATH:-}"
fi

# --- Log summary
echo "📁 Repo root:        ${repo_root}"
echo "🏗️  Build dir:        ${build_dir}"
echo "🧰 Generator:         ${generator}"
echo "🔧 Build type:        ${CMAKE_BUILD_TYPE}"
echo "🧪 Build tests:       ${BUILD_TESTS}"
echo "🧵 Parallel jobs:     ${NPROCS}"
echo "🧵 OMP_NUM_THREADS:   ${OMP_NUM_THREADS}"
echo "📍 OMP_PLACES:        ${OMP_PLACES}"
echo "🧲 OMP_PROC_BIND:     ${OMP_PROC_BIND}"
echo "🌀 OMP_DYNAMIC:       ${OMP_DYNAMIC}"
[[ -n "${MPIEXEC_PREFLAGS:-}" ]] && echo "🚀 MPI pref.:          ${MPIEXEC_PREFLAGS}"
[[ -n "${EXTRA_CMAKE_ARGS:-}" ]] && echo "➕ Extra CMake args:   ${EXTRA_CMAKE_ARGS}"
[[ "${OFFLINE:-0}" == "1" || "$all_offline" == true ]] && echo "📦 FetchContent:       disconnected"
echo "🐾 PETSC_DIR:          ${PETSC_DIR:-<unset>}"
echo "🐾 PKG_CONFIG_PATH:    ${PKG_CONFIG_PATH:-<unset>}"

petsc_doctor() {
  # Only run if pkg-config can find PETSc (system or vendored staged)
  if ! pkg-config --exists petsc 2>/dev/null; then
    echo "ℹ️  PETSc not on pkg-config yet (vendored build may add it later)."
    return 0
  fi
  local ver pref req mpiexec cflags xflags
  ver="$(pkg-config --modversion petsc 2>/dev/null || true)"
  pref="$(pkg-config --variable=prefix petsc 2>/dev/null || true)"
  req="$(pkg-config --variable=requires_private petsc 2>/dev/null || true)"
  mpiexec="$(pkg-config --variable=mpiexec petsc 2>/dev/null || true)"
  cflags="$(pkg-config --variable=cflags_extra petsc 2>/dev/null || true)"
  xflags="$(pkg-config --variable=cxxflags_extra petsc 2>/dev/null || true)"
  echo "🔎 PETSc doctor:"
  echo "   • version:   ${ver}"
  echo "   • prefix:    ${pref}"
  echo "   • requires:  ${req}"
  echo "   • mpiexec:   ${mpiexec}"
  echo "   • cflags:    ${cflags}"
  echo "   • cxxflags:  ${xflags}"

  # Hard failures when the selection contradicts the build mode
  if [[ "${mpiexec}" == *"petsc-mpiexec.uni"* ]]; then
    echo "❌ PETSc is serial (mpiuni)."
    echo "   Fix by exporting PETSC_DIR (MPI build) or rebuilding vendored PETSc with --with-mpi=1."
    exit 2
  fi
  if [[ "${CMAKE_BUILD_TYPE:-Release}" =~ ^[Rr]elease$ ]] && [[ "${cflags}${xflags}" == *"-O0"* ]]; then
    echo "❌ PETSc was compiled with -O0 (debug) but you requested a Release build."
    echo "   Fix by using an optimized system PETSc, or pass optimized configure options to vendored PETSc."
    exit 3
  fi
}

# --- Configure & build
echo -e "\n⚙️  Configuring CMake…"
petsc_doctor
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
