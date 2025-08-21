#!/usr/bin/env bash
# Run MPI tests via CTest with a laptop/cluster-friendly environment.
# Requires: build.sh in repo root, mpi_env.sh in repo root.
#
# Examples:
#   ./run_mpi_tests.sh                         # auto mode, default dirs
#   ./run_mpi_tests.sh --mode emulate --np 4   # laptop OpenMPI dev
#   ./run_mpi_tests.sh --mode cluster --label mpi
#   BUILD_DIR=build-mpi ./run_mpi_tests.sh --np 8

set -Eeuo pipefail

usage() {
  cat <<'USAGE'
Usage: run_mpi_tests.sh [options] [-- [extra ctest args]]
Options:
  --mode {auto|emulate|cluster}  MPI environment mode (default: auto)
  --np N                         Default MPI ranks for tests that consult NP (default: 2)
  --build-dir DIR                Build directory (default: build-mpi, or $BUILD_DIR)
  --type {Debug|Release|RelWithDebInfo|MinSizeRel}  CMAKE_BUILD_TYPE (default: Debug)
  --label REGEX                  ctest --label-regex (default: mpi)
  --keep-going                   pass -K/--no-stop-on-failure to ctest
  -h, --help                     Show this help
Everything after "--" is passed through to ctest.
USAGE
}

# --- defaults ---
MODE="${MODE:-auto}"
NP="${NP:-2}"
BUILD_DIR="${BUILD_DIR:-build-mpi}"
CMAKE_BUILD_TYPE="${CMAKE_BUILD_TYPE:-Debug}"
LABEL_RE="mpi"
CTEST_PASSTHRU=()
KEEP_GOING=0

# --- parse args ---
while [[ $# -gt 0 ]]; do
  case "$1" in
    --mode) MODE="${2:-}"; shift ;;
    --np) NP="${2:-}"; shift ;;
    --build-dir) BUILD_DIR="${2:-}"; shift ;;
    --type) CMAKE_BUILD_TYPE="${2:-}"; shift ;;
    --label) LABEL_RE="${2:-}"; shift ;;
    --keep-going) KEEP_GOING=1 ;;
    -h|--help) usage; exit 0 ;;
    --) shift; CTEST_PASSTHRU+=("$@"); break ;;
    *) CTEST_PASSTHRU+=("$1") ;;
  esac
  shift || true
done

# --- repo root ---
_script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
repo_root="$_script_dir"
while [[ ! -f "$repo_root/scripts/mpi_env.sh" && "$repo_root" != "/" ]]; do
  repo_root="$(dirname "$repo_root")"
done
if [[ ! -f "$repo_root/scripts/mpi_env.sh" ]]; then
  echo "❌ mpi_env.sh not found by walking up from $_script_dir"
  echo "   Looked for <repo>/scripts/mpi_env.sh"
  exit 1
fi


# --- source MPI env (exports MPIEXEC_* variables) ---
if [[ ! -f "$repo_root/scripts/mpi_env.sh" ]]; then
  echo "❌ mpi_env.sh not found at repo root ($repo_root)."; exit 1
fi
# shellcheck disable=SC1091
source "$repo_root/scripts/mpi_env.sh" "$MODE"

# Fallback if launcher isn't found
if [[ -z "${MPIEXEC_EXECUTABLE:-}" ]]; then
  if command -v srun >/dev/null && [[ -n "${SLURM_JOB_ID:-}" ]]; then
    export MPIEXEC_EXECUTABLE="srun"; export MPIEXEC_NUMPROC_FLAG="-n"
  elif command -v mpirun >/dev/null; then
    export MPIEXEC_EXECUTABLE="mpirun"; export MPIEXEC_NUMPROC_FLAG="-np"
  elif command -v mpiexec >/dev/null; then
    export MPIEXEC_EXECUTABLE="mpiexec"; export MPIEXEC_NUMPROC_FLAG="-np"
  else
    echo "⚠️  No MPI launcher (srun/mpirun/mpiexec) found; skipping MPI tests."
    # Leave a placeholder JUnit so CI artifact upload succeeds
    report_dir="${BUILD_DIR}/test-reports/mpi"
    mkdir -p "$report_dir"
    printf '<testsuite name="mpi" tests="0" failures="0" skipped="0"/>\n' > "${report_dir}/ctest-mpi.xml"
    echo "📄 JUnit: ${report_dir}/ctest-mpi.xml"
    exit 0
  fi
fi

# --- configure/build with MPI and pass launcher hints to CMake/CTest ---
echo "🔧 Configuring & building in ${BUILD_DIR} (type=${CMAKE_BUILD_TYPE})"
ENABLE_MPI=ON BUILD_TESTS=ON \
ENABLE_CUDA="${ENABLE_CUDA:-OFF}" \
USE_CUDA_UM="${USE_CUDA_UM:-OFF}" \
MPIEXEC_EXECUTABLE="${MPIEXEC_EXECUTABLE:-}" \
MPIEXEC_NUMPROC_FLAG="${MPIEXEC_NUMPROC_FLAG:-}" \
MPIEXEC_PREFLAGS="${MPIEXEC_PREFLAGS:-}" \
MPIEXEC_POSTFLAGS="${MPIEXEC_POSTFLAGS:-}" \
CMAKE_BUILD_TYPE="${CMAKE_BUILD_TYPE}" \
BUILD_DIR="${BUILD_DIR}" \
  EXTRA_CMAKE_ARGS="-DENABLE_TESTS_UNIT=OFF -DENABLE_TESTS_MPI=ON -DENABLE_TESTS_PERF=OFF -DENABLE_TESTS_REGRESSION=OFF" \
$repo_root/scripts/build.sh

# --- test report path (JUnit if supported) ---
now="$(date +%Y%m%d-%H%M%S)"
report_dir="${BUILD_DIR}/test-reports/mpi"
mkdir -p "$report_dir"
junit_xml="${report_dir}/ctest-mpi-${now}.xml"

# Detect if this ctest supports --output-junit; no pre-run placeholder
if ctest --help 2>/dev/null | grep -q -- '--output-junit'; then
  CTEST_JUNIT_OPTS=( --output-junit "$junit_xml" )
else
  CTEST_JUNIT_OPTS=()
fi

# Keep-going flag
if [[ $KEEP_GOING -eq 1 ]]; then
  CTEST_K_FLAG=( --no-stop-on-failure )
else
  CTEST_K_FLAG=()
fi

# --- run tests ---
echo "🧪 Running MPI tests (label-regex=${LABEL_RE})"
(
  cd "$BUILD_DIR"
  # Allow users to hint CTest parallelism for independent tests (not MPI ranks)
  : "${CTEST_PARALLEL_LEVEL:=1}"
  export CTEST_PARALLEL_LEVEL
  set -x
  ctest --label-regex "${LABEL_RE}" \
        --output-on-failure \
        "${CTEST_JUNIT_OPTS[@]}" \
        "${CTEST_K_FLAG[@]}" \
        "${CTEST_PASSTHRU[@]}"
)
rc=$?

# Post-run fallback: only if no XML was produced by CTest or test binaries
shopt -s nullglob
mpi_xmls=("${report_dir}"/*.xml)
if [[ ${#mpi_xmls[@]} -eq 0 ]]; then
  printf '<testsuite name="mpi" tests="0" failures="0" skipped="0"/>\n' > "$junit_xml"
  echo "📄 JUnit (fallback): ${junit_xml}"
elif [[ ${#CTEST_JUNIT_OPTS[@]} -gt 0 ]]; then
  # CTest wrote the summary XML we asked for
  echo "📄 JUnit: ${junit_xml}"
fi

exit "$rc"
