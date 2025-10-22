#!/usr/bin/env bash
set -euo pipefail
# Build (optional) and run *regression/integration* tests using build.sh
#
# Env (all optional, same model as other runners):
#   BUILD_DIR           default: build-regression
#   SKIP_BUILD          0|1 (default 0)
#   CMAKE_BUILD_TYPE    Debug|Release|RelWithDebInfo|MinSizeRel (default Release)
#   ENABLE_CUDA         ON|OFF (default OFF)       # forwarded; usually OFF here
#   ENABLE_MPI          ON|OFF (default OFF)       # forwarded; usually OFF here
#   CTEST_PARALLEL_LEVEL  integer (default: nproc or 2)
#   CTEST_TIMEOUT       seconds (default 900)
#   REPORT_DIR          default: $BUILD_DIR/test-reports/regression
#   LABEL_RE            ctest --label-regex (default "regression")
#   CTEST_NAME_REGEX    if set, use ctest -R <regex> instead of -L (default unset)
#   EXTRA_CMAKE_ARGS_USER  extra args appended to CMake (optional)
#
# Examples:
#   ./scripts/run_regression_tests.sh
#   SKIP_BUILD=1 ./scripts/run_regression_tests.sh
#   LABEL_RE='regression|integ' ./scripts/run_regression_tests.sh
#   CTEST_NAME_REGEX='^integ::' ./scripts/run_regression_tests.sh

# Set MPI env
if [[ -f "scripts/mpi_env.sh" ]]; then
  # Override with MPI_MODE=cluster|emulate when needed.
  # shellcheck source=/dev/null
  source scripts/mpi_env.sh "${MPI_MODE:-auto}"
  # Perf-critical: turn on strict PE by default (ranks × PE must fit cores)
  export MPI_STRICT_PE="${MPI_STRICT_PE:-1}"
  export OMP_NUM_THREADS=1 # limit thread number for small runs (32^3)
  # If a launcher is available, default to building PETSc with MPI.
  if [[ -n "${MPIEXEC_EXECUTABLE:-}" ]]; then
    : "${ENABLE_MPI:=On}"
  fi
fi

BUILD_DIR=${BUILD_DIR:-build-regression}
SKIP_BUILD=${SKIP_BUILD:-0}
CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE:-Release}
CTEST_PARALLEL_LEVEL=${CTEST_PARALLEL_LEVEL:-$(command -v nproc >/dev/null && nproc || echo 2)}
CTEST_TIMEOUT=${CTEST_TIMEOUT:-3000}
REPORT_DIR=${REPORT_DIR:-"${BUILD_DIR}/test-reports/regression"}
LABEL_RE=${LABEL_RE:-regression}
CTEST_NAME_REGEX=${CTEST_NAME_REGEX:-}

PETSC_OPTIONS="-ksp_converged_reason \
-ksp_error_if_not_converged true \
-ksp_view \
-ksp_monitor_short \
-mg_levels_ksp_view \
-mg_levels_ksp_monitor"

# Uncomment for PETSC linear solve logging
# export PETSC_OPTIONS

if [[ "${SKIP_BUILD}" != "1" ]]; then
  # Suggested PETSc configure for vendored builds (ignored if USE_SYSTEM_PETSC=ON)
 : "${OPT_LEVEL:=3}"
 : "${ARCH_FLAGS:=-march=native}"
 : "${PETSC_CONFIGURE_OPTS:=--with-debugging=0 --with-mpi=1 \
   COPTFLAGS='-O'${OPT_LEVEL}' '${ARCH_FLAGS} \
   CXXOPTFLAGS='-O'${OPT_LEVEL}' '${ARCH_FLAGS} \
   FOPTFLAGS='-O'${OPT_LEVEL}' '${ARCH_FLAGS} }"
  export PETSC_CONFIGURE_OPTS
  BUILD_DIR="${BUILD_DIR}" \
  CMAKE_BUILD_TYPE="${CMAKE_BUILD_TYPE:-Release}" \
  ENABLE_MPI="${ENABLE_MPI}" \
  ENABLE_CUDA="${ENABLE_CUDA:-OFF}" \
  USE_CUDA_UM="${USE_CUDA_UM:-OFF}" \
  BUILD_TESTS=ON \
  TEST_TOGGLES="-DENABLE_TESTS_UNIT=OFF -DENABLE_TESTS_PERF=OFF -DENABLE_TESTS_REGRESSION=ON" \
  EXTRA_CMAKE_ARGS="${EXTRA_CMAKE_ARGS:+$EXTRA_CMAKE_ARGS }${TEST_TOGGLES}${EXTRA_CMAKE_ARGS_USER:+ ${EXTRA_CMAKE_ARGS_USER}}" \
  scripts/build.sh
else
  [[ -d "${BUILD_DIR}" ]] || { echo "❌ BUILD_DIR ${BUILD_DIR} not found; remove SKIP_BUILD=1"; exit 1; }
fi

mkdir -p "${REPORT_DIR}"

# Ask ctest to emit JUnit if it supports it
if ctest --help 2>/dev/null | grep -q -- '--output-junit'; then
  CTEST_JUNIT_OPTS=( --output-junit "${REPORT_DIR}/ctest-regression.xml" )
else
  CTEST_JUNIT_OPTS=()
fi

# Choose filter mode: by label (default) or by test name regex (if CTEST_NAME_REGEX set)
CTEST_FILTER_OPTS=()
if [[ -n "${CTEST_NAME_REGEX}" ]]; then
  CTEST_FILTER_OPTS=( -R "${CTEST_NAME_REGEX}" )
else
  CTEST_FILTER_OPTS=( --label-regex "${LABEL_RE}" )
fi

ctest -V --test-dir "${BUILD_DIR}" \
      "${CTEST_FILTER_OPTS[@]}" \
      -j "${CTEST_PARALLEL_LEVEL}" \
      --timeout "${CTEST_TIMEOUT}" \
      --no-tests=error \
      --output-on-failure \
      "${CTEST_JUNIT_OPTS[@]}"

# Post-run fallback: if nothing produced JUnit, create a tiny placeholder
shopt -s nullglob
reg_xmls=("${REPORT_DIR}"/*.xml)
if [[ ${#reg_xmls[@]} -eq 0 ]]; then
  printf '<testsuite name="regression" tests="0" failures="0" skipped="0"/>\n' > "${REPORT_DIR}/ctest-regression.xml"
  echo "📄 JUnit (fallback): ${REPORT_DIR}/ctest-regression.xml"
fi
