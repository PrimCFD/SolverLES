#!/usr/bin/env bash
set -euo pipefail
# Build (optional) and run *perf* tests using build.sh


# Set MPI env

if [[ -f "scripts/mpi_env.sh" ]]; then
  # Override with MPI_MODE=cluster|emulate when needed.
  # shellcheck source=/dev/null
  source scripts/mpi_env.sh "${MPI_MODE:-auto}"
  # Perf-critical: turn on strict PE by default (ranks × PE must fit cores)
  export MPI_STRICT_PE="${MPI_STRICT_PE:-1}"
  export OMP_NUM_THREADS=2 # limit thread number/MPI
  # If a launcher is available, default to building PETSc with MPI.
  if [[ -n "${MPIEXEC_EXECUTABLE:-}" ]]; then
    : "${ENABLE_MPI:=On}"
  fi
fi

PETSC_OPTIONS="-ksp_converged_reason \
-ksp_monitor_short"

# Uncomment for PETSC linear solve logging
# export PETSC_OPTIONS

BUILD_DIR=${BUILD_DIR:-build-perf}
SKIP_BUILD=${SKIP_BUILD:-0}
CTEST_PARALLEL_LEVEL=${CTEST_PARALLEL_LEVEL:-$(command -v nproc >/dev/null && nproc || echo 2)}
CTEST_TIMEOUT=${CTEST_TIMEOUT:-900}
REPORT_DIR=${REPORT_DIR:-"${BUILD_DIR}/test-reports/perf"}

if [[ "${SKIP_BUILD}" != "1" ]]; then
  # Suggested PETSc configure for vendored builds (ignored if USE_SYSTEM_PETSC=ON)
 : "${OPT_LEVEL:=3}"
 : "${ARCH_FLAGS:=-march=native}"
 : "${PETSC_CONFIGURE_OPTS:=--with-debugging=0 --with-mpi=1 \
   COPTFLAGS='-O'${OPT_LEVEL}' '${ARCH_FLAGS} \
   CXXOPTFLAGS='-O'${OPT_LEVEL}' '${ARCH_FLAGS} \
   FOPTFLAGS='-O'${OPT_LEVEL}' '${ARCH_FLAGS} }"
 export PETSC_CONFIGURE_OPTS
  # Auto-detect CUDA unless user forced it
  if [[ -z "${ENABLE_CUDA:-}" || "${ENABLE_CUDA}" == "AUTO" ]]; then
    if command -v nvcc >/dev/null || [[ -d "${CUDAToolkit_ROOT:-/usr/local/cuda}" ]]; then
      ENABLE_CUDA=ON
      : "${USE_CUDA_UM:=ON}"
    else
      ENABLE_CUDA=OFF
    fi
  fi

  BUILD_DIR="${BUILD_DIR}" \
  CMAKE_BUILD_TYPE="${CMAKE_BUILD_TYPE:-Release}" \
  ENABLE_MPI="${ENABLE_MPI}" \
  ENABLE_CUDA="${ENABLE_CUDA}" \
  USE_CUDA_UM="${USE_CUDA_UM:-OFF}" \
  BUILD_TESTS=ON \
  TEST_TOGGLES="-DENABLE_TESTS_UNIT=OFF -DENABLE_TESTS_PERF=ON -DENABLE_TESTS_REGRESSION=OFF" \
  EXTRA_CMAKE_ARGS="${EXTRA_CMAKE_ARGS:+$EXTRA_CMAKE_ARGS }${TEST_TOGGLES}${EXTRA_CMAKE_ARGS_USER:+ ${EXTRA_CMAKE_ARGS_USER}}" \
  scripts/build.sh
else
  [[ -d "${BUILD_DIR}" ]] || { echo "❌ BUILD_DIR ${BUILD_DIR} not found; remove SKIP_BUILD=1"; exit 1; }
fi


mkdir -p "${REPORT_DIR}"

if ctest --help 2>/dev/null | grep -q -- '--output-junit'; then
  CTEST_JUNIT_OPTS=( --output-junit "${REPORT_DIR}/ctest-perf.xml" )
else
  CTEST_JUNIT_OPTS=()
fi

ctest -V --test-dir "${BUILD_DIR}" \
      -L perf \
      -j "${CTEST_PARALLEL_LEVEL}" \
      --timeout "${CTEST_TIMEOUT}" \
      --no-tests=error \
      --output-on-failure \
      "${CTEST_JUNIT_OPTS[@]}" \

# Post-run fallback: only if no XML was produced by CTest/Catch2
shopt -s nullglob
perf_xmls=("${REPORT_DIR}"/*.xml)
if [[ ${#perf_xmls[@]} -eq 0 ]]; then
  printf '<testsuite name="perf" tests="0" failures="0" skipped="0"/>\n' > "${REPORT_DIR}/ctest-perf.xml"
  echo "📄 JUnit (fallback): ${REPORT_DIR}/ctest-perf.xml"
fi
