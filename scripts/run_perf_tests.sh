#!/usr/bin/env bash
set -euo pipefail
# Build (optional) and run *perf* tests using build.sh

BUILD_DIR=${BUILD_DIR:-build-perf}
SKIP_BUILD=${SKIP_BUILD:-0}
CTEST_PARALLEL_LEVEL=${CTEST_PARALLEL_LEVEL:-$(command -v nproc >/dev/null && nproc || echo 2)}
CTEST_TIMEOUT=${CTEST_TIMEOUT:-900}
REPORT_DIR=${REPORT_DIR:-"${BUILD_DIR}/test-reports/perf"}

if [[ "${SKIP_BUILD}" != "1" ]]; then
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
  ENABLE_MPI=OFF \
  ENABLE_CUDA="${ENABLE_CUDA}" \
  USE_CUDA_UM="${USE_CUDA_UM:-OFF}" \
  BUILD_TESTS=ON \
  EXTRA_CMAKE_ARGS="-DENABLE_TESTS_UNIT=OFF -DENABLE_TESTS_MPI=OFF -DENABLE_TESTS_PERF=ON -DENABLE_TESTS_REGRESSION=OFF ${EXTRA_CMAKE_ARGS_USER:-}" \
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

ctest --test-dir "${BUILD_DIR}" \
      -L perf \
      -j "${CTEST_PARALLEL_LEVEL}" \
      --timeout "${CTEST_TIMEOUT}" \
      --no-tests=error \
      --output-on-failure \
      "${CTEST_JUNIT_OPTS[@]}"

# Post-run fallback: only if no XML was produced by CTest/Catch2
shopt -s nullglob
perf_xmls=("${REPORT_DIR}"/*.xml)
if [[ ${#perf_xmls[@]} -eq 0 ]]; then
  printf '<testsuite name="perf" tests="0" failures="0" skipped="0"/>\n' > "${REPORT_DIR}/ctest-perf.xml"
  echo "📄 JUnit (fallback): ${REPORT_DIR}/ctest-perf.xml"
fi
