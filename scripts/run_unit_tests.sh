#!/usr/bin/env bash
set -euo pipefail
# Build (optional) and run *unit* tests using build.sh

BUILD_DIR=${BUILD_DIR:-build-unit}
SKIP_BUILD=${SKIP_BUILD:-0}
CTEST_PARALLEL_LEVEL=${CTEST_PARALLEL_LEVEL:-$(command -v nproc >/dev/null && nproc || echo 2)}
CTEST_TIMEOUT=${CTEST_TIMEOUT:-900}
REPORT_DIR=${REPORT_DIR:-"${BUILD_DIR}/test-reports/unit"}

if [[ "${SKIP_BUILD}" != "1" ]]; then
  BUILD_DIR="${BUILD_DIR}" \
  CMAKE_BUILD_TYPE="${CMAKE_BUILD_TYPE:-Debug}" \
  ENABLE_CUDA="${ENABLE_CUDA:-OFF}" \
  USE_CUDA_UM="${USE_CUDA_UM:-OFF}" \
  BUILD_TESTS=ON \
  TEST_TOGGLES="-DENABLE_TESTS_UNIT=ON -DENABLE_TESTS_PERF=OFF -DENABLE_TESTS_REGRESSION=OFF" \
  EXTRA_CMAKE_ARGS="${EXTRA_CMAKE_ARGS:+$EXTRA_CMAKE_ARGS }${TEST_TOGGLES}${EXTRA_CMAKE_ARGS_USER:+ ${EXTRA_CMAKE_ARGS_USER}}" \
  scripts/build.sh
else
  [[ -d "${BUILD_DIR}" ]] || { echo "❌ BUILD_DIR ${BUILD_DIR} not found; remove SKIP_BUILD=1"; exit 1; }
fi


mkdir -p "${REPORT_DIR}"

if ctest --help 2>/dev/null | grep -q -- '--output-junit'; then
  CTEST_JUNIT_OPTS=( --output-junit "${REPORT_DIR}/ctest-unit.xml" )
else
  CTEST_JUNIT_OPTS=()
fi

ctest --test-dir "${BUILD_DIR}" \
      -L unit \
      -j "${CTEST_PARALLEL_LEVEL}" \
      --timeout "${CTEST_TIMEOUT}" \
      --no-tests=error \
      --output-on-failure \
      "${CTEST_JUNIT_OPTS[@]}"

# Post-run fallback: only if no XML was produced by CTest/Catch2
shopt -s nullglob
unit_xmls=("${REPORT_DIR}"/*.xml)
if [[ ${#unit_xmls[@]} -eq 0 ]]; then
  printf '<testsuite name="unit" tests="0" failures="0" skipped="0"/>\n' > "${REPORT_DIR}/ctest-unit.xml"
  echo "📄 JUnit (fallback): ${REPORT_DIR}/ctest-unit.xml"
fi
