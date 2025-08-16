#!/usr/bin/env bash
set -euo pipefail
# Build (optional) and run *unit* tests using build.sh

BUILD_DIR=${BUILD_DIR:-build-unit}
SKIP_BUILD=${SKIP_BUILD:-0}
CTEST_PARALLEL_LEVEL=${CTEST_PARALLEL_LEVEL:-$(command -v nproc >/dev/null && nproc || echo 2)}
CTEST_TIMEOUT=${CTEST_TIMEOUT:-900}
REPORT_DIR=${REPORT_DIR:-"${BUILD_DIR}/test-reports/unit"}

if [[ "${SKIP_BUILD}" != "1" ]]; then
  # Ensure dependencies & tests are configured consistently via build.sh
  BUILD_DIR="${BUILD_DIR}" \
  CMAKE_BUILD_TYPE="${CMAKE_BUILD_TYPE:-Debug}" \
  ENABLE_MPI=OFF \
  BUILD_TESTS=ON \
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
