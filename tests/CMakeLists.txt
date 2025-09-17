list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
include(Testing)

option(ENABLE_TESTS_UNIT "Build unit tests" ON)
option(ENABLE_TESTS_MPI "Build MPI tests (under unit/core)" ON)
option(ENABLE_TESTS_PERF "Build performance tests" ON)
option(ENABLE_TESTS_REGRESSION "Build regression tests" ON)

# Add suites
if(ENABLE_TESTS_UNIT OR ENABLE_TESTS_MPI)
  add_subdirectory(unit)
endif()

if(ENABLE_TESTS_PERF)
  add_subdirectory(performance)
endif()

if(ENABLE_TESTS_REGRESSION)
  add_subdirectory(regression)
endif()
