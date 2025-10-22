list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
include(Testing)

option(ENABLE_TESTS_UNIT "Build unit tests" ON)
option(ENABLE_TESTS_PERF "Build performance tests" ON)
option(ENABLE_TESTS_REGRESSION "Build regression tests" ON)

# Add suites
if(ENABLE_TESTS_UNIT)
  add_subdirectory(unit)
endif()

if(ENABLE_TESTS_PERF)
  add_subdirectory(performance)
endif()

if(ENABLE_TESTS_REGRESSION)
  add_subdirectory(regression)
endif()
