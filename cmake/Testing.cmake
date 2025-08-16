include(CTest)
include(CMakeParseArguments)

# Default mpiexec flags (OpenMPI on CI needs --oversubscribe)
if(ENABLE_MPI AND NOT DEFINED MPIEXEC_PREFLAGS)
  set(MPIEXEC_PREFLAGS
      "--oversubscribe"
      CACHE STRING "Extra flags passed to mpiexec")
endif()

# Robust MPI test helper: becomes a no-op if MPI is disabled. Usage: Usage:
function(add_mpi_test name target np)
  add_test(NAME ${name} COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG}
                                ${np} $<TARGET_FILE:${target}>)
  set_tests_properties(${name} PROPERTIES LABELS mpi)
endfunction()

# Catch2 discovery + labeling + JUnit output in one place. Usage:
# register_catch_tests(<target> LABEL <lbl> TEST_PREFIX <pfx> OUTPUT_DIR <dir>
# [EXTRA_ARGS ...])
function(register_catch_tests target)
  set(options)
  set(oneValueArgs LABEL TEST_PREFIX OUTPUT_DIR)
  set(multiValueArgs EXTRA_ARGS)
  cmake_parse_arguments(R "${options}" "${oneValueArgs}" "${multiValueArgs}"
                        ${ARGN})

  # Sensible defaults
  if(NOT R_OUTPUT_DIR)
    set(R_OUTPUT_DIR "${CMAKE_BINARY_DIR}/test-reports")
  endif()
  file(MAKE_DIRECTORY "${R_OUTPUT_DIR}")

  if(USE_SYSTEM_CATCH2)
    find_package(Catch2 3 REQUIRED)
  else()
    if(NOT TARGET Catch2::Catch2WithMain)
      message(
        FATAL_ERROR
          "Catch2::Catch2WithMain not found.\n"
          "Either build with -DUSE_SYSTEM_CATCH2=ON and set Catch2_DIR (or CMAKE_PREFIX_PATH),\n"
          "or ensure the root included cmake/FetchCatch2.cmake when BUILD_TESTS=ON."
      )
    endif()
  endif()

  include(Catch) # ships with Catch2
  catch_discover_tests(
    ${target}
    TEST_PREFIX
    "${R_TEST_PREFIX}"
    WORKING_DIRECTORY
    ${CMAKE_CURRENT_BINARY_DIR}
    EXTRA_ARGS
    ${R_EXTRA_ARGS}
    REPORTER
    junit
    OUTPUT_DIR
    "${R_OUTPUT_DIR}"
    OUTPUT_SUFFIX
    .xml
    PROPERTIES
    LABELS
    "${R_LABEL}")
endfunction()

# Performance tests Register a perf test with sane defaults: - Label "perf" -
# Run serially (to reduce noise) - Treat lines containing "skipped" as a skipped
# test (not a failure) - Allow extra args after the target name
function(add_perf_test name target)
  add_test(NAME ${name} COMMAND $<TARGET_FILE:${target}> ${ARGN})
  set_tests_properties(
    ${name}
    PROPERTIES LABELS
               "perf"
               RUN_SERIAL
               TRUE
               TIMEOUT
               600
               SKIP_REGULAR_EXPRESSION
               "skipped|SKIPPED")
endfunction()
