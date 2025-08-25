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
  if(NOT ENABLE_MPI OR NOT MPIEXEC_EXECUTABLE)
    message(STATUS "MPI disabled or launcher missing; skipping ${name}")
    return()
  endif()

  # Normalize pre/post flags (env/CLI may provide them as a single
  # space-separated string)
  set(_mpipre "${MPIEXEC_PREFLAGS}")
  set(_mpipost "${MPIEXEC_POSTFLAGS}")
  if(_mpipre)
    separate_arguments(_mpipre NATIVE_COMMAND)
  endif()
  if(_mpipost)
    separate_arguments(_mpipost NATIVE_COMMAND)
  endif()

  # JUnit target directory for MPI
  set(_junit_dir "${CMAKE_BINARY_DIR}/test-reports/mpi")
  file(MAKE_DIRECTORY "${_junit_dir}")

  if(TARGET Catch2::Catch2WithMain OR TARGET Catch2)
    add_test(
      NAME ${name}
      COMMAND
        ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np} ${_mpipre}
        $<TARGET_FILE:${target}> ${_mpipost} --reporter junit --out
        "${_junit_dir}/${name}.xml")
  else()
    add_test(NAME ${name}
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
                     ${_mpipre} $<TARGET_FILE:${target}> ${_mpipost})
  endif()

  set_property(
    TEST ${name}
    APPEND
    PROPERTY ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")

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

  # Default label if none provided
  if(NOT R_LABEL)
    set(R_LABEL unit)
  endif()

  # Skip entirely when the suite is disabled
  if(R_LABEL STREQUAL "unit"
     AND DEFINED ENABLE_TESTS_UNIT
     AND NOT ENABLE_TESTS_UNIT)
    return()
  endif()
  if(R_LABEL STREQUAL "mpi"
     AND DEFINED ENABLE_TESTS_MPI
     AND NOT ENABLE_TESTS_MPI)
    return()
  endif()
  if(R_LABEL STREQUAL "perf"
     AND DEFINED ENABLE_TESTS_PERF
     AND NOT ENABLE_TESTS_PERF)
    return()
  endif()
  if(R_LABEL STREQUAL "regression"
     AND DEFINED ENABLE_TESTS_REGRESSION
     AND NOT ENABLE_TESTS_REGRESSION)
    return()
  endif()

  # Smarter default report location (per label)
  if(NOT R_OUTPUT_DIR)
    set(R_OUTPUT_DIR "${CMAKE_BINARY_DIR}/test-reports/${R_LABEL}")
  endif()

  # Create the directory only if we're actually registering this suite
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
