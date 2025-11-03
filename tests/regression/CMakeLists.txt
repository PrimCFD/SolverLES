# --- build the decay checker
add_executable(check_tgv_decay ${CMAKE_CURRENT_SOURCE_DIR}/check_tgv_decay.cpp)
target_compile_features(check_tgv_decay PRIVATE cxx_std_17)
target_link_libraries(check_tgv_decay PRIVATE core CGNS::cgns)

# pick up MPI launcher settings from env
set(_MPIEXEC_EXECUTABLE "$ENV{MPIEXEC_EXECUTABLE}")
set(_MPIEXEC_NUMPROC_FLAG "$ENV{MPIEXEC_NUMPROC_FLAG}")
if(_MPIEXEC_NUMPROC_FLAG STREQUAL "")
  set(_MPIEXEC_NUMPROC_FLAG "-np") # sane default
endif()

# Split pre/post flags into proper CMake lists (handles spaces/quotes)
set(_MPIEXEC_PREFLAGS_STR "$ENV{MPIEXEC_PREFLAGS}")
set(_MPIEXEC_POSTFLAGS_STR "$ENV{MPIEXEC_POSTFLAGS}")
separate_arguments(_MPIEXEC_PREFLAGS NATIVE_COMMAND "${_MPIEXEC_PREFLAGS_STR}")
separate_arguments(_MPIEXEC_POSTFLAGS NATIVE_COMMAND
                   "${_MPIEXEC_POSTFLAGS_STR}")

# Helper to add a TG decay assertion for a case produced by add_fluids_case()
function(
  add_tgv_decay_check
  name
  nu
  Lx
  Ly
  Lz
  dt
  t_end
  U0
  tol_rel)
  # Mirrors how add_fluids_case lays out the case
  set(case_bin "${CMAKE_CURRENT_BINARY_DIR}/${name}")
  add_test(NAME regression::fluids::${name}::decay
           COMMAND $<TARGET_FILE:check_tgv_decay> "${case_bin}/out" ${nu} ${Lx}
                   ${Ly} ${Lz} ${dt} ${t_end} ${U0} ${tol_rel})
  set_tests_properties(
    regression::fluids::${name}::decay
    PROPERTIES LABELS "regression;fluids" DEPENDS
               regression::fluids::${name}::run)
endfunction()

# Helper to add a fluids regression case Usage: default (MODE auto; NP from env
# or 1) add_fluids_case(name) force sequential run (no mpiexec)
# add_fluids_case(name MODE seq) force MPI run with 2 ranks add_fluids_case(name
# MODE mpi NP 2)
function(add_fluids_case name)
  # Parse optional keywords: MODE <seq|mpi>, NP <int>
  cmake_parse_arguments(CASE "" "MODE;NP" "" ${ARGN})
  string(TOLOWER "${CASE_MODE}" _mode)
  if(_mode STREQUAL "")
    set(_mode "auto")
  elseif(NOT _mode STREQUAL "seq" AND NOT _mode STREQUAL "mpi")
    message(
      FATAL_ERROR
        "add_fluids_case(${name} ...): MODE must be 'seq' or 'mpi' (got '${CASE_MODE}')"
    )
  endif()
  set(case_src "${CMAKE_CURRENT_SOURCE_DIR}/${name}/config.yml")
  if(NOT EXISTS "${case_src}")
    message(FATAL_ERROR "Missing ${case_src}")
  endif()

  # Resolve desired number of ranks
  if(_mode STREQUAL "seq")
    set(np "1")
  else()
    if(DEFINED CASE_NP AND NOT "${CASE_NP}" STREQUAL "")
      set(np "${CASE_NP}")
    elseif(DEFINED ENV{MPI_TEST_NP} AND NOT "$ENV{MPI_TEST_NP}" STREQUAL "")
      set(np "$ENV{MPI_TEST_NP}")
    else()
      set(np "1")
    endif()
  endif()

  set(case_bin "${CMAKE_CURRENT_BINARY_DIR}/${name}")
  file(MAKE_DIRECTORY "${case_bin}")

  # Inject an absolute output directory into the YAML and ensure it exists
  set(OUT_DIR "${case_bin}/out")
  file(MAKE_DIRECTORY "${OUT_DIR}")
  configure_file("${case_src}" "${case_bin}/config.yml" @ONLY)

  # --- plugin search path for dlopen()
  if(APPLE)
    set(_PLUGIN_ENV
        "DYLD_LIBRARY_PATH=$<TARGET_FILE_DIR:physics_fluids>:$ENV{DYLD_LIBRARY_PATH}"
    )
  elseif(WIN32)
    set(_PLUGIN_ENV "PATH=$<TARGET_FILE_DIR:physics_fluids>;$ENV{PATH}")
  else()
    set(_PLUGIN_ENV
        "LD_LIBRARY_PATH=$<TARGET_FILE_DIR:physics_fluids>:$ENV{LD_LIBRARY_PATH}"
    )
  endif()

  # --- Run solver in the case dir Launch rule: MODE seq           -> run
  # directly (no mpiexec) MODE mpi           -> run under mpiexec (requires
  # launcher) with np ranks MODE auto (default)-> run under mpiexec iff launcher
  # is found AND np>1
  set(_launch_direct 0)
  if(_mode STREQUAL "seq")
    set(_launch_direct 1)
  elseif(_mode STREQUAL "mpi")
    if(NOT _MPIEXEC_EXECUTABLE)
      message(
        WARNING
          "add_fluids_case(${name}): MODE=mpi but no MPI launcher found; falling back to direct run"
      )
      set(_launch_direct 1)
    endif()
  endif()

  if((_mode STREQUAL "mpi" AND NOT _launch_direct)
     OR (_mode STREQUAL "auto"
         AND _MPIEXEC_EXECUTABLE
         AND np GREATER 1))
    add_test(
      NAME regression::fluids::${name}::run
      COMMAND
        ${_MPIEXEC_EXECUTABLE} ${_MPIEXEC_NUMPROC_FLAG} ${np}
        ${_MPIEXEC_PREFLAGS} $<TARGET_FILE:solver> config.yml
        ${_MPIEXEC_POSTFLAGS}
      WORKING_DIRECTORY "${case_bin}")
  else()
    add_test(
      NAME regression::fluids::${name}::run
      COMMAND $<TARGET_FILE:solver> config.yml
      WORKING_DIRECTORY "${case_bin}")
  endif()

  set_tests_properties(
    regression::fluids::${name}::run
    PROPERTIES LABELS "regression;fluids" ENVIRONMENT "${_PLUGIN_ENV}"
               FAIL_REGULAR_EXPRESSION "error;ERROR;Error;terminate called")

endfunction()

# Cases: - Keep cavity fast: force sequential launch (2 ranks) - Force MPI for
# TG to exercise single rank plumbing (1 rank)
add_fluids_case(cavity32 MODE mpi NP 8)
add_fluids_case(tgv32 MODE seq)

# Assert TG decay at final time (2% rel tol)
add_tgv_decay_check(
  tgv32
  0.1
  6.283185307179586
  6.283185307179586
  6.283185307179586
  0.05
  5.0
  1.0
  0.02)
