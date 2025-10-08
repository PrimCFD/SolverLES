# Checker (parses Dimensions="..." from the single XMF/CGNS in out/)
add_executable(check_dims ${CMAKE_CURRENT_SOURCE_DIR}/check_dims.cpp)
target_compile_features(check_dims PRIVATE cxx_std_17)
target_link_libraries(check_dims PRIVATE CGNS::cgns)

# --- build the decay checker
add_executable(check_tgv_decay ${CMAKE_CURRENT_SOURCE_DIR}/check_tgv_decay.cpp)
target_compile_features(check_tgv_decay PRIVATE cxx_std_17)
target_link_libraries(check_tgv_decay PRIVATE CGNS::cgns)

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

# Helper to add a fluids regression case Usage: add_fluids_case(name nx ny nz) #
# uses $ENV{MPI_TEST_NP} if set, else 1 add_fluids_case(name nx ny nz np) #
# explicit mpi np
function(add_fluids_case name nx ny nz)
  set(case_src "${CMAKE_CURRENT_SOURCE_DIR}/${name}/config.yml")
  if(NOT EXISTS "${case_src}")
    message(FATAL_ERROR "Missing ${case_src}")
  endif()

  # Optional 5th argument: desired MPI ranks for this test
  if(ARGC GREATER 4)
    set(np "${ARGV4}")
  else()
    # Default from environment, falling back to 1
    if(DEFINED ENV{MPI_TEST_NP} AND NOT "$ENV{MPI_TEST_NP}" STREQUAL "")
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

  # --- Run solver in the case dir
  if(_MPIEXEC_EXECUTABLE AND np GREATER 1)
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

  # --- Validate output directory (checker finds the single *.xmf)
  add_test(NAME regression::fluids::${name}::check
           COMMAND $<TARGET_FILE:check_dims> "${case_bin}/out" ${nx} ${ny}
                   ${nz})
  set_tests_properties(
    regression::fluids::${name}::check
    PROPERTIES LABELS "regression;fluids" DEPENDS
               regression::fluids::${name}::run)
endfunction()

# Cases - Keep cavity fast (NP=1) - Force MPI for TG to exercise plumbing (NP=2)
add_fluids_case(cavity32 32 32 32 1)
add_fluids_case(tgv32 32 32 32 2)

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
