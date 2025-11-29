# Single Catch2 runner for fluids unit tests
add_executable(
  unit_fluids
  test_sgs.cpp
  test_grad_div.cpp
  test_diffuse.cpp
  test_program_plan.cpp
  test_correct_velocity.cpp
  test_poisson.cpp
  test_diffuse_be.cpp
  test_poisson_varcoef.cpp
  test_correct_velocity_varrho.cpp)

# Link: plugin + core + Catch2
target_link_libraries(unit_fluids PRIVATE physics_fluids core numerics
                                          Catch2::Catch2WithMain PETSC::petsc)

# Include plugin headers (for kernels_fluids.h) and public core headers
target_include_directories(
  unit_fluids PRIVATE ${CMAKE_SOURCE_DIR}/src/physics/fluids/include
                      ${CMAKE_SOURCE_DIR}/src/core/include)

# Make PETSc shared libs resolvable at runtime
set(CMAKE_BUILD_RPATH_USE_ORIGIN ON)
set(_origin "\$ORIGIN")
set_target_properties(
  unit_fluids
  PROPERTIES
    BUILD_RPATH
    "${_origin};${_origin}/../../core;${_origin}/../../physics/fluids;${_origin}/../../physics/em;$<TARGET_FILE_DIR:PETSC::petsc_lib>"
    INSTALL_RPATH
    "${_origin}/../lib;${_origin}/../lib/plugins;$<TARGET_FILE_DIR:PETSC::petsc_lib>"
)

# Auto-register each TEST_CASE as a ctest (label "unit"), JUnit to
# test-reports/unit
register_catch_tests(
  unit_fluids
  LABEL
  unit
  TEST_PREFIX
  "fluids::"
  REPORTER
  junit
  OUTPUT_DIR
  "${CMAKE_BINARY_DIR}/test-reports/unit"
  OUTPUT_PREFIX
  "catch-"
  OUTPUT_SUFFIX
  ".xml")
