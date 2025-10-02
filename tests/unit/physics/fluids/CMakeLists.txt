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
  test_divergence_rhie_chow.cpp
  test_poisson_varcoef.cpp
  test_correct_velocity_varrho.cpp)

# Link: plugin + core + Catch2
target_link_libraries(unit_fluids PRIVATE physics_fluids core
                                          Catch2::Catch2WithMain)

# Include plugin headers (for kernels_fluids.h) and public core headers
target_include_directories(
  unit_fluids PRIVATE ${CMAKE_SOURCE_DIR}/src/physics/fluids/include
                      ${CMAKE_SOURCE_DIR}/src/core/include)

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

if(ENABLE_MPI)
  add_subdirectory(mpi)
endif()
