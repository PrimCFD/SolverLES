add_executable(unit_core test_alignment.cpp test_registry.cpp
                         test_field_index.cpp)

target_link_libraries(unit_core PRIVATE core Catch2::Catch2WithMain)

# Auto-register every TEST_CASE as a CTest test, all labeled "unit"
register_catch_tests(
  unit_core
  LABEL
  unit
  TEST_PREFIX
  "core::"
  REPORTER
  junit
  OUTPUT_DIR
  "${CMAKE_BINARY_DIR}/test-reports/unit"
  OUTPUT_PREFIX
  "catch-"
  OUTPUT_SUFFIX
  ".xml")

add_executable(test_halo test_halo.cpp mpi_main.cpp)
target_link_libraries(test_halo PRIVATE core Catch2::Catch2 MPI::MPI_CXX)

# MPI tests for this suite (registered only if ENABLE_MPI=ON)
add_mpi_test(core_halo_2rank test_halo 2)
add_mpi_test(core_halo_4rank test_halo 4)
