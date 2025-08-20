add_executable(
  unit_core test_alignment.cpp test_registry.cpp test_field_index.cpp
            test_copy_roundtrip.cpp test_device_mirror.cpp)

target_link_libraries(unit_core PRIVATE core Catch2::Catch2WithMain)

if(ENABLE_CUDA)
  target_link_libraries(unit_core PRIVATE CUDA::cudart)
  target_compile_definitions(unit_core PRIVATE HAVE_CUDA)
endif()

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

if(ENABLE_MPI)
  find_package(MPI REQUIRED) # in the MPI job only
  add_executable(test_halo mpi_main.cpp test_halo.cpp test_halo_ng2.cpp)
  target_link_libraries(test_halo PRIVATE core Catch2::Catch2 MPI::MPI_CXX)

  add_mpi_test(core_halo_2rank test_halo 2)
  add_mpi_test(core_halo_4rank test_halo 4)
endif()
