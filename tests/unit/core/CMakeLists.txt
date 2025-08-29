# General
add_executable(
  unit_core
  test_alignment.cpp
  test_registry_memory.cpp
  test_field_index.cpp
  test_copy_roundtrip.cpp
  test_device_mirror.cpp
  test_boundary.cpp
  test_fieldcatalog.cpp
  test_fieldcatalog_more.cpp
  test_scheduler_noop.cpp
  test_scheduler_calls.cpp
  test_plugin_register.cpp
  test_registry_error_master.cpp)

# Seam lib
add_library(mesh_seams STATIC mesh_seams.cpp)
target_include_directories(mesh_seams PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(unit_core PRIVATE core Catch2::Catch2WithMain mesh_seams)

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

if(ENABLE_CUDA)
  add_executable(unit_gpu gpu/test_um_fields.cpp)
  target_link_libraries(unit_gpu PRIVATE CUDA::cudart core
                                         Catch2::Catch2WithMain)
  target_compile_definitions(unit_gpu PRIVATE HAVE_CUDA)
endif()

if(ENABLE_MPI)
  find_package(MPI REQUIRED) # in the MPI job only

  # Halo exchange tests
  add_executable(test_halo mpi/mpi_main.cpp mpi/test_halo.cpp
                           mpi/test_halo_ng2.cpp)

  target_link_libraries(test_halo PRIVATE core Catch2::Catch2 MPI::MPI_CXX)

  add_mpi_test(core_halo_2rank test_halo 2)
  add_mpi_test(core_halo_4rank test_halo 4)

  # Boundary condition tests
  add_executable(test_boundary_mpi mpi/mpi_main.cpp mpi/test_boundary_mpi.cpp)

  target_link_libraries(test_boundary_mpi PRIVATE core Catch2::Catch2
                                                  MPI::MPI_CXX)

  # Run on a single rank to sanity-check "physical face" behavior
  add_mpi_test(core_boundary_mpi_1rank test_boundary_mpi 1)

endif()
