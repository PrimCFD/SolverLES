add_executable(test_rbgs_orchestration mpi_main.cpp
                                       test_rbgs_orchestration_mpi.cpp)
target_link_libraries(test_rbgs_orchestration
                      PRIVATE physics_fluids core Catch2::Catch2 MPI::MPI_CXX)
target_compile_definitions(test_rbgs_orchestration PRIVATE HAVE_MPI)
target_include_directories(
  test_rbgs_orchestration PRIVATE ${CMAKE_SOURCE_DIR}/src/physics/fluids/include
                                  ${CMAKE_SOURCE_DIR}/src/core/include)

add_executable(test_poisson_orchestration mpi_main.cpp
                                          test_poisson_orchestration_mpi.cpp)
target_link_libraries(test_poisson_orchestration
                      PRIVATE physics_fluids core Catch2::Catch2 MPI::MPI_CXX)
target_compile_definitions(test_poisson_orchestration PRIVATE HAVE_MPI)
target_include_directories(
  test_poisson_orchestration
  PRIVATE ${CMAKE_SOURCE_DIR}/src/physics/fluids/include
          ${CMAKE_SOURCE_DIR}/src/core/include)

add_executable(test_projection_loop mpi_main.cpp test_projection_loop_mpi.cpp)
target_link_libraries(test_projection_loop PRIVATE physics_fluids core
                                                   Catch2::Catch2 MPI::MPI_CXX)
target_compile_definitions(test_projection_loop PRIVATE HAVE_MPI)
target_include_directories(
  test_projection_loop PRIVATE ${CMAKE_SOURCE_DIR}/src/physics/fluids/include
                               ${CMAKE_SOURCE_DIR}/src/core/include)

# Run them under MPI; keep label consistent with your existing suite.
add_mpi_test(fluids_orch_rbgs_2 test_rbgs_orchestration 2 LABEL mpi)
add_mpi_test(fluids_orch_poisson_2 test_poisson_orchestration 2 LABEL mpi)
add_mpi_test(fluids_orch_projection_2 test_projection_loop 2 LABEL mpi)

# Feel free to enable 4-rank runs as well:
add_mpi_test(fluids_orch_rbgs_4 test_rbgs_orchestration 4 LABEL mpi)
add_mpi_test(fluids_orch_poisson_4 test_poisson_orchestration 4 LABEL mpi)
add_mpi_test(fluids_orch_projection_4 test_projection_loop 4 LABEL mpi)
