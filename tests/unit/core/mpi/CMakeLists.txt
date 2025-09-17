# Halo exchange tests
add_executable(test_halo mpi_main.cpp test_halo.cpp test_halo_ng2.cpp)

target_link_libraries(test_halo PRIVATE core Catch2::Catch2 MPI::MPI_CXX)

add_mpi_test(core_halo_2rank test_halo 2)
add_mpi_test(core_halo_4rank test_halo 4)

# Boundary condition tests
add_executable(test_boundary_mpi mpi_main.cpp test_boundary_mpi.cpp)

target_link_libraries(test_boundary_mpi PRIVATE core Catch2::Catch2
                                                MPI::MPI_CXX)

# Run on a single rank to sanity-check "physical face" behavior
add_mpi_test(core_boundary_mpi_1rank test_boundary_mpi 1)

# Parralel io tests
add_executable(test_io mpi_main.cpp test_io_parallel.cpp)

target_link_libraries(test_io PRIVATE core Catch2::Catch2 MPI::MPI_CXX)

add_mpi_test(test_io_2rank test_io 2)
