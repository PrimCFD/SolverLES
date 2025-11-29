# Performance benches for fluids kernels

# Helper: add a perf test that runs under MPI if a launcher is available
function(add_mpi_perf_test name exe)
  # default np (can override via env: MPI_TEST_NP)
  if(DEFINED ENV{MPI_TEST_NP})
    set(_np $ENV{MPI_TEST_NP})
  else()
    set(_np 2)
  endif()

  if(DEFINED ENV{MPIEXEC_EXECUTABLE})
    # Convert space-separated env strings into proper CMake arg lists
    set(_mpi_launcher "$ENV{MPIEXEC_EXECUTABLE}")
    set(_pre "$ENV{MPIEXEC_PREFLAGS}")
    set(_post "$ENV{MPIEXEC_POSTFLAGS}")
    if(_pre)
      separate_arguments(_pre) # -> list
    endif()
    if(_post)
      separate_arguments(_post) # -> list
    endif()

    # Many mpiexec/mpirun variants accept options in any order, but keep:
    # launcher, preflags, -np, N, postflags, executable
    add_test(NAME ${name}
             COMMAND ${_mpi_launcher} ${_pre} $ENV{MPIEXEC_NUMPROC_FLAG} ${_np}
                     ${_post} $<TARGET_FILE:${exe}>)
  else()
    message(WARNING "No MPI launcher detected; ${name} will run single-rank.")
    add_test(NAME ${name} COMMAND $<TARGET_FILE:${exe}>)
  endif()
  set_tests_properties(${name} PROPERTIES LABELS "perf")
endfunction()

# Header-only micro-bench helper
set(SIMPLE_BENCH_HDR ${CMAKE_SOURCE_DIR}/tests/performance/simple_bench.hpp)

# SGS (Smagorinsky)
add_executable(bench_sgs_fluids bench_sgs.cpp ${SIMPLE_BENCH_HDR})
target_link_libraries(bench_sgs_fluids PRIVATE physics_fluids core numerics)
target_include_directories(
  bench_sgs_fluids
  PRIVATE ${CMAKE_SOURCE_DIR}/src/physics/fluids/include
          ${CMAKE_SOURCE_DIR}/src/core/include
          ${CMAKE_SOURCE_DIR}/tests/performance)
add_perf_test(perf_sgs_fluids bench_sgs_fluids)

# Diffusion - Forward Euler (FE)
add_executable(bench_diffuse_fluids bench_diffuse.cpp ${SIMPLE_BENCH_HDR})
target_link_libraries(bench_diffuse_fluids PRIVATE physics_fluids core numerics)
target_include_directories(
  bench_diffuse_fluids
  PRIVATE ${CMAKE_SOURCE_DIR}/src/physics/fluids/include
          ${CMAKE_SOURCE_DIR}/src/core/include
          ${CMAKE_SOURCE_DIR}/tests/performance)
add_perf_test(perf_diffuse_fluids bench_diffuse_fluids)

# Diffusion - Backward Euler (BE) - Jacobi sweep (ping-pong) - Red/Black
# Gauss–Seidel (in-place, red+black)
add_executable(bench_diffuse_be_fluids bench_diffuse_be.cpp ${SIMPLE_BENCH_HDR})
target_link_libraries(bench_diffuse_be_fluids PRIVATE physics_fluids core numerics)
target_include_directories(
  bench_diffuse_be_fluids
  PRIVATE ${CMAKE_SOURCE_DIR}/src/physics/fluids/include
          ${CMAKE_SOURCE_DIR}/src/core/include
          ${CMAKE_SOURCE_DIR}/tests/performance)
add_perf_test(perf_be_jacobi_fluids bench_diffuse_be_fluids)

# Poisson - PETSc MG
add_executable(bench_poisson_fluids bench_poisson.cpp ${SIMPLE_BENCH_HDR})
target_link_libraries(bench_poisson_fluids PRIVATE physics_fluids core numerics
                                                   PETSC::petsc)
target_include_directories(
  bench_poisson_fluids
  PRIVATE ${CMAKE_SOURCE_DIR}/src/physics/fluids/include
          ${CMAKE_SOURCE_DIR}/src/core/include
          ${CMAKE_SOURCE_DIR}/tests/performance)
# Make PETSc shared libs resolvable at runtime
set(CMAKE_BUILD_RPATH_USE_ORIGIN ON)
set(_origin "\$ORIGIN")
set_target_properties(
  bench_poisson_fluids
  PROPERTIES
    BUILD_RPATH
    "${_origin};${_origin}/../../core;${_origin}/../../physics/fluids;${_origin}/../../physics/em;$<TARGET_FILE_DIR:PETSC::petsc_lib>"
    INSTALL_RPATH
    "${_origin}/../lib;${_origin}/../lib/plugins;$<TARGET_FILE_DIR:PETSC::petsc_lib>"
)
add_mpi_perf_test(perf_poisson_petsc_fluids bench_poisson_fluids)

# Gradient / Divergence / Corrector
add_executable(bench_grad_div_correct_fluids bench_grad_div_correct.cpp
                                             ${SIMPLE_BENCH_HDR})
target_link_libraries(bench_grad_div_correct_fluids PRIVATE physics_fluids core numerics)
target_include_directories(
  bench_grad_div_correct_fluids
  PRIVATE ${CMAKE_SOURCE_DIR}/src/physics/fluids/include
          ${CMAKE_SOURCE_DIR}/src/core/include
          ${CMAKE_SOURCE_DIR}/tests/performance)
add_perf_test(perf_grad_div_correct_fluids bench_grad_div_correct_fluids)
