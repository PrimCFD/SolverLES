# Performance benches for fluids kernels

# Header-only micro-bench helper
set(SIMPLE_BENCH_HDR ${CMAKE_SOURCE_DIR}/tests/performance/simple_bench.hpp)

# SGS (Smagorinsky)
add_executable(bench_sgs_fluids bench_sgs.cpp ${SIMPLE_BENCH_HDR})
target_link_libraries(bench_sgs_fluids PRIVATE physics_fluids core)
target_include_directories(
  bench_sgs_fluids
  PRIVATE ${CMAKE_SOURCE_DIR}/src/physics/fluids/include
          ${CMAKE_SOURCE_DIR}/src/core/include
          ${CMAKE_SOURCE_DIR}/tests/performance)
add_perf_test(perf_sgs_fluids bench_sgs_fluids)

# Diffusion - Forward Euler (FE)
add_executable(bench_diffuse_fluids bench_diffuse.cpp ${SIMPLE_BENCH_HDR})
target_link_libraries(bench_diffuse_fluids PRIVATE physics_fluids core)
target_include_directories(
  bench_diffuse_fluids
  PRIVATE ${CMAKE_SOURCE_DIR}/src/physics/fluids/include
          ${CMAKE_SOURCE_DIR}/src/core/include
          ${CMAKE_SOURCE_DIR}/tests/performance)
add_perf_test(perf_diffuse_fluids bench_diffuse_fluids)

# Diffusion - Backward Euler (BE) - Jacobi sweep (ping-pong) - Red/Black
# Gauss–Seidel (in-place, red+black)
add_executable(bench_diffuse_be_fluids bench_diffuse_be.cpp ${SIMPLE_BENCH_HDR})
target_link_libraries(bench_diffuse_be_fluids PRIVATE physics_fluids core)
target_include_directories(
  bench_diffuse_be_fluids
  PRIVATE ${CMAKE_SOURCE_DIR}/src/physics/fluids/include
          ${CMAKE_SOURCE_DIR}/src/core/include
          ${CMAKE_SOURCE_DIR}/tests/performance)
add_perf_test(perf_be_jacobi_fluids bench_diffuse_be_fluids)

# Poisson - PETSc MG
add_executable(bench_poisson_fluids bench_poisson.cpp ${SIMPLE_BENCH_HDR})
target_link_libraries(bench_poisson_fluids PRIVATE physics_fluids core
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
add_perf_test(perf_poisson_petsc_fluids bench_poisson_fluids)

# Gradient / Divergence / Corrector
add_executable(bench_grad_div_correct_fluids bench_grad_div_correct.cpp
                                             ${SIMPLE_BENCH_HDR})
target_link_libraries(bench_grad_div_correct_fluids PRIVATE physics_fluids core)
target_include_directories(
  bench_grad_div_correct_fluids
  PRIVATE ${CMAKE_SOURCE_DIR}/src/physics/fluids/include
          ${CMAKE_SOURCE_DIR}/src/core/include
          ${CMAKE_SOURCE_DIR}/tests/performance)
add_perf_test(perf_grad_div_correct_fluids bench_grad_div_correct_fluids)
