cmake_minimum_required(VERSION 3.22)
project(
  SolverLES
  VERSION 0.1
  DESCRIPTION "Modular LES Plasma Solver (C++/Fortran)"
  LANGUAGES C CXX Fortran)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Enable verbose output for debugging
option(VERBOSE_MAKEFILE "Enable verbose makefile output" ON)

# Enable OpenMP if available
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  message(STATUS "OpenMP enabled")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
endif()

include(cmake/FetchCGNS.cmake)
include(cmake/FetchPETSc.cmake)

include(cmake/CompileOptions.cmake)
add_subdirectory(src/core)
add_subdirectory(src/kernels)
add_subdirectory(src/bindings)
add_subdirectory(src/plugins/flux)
add_subdirectory(src/plugins/sgs)
add_subdirectory(src/plugins/time)
add_subdirectory(src/gui)

option(BUILD_TESTS "Enable unit and regression tests" ON)
option(BUILD_EXAMPLES "Enable example cases" ON)

if(BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests/ci_inputs)
  add_subdirectory(tests/performance)
  add_subdirectory(tests/regression)
  add_subdirectory(tests/unit)
endif()

if(BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

install(
  TARGETS solver
  EXPORT SolverTargets
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

option(BUILD_DOCS "Build developer documentation" OFF)
if(BUILD_DOCS)
  find_package(Sphinx REQUIRED)
  add_subdirectory(docs)
endif()
