# Root listfile – only global policies, dependencies, and one add_subdirectory

cmake_minimum_required(VERSION 3.24)

# Prefetch mode (nothing else gets configured)
if(PREFETCH_THIRD_PARTY)
  message(STATUS "Prefetch mode: downloading archives and exiting")
  include(cmake/PrefetchDependencies.cmake)
  return()
endif()

project(
  SolverLES
  VERSION 0.1
  DESCRIPTION "Modular LES Plasma Solver (C++/Fortran)"
  LANGUAGES C CXX Fortran)

# Global build options
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON) # shared libs & plugins

option(VERBOSE_MAKEFILE "Enable verbose makefile output" ON)
option(USE_CUDA_UM "Enable CUDA Unified Memory path" OFF)
option(BUILD_TESTS "Enable unit and regression tests" ON)
option(BUILD_EXAMPLES "Build example input cases" ON)
option(BUILD_GUI "Build GUI front-end + launcher" OFF)
option(BUILD_PLUGINS "Build run-time plug-ins" ON)

# Third‑party toggles
option(USE_SYSTEM_HDF5 "Use HDF5 provided by the system" OFF)
option(USE_SYSTEM_CGNS "Use CGNS provided by the system" OFF)
option(USE_SYSTEM_PETSC "Use PETSc provided by the system" OFF)
option(USE_SYSTEM_BLAS "Use BLAS/LAPACK provided by the system" OFF)
option(USE_SYSTEM_YAML_CPP "Use yaml-cpp provided by the system" OFF)
option(USE_SYSTEM_CATCH2 "Use Catch2 provided by the system" OFF)

# Make it easy to pre‑download archives (& exit)
option(PREFETCH_THIRD_PARTY "Download archives into extern/ and exit" OFF)

set(FETCHCONTENT_FULLY_DISCONNECTED
    OFF
    CACHE BOOL "Never attempt to reach the network even if URLs are missing")

# Third–party dependencies that every target needs
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  message(STATUS "OpenMP enabled")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
endif()

option(ENABLE_MPI "Build with MPI support" ON)
if(ENABLE_MPI)
  find_package(MPI REQUIRED)
  add_compile_definitions(HAVE_MPI)
endif()

# CUDA toggle
option(ENABLE_CUDA "Enable CUDA features (HAVE_CUDA, link cudart)" OFF)

if(ENABLE_CUDA)
  find_package(CUDAToolkit REQUIRED)
  add_compile_definitions(HAVE_CUDA)

  if(USE_CUDA_UM)
    add_compile_definitions(USE_CUDA_UM)
  endif()
endif()

# Helper modules containing reused CMake logic (warnings, sanitizers…)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Normal mode
include(CompileOptions)

# BLAS/LAPACK
if(USE_SYSTEM_BLAS)
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)
else()
  include(FetchOpenBLAS)
  if(NOT TARGET BLAS::BLAS)
    add_library(BLAS::BLAS INTERFACE IMPORTED)
    target_link_libraries(BLAS::BLAS INTERFACE OpenBLAS::openblas)
  endif()
  if(NOT TARGET LAPACK::LAPACK)
    add_library(LAPACK::LAPACK ALIAS BLAS::BLAS)
  endif()
endif()

# PETSc
if(USE_SYSTEM_PETSC)
  find_package(PETSc CONFIG REQUIRED)
else()
  include(FetchPETSc)

  file(MAKE_DIRECTORY "${_petsc_install}/include" "${_petsc_install}/lib")
  add_library(PETSC::petsc INTERFACE IMPORTED GLOBAL)

  set_target_properties(PETSC::petsc PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
                                                "${_petsc_install}/include")

  if(WIN32)
    set(_petsc_lib "${_petsc_install}/lib/petsc.lib")
  elseif(APPLE)
    set(_petsc_lib "${_petsc_install}/lib/libpetsc.dylib")
  else()
    set(_petsc_lib
        "${_petsc_install}/lib/libpetsc${CMAKE_SHARED_LIBRARY_SUFFIX}")
  endif()

  set_target_properties(PETSC::petsc PROPERTIES IMPORTED_LOCATION
                                                "${_petsc_lib}")

  if(TARGET petsc)
    add_dependencies(PETSC::petsc petsc)
  endif()

endif()

# HDF5
if(USE_SYSTEM_HDF5)
  find_package(HDF5 CONFIG REQUIRED)
else()
  include(FetchHDF5)
  if(NOT TARGET HDF5::hdf5 AND TARGET hdf5)
    add_library(HDF5::hdf5 ALIAS hdf5)
  endif()
endif()

# CGNS
if(USE_SYSTEM_CGNS)
  find_package(CGNS CONFIG REQUIRED)
else()
  include(FetchCGNS)
  if(NOT TARGET CGNS::cgns)
    if(TARGET cgns_shared)
      add_library(CGNS::cgns ALIAS cgns_shared)
    elseif(TARGET cgns_static)
      add_library(CGNS::cgns ALIAS cgns_static)
    endif()
  endif()
endif()

# yaml-cpp
if(USE_SYSTEM_YAML_CPP)
  find_package(yaml-cpp CONFIG REQUIRED)
  # Normalize to the same alias we use everywhere
  if(NOT TARGET YAML_CPP::yaml-cpp)
    if(TARGET yaml-cpp)
      add_library(YAML_CPP::yaml-cpp ALIAS yaml-cpp)
    elseif(TARGET yaml-cpp::yaml-cpp)
      add_library(YAML_CPP::yaml-cpp ALIAS yaml-cpp::yaml-cpp)
    endif()
  endif()
else()
  include(FetchYamlCpp)
  # The Fetch module defines a target (yaml-cpp or yaml-cpp::yaml-cpp); expose a
  # consistent alias for consumers.
  if(NOT TARGET YAML_CPP::yaml-cpp)
    if(TARGET yaml-cpp)
      add_library(YAML_CPP::yaml-cpp ALIAS yaml-cpp)
    elseif(TARGET yaml-cpp::yaml-cpp)
      add_library(YAML_CPP::yaml-cpp ALIAS yaml-cpp::yaml-cpp)
    endif()
  endif()
endif()

# Only pull Catch2 when unit tests are enabled
if(BUILD_TESTS AND NOT USE_SYSTEM_CATCH2)
  include(FetchCatch2)
  if(TARGET Catch2WithMain AND NOT TARGET Catch2::Catch2WithMain)
    add_library(Catch2::Catch2WithMain ALIAS Catch2WithMain)
  endif()
endif()

# Source tree
add_subdirectory(src)

# Optional test & example trees (live *outside* src/)
if(BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
  add_custom_target(
    test_unit
    COMMAND ${CMAKE_CTEST_COMMAND} -L unit --output-on-failure
    COMMENT "Run Catch2 unit tests")

  add_custom_target(
    test_mpi
    COMMAND ${CMAKE_CTEST_COMMAND} -L mpi --output-on-failure
    COMMENT "Run MPI functional tests")

  add_custom_target(
    test_perf
    COMMAND ${CMAKE_CTEST_COMMAND} -L perf --output-on-failure
    COMMENT "Run micro-benchmarks (perf)")
endif()

if(BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

# Developer documentation toggle
option(BUILD_DOCS "Build developer documentation" OFF)
if(BUILD_DOCS)
  find_package(Doxygen REQUIRED)
  find_program(SPHINX_BUILD sphinx-build REQUIRED)

  set(DOCS_SRC ${CMAKE_SOURCE_DIR}/docs)
  set(DOCS_BIN ${CMAKE_BINARY_DIR}/docs)

  configure_file(${DOCS_SRC}/Doxyfile.in ${DOCS_BIN}/Doxyfile @ONLY)

  add_custom_target(
    docs ALL
    COMMAND ${SPHINX_BUILD} -b html ${DOCS_SRC} ${DOCS_BIN}/html
    DEPENDS ${DOCS_BIN}/Doxyfile
    COMMENT "Building Sphinx docs")
endif()

# Attach warnings to everything under src/ and tests/ (and nowhere else)
project_enable_warnings_under("${CMAKE_SOURCE_DIR}/src"
                              "${CMAKE_SOURCE_DIR}/tests")
