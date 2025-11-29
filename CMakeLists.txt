# Root listfile – only global policies, dependencies, and one add_subdirectory

cmake_minimum_required(VERSION 3.24)

# Prefetch mode (nothing else gets configured)
if(PREFETCH_THIRD_PARTY)
  message(STATUS "Prefetch mode: downloading archives and exiting")
  include(cmake/PrefetchDependencies.cmake)
  return()
endif()

project(
  KolmoPlas
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

find_package(MPI REQUIRED)

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

# Third–party dependencies that every target needs

# OpenMP System-only OpenMP (no vendoring)
find_package(OpenMP REQUIRED) # expects libgomp/libomp present on the node

# Prove it actually links
include(CheckCXXSourceCompiles)
set(CMAKE_REQUIRED_LIBRARIES OpenMP::OpenMP_CXX)
check_cxx_source_compiles(
  "
  #include <omp.h>
  int main(){
    int n=0;
    #pragma omp parallel reduction(+:n)
    for(int i=0;i<32;++i) n+=i;
    return n==496 ? 0 : 1;
  }"
  OMP_LINKS_OK)

if(NOT OMP_LINKS_OK)
  message(
    FATAL_ERROR
      "OpenMP was detected but failed to link. On this cluster, load a toolchain module, e.g.:
       module load gcc/12   # libgomp
     or
       module load llvm/18  # libomp
     Then re-configure.")
endif()

# Wire flags like your other deps
if(OpenMP_CXX_FOUND)
  message(STATUS "OpenMP enabled (system runtime)")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()
if(OpenMP_Fortran_FOUND)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
endif()

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

  # Create a real imported library target for libpetsc.so so we can use
  # TARGET_FILE_DIR
  add_library(PETSC::petsc_lib SHARED IMPORTED GLOBAL)
  set_target_properties(PETSC::petsc_lib PROPERTIES IMPORTED_LOCATION
                                                    "${_petsc_lib}")

  # Publish PETSc’s common deps via the INTERFACE target you already use
  set(_petsc_public_deps "")
  list(APPEND _petsc_public_deps MPI::MPI_C MPI::MPI_CXX)
  list(APPEND _petsc_public_deps BLAS::BLAS LAPACK::LAPACK)

  set_target_properties(
    PETSC::petsc PROPERTIES INTERFACE_LINK_LIBRARIES
                            "PETSC::petsc_lib;${_petsc_public_deps}")

  if(TARGET petsc)
    add_dependencies(PETSC::petsc petsc)
  endif()

endif()

# HDF5
if(USE_SYSTEM_HDF5)
  find_package(HDF5 CONFIG REQUIRED)
else()
  include(FetchHDF5)

  if(DEFINED _hdf5_install)
    file(MAKE_DIRECTORY "${_hdf5_install}/include" "${_hdf5_install}/lib")
  endif()

  if(NOT TARGET HDF5::hdf5 AND TARGET hdf5)
    add_library(HDF5::hdf5 ALIAS hdf5)
  endif()
endif()

# CGNS
if(USE_SYSTEM_CGNS)
  find_package(CGNS CONFIG REQUIRED)
else()
  include(FetchCGNS)

  if(DEFINED _cgns_install)
    file(MAKE_DIRECTORY "${_cgns_install}/include" "${_cgns_install}/lib")
  endif()

  if(NOT TARGET CGNS::cgns)
    if(TARGET cgns_shared)
      add_library(CGNS::cgns ALIAS cgns_shared)
    elseif(TARGET cgns_static)
      add_library(CGNS::cgns ALIAS cgns_static)
    endif()
  endif()

  # If the fetch module defines an EP target called 'cgns', make the alias
  # depend on it.
  if(TARGET cgns AND TARGET CGNS::cgns)
    add_dependencies(CGNS::cgns cgns)
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
