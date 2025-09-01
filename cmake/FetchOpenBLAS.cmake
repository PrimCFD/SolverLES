include(FetchContent)
include(ExternalProject)
include(${CMAKE_SOURCE_DIR}/cmake/FetchUtils.cmake)

find_program(MAKE_EXECUTABLE NAMES make gmake REQUIRED)

# Version and URL (unchanged)
if(NOT DEFINED OPENBLAS_VERSION)
  set(OPENBLAS_VERSION 0.3.30)
endif()
set(OPENBLAS_URL
    "https://github.com/xianyi/OpenBLAS/archive/refs/tags/v${OPENBLAS_VERSION}.tar.gz"
    CACHE STRING "OpenBLAS source URL")

fu_dep_base(_DEPS_BASE)
set(_openblas_prefix "${_DEPS_BASE}/openblas")
set(_openblas_install "${_DEPS_BASE}/openblas-install")

# Prefer extern tarball, fall back to online
fu_url_args("openblas" "${OPENBLAS_URL}" _OB_URL_ARGS)

# Prefetch-only: just populate sources and exit
if(PREFETCH_THIRD_PARTY)
  FetchContent_Declare(openblas ${_OB_URL_ARGS})
  FetchContent_Populate(openblas)
  add_custom_target(openblas_src COMMENT "OpenBLAS sources are populated")
  _fu_register_prefetch_target(openblas_src)
  return()
endif()

# Makefile build knobs - Build SHARED only (NO_STATIC=1) so downstream can
# reliably link a .so/.dylib/.dll - Do NOT build or run tests/benchmarks -
# DYNAMIC_ARCH=1 gives broad CPU coverage on x86_64
set(_ob_make_opts)
list(
  APPEND
  _ob_make_opts
  "USE_OPENMP=0"
  "DYNAMIC_ARCH=1"
  "NO_FORTRAN=0" # needed to also get LAPACK
  "NO_STATIC=1" # we want a shared library artifact
  "NO_SHARED=0"
  "NO_TEST=1"
  "NO_CBLAS_TEST=1"
  "NO_LAPACK_TEST=1"
  "NO_BENCHMARK=1"
  "BINARY=64")

# Build only what we need: libraries + the shared library target. No tests will
# be invoked.
ExternalProject_Add(
  openblas_ep
  ${_OB_URL_ARGS}
  PREFIX "${_openblas_prefix}"
  UPDATE_COMMAND ""
  CONFIGURE_COMMAND ""
  BUILD_IN_SOURCE 1
  BUILD_COMMAND ${MAKE_EXECUTABLE} -s ${_ob_make_opts} libs shared
  INSTALL_COMMAND ${MAKE_EXECUTABLE} -s ${_ob_make_opts} PREFIX=<INSTALL_DIR>
                  install
  LOG_DOWNLOAD 1
  LOG_UPDATE 1
  LOG_CONFIGURE 1
  LOG_BUILD 1
  LOG_INSTALL 1
  INSTALL_DIR "${_openblas_install}"
  UPDATE_DISCONNECTED ON
  BUILD_BYPRODUCTS
    "${_openblas_install}/lib/libopenblas${CMAKE_SHARED_LIBRARY_SUFFIX}")

# Imported target: prefer the canonical soname (libopenblas.so / .dylib / .dll)
# We install to <prefix>/lib on Linux/macOS. (Windows Makefile builds are
# uncommon; adjust if needed.)
set(_ob_lib
    "${_openblas_install}/lib/libopenblas${CMAKE_SHARED_LIBRARY_SUFFIX}")

# If a distro emits only a fully versioned filename, fall back to first matching
# .so* / .dylib* / .dll*
if(NOT EXISTS "${_ob_lib}")
  if(WIN32)
    file(GLOB _ob_candidates "${_openblas_install}/lib/libopenblas*.dll")
  elseif(APPLE)
    file(GLOB _ob_candidates "${_openblas_install}/lib/libopenblas*.dylib*")
  else()
    file(GLOB _ob_candidates "${_openblas_install}/lib/libopenblas*.so*")
  endif()
  list(LENGTH _ob_candidates _n)
  if(_n GREATER 0)
    list(GET _ob_candidates 0 _ob_lib)
  endif()
endif()

# Final safety: error early if we still haven't got a shared lib path
if(NOT _ob_lib)
  message(
    FATAL_ERROR
      "OpenBLAS shared library not found under ${_openblas_install}/lib after build."
  )
endif()

# Define the imported target (supports .so/.dylib/.dll)
add_library(OpenBLAS::OpenBLAS UNKNOWN IMPORTED GLOBAL)
set_target_properties(
  OpenBLAS::OpenBLAS
  PROPERTIES IMPORTED_LOCATION "${_ob_lib}" INTERFACE_INCLUDE_DIRECTORIES
                                            "${_openblas_install}/include")
add_dependencies(OpenBLAS::OpenBLAS openblas_ep)

# Back-compat alias some code expects
add_library(OpenBLAS::openblas ALIAS OpenBLAS::OpenBLAS)

# Provide BLAS/LAPACK interface targets if not already defined
if(NOT TARGET BLAS::BLAS)
  add_library(BLAS::BLAS INTERFACE IMPORTED)
  target_link_libraries(BLAS::BLAS INTERFACE OpenBLAS::OpenBLAS)
endif()
if(NOT TARGET LAPACK::LAPACK)
  add_library(LAPACK::LAPACK INTERFACE IMPORTED)
  target_link_libraries(LAPACK::LAPACK INTERFACE OpenBLAS::OpenBLAS)
endif()

# Convenience cache vars used elsewhere
set(OpenBLAS_INCLUDE_DIRS
    "${_openblas_install}/include"
    CACHE PATH "" FORCE)
set(OpenBLAS_LIBRARIES
    "${_ob_lib}"
    CACHE FILEPATH "" FORCE)
list(APPEND CMAKE_PREFIX_PATH "${_openblas_install}")
