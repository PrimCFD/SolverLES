cmake_minimum_required(VERSION 3.24)

include(FetchContent)
include(ExternalProject)
include(${CMAKE_SOURCE_DIR}/cmake/FetchUtils.cmake)

find_program(MAKE_EXECUTABLE NAMES make gmake REQUIRED)

fu_dep_base(_DEPS_BASE)
set(_petsc_prefix "${_DEPS_BASE}/petsc")
set(_petsc_install
    "${_DEPS_BASE}/petsc-install"
    CACHE PATH "")

# Upstream tarball (unchanged)
set(_petsc_url
    "https://gitlab.com/petsc/petsc/-/archive/v3.20.0/petsc-v3.20.0.tar.gz")
fu_url_args("petsc" "${_petsc_url}" _P_URL_ARGS)

# Prefetch-only
if(PREFETCH_THIRD_PARTY)
  FetchContent_Declare(petsc ${_P_URL_ARGS})
  FetchContent_Populate(petsc)
  add_custom_target(petsc_src COMMENT "PETSc sources are populated")
  _fu_register_prefetch_target(petsc_src)
  return()
endif()

# BLAS/LAPACK: prefer built OpenBLAS (module) unless USE_SYSTEM_BLAS is set by
# the root
set(_petsc_deps)
if(TARGET OpenBLAS::openblas)
  set(_blaslib OpenBLAS::openblas)
  list(APPEND _petsc_deps OpenBLAS::openblas)
elseif(TARGET BLAS::BLAS AND TARGET LAPACK::LAPACK)
  # System BLAS/LAPACK already discovered by the root listfile
  set(_blaslib BLAS::BLAS)
else()
  message(
    FATAL_ERROR
      "No BLAS/LAPACK available. Include FetchOpenBLAS.cmake or set USE_SYSTEM_BLAS=ON."
  )
endif()

# Figure out the actual library path to pass PETSc's configure
get_target_property(_blas_libloc "${_blaslib}" IMPORTED_LOCATION)
if(NOT _blas_libloc)
  # BLAS::BLAS/LAPACK::LAPACK are INTERFACE; try common locations from
  # CMAKE_PREFIX_PATH (best-effort; usually OpenBLAS::openblas gives us the
  # path)
  set(_blas_libloc "")
endif()

# MPI toggle
if(DEFINED ENABLE_MPI AND ENABLE_MPI)
  set(_with_mpi "--with-mpi=1")
else()
  set(_with_mpi "--with-mpi=0")
endif()

# Construct configure command
set(_cfg_cmd ./configure "--prefix=${_petsc_install}" "${_with_mpi}")
if(_blas_libloc)
  list(APPEND _cfg_cmd "--with-blaslapack-lib=${_blas_libloc}")
else()
  # Let PETSc do its own BLAS/LAPACK discovery via pkg-config / system search
  list(APPEND _cfg_cmd "--with-blaslapack=1")
endif()

# Build & install
set(_build_cmd ${MAKE_EXECUTABLE})
set(_install_cmd ${MAKE_EXECUTABLE} install)

ExternalProject_Add(
  petsc
  ${_P_URL_ARGS}
  PREFIX "${_petsc_prefix}"
  SOURCE_SUBDIR "" # top-level configure
  UPDATE_COMMAND ""
  CONFIGURE_COMMAND ${_cfg_cmd}
  BUILD_COMMAND ${_build_cmd}
  INSTALL_COMMAND ${_install_cmd}
  BUILD_IN_SOURCE 1
  UPDATE_DISCONNECTED ON
  DEPENDS ${_petsc_deps}
  BUILD_BYPRODUCTS
    "${_petsc_install}/lib/libpetsc${CMAKE_SHARED_LIBRARY_SUFFIX}")

# Expose for find_package(PETSc CONFIG) if a config is generated
set(PETSC_DIR
    "${_petsc_install}"
    CACHE PATH "" FORCE)
list(APPEND CMAKE_PREFIX_PATH "${_petsc_install}")
