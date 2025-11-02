include(FetchContent)
include(ExternalProject)
include(${CMAKE_SOURCE_DIR}/cmake/FetchUtils.cmake)

find_program(MPIEXEC_EXECUTABLE NAMES mpiexec mpirun srun)
find_program(MAKE_EXECUTABLE NAMES make gmake REQUIRED)

fu_dep_base(_DEPS_BASE)
set(_petsc_prefix "${_DEPS_BASE}/petsc")
set(_petsc_install
    "${_DEPS_BASE}/petsc-install"
    CACHE PATH "")

# Upstream tarball (unchanged)
set(_petsc_url
    "https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-3.24.0.tar.gz"
)
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

# --- Derive PETSc debugging toggle & optimization flags from the superbuild
string(TOLOWER "${CMAKE_BUILD_TYPE}" _bt)
if(_bt STREQUAL "debug")
  set(_petsc_debug 1)
else()
  set(_petsc_debug 0)
endif()

# Pull Release flags from the toolchain so PETSc pkg-config reflects
# -O3/-DNDEBUG (fallbacks keep things sensible if variables are empty)
set(_copt "${CMAKE_C_FLAGS_RELEASE}")
if(_copt STREQUAL "")
  set(_copt "-O3 -DNDEBUG")
endif()
set(_cxxopt "${CMAKE_CXX_FLAGS_RELEASE}")
if(_cxxopt STREQUAL "")
  set(_cxxopt "-O3 -DNDEBUG")
endif()
set(_fopt "${CMAKE_Fortran_FLAGS_RELEASE}")
if(_fopt STREQUAL "")
  set(_fopt "-O3 -DNDEBUG")
endif()

# Location where prefetch caches PETSc third-party tarballs (used offline)
set(PETSC_PKG_CACHE
    "${CMAKE_SOURCE_DIR}/extern/.petsc-downloads"
    CACHE PATH "Local cache of PETSc package tarballs (prefetched)")

# Construct configure command
set(_cfg_cmd ./configure "--prefix=${_petsc_install}" "--with-mpi=1"
             "--with-debugging=${_petsc_debug}")

# If we can find an MPI launcher, hint PETSc so it doesn't fall back to MPIUNI
if(MPIEXEC_EXECUTABLE)
  list(APPEND _cfg_cmd "--with-mpiexec=${MPIEXEC_EXECUTABLE}")
endif()

# Always tell PETSc where to find locally prefetched tarballs (works
# online/offline)
list(APPEND _cfg_cmd "--with-packages-download-dir=${PETSC_PKG_CACHE}")

if(_blas_libloc)
  list(APPEND _cfg_cmd "--with-blaslapack-lib=${_blas_libloc}")
else()
  # Let PETSc do its own BLAS/LAPACK discovery via pkg-config / system search
  list(APPEND _cfg_cmd "--with-blaslapack=1")
endif()

# --- Enable parallel LU (MUMPS + deps) downloads when needed PETSc will
# fetch/build these if system libs are not provided. See PETSc docs: install +
# MATSOLVERMUMPS.  (Runtime: -pc_type lu -pc_factor_mat_solver_type mumps)
# https://petsc.org/main/install/install/
# https://petsc.org/release/manualpages/Mat/MATSOLVERMUMPS/
list(
  APPEND
  _cfg_cmd
  "--download-mumps"
  "--download-scalapack"
  "--download-blacs"
  "--download-parmetis"
  "--download-metis")

# Optimization flags (propagate superbuild Release flags to PETSc)
list(APPEND _cfg_cmd "COPTFLAGS=${_copt}" "CXXOPTFLAGS=${_cxxopt}"
     "FOPTFLAGS=${_fopt}")

# Vendor-agnostic: only request MPICH if explicitly asked, or if we already
# *have* an mpich tarball in the local cache.
#
# Usage examples: - Default (vendor-agnostic): no -DPETSC_VENDOR_MPI → no
# --download-mpich - Force MPICH:               -DPETSC_VENDOR_MPI=mpich
#
# Also auto-enable if a tarball is available offline to avoid network.
set(PETSC_VENDOR_MPI
    "${PETSC_VENDOR_MPI}"
    CACHE STRING "MPI vendor for PETSc offline configure (empty or 'mpich')")

# Only request MPICH if explicitly asked for.
set(PETSC_VENDOR_MPI
    "${PETSC_VENDOR_MPI}"
    CACHE STRING "MPI vendor for PETSc offline configure (empty or 'mpich')")
if(PETSC_VENDOR_MPI STREQUAL "mpich")
  list(APPEND _cfg_cmd "--download-mpich")
endif()

# --- Allow scripts to inject extra PETSc configure options e.g. export
# PETSC_CONFIGURE_OPTS="--with-debugging=0 --with-mpi=1 COPTFLAGS='-O3
# -march=native' ..."
if(DEFINED ENV{PETSC_CONFIGURE_OPTS} AND NOT "$ENV{PETSC_CONFIGURE_OPTS}"
                                         STREQUAL "")
  separate_arguments(_petsc_cfg_envopts NATIVE_COMMAND
                     "$ENV{PETSC_CONFIGURE_OPTS}")
  list(APPEND _cfg_cmd ${_petsc_cfg_envopts})
endif()

# --- Force MPI wrapper compilers into PETSc configure
set(_cc "$ENV{MPI_C_COMPILER}")
set(_cxx "$ENV{MPI_CXX_COMPILER}")
set(_fc "$ENV{MPI_Fortran_COMPILER}")

# --- Pass wrapper compilers to PETSc via configure args (env is ignored) Only
# append when set, so users can still override in PETSC_CONFIGURE_OPTS.
if(NOT "${_cc}" STREQUAL "")
  list(APPEND _cfg_cmd "CC=${_cc}")
endif()
if(NOT "${_cxx}" STREQUAL "")
  list(APPEND _cfg_cmd "CXX=${_cxx}")
endif()
if(NOT "${_fc}" STREQUAL "")
  list(APPEND _cfg_cmd "FC=${_fc}")
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
  CONFIGURE_COMMAND
  CONFIGURE_COMMAND ${_cfg_cmd}
  BUILD_COMMAND ${_build_cmd}
  INSTALL_COMMAND ${_install_cmd}
  BUILD_IN_SOURCE 1
  UPDATE_DISCONNECTED ON
  DEPENDS ${_petsc_deps}
  BUILD_BYPRODUCTS
    "${_petsc_install}/lib/libpetsc${CMAKE_SHARED_LIBRARY_SUFFIX}")

# Expose cache path for downstream tools
set(PETSC_PKG_CACHE
    "${PETSC_PKG_CACHE}"
    CACHE PATH "" FORCE)

# Expose for find_package(PETSc CONFIG) if a config is generated
set(PETSC_DIR
    "${_petsc_install}"
    CACHE PATH "" FORCE)
list(APPEND CMAKE_PREFIX_PATH "${_petsc_install}")
