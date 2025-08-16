include(FetchContent)
include(FetchUtils)

if(NOT DEFINED FETCHCONTENT_UPDATES_DISCONNECTED)
  set(FETCHCONTENT_UPDATES_DISCONNECTED
      ON
      CACHE BOOL "Disable URL re‑checks")
endif()

set(PETSC_VERSION 3.23.5)
set(PETSC_TAR petsc-${PETSC_VERSION}.tar.gz)

set(_petsc_primary_url
    "https://web.cels.anl.gov/projects/petsc/download/release-snapshots/${PETSC_TAR}"
)
resolve_url(${PETSC_TAR} _petsc_url ${_petsc_primary_url})

# network fallback
if(_petsc_url STREQUAL ${_petsc_primary_url}
   AND NOT FETCHCONTENT_FULLY_DISCONNECTED)
  set(_petsc_fallback_url "https://fossies.org/linux/misc/${PETSC_TAR}")
endif()

set(_petsc_sha256
    b0bb614dfbf36c286c8cad30912fe77359dccbf6b65a5edd1dde82af293f21fc)
set(_petsc_install ${CMAKE_BINARY_DIR}/_deps/petsc-install)

FetchContent_Declare(
  petsc
  URL ${_petsc_url}
  URL_HASH SHA256=${_petsc_sha256}
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${_petsc_install})

if(PREFETCH_THIRD_PARTY)
  FetchContent_Populate(petsc)
  return()
endif()

FetchContent_MakeAvailable(petsc)

list(APPEND CMAKE_PREFIX_PATH "${_petsc_install}")
