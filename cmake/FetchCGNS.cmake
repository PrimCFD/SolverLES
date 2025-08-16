include(FetchContent)
include(FetchUtils)

if(NOT DEFINED FETCHCONTENT_UPDATES_DISCONNECTED)
  set(FETCHCONTENT_UPDATES_DISCONNECTED
      ON
      CACHE BOOL "Disable URL re-checks for third-party archives")
endif()

set(CGNS_VERSION 4.4.0)
set(CGNS_TARBALL "cgns-${CGNS_VERSION}.tar.gz") # logical name for resolve_url
set(_cgns_fallback_url
    "https://github.com/CGNS/CGNS/archive/refs/tags/v${CGNS_VERSION}.tar.gz")

resolve_url(${CGNS_TARBALL} _cgns_url ${_cgns_fallback_url})

set(_cgns_sd_arg)
if(DEFINED FETCHCONTENT_SOURCE_DIR_CGNS AND FETCHCONTENT_SOURCE_DIR_CGNS)
  set(_cgns_sd_arg SOURCE_DIR "${FETCHCONTENT_SOURCE_DIR_CGNS}")
endif()

set(_cgns_install "${FETCHCONTENT_BASE_DIR}/cgns-install")

FetchContent_Declare(
  cgns
  URL "${_cgns_url}"
      ${_cgns_sd_arg}
      CMAKE_ARGS
      -DCMAKE_INSTALL_PREFIX=${_cgns_install}
      -DBUILD_SHARED_LIBS=ON
      -DBUILD_TESTING=OFF
      -DCGNS_BUILD_TESTING=OFF
      -DCGNS_BUILD_TOOLS=OFF
      -DCGNS_BUILD_CGNSTOOLS=OFF
      -DCGNS_BUILD_EXAMPLES=OFF
      -DCGNS_ENABLE_HDF5=ON
      -DCGNS_ENABLE_HDF5_HL=ON
      -DCGNS_ENABLE_PARALLEL=${ENABLE_MPI}
      -DHDF5_DIR=${HDF5_DIR}
      -DCMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH})

if(PREFETCH_THIRD_PARTY)
  FetchContent_Populate(cgns)
  return()
endif()

FetchContent_MakeAvailable(cgns)

list(APPEND CMAKE_PREFIX_PATH "${_cgns_install}")
set(_cgns_install
    "${_cgns_install}"
    PARENT_SCOPE)
