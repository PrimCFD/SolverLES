include(FetchContent)
include(FetchUtils)

if(NOT DEFINED FETCHCONTENT_UPDATES_DISCONNECTED)
  set(FETCHCONTENT_UPDATES_DISCONNECTED
      ON
      CACHE BOOL "Disable URL re-checks for third-party archives")
endif()

set(HDF5_VERSION 1_14_3)
set(HDF5_TARBALL "hdf5-${HDF5_VERSION}.tar.gz")
set(_hdf5_fallback_url
    "https://github.com/HDFGroup/hdf5/archive/refs/tags/${HDF5_TARBALL}")

resolve_url(${HDF5_TARBALL} _hdf5_url ${_hdf5_fallback_url})

set(_hdf5_sd_arg)
if(DEFINED FETCHCONTENT_SOURCE_DIR_HDF5 AND FETCHCONTENT_SOURCE_DIR_HDF5)
  set(_hdf5_sd_arg SOURCE_DIR "${FETCHCONTENT_SOURCE_DIR_HDF5}")
endif()

set(_hdf5_install "${FETCHCONTENT_BASE_DIR}/hdf5-install")
set(_hdf5_ext_build "${FETCHCONTENT_BASE_DIR}/hdf5-ext-build")

FetchContent_Declare(hdf5 URL "${_hdf5_url}" ${_hdf5_sd_arg})

if(PREFETCH_THIRD_PARTY)
  FetchContent_Populate(hdf5)
  return()
endif()

FetchContent_Populate(hdf5)
FetchContent_GetProperties(hdf5)

set(_hdf5_config_path "${_hdf5_install}/lib/cmake/hdf5/hdf5-config.cmake")
if(NOT EXISTS "${_hdf5_config_path}")
  file(MAKE_DIRECTORY "${_hdf5_ext_build}")

  set(_hdf5_args
      -DCMAKE_INSTALL_PREFIX=${_hdf5_install}
      -DBUILD_SHARED_LIBS=ON
      -DHDF5_BUILD_HL_LIB=ON
      -DHDF5_BUILD_TOOLS=ON
      -DHDF5_BUILD_EXAMPLES=OFF
      -DHDF5_ENABLE_Z_LIB_SUPPORT=OFF
      -DHDF5_ENABLE_SZIP_SUPPORT=OFF
      -DHDF5_BUILD_TESTING=OFF
      -DBUILD_TESTING=OFF
      -DHDF5_ENABLE_PARALLEL=${ENABLE_MPI})

  if(CMAKE_TOOLCHAIN_FILE)
    list(APPEND _hdf5_args -DCMAKE_TOOLCHAIN_FILE:PATH=${CMAKE_TOOLCHAIN_FILE})
  endif()
  if(CMAKE_C_COMPILER)
    list(APPEND _hdf5_args -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER})
  endif()
  if(CMAKE_CXX_COMPILER)
    list(APPEND _hdf5_args -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER})
  endif()
  if(CMAKE_Fortran_COMPILER)
    list(APPEND _hdf5_args
         -DCMAKE_Fortran_COMPILER:FILEPATH=${CMAKE_Fortran_COMPILER})
  endif()
  if(CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    list(APPEND _hdf5_args -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE})
  endif()

  execute_process(
    COMMAND ${CMAKE_COMMAND} -S ${hdf5_SOURCE_DIR} -B ${_hdf5_ext_build}
            ${_hdf5_args} RESULT_VARIABLE _hdf5_cfg_rv)
  if(NOT _hdf5_cfg_rv EQUAL 0)
    message(FATAL_ERROR "Configuring HDF5 failed (rv=${_hdf5_cfg_rv})")
  endif()

  set(_hdf5_build_cmd ${CMAKE_COMMAND} --build ${_hdf5_ext_build} --target
                      install)
  if(CMAKE_CONFIGURATION_TYPES)
    list(APPEND _hdf5_build_cmd --config Release)
  endif()
  execute_process(COMMAND ${_hdf5_build_cmd} RESULT_VARIABLE _hdf5_bld_rv)
  if(NOT _hdf5_bld_rv EQUAL 0)
    message(FATAL_ERROR "Building/Installing HDF5 failed (rv=${_hdf5_bld_rv})")
  endif()
endif()

set(HDF5_DIR
    "${_hdf5_install}/lib/cmake/hdf5"
    CACHE PATH "" FORCE)
set(HDF5_ROOT
    "${_hdf5_install}"
    CACHE PATH "" FORCE)
list(APPEND CMAKE_PREFIX_PATH "${_hdf5_install}")

set(_hdf5_install
    "${_hdf5_install}"
    PARENT_SCOPE)
set(_hdf5_ext_build
    "${_hdf5_ext_build}"
    PARENT_SCOPE)
