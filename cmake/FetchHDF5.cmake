cmake_minimum_required(VERSION 3.24)

include(FetchContent)
include(ExternalProject)
include(${CMAKE_SOURCE_DIR}/cmake/FetchUtils.cmake)

fu_dep_base(_DEPS_BASE)

# Layout
set(_hdf5_prefix "${_DEPS_BASE}/hdf5") # EP root
set(_hdf5_build "${_hdf5_prefix}/build") # out-of-source build
set(_hdf5_install "${_DEPS_BASE}/hdf5-install") # install prefix

# Source (prefer extern tarball)
set(_hdf5_url
    "https://github.com/HDFGroup/hdf5/archive/refs/tags/hdf5-1_14_3.tar.gz")
fu_url_args("hdf5" "${_hdf5_url}" _H_URL_ARGS)

# Prefetch-only mode: populate sources and stop
if(PREFETCH_THIRD_PARTY)
  FetchContent_Declare(hdf5 ${_H_URL_ARGS})
  FetchContent_Populate(hdf5)
  add_custom_target(hdf5_src COMMENT "HDF5 sources are populated")
  _fu_register_prefetch_target(hdf5_src)
  return()
endif()

# Configure options – lean build, no tests/tools/examples
set(_hdf5_cmake_opts
    -DBUILD_TESTING=OFF
    -DHDF5_BUILD_EXAMPLES=OFF
    -DHDF5_BUILD_TOOLS=OFF
    -DHDF5_BUILD_CPP_LIB=OFF
    -DHDF5_ENABLE_Z_LIB_SUPPORT=ON
    -DHDF5_ENABLE_SZIP_SUPPORT=OFF
    -DBUILD_SHARED_LIBS=ON
    -DCMAKE_MESSAGE_LOG_LEVEL=WARNING
    -Wno-dev
    -DCMAKE_RULE_MESSAGES=OFF)

ExternalProject_Add(
  hdf5_ep
  ${_H_URL_ARGS}
  PREFIX "${_hdf5_prefix}"
  UPDATE_COMMAND ""
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${_hdf5_install}
             -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} ${_hdf5_cmake_opts}
  BINARY_DIR "${_hdf5_build}"
  INSTALL_DIR "${_hdf5_install}"
  UPDATE_DISCONNECTED ON
  BUILD_BYPRODUCTS "${_hdf5_install}/lib/libhdf5${CMAKE_SHARED_LIBRARY_SUFFIX}"
  LOG_DOWNLOAD 1
  LOG_UPDATE 1
  LOG_CONFIGURE 1
  LOG_BUILD 1
  LOG_INSTALL 1)

# Create an imported target HDF5::hdf5
add_library(HDF5::hdf5 SHARED IMPORTED GLOBAL)
if(WIN32)
  set(_hdf5_lib "${_hdf5_install}/bin/hdf5.dll")
elseif(APPLE)
  set(_hdf5_lib "${_hdf5_install}/lib/libhdf5.dylib")
else()
  set(_hdf5_lib "${_hdf5_install}/lib/libhdf5${CMAKE_SHARED_LIBRARY_SUFFIX}")
endif()
set_target_properties(
  HDF5::hdf5
  PROPERTIES IMPORTED_LOCATION "${_hdf5_lib}" INTERFACE_INCLUDE_DIRECTORIES
                                              "${_hdf5_install}/include")
add_dependencies(HDF5::hdf5 hdf5_ep)

# Make package config discoverable for consumers: find_package(HDF5 CONFIG)
set(HDF5_DIR
    "${_hdf5_install}/cmake"
    CACHE PATH "" FORCE)
list(APPEND CMAKE_PREFIX_PATH "${_hdf5_install}")
