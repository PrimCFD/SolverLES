include(FetchContent)
include(ExternalProject)
include(${CMAKE_SOURCE_DIR}/cmake/FetchUtils.cmake)

fu_dep_base(_DEPS_BASE)
set(_cgns_prefix "${_DEPS_BASE}/cgns")
set(_cgns_build "${_cgns_prefix}/build")
set(_cgns_install "${_DEPS_BASE}/cgns-install")

# Prefer extern tarball
set(_cgns_url "https://github.com/CGNS/CGNS/archive/refs/tags/v4.4.0.tar.gz")
fu_url_args("cgns" "${_cgns_url}" _C_URL_ARGS)

# Prefetch-only
if(PREFETCH_THIRD_PARTY)
  FetchContent_Declare(cgns ${_C_URL_ARGS})
  FetchContent_Populate(cgns)
  add_custom_target(cgns_src COMMENT "CGNS sources are populated")
  _fu_register_prefetch_target(cgns_src)
  return()
endif()

# Require HDF5::hdf5 from our HDF5 module (or system)
if(NOT TARGET HDF5::hdf5)
  find_package(HDF5 REQUIRED)
endif()

# Pass HDF5 into CGNS configure and ensure consistent detection
set(_cgns_cmake_opts
    -DBUILD_SHARED_LIBS=ON
    -DCGNS_ENABLE_HDF5=ON
    -DHDF5_DIR=${HDF5_DIR}
    -DHDF5_ROOT=${HDF5_DIR}
    -DCGNS_BUILD_SHARED=ON
    -DCGNS_ENABLE_FORTRAN=OFF
    -DCGNS_ENABLE_TESTS=OFF
    -DCGNS_BUILD_EXAMPLES=OFF)

ExternalProject_Add(
  cgns_ep
  ${_C_URL_ARGS}
  PREFIX "${_cgns_prefix}"
  UPDATE_COMMAND ""
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${_cgns_install}
             -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} ${_cgns_cmake_opts}
  BINARY_DIR "${_cgns_build}"
  INSTALL_DIR "${_cgns_install}"
  UPDATE_DISCONNECTED ON
  DEPENDS HDF5::hdf5
  BUILD_BYPRODUCTS "${_cgns_install}/lib/libcgns${CMAKE_SHARED_LIBRARY_SUFFIX}")

# Imported target CGNS::cgns
add_library(CGNS::cgns SHARED IMPORTED GLOBAL)
if(WIN32)
  set(_cgns_lib "${_cgns_install}/bin/cgns.dll")
elseif(APPLE)
  set(_cgns_lib "${_cgns_install}/lib/libcgns.dylib")
else()
  set(_cgns_lib "${_cgns_install}/lib/libcgns${CMAKE_SHARED_LIBRARY_SUFFIX}")
endif()
set_target_properties(
  CGNS::cgns
  PROPERTIES IMPORTED_LOCATION "${_cgns_lib}" INTERFACE_INCLUDE_DIRECTORIES
                                              "${_cgns_install}/include")
add_dependencies(CGNS::cgns cgns_ep)

# Help find_package(CGNS CONFIG) downstream if desired
set(CGNS_DIR
    "${_cgns_install}/lib/cmake/CGNS"
    CACHE PATH "" FORCE)
list(APPEND CMAKE_PREFIX_PATH "${_cgns_install}")
