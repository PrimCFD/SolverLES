cmake_minimum_required(VERSION 3.24)
project(Prefetch NONE)

get_filename_component(_prefetch_root "${CMAKE_CURRENT_LIST_DIR}" ABSOLUTE)
list(APPEND CMAKE_MODULE_PATH "${_prefetch_root}")

# Network hygiene — disconnected by default
if(NOT DEFINED FETCHCONTENT_UPDATES_DISCONNECTED)
  set(FETCHCONTENT_UPDATES_DISCONNECTED
      ON
      CACHE BOOL "Disable URL re‑checks")
endif()

# Bring in every third‑party component
include(${CMAKE_SOURCE_DIR}/cmake/FetchHDF5.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/FetchCGNS.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/FetchPETSc.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/FetchCatch2.cmake)

add_custom_target(prefetch-archives)
foreach(pkg catch2 cgns hdf5 petsc)
  if(TARGET ${pkg}-populate)
    add_dependencies(prefetch-archives ${pkg}-populate)
  endif()
endforeach()
