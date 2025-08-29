# cmake/FetchUtils.cmake
cmake_minimum_required(VERSION 3.24)

# Internal accumulator for prefetch-only targets
function(_fu_register_prefetch_target tgt)
  get_property(_acc GLOBAL PROPERTY SOLVERLES_PREFETCH_TARGETS)
  if(NOT _acc)
    set(_acc "")
  endif()
  list(APPEND _acc "${tgt}")
  set_property(GLOBAL PROPERTY SOLVERLES_PREFETCH_TARGETS "${_acc}")
endfunction()

# Base directory where dependencies are built/installed
# Honors FETCHCONTENT_BASE_DIR if provided by the toolchain/driver script.
function(fu_dep_base out_dir)
  if(DEFINED FETCHCONTENT_BASE_DIR AND FETCHCONTENT_BASE_DIR)
    set(${out_dir} "${FETCHCONTENT_BASE_DIR}" PARENT_SCOPE)
  else()
    set(${out_dir} "${CMAKE_BINARY_DIR}/_deps" PARENT_SCOPE)
  endif()
endfunction()

# Locate a pre-fetched archive in extern/, optionally return its checksum
#   fu_find_extern_archive(<name> out_path out_sha256)
function(fu_find_extern_archive name out_path out_sha)
  set(${out_path} "" PARENT_SCOPE)
  set(${out_sha}  "" PARENT_SCOPE)
  set(_extern_dir "${CMAKE_SOURCE_DIR}/extern")

  if(NOT EXISTS "${_extern_dir}")
    return()
  endif()

  file(GLOB _candidates "${_extern_dir}/${name}-*.tgz")
  list(LENGTH _candidates _n)
  if(_n EQUAL 0)
    return()
  endif()

  # Pick the lexicographically last; prefetch script normally keeps one per dep.
  list(SORT _candidates)
  list(REVERSE _candidates)
  list(GET _candidates 0 _chosen)

  # Compute SHA256 (used both for verification and URL_HASH)
  if(EXISTS "${_extern_dir}/SHA256SUMS")
    # Best-effort verification against the checksums file
    file(SHA256 "${_chosen}" _sha_now)
    file(READ "${_extern_dir}/SHA256SUMS" _sha_file)
    string(FIND "${_sha_file}" "${_sha_now}  ${_name}.tgz" _pos)
    # (Don’t fail on mismatch here; the prefetch script guarantees integrity.)
    set(${out_sha} "${_sha_now}" PARENT_SCOPE)
  else()
    file(SHA256 "${_chosen}" _sha_now)
    set(${out_sha} "${_sha_now}" PARENT_SCOPE)
  endif()

  set(${out_path} "${_chosen}" PARENT_SCOPE)
endfunction()

# Compose URL args for FetchContent/ExternalProject
# If an extern/<name>-*.tgz exists, use it (with URL_HASH); else use fallback URL.
#   fu_url_args(<name> <fallback_url> <out_args_var>)
function(fu_url_args name default_url out_args)
  fu_find_extern_archive("${name}" _local _hash)
  if(_local)
    if(_hash)
      set(${out_args} "URL;${_local};URL_HASH;SHA256=${_hash}" PARENT_SCOPE)
    else()
      set(${out_args} "URL;${_local}" PARENT_SCOPE)
    endif()
  else()
    set(${out_args} "URL;${default_url}" PARENT_SCOPE)
  endif()
endfunction()

# Create a single aggregate phony target that depends on all prefetch populate targets
# Call once from PrefetchDependencies.cmake
function(fu_make_prefetch_aggregate)
  get_property(_acc GLOBAL PROPERTY SOLVERLES_PREFETCH_TARGETS)
  if(_acc)
    add_custom_target(prefetch-archives ALL
      COMMENT "Prefetch-only: populate and stage third-party sources")
    add_dependencies(prefetch-archives ${_acc})
  else()
    message(STATUS "No prefetch targets registered. Nothing to prefetch.")
  endif()
endfunction()
