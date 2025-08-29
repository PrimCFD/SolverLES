cmake_minimum_required(VERSION 3.24)

include(FetchContent)
include(${CMAKE_SOURCE_DIR}/cmake/FetchUtils.cmake)

# Prefer extern/catch2-*.tgz; else fall back online.
set(_catch2_url
    "https://github.com/catchorg/Catch2/archive/refs/tags/v3.5.3.tar.gz")
fu_url_args("catch2" "${_catch2_url}" _C2_URL_ARGS)

if(PREFETCH_THIRD_PARTY)
  FetchContent_Declare(Catch2 ${_C2_URL_ARGS})
  FetchContent_Populate(Catch2)
  add_custom_target(catch2_src COMMENT "Catch2 sources are populated")
  _fu_register_prefetch_target(catch2_src)
  return()
endif()

FetchContent_Declare(Catch2 ${_C2_URL_ARGS})
FetchContent_MakeAvailable(Catch2) # Defines Catch2::Catch2 and
                                   # Catch2::Catch2WithMain
