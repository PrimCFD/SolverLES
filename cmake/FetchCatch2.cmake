include(FetchContent)
include(FetchUtils)

if(NOT DEFINED FETCHCONTENT_UPDATES_DISCONNECTED)
  set(FETCHCONTENT_UPDATES_DISCONNECTED
      ON
      CACHE BOOL "Disable URL re‑checks")
endif()

set(CATCH2_VERSION 3.6.0)
set(CATCH2_TAR catch2-${CATCH2_VERSION}.tar.gz)
set(_catch2_default_url
    "https://github.com/catchorg/Catch2/archive/refs/tags/v${CATCH2_VERSION}.tar.gz"
)

resolve_url(${CATCH2_TAR} _catch2_url ${_catch2_default_url})
set(_catch2_sha256
    485932259a75c7c6b72d4b874242c489ea5155d17efa345eb8cc72159f49f356)
set(_catch2_install ${CMAKE_BINARY_DIR}/_deps/catch2-install)

FetchContent_Declare(
  catch2
  URL ${_catch2_url}
  URL_HASH SHA256=${_catch2_sha256}
  CMAKE_ARGS
  -DCMAKE_INSTALL_PREFIX=${_catch2_install}
  -DCATCH_BUILD_TESTING=OFF
  -DCATCH_BUILD_EXAMPLES=OFF
  -DCATCH_BUILD_EXTRA_TESTS=OFF
  -DCATCH_INSTALL_DOCS=OFF
  -DCATCH_INSTALL_EXTRAS=OFF)

if(PREFETCH_THIRD_PARTY)
  FetchContent_Populate(catch2)
  return()
endif()

FetchContent_MakeAvailable(catch2)

list(APPEND CMAKE_PREFIX_PATH "${_catch2_install}")
