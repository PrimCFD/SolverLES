# Ensure the cache has a project name for scripts that read CMakeCache.txt
# (e.g., scripts/prefetch_third_party.sh). Using LANGUAGES NONE avoids finding
# compilers during prefetch-only builds.
if(NOT CMAKE_PROJECT_NAME)
  project(KolmoPlas LANGUAGES NONE)
endif()

# Turn on "prefetch-only" behavior inside each Fetch*.cmake
set(PREFETCH_THIRD_PARTY
    ON
    CACHE BOOL "" FORCE)

include(${CMAKE_SOURCE_DIR}/cmake/FetchUtils.cmake)

# Populate but do not build. Each file registers a *_src target.
include(${CMAKE_SOURCE_DIR}/cmake/FetchHDF5.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/FetchCGNS.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/FetchOpenBLAS.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/FetchPETSc.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/FetchCatch2.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/FetchYamlCpp.cmake)

# Aggregate target used by scripts/prefetch_third_party.sh
fu_make_prefetch_aggregate()
