cmake_minimum_required(VERSION 3.24)

# Public headers (everything under include/)
file(GLOB_RECURSE CORE_PUBLIC_HEADERS CONFIGURE_DEPENDS
     "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp")

# Private sources (all .cpp in src/)
file(GLOB_RECURSE  CORE_SOURCES CONFIGURE_DEPENDS
     "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")

add_library(core STATIC ${CORE_SOURCES} ${CORE_PUBLIC_HEADERS})

# Public include path (build tree and install tree)
target_include_directories(
  core PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
              $<INSTALL_INTERFACE:include>)

# Require C++23 for all consumers
set_target_properties(core PROPERTIES CXX_STANDARD 23 CXX_STANDARD_REQUIRED YES)

# Optional dependencies ------------------------------------------------------
if(MPI_FOUND)
  target_link_libraries(core PUBLIC MPI::MPI_CXX)
endif()

if(ENABLE_CUDA)
  target_link_libraries(core PUBLIC CUDA::cudart)
  target_compile_definitions(core PUBLIC HAVE_CUDA)
  # Optional UM default:
  if(USE_CUDA_UM)
    target_compile_definitions(core PUBLIC USE_CUDA_UM)
  endif()
endif()

# Needed for PluginHost on Linux
if(UNIX AND NOT APPLE)
  target_link_libraries(core PUBLIC dl)
endif()

# Propagate global compile definitions (USE_CUDA_UM, HAVE_MPI, HAVE_CUDA) They
# are already set at directory scope by the root listfile, nothing to do.

# Installation helpers (optional but recommended) ----------------------------
install(
  TARGETS core
  EXPORT coreTargets
  ARCHIVE DESTINATION lib)
install(DIRECTORY include/ DESTINATION include)
