cmake_minimum_required(VERSION 3.24)

# Gather headers and sources
file(GLOB CORE_HEADERS CONFIGURE_DEPENDS include/*.hpp)
set(CORE_SOURCES memory/MemoryManager.cpp)

add_library(core STATIC ${CORE_SOURCES} ${CORE_HEADERS})

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

if(CUDAToolkit_FOUND)
  target_link_libraries(core PUBLIC CUDA::cudart)
  # Build any .cu files with correct arch; adjust as needed
  set_target_properties(core PROPERTIES CUDA_ARCHITECTURES 70)
endif()

# Propagate global compile definitions (USE_CUDA_UM, HAVE_MPI, HAVE_CUDA) They
# are already set at directory scope by the root listfile, nothing to do.

# Installation helpers (optional but recommended) ----------------------------
install(
  TARGETS core
  EXPORT coreTargets
  ARCHIVE DESTINATION lib)
install(DIRECTORY include/ DESTINATION include)
