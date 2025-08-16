list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
include(Testing)

# Add suites
add_subdirectory(unit/core)

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/performance/core/CMakeLists.txt")
  add_subdirectory(performance/core)
endif()

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/regression/CMakeLists.txt")
  add_subdirectory(regression)
endif()
