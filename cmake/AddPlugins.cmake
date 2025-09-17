# Simple helper to declare a physics plugin.
#
# add_physics_plugin( TARGET   physics_fluids SOURCES  src/a.cpp;src/b.cpp
# KERNELS  kernels/x.f90;kernels/y.f90 INCLUDE  include LINK core;PETSC::petsc )
function(add_physics_plugin)
  cmake_parse_arguments(APP "" "TARGET;INCLUDE" "SOURCES;KERNELS;LINK" ${ARGN})

  if(NOT APP_TARGET)
    message(FATAL_ERROR "add_physics_plugin: TARGET is required")
  endif()

  set(_all_sources ${APP_SOURCES} ${APP_KERNELS})
  add_library(${APP_TARGET} SHARED ${_all_sources})

  if(APP_INCLUDE)
    target_include_directories(
      ${APP_TARGET} PRIVATE ${APP_INCLUDE} ${CMAKE_SOURCE_DIR}/src/core/include)
  else()
    target_include_directories(${APP_TARGET}
                               PRIVATE ${CMAKE_SOURCE_DIR}/src/core/include)
  endif()

  if(APP_LINK)
    target_link_libraries(${APP_TARGET} PRIVATE ${APP_LINK})
  endif()

  set_target_properties(${APP_TARGET} PROPERTIES OUTPUT_NAME ${APP_TARGET}
                                                 POSITION_INDEPENDENT_CODE ON)
endfunction()
