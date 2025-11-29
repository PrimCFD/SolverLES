# Simple helper to declare a physics plugin.
#
# New usage: add_physics_plugin( TARGET   physics_fluids SOURCES
# src/a.cpp;src/b.cpp INCLUDE  include LINK     core;numerics;PETSC::petsc)
function(add_physics_plugin)
  cmake_parse_arguments(APP "" "TARGET;INCLUDE" "SOURCES;LINK" ${ARGN})

  if(NOT APP_TARGET)
    message(FATAL_ERROR "add_physics_plugin: TARGET is required")
  endif()

  if(NOT APP_SOURCES)
    message(FATAL_ERROR "add_physics_plugin: SOURCES is required")
  endif()

  # Plugins are built as shared libraries from C++ sources only.
  add_library(${APP_TARGET} SHARED ${APP_SOURCES})

  # Always expose core headers; plugin-specific headers are optional.
  if(APP_INCLUDE)
    target_include_directories(
      ${APP_TARGET} PRIVATE ${APP_INCLUDE} ${CMAKE_SOURCE_DIR}/src/core/include)
  else()
    target_include_directories(${APP_TARGET}
                               PRIVATE ${CMAKE_SOURCE_DIR}/src/core/include)
  endif()

  # Link against requested libs (e.g. core, numerics, PETSC::petsc)
  if(APP_LINK)
    target_link_libraries(${APP_TARGET} PRIVATE ${APP_LINK})
  endif()

  set_target_properties(${APP_TARGET} PROPERTIES OUTPUT_NAME ${APP_TARGET}
                                                 POSITION_INDEPENDENT_CODE ON)
endfunction()
