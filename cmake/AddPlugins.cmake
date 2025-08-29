# Helper to build runtime plugins cleanly, with correct PETSc ordering.

include_guard(GLOBAL)

function(add_plugin target)
  set(options)
  set(oneValueArgs)
  set(multiValueArgs SOURCES LINK_LIBS INCLUDE_DIRS DEFINES)
  cmake_parse_arguments(AP "${options}" "${oneValueArgs}" "${multiValueArgs}"
                        ${ARGN})

  if(NOT AP_SOURCES)
    message(FATAL_ERROR "add_plugin(${target}) requires: SOURCES <...>")
  endif()

  add_library(${target} SHARED ${AP_SOURCES})
  target_compile_features(${target} PRIVATE cxx_std_23)

  if(AP_INCLUDE_DIRS)
    target_include_directories(${target} PRIVATE ${AP_INCLUDE_DIRS})
  endif()

  if(AP_DEFINES)
    target_compile_definitions(${target} PRIVATE ${AP_DEFINES})
  endif()

  if(AP_LINK_LIBS)
    target_link_libraries(${target} PRIVATE ${AP_LINK_LIBS})
  endif()

  # Project-wide warnings/options hook (if you have it)
  if(COMMAND project_apply_common_warnings)
    project_apply_common_warnings(${target})
  endif()

  # --- PETSc ordering: if this plugin links PETSC::petsc, depend on it.
  if(TARGET PETSC::petsc)
    get_target_property(_ll ${target} LINK_LIBRARIES)
    if(_ll)
      string(REPLACE ";" ";" _ll_sc ";${_ll};")
      if(_ll_sc MATCHES ";PETSC::petsc;")
        add_dependencies(${target} PETSC::petsc)
      endif()
    endif()
  endif()
endfunction()
