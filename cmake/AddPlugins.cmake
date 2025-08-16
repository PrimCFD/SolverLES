function(add_plugin name)
  add_library(${name} SHARED ${ARGN})
  target_link_libraries(${name} PUBLIC core)
  set_target_properties(${name} PROPERTIES LIBRARY_OUTPUT_DIRECTORY
                                           ${PROJECT_BINARY_DIR}/lib/plugins)
  install(TARGETS ${name} LIBRARY DESTINATION lib/plugins)
endfunction()
