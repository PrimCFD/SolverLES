include(AddPlugins) # helper that defines the function below

# iterate over every immediate child dir (flux/ sgs/ time/)
file(
  GLOB children
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  *)
foreach(subdir ${children})
  if(IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${subdir}")
    add_subdirectory(${subdir})
  endif()
endforeach()
