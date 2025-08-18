add_subdirectory(core)
add_subdirectory(kernels)

if(BUILD_PLUGINS)
  add_subdirectory(plugins)
endif()

add_subdirectory(apps/solver) # solver

if(BUILD_GUI)
  add_subdirectory(ipc)
  add_subdirectory(gui) # shared lib
  add_subdirectory(apps/gui_launcher)
endif()

target_link_libraries(core PUBLIC CGNS::cgns PETSC::petsc)
