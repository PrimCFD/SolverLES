add_subdirectory(core)
add_subdirectory(kernels)

# Physics stacks (each may provide plugins and helper libs)
# add_subdirectory(physics)

add_subdirectory(apps/solver) # solver

if(BUILD_GUI)
  add_subdirectory(ipc)
  add_subdirectory(gui) # shared lib
  add_subdirectory(apps/gui_launcher)
endif()
