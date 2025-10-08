add_executable(solver main.cpp)
target_link_libraries(solver PRIVATE core dl PETSC::petsc)

# Use $ORIGIN so the binary is relocatable
set(CMAKE_BUILD_RPATH_USE_ORIGIN ON)
set(_origin "\$ORIGIN")

# Search paths at runtime
set_target_properties(
  solver
  PROPERTIES
    # Run from the build tree
    BUILD_RPATH
    "${_origin};${_origin}/../../core;${_origin}/../../physics/fluids;${_origin}/../../physics/em;$<TARGET_FILE_DIR:PETSC::petsc_lib>"
    # Run after `cmake --install` (assuming bin/, lib/, lib/plugins install
    # layout)
    INSTALL_RPATH
    "${_origin}/../lib;${_origin}/../lib/plugins;$<TARGET_FILE_DIR:PETSC::petsc_lib>"
)

install(TARGETS solver RUNTIME DESTINATION bin)
