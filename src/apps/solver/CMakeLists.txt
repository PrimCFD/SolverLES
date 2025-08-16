add_executable(solver main.cpp)
target_link_libraries(solver PRIVATE core kernels ipc)

# RPATH so it finds lib/plugins at runtime
set_target_properties(solver PROPERTIES INSTALL_RPATH "$ORIGIN/../lib/plugins")

install(TARGETS solver RUNTIME DESTINATION bin)
