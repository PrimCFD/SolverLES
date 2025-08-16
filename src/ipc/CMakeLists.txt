# Static library for IPC helpers (ZeroMQ / gRPC).

add_library(ipc STATIC stub_ipc.cpp)

target_link_libraries(ipc PUBLIC core)

install(TARGETS ipc ARCHIVE DESTINATION lib)
