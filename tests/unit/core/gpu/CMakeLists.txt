add_executable(unit_gpu test_um_fields.cpp test_io_residency.cpp)
target_link_libraries(unit_gpu PRIVATE CUDA::cudart core Catch2::Catch2WithMain)
target_compile_definitions(unit_gpu PRIVATE HAVE_CUDA)
