# Static library for compute kernels (Fortran/CUDA).

add_library(kernels STATIC stub_kernel.f90)

target_link_libraries(kernels PUBLIC core)

install(TARGETS kernels ARCHIVE DESTINATION lib)
