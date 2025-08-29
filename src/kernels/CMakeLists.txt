# Static library for compute kernels (Fortran/CUDA).

add_library(kernels OBJECT flux_kernels.f90)
target_link_libraries(kernels PUBLIC OpenMP::OpenMP_Fortran) 
