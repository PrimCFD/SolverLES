project(numerics LANGUAGES CXX Fortran)

# Public headers (exported to plugins)
file(GLOB_RECURSE NUMERICS_HEADERS CONFIGURE_DEPENDS
     "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp")

# C++ sources
set(NUMERICS_SOURCES src/MacOps.cpp src/PoissonScalar.cpp)

# Fortran kernels
set(NUMERICS_KERNELS kernels/kernels.f90)

add_library(numerics STATIC ${NUMERICS_SOURCES} ${NUMERICS_HEADERS}
                            ${NUMERICS_KERNELS})

# Include paths: install as "numerics/..."
target_include_directories(
  numerics PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                  $<INSTALL_INTERFACE:include>)

# Depends on core's types (Mesh, FieldCatalog, views, etc.)
target_link_libraries(
  numerics PUBLIC core # reuse Mesh, RunContext, FieldCatalog, etc.
                  PETSC::petsc # for ScalarPoisson
)

find_package(OpenMP REQUIRED)

target_link_libraries(numerics PRIVATE OpenMP::OpenMP_Fortran
                                       OpenMP::OpenMP_CXX)

# Generic high opt for both C++ and Fortran in numerics
target_compile_options(
  numerics PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-O3 -march=native>
                   $<$<COMPILE_LANGUAGE:Fortran>:-O3 -march=native>)

# Extra SIMD/bandwidth-friendly switches per compiler (Fortran only)
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  target_compile_options(
    numerics PRIVATE $<$<COMPILE_LANGUAGE:Fortran>: -fopenmp-simd -ffast-math
                     -fno-math-errno -fno-trapping-math -funroll-loops>)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM")
  target_compile_options(numerics PRIVATE $<$<COMPILE_LANGUAGE:Fortran>: -Ofast
                                          -qopt-zmm-usage=high>)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  target_compile_options(numerics PRIVATE $<$<COMPILE_LANGUAGE:Fortran>: -O3
                                          -xHost>)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Flang|LLVMFlang")
  target_compile_options(numerics PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:
                                          -ffast-math>)
endif()

# (Optional) vectorization reports just for numerics
option(NUMERICS_VEC_REPORT "Enable vectorization reports for kernels" OFF)
if(NUMERICS_VEC_REPORT)
  if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    target_compile_options(
      numerics PRIVATE $<$<COMPILE_LANGUAGE:Fortran>: -fopt-info-vec-missed
                       -fopt-info-vec-optimized>)
  elseif(CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM")
    target_compile_options(
      numerics PRIVATE $<$<COMPILE_LANGUAGE:Fortran>: -qopt-report=5
                       -qopt-report-phase=vec>)
  endif()
endif()

# Install alongside core
install(
  TARGETS numerics
  EXPORT numericsTargets
  ARCHIVE DESTINATION lib)

install(DIRECTORY include/ DESTINATION include)
