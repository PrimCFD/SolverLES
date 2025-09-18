project(physics_fluids LANGUAGES CXX Fortran)

include(${CMAKE_SOURCE_DIR}/cmake/AddPlugins.cmake)

set(SRC
    src/Register.cpp
    src/ApplyBCs.cpp
    src/Program.cpp
    src/ProjectionLoop.cpp
    src/InitTG.cpp
    src/SGS.cpp
    src/MomentumPredictor.cpp
    src/PressurePoisson.cpp
    src/VelocityCorrector.cpp)

set(KRN kernels/fluids.f90)

add_physics_plugin(
  TARGET
  physics_fluids
  SOURCES
  "${SRC}"
  KERNELS
  "${KRN}"
  INCLUDE
  "${CMAKE_CURRENT_SOURCE_DIR}/include"
  LINK
  core)

# --- OpenMP + per-target tuning ----
find_package(OpenMP REQUIRED)

# Propagate the correct compile/link flags for both languages
target_link_libraries(physics_fluids PRIVATE OpenMP::OpenMP_Fortran
                                             OpenMP::OpenMP_CXX)

# Generic high opt for both C++ and Fortran
target_compile_options(
  physics_fluids PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-O3 -march=native>
                         $<$<COMPILE_LANGUAGE:Fortran>:-O3 -march=native>)

# Extra SIMD/bandwidth-friendly switches per compiler (Fortran)
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  target_compile_options(
    physics_fluids
    PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-fopenmp-simd -ffast-math
            -fno-math-errno -fno-trapping-math -funroll-loops>)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM")
  # ifx / Intel oneAPI (LLVM-based)
  target_compile_options(
    physics_fluids PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-Ofast
                           -qopt-zmm-usage=high>)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  # classic ifort
  target_compile_options(physics_fluids
                         PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-O3 -xHost>)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Flang|LLVMFlang")
  target_compile_options(physics_fluids
                         PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-ffast-math>)
endif()

# (Optional) vectorization reports: enable when investigating perf
option(PHYS_FLUIDS_VEC_REPORT "Enable vectorization reports" OFF)
if(PHYS_FLUIDS_VEC_REPORT)
  if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    target_compile_options(
      physics_fluids PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-fopt-info-vec-missed
                             -fopt-info-vec-optimized>)
  elseif(CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM")
    target_compile_options(
      physics_fluids PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-qopt-report=5
                             -qopt-report-phase=vec>)
  endif()
endif()
