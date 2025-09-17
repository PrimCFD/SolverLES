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
