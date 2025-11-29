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

add_physics_plugin(
  TARGET  physics_fluids
  SOURCES "${SRC}"
  INCLUDE "${CMAKE_CURRENT_SOURCE_DIR}/include"
  LINK    core;numerics;PETSC::petsc)

# --- OpenMP + per-target tuning ----
find_package(OpenMP REQUIRED)

# Propagate the correct compile/link flags for both languages
target_link_libraries(physics_fluids PRIVATE OpenMP::OpenMP_CXX)

# Generic high opt for both C++ and Fortran
target_compile_options(
  physics_fluids PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-O3 -march=native>)

set_target_properties(
  physics_fluids
  PROPERTIES
    BUILD_RPATH
    "$ORIGIN;$ORIGIN/..;$ORIGIN/../lib;$<TARGET_FILE_DIR:PETSC::petsc_lib>"
    INSTALL_RPATH
    "$ORIGIN;$ORIGIN/..;$ORIGIN/../lib;$<TARGET_FILE_DIR:PETSC::petsc_lib>")

target_link_options(physics_fluids PRIVATE $<$<PLATFORM_ID:Linux>:-Wl,-z,defs>)
