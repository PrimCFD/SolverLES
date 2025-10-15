# General
add_executable(
  unit_core
  test_alignment.cpp
  test_registry_memory.cpp
  test_field_index.cpp
  test_copy_roundtrip.cpp
  test_device_mirror.cpp
  test_boundary.cpp
  test_fieldcatalog.cpp
  test_fieldcatalog_more.cpp
  test_scheduler_noop.cpp
  test_scheduler_calls.cpp
  test_plugin_register.cpp
  test_registry_error_master.cpp
  test_async_backpressure.cpp
  test_cgns_roundtrip.cpp
  test_xdmf_roundtrip.cpp
  test_preflight.cpp
  test_writeplan_pack.cpp
  test_config_yaml.cpp)

target_link_libraries(unit_core PRIVATE core Catch2::Catch2WithMain)

# Auto-register every TEST_CASE as a CTest test, all labeled "unit"
register_catch_tests(
  unit_core
  LABEL
  unit
  TEST_PREFIX
  "core::"
  REPORTER
  junit
  OUTPUT_DIR
  "${CMAKE_BINARY_DIR}/test-reports/unit"
  OUTPUT_PREFIX
  "catch-"
  OUTPUT_SUFFIX
  ".xml")

if(ENABLE_CUDA)
  add_subdirectory(gpu)
endif()
