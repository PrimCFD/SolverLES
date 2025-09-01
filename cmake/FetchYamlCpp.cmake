include(FetchContent)
include(${CMAKE_SOURCE_DIR}/cmake/FetchUtils.cmake)

# Prefer extern/yamlcpp-*.tgz; else fall back online.
set(_yamlcpp_url
    "https://github.com/jbeder/yaml-cpp/archive/refs/tags/0.8.0.tar.gz")
fu_url_args("yaml-cpp" "${_yamlcpp_url}" _yamlcpp_URL_ARGS)

if(PREFETCH_THIRD_PARTY)
  FetchContent_Declare(yaml-cpp ${_yamlcpp_URL_ARGS})
  FetchContent_Populate(yaml-cpp)
  add_custom_target(yaml-cpp_src COMMENT "yaml-cpp sources are populated")
  _fu_register_prefetch_target(yaml-cpp_src)
  return()
endif()

FetchContent_Declare(yaml-cpp ${_yamlcpp_URL_ARGS})
FetchContent_MakeAvailable(yaml-cpp) # Defines YAML_CPP::yaml-cpp
