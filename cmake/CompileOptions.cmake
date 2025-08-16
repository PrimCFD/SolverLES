# Global warnings & sanitizers — GCC/Clang only
if(TARGET project_warnings)
  return()
endif()

option(WARNINGS_AS_ERRORS "Treat warnings as errors for project code" OFF)
option(STRICT_WARNINGS "Turn on extra noisy warnings (dev/CI)" OFF)

add_library(project_warnings INTERFACE)

# Base: safe, low-noise
target_compile_options(
  project_warnings
  INTERFACE
    $<$<COMPILE_LANG_AND_ID:C,GNU,Clang,AppleClang>:
    -Wall
    -Wextra
    -Wpedantic
    -Wformat=2
    -Wcast-align
    -Wunused
    >
    $<$<COMPILE_LANG_AND_ID:CXX,GNU,Clang,AppleClang>:
    -Wall
    -Wextra
    -Wpedantic
    -Wformat=2
    -Wcast-align
    -Wunused
    -Woverloaded-virtual
    -Wnon-virtual-dtor
    >
    $<$<AND:$<BOOL:${WARNINGS_AS_ERRORS}>,$<COMPILE_LANG_AND_ID:C,GNU,Clang,AppleClang>>:-Werror>
    $<$<AND:$<BOOL:${WARNINGS_AS_ERRORS}>,$<COMPILE_LANG_AND_ID:CXX,GNU,Clang,AppleClang>>:-Werror>
)

# Optional: chase conversions/shadows, etc.
if(STRICT_WARNINGS)
  target_compile_options(
    project_warnings
    INTERFACE $<$<COMPILE_LANG_AND_ID:C,GNU,Clang,AppleClang>:
              -Wshadow
              -Wconversion
              -Wdouble-promotion
              -Wnull-dereference
              -Wundef
              >
              $<$<COMPILE_LANG_AND_ID:CXX,GNU,Clang,AppleClang>:
              -Wshadow
              -Wconversion
              -Wsign-conversion
              -Wold-style-cast
              -Wnull-dereference
              -Wundef
              >)
endif()

option(SANITIZE "Build with ASan/UBSan" OFF)
if(SANITIZE)
  add_compile_options(-fsanitize=address,undefined)
  add_link_options(-fsanitize=address,undefined)
endif()

# Auto-attach warnings to every non-imported target whose SOURCE_DIR is inside
# any of the given roots (e.g. ${CMAKE_SOURCE_DIR}/src, tests).
function(project_enable_warnings_under)
  if(NOT ARGN)
    message(
      FATAL_ERROR
        "project_enable_warnings_under(<dir> [...]) requires at least one directory"
    )
  endif()

  # Normalize roots to absolute paths
  set(_roots_abs "")
  foreach(_r IN LISTS ARGN)
    if(NOT IS_ABSOLUTE "${_r}")
      set(_r "${CMAKE_SOURCE_DIR}/${_r}")
    endif()
    file(REAL_PATH "${_r}" _r_abs)
    list(APPEND _roots_abs "${_r_abs}")
  endforeach()

  # Scan all targets created so far
  get_property(_all GLOBAL PROPERTY TARGETS)
  set(_enabled "")
  foreach(tgt IN LISTS _all)
    # skip imported/alias/unknown
    get_target_property(_imp "${tgt}" IMPORTED)
    if(_imp)
      continue()
    endif()
    get_target_property(_srcdir "${tgt}" SOURCE_DIR)
    if(NOT _srcdir)
      continue()
    endif()
    file(REAL_PATH "${_srcdir}" _src_abs)

    # If the target's source dir is inside any root, attach warnings
    set(_hit FALSE)
    foreach(_root IN LISTS _roots_abs)
      file(RELATIVE_PATH _rel "${_root}" "${_src_abs}")
      if(NOT _rel MATCHES "^[.][.]") # not starting with '..' → inside
        set(_hit TRUE)
        break()
      endif()
    endforeach()

    if(_hit)
      target_link_libraries("${tgt}" PRIVATE project_warnings)
      list(APPEND _enabled "${tgt}")
    endif()
  endforeach()

  list(REMOVE_DUPLICATES _enabled)
  if(_enabled)
    message(STATUS "project_warnings → ${_enabled}")
  endif()
endfunction()
