if(COMMAND resolve_url) # guard so no redefine on repeat include()
  return()
endif()

function(resolve_url _tarball _outvar _remote_url)
  # Prefer a pre-downloaded archive under <source-root>/extern
  if(EXISTS ${CMAKE_SOURCE_DIR}/extern/${_tarball})
    set(${_outvar}
        "file://${CMAKE_SOURCE_DIR}/extern/${_tarball}"
        PARENT_SCOPE)
  else()
    set(${_outvar}
        "${_remote_url}"
        PARENT_SCOPE)
  endif()
endfunction()
