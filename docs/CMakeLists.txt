find_package(Sphinx REQUIRED)

set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR})
set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/html)

add_custom_target(docs ALL
    COMMAND ${Sphinx_EXECUTABLE} -b html ${SPHINX_SOURCE} ${SPHINX_BUILD}
    COMMENT "Generating documentation with Sphinx"
)