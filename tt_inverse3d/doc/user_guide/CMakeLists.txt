include(${CMAKE_SOURCE_DIR}/util/cmake/UseLATEX.cmake)

if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(WARN "    Only out of source build supported for latex file generation.")
else ()
  set(LATEX_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR})
  add_latex_document(
    user.tex
    IMAGE_DIRS fig
    DEFAULT_PDF
    MANGLE_TARGET_NAMES
    )
endif ()
