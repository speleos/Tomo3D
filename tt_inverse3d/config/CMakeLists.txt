try_compile( 
  AZUR_HAVE_PLAIN_RESTRICT 
  ${CMAKE_BINARY_DIR}/prestrict
  ${PROJECT_SOURCE_DIR}/config/has_plain_restrict.cpp 
  )

try_compile( 
  AZUR_HAVE_UNDERSCORED_RESTRICT 
  ${CMAKE_BINARY_DIR}/urestrict
  ${PROJECT_SOURCE_DIR}/config/has_underscored_restrict.cpp 
  OUTPUT_VARIABLE OUTPUT
  )

if(AZUR_HAVE_UNDERSCORED_RESTRICT)
  message( STATUS "The restrict keyword is known as __restrict__.")
else(AZUR_HAVE_UNDERSCORED_RESTRICT)
  message(${OUTPUT})
endif(AZUR_HAVE_UNDERSCORED_RESTRICT)

try_compile( 
  AZUR_CAN_INHERIT_CTOR
  ${CMAKE_BINARY_DIR}/cxx11ctor
  ${PROJECT_SOURCE_DIR}/config/can_inherit_ctor.cpp
  COMPILE_DEFINITIONS "${CMAKE_CXX_FLAGS}"
  OUTPUT_VARIABLE CXX11_CTOR_OUTPUT
  )

if (NOT AZUR_CAN_INHERIT_CTOR)
  message( FATAL_ERROR "Cannot inherit ctor (required), did you enable C++11 ?\n ${CXX11_CTOR_OUTPUT}")
else()
  message( STATUS "Ctor inheritance supported.")
endif()

configure_file(
  ${CMAKE_SOURCE_DIR}/config/config.hpp.in
  ${CMAKE_BINARY_DIR}/include/config.hpp
  )

