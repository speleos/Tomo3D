cmake_minimum_required(VERSION 2.6)

project (TT_Inverse3D CXX C)

include(CTest)

enable_testing()

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-std=c++11" COMPILER_CXX11_OPTION)
if(COMPILER_CXX11_OPTION)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
  message( STATUS "adding \"-std=c++11\" for C++11 support.")
else()
  message( STATUS "No C++11 specific option found, either not supported or enable by default.")
  message( STATUS "(short term) future will tell.")
endif()

set(Boost_USE_MULTITHREADED      ON) 
set(Boost_USE_STATIC_LIBS        ON) 
set(Boost_DEBUG      OFF) 
#find_package( Boost 1.59.0 COMPONENTS mpi serialization thread system program_options)
find_package( Boost 1.78.0 COMPONENTS mpi serialization thread system program_options)

message( "Boost root set in ${BOOST_ROOT}")
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
  message( STATUS "Boost ${Boost_VERSION} found in ${BOOST_ROOT}")
else (Boost_FOUND)
  message( FATAL_ERROR "${Boost_ERROR_REASON}")
endif(Boost_FOUND)

find_library(MKL mkl_intel_lp64)
if(${MKL} STREQUAL "MKL-NOTFOUND")
  message( STATUS "PARDISO not found.")
  set(TOMO3D_PARDISO_AVAILABLE 1)
  set(PARDISO_LIBS "")
  set(PARDISO_SRC pardiso_solver_stub.cpp)
#else(${MKL} EQUAL "MKL-NOTFOUND")   #EQUAL compares numbers, not strings!
else(${MKL} STREQUAL "MKL")  #If library is found find_library returns <VAR>, in this case MKL
  message( STATUS "PARDISO found. (in ${MKL})")
  set(TOMO3D_PARDISO_AVAILABLE 1)
  set(PARDISO_LIBS mkl_intel_lp64 mkl_intel_thread mkl_core iomp5)
  set(PARDISO_SRC pardiso_solver.cpp)
endif()

if (NOT (MPI_COMPILE_FLAGS OR MPI_LINK_FLAGS) )
  find_package( MPI )

  if(MPI_FOUND)
    message( STATUS "Found MPI")
    message( STATUS "   compile with flags \"${MPI_C_COMPILE_FLAGS}\" and include dir \"${MPI_INCLUDE_PATH}\"")
    message( STATUS "   and link with \"${MPI_C_LINK_FLAGS}\"  and \"${MPI_C_LIBRARIES}\"")
    include_directories(${MPI_INCLUDE_PATH})
    set( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${MPI_LINK_FLAGS}" )
  else(MPI_FOUND)
    message( FATAL_ERROR "Could not find MPI")
  endif(MPI_FOUND)
else(NOT (MPI_COMPILE_FLAGS OR MPI_LINK_FLAGS) )
  set(MPI_FOUND True)
#  set({MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2")
endif(NOT (MPI_COMPILE_FLAGS OR MPI_LINK_FLAGS) )

find_package(OpenMP)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} ${MPI_COMPILE_FLAGS}")

include_directories(${CMAKE_BINARY_DIR}/include)
add_subdirectory(config)

add_library(tomo3d
  ani3d.cpp
  axb_solver.cpp
  in_house_solver.cpp
  in_house_omp_solver.cpp
  ${PARDISO_SRC}
  inverse3d.cpp
  inverse3d_solve.cpp
  sparse_rect.cpp
  mesh.cpp
  smesh3d.cpp
  graph3d.cpp 
  traveltime3d.cpp
  betaspline3d.cpp 
  bend3d.cpp
  bend_mnbrak.cpp
  bend_brent.cpp 
  interface3d.cpp
  corrlen3d.cpp 
  geom3d.cpp
  error.cpp
  d_jacobi.cpp
  d_choldc.cpp 
  nrutil.cpp
  d_realft.cpp)
target_link_libraries( tomo3d ${PARDISO_LIBS} ${Boost_LIBRARIES} ${MPI_LIBRARIES} pthread)

add_executable(tt_inverse3d tt_inverse3d.cpp)
set_target_properties( tt_inverse3d
  PROPERTIES
  LINK_FLAGS "${OpenMP_CXX_FLAGS} ${MPI_LINK_FLAGS}")
target_link_libraries( tt_inverse3d tomo3d)

add_executable(offline_solver offline_solver.cpp)
set_target_properties( offline_solver
  PROPERTIES
  LINK_FLAGS "${OpenMP_CXX_FLAGS} ${MPI_LINK_FLAGS}")
target_link_libraries( offline_solver tomo3d)

add_subdirectory(test)

if (EXISTS ${CMAKE_SOURCE_DIR}/doc
    AND ( NOT EXISTS ${CMAKE_SOURCE_DIR}/SKIP_DOC))
  add_subdirectory(doc)
endif()
