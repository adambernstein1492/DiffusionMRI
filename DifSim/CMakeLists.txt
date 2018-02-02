PROJECT(DIFSIM)

cmake_minimum_required(VERSION 2.8)
SET(CMAKE_BUILD_TYPE Release)

set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "")
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake-modules")

if( "$ENV{TARGET}" STREQUAL "" )
  set(TARGET "x86_64")
else()
  set(TARGET "$ENV{TARGET}" )
endif()
if( "$ENV{CC}" STREQUAL "" )
  set(ENV{CC} "gcc")
  set(ENV{CXX} "g++")
endif()
EXECUTE_PROCESS(
  COMMAND $ENV{CC} -v 
  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
  OUTPUT_QUIET
  ERROR_VARIABLE "ICC"
  OUTPUT_STRIP_TRAILING_WHITESPACE)
string( FIND ${ICC} "icc" POS )
if( "${POS}" STREQUAL "-1" )
  set(ENV{CC} "gcc")
  set(ENV{CXX} "g++")
else()
  set(ENV{CC} "icc")
  set(ENV{CXX} "icpc")
endif()
MESSAGE( STATUS "ICC: ${POS} : $ENV{CC}" )

if( ${TARGET} STREQUAL "mic" )

  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mmic -mkl -std=c++0x")
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mmic -mkl")

  include(ExternalProject)

  ExternalProject_Add(mcell-3
  PREFIX ${CMAKE_SOURCE_DIR}/mcell-3
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/mcell-3
  TMP_DIR ${CMAKE_SOURCE_DIR}/mcell-3
  STAMP_DIR ${CMAKE_SOURCE_DIR}/mcell-3
  BINARY_DIR ${CMAKE_SOURCE_DIR}/mcell-3
  CONFIGURE_COMMAND ${CMAKE_SOURCE_DIR}/mcell_configure.sh ${CMAKE_SOURCE_DIR} ${TARGET}
  )
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/mcell-3 )
  set(MCELL_LIBRARY "${CMAKE_SOURCE_DIR}/mcell-3/libmcell.a")

  FIND_PATH(FFTW_INCLUDE_DIR
  NAMES fftw3.h
  PATHS $ENV{MKLROOT}/include/fftw
  NO_DEFAULT_PATH
  DOC "FFTW header include dir"
  )
  SET(FFTW_LIBRARIES '')

  message("-- Found FFTW:")
  message("--   " ${FFTW_INCLUDE_DIR})

  add_definitions(-DMRI_DIFF_SIM)

  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/mcell-3 ${FFTW_INCLUDE_DIR})

  SET(MPI_DIR $ENV{I_MPI_ROOT})

  if (MPI_DIR)
    message("-- Found MPI:")
    message("--   " ${MPI_DIR})
    add_definitions(-DDO_MPI -DMPICH_SKIP_MPICXX)
    INCLUDE_DIRECTORIES(${CMAKE_INCLUDE_PATH} ${MPI_DIR}/include)
    SET(CMAKE_C_COMPILER ${MPI_DIR}/bin64/mpiicc)
    SET(CMAKE_CXX_COMPILER ${MPI_DIR}/bin64/mpiicpc)
    SET(CMAKE_LINKER ${MPI_DIR}/bin64/mpiicpc)
    message("--  " ${CMAKE_LINKER})
    message("--  " ${CMAKE_CXX_COMPILER})
  else ()
    message("-- Not found MPI")
  endif ()
  
  set(difsim_kernel_SRCS difsim_kernel.cc gtb_histogram.cc sphtes.cc ConfigFile.cc )
  add_executable(difsim ${difsim_kernel_SRCS})
  target_link_libraries (difsim ${MCELL_LIBRARY})
  add_dependencies(difsim mcell-3)

else ()
  if ( $ENV{CC} STREQUAL "icc" )

  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mkl -std=c++0x")
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mkl")

  include(ExternalProject)

  ExternalProject_Add(mcell-3
  PREFIX ${CMAKE_SOURCE_DIR}/mcell-3
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/mcell-3
  TMP_DIR ${CMAKE_SOURCE_DIR}/mcell-3
  STAMP_DIR ${CMAKE_SOURCE_DIR}/mcell-3
  BINARY_DIR ${CMAKE_SOURCE_DIR}/mcell-3
  CONFIGURE_COMMAND ${CMAKE_SOURCE_DIR}/mcell_configure.sh ${CMAKE_SOURCE_DIR} ${TARGET}
  )
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/mcell-3 )
  set(MCELL_LIBRARY "${CMAKE_SOURCE_DIR}/mcell-3/libmcell.a")

  FIND_PATH(FFTW_INCLUDE_DIR
  NAMES fftw3.h
  PATHS $ENV{MKLROOT}/include/fftw
  NO_DEFAULT_PATH
  DOC "FFTW header include dir"
  )
  SET(FFTW_LIBRARIES '')

  message("-- Found FFTW:")
  message("--   " ${FFTW_INCLUDE_DIR})

  add_definitions(-DMRI_DIFF_SIM)

  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/mcell-3 ${FFTW_INCLUDE_DIR})

  SET(MPI_DIR $ENV{I_MPI_ROOT})

  if (MPI_DIR)
    message("-- Found MPI:")
    message("--   " ${MPI_DIR})
    add_definitions(-DDO_MPI -DMPICH_SKIP_MPICXX)
    INCLUDE_DIRECTORIES(${CMAKE_INCLUDE_PATH} ${MPI_DIR}/include)
    SET(CMAKE_C_COMPILER ${MPI_DIR}/bin64/mpiicc)
    SET(CMAKE_CXX_COMPILER ${MPI_DIR}/bin64/mpiicpc)
    SET(CMAKE_LINKER ${MPI_DIR}/bin64/mpiicpc)
    message("--  " ${CMAKE_LINKER})
    message("--  " ${CMAKE_CXX_COMPILER})
  else ()
    message("-- Not found MPI")
  endif ()
  
  set(difsim_kernel_SRCS difsim_kernel.cc gtb_histogram.cc sphtes.cc ConfigFile.cc )
  add_executable(difsim ${difsim_kernel_SRCS})
  target_link_libraries (difsim ${MCELL_LIBRARY})
  add_dependencies(difsim mcell-3)

  else ()

  execute_process(COMMAND ${CMAKE_C_COMPILER} -dumpversion
                  OUTPUT_VARIABLE GCC_VERSION)
  if (GCC_VERSION VERSION_GREATER 4.0 OR GCC_VERSION VERSION_EQUAL 4.0)
    message(STATUS "GCC Version >= 4.0: Ok")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
  else(GCC_VERSION VERSION_GREATER 4.0 OR GCC_VERSION VERSION_EQUAL 4.0)
    message(STATUS "GCC Version < 4.0: do not use -std=c++0x flag")
  endif(GCC_VERSION VERSION_GREATER 4.0 OR GCC_VERSION VERSION_EQUAL 4.0)

  include(ExternalProject)

  ExternalProject_Add(mcell-3
  PREFIX ${CMAKE_SOURCE_DIR}/mcell-3
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/mcell-3
  TMP_DIR ${CMAKE_SOURCE_DIR}/mcell-3
  STAMP_DIR ${CMAKE_SOURCE_DIR}/mcell-3
  BINARY_DIR ${CMAKE_SOURCE_DIR}/mcell-3
  CONFIGURE_COMMAND ${CMAKE_SOURCE_DIR}/mcell_configure.sh ${CMAKE_SOURCE_DIR} ${TARGET}
)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/mcell-3 )
  set(MCELL_LIBRARY "${CMAKE_SOURCE_DIR}/mcell-3/libmcell.a")

  FIND_PATH(FFTW_INCLUDE_DIR
  NAMES fftw3.h
  PATHS ${FFTW_POSSIBLE_ROOT_DIRS}
  PATH_SUFFIXES include
  DOC "FFTW header include dir"
  )

  FIND_LIBRARY(FFTW_LIBRARY
  NAMES fftw3 libfftw3
  PATHS  ${FFTW_POSSIBLE_ROOT_DIRS}
  PATH_SUFFIXES lib lib64
  DOC "FFTW library dir" )

  SET(FFTW_LIBRARIES ${FFTW_LIBRARY})

  message("-- Found FFTW libraries:")
  message("--   " ${FFTW_LIBRARY})

  add_definitions(-DMRI_DIFF_SIM)

  IF (NOT APPLE)
    FIND_PATH(MPI_DIR
    NAMES include/mpi.h
    PATHS ${MPI_POSSIBLE_ROOT_DIRS}
    DOC "MPI dir"
    )
  else ()
    FIND_PATH(MPI_DIR
    NAMES include/openmpi-clang/mpi.h
    PATHS ${MPI_POSSIBLE_ROOT_DIRS}
    DOC "MPI dir"
    )
  endif ()
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/mcell-3 ${FFTW_INCLUDE_DIR})

  if (MPI_DIR)
    message("-- Found MPI:")
    message("--   " ${MPI_DIR})
    add_definitions(-DDO_MPI -DMPICH_SKIP_MPICXX)
    INCLUDE_DIRECTORIES(${CMAKE_INCLUDE_PATH} ${MPI_DIR}/include)
    SET(CMAKE_C_COMPILER ${MPI_DIR}/bin/mpicc)
    IF (NOT APPLE)
      SET(CMAKE_CXX_COMPILER ${MPI_DIR}/bin/mpic++)
      SET(CMAKE_LINKER ${MPI_DIR}/bin/mpic++)
    else ()
      SET(CMAKE_CXX_COMPILER ${MPI_DIR}/bin/mpicxx)
      SET(CMAKE_LINKER ${MPI_DIR}/bin/mpicxx)
    endif ()
    message("--  " ${CMAKE_LINKER})
    message("--  " ${CMAKE_CXX_COMPILER})
  else ()
    message("-- Not found MPI")
  endif ()

  set(difsim_kernel_SRCS difsim_kernel.cc gtb_histogram.cc sphtes.cc ConfigFile.cc )
  add_executable(difsim ${difsim_kernel_SRCS})
  target_link_libraries (difsim ${MCELL_LIBRARY} ${FFTW_LIBRARIES})
  add_dependencies(difsim mcell-3)
  endif ()

endif ()

set(small_test_SRCS small_test.cc ConfigFile.cc mcell-3/isaac64.c)
set_source_files_properties(mcell-3/isaac64.c PROPERTIES LANGUAGE CXX )
add_executable(small_test ${small_test_SRCS})
#target_link_libraries (small_test ${MCELL_LIBRARY} )
add_dependencies(small_test mcell-3)

IF (UNIX)
  # since it's unix-specific we will use bash
  ADD_CUSTOM_TARGET (distclean echo cleaning ${CMAKE_SOURCE_DIR} for source distribution)

  ADD_CUSTOM_COMMAND(TARGET distclean
    COMMAND make ARGS clean
    COMMAND find ARGS ${CMAKE_SOURCE_DIR} -depth -type d -iname mcell -exec rm -rf {} \;
    COMMAND find ARGS ${CMAKE_SOURCE_DIR} -iname mcell-3-* -delete
    COMMAND find ARGS ${CMAKE_SOURCE_DIR} -iname mcell -delete
    COMMAND find ARGS ${CMAKE_SOURCE_DIR} -iname ar-lib -delete
    COMMAND find ARGS ${CMAKE_SOURCE_DIR} -iname depcomp -delete
    COMMAND find ARGS ${CMAKE_SOURCE_DIR} -iname ylwrap -delete
    COMMAND find ARGS ${CMAKE_SOURCE_DIR} -iname missing -delete
    COMMAND find ARGS ${CMAKE_SOURCE_DIR} -iname config.log -delete
    COMMAND find ARGS ${CMAKE_SOURCE_DIR} -iname config.h.in -delete
    COMMAND find ARGS ${CMAKE_SOURCE_DIR} -iname aclocal.m4 -delete
    COMMAND find ARGS ${CMAKE_SOURCE_DIR} -iname config.status -delete
    COMMAND find ARGS ${CMAKE_SOURCE_DIR} -iname CMakeCache.txt -delete
    COMMAND find ARGS ${CMAKE_SOURCE_DIR} -iname CTestConfig.cmake -delete
    COMMAND find ARGS ${CMAKE_SOURCE_DIR} -iname CTestTestfile.cmake -delete
    COMMAND find ARGS ${CMAKE_SOURCE_DIR} -iname CPackConfig.cmake -delete
    COMMAND find ARGS ${CMAKE_SOURCE_DIR} -iname CPackSourceConfig.cmake -delete
    COMMAND find ARGS ${CMAKE_SOURCE_DIR} -iname cmake_install.cmake -delete
    COMMAND find ARGS ${CMAKE_SOURCE_DIR} -iname install_manifest.txt -delete
    COMMAND find ARGS ${CMAKE_SOURCE_DIR} -iname Makefile -delete
    COMMAND find ARGS ${CMAKE_SOURCE_DIR} -iname Makefile.in -delete
    COMMAND find ARGS ${CMAKE_SOURCE_DIR} -iname shapeVersion.h -delete
    COMMAND find ARGS ${CMAKE_SOURCE_DIR} -iname shape.doxyfile -delete
    COMMAND find ARGS ${CMAKE_SOURCE_DIR} -depth -type d -iname CMakeFiles -exec rm -rf {} \;
    COMMAND find ARGS ${CMAKE_SOURCE_DIR} -depth -type d -iname autom4te.cache -exec rm -rf {} \;
    COMMAND find ARGS ${CMAKE_SOURCE_DIR} -depth -type d -iname bin -exec rm -rf {} \;
    COMMAND find ARGS ${CMAKE_SOURCE_DIR} -depth -type d -iname lib -exec rm -rf {} \;
    COMMAND find ARGS ${CMAKE_SOURCE_DIR} -depth -type d -iname .deps -exec rm -rf {} \;
    COMMAND find ARGS ${CMAKE_SOURCE_DIR} -depth -type f -iname *.o -exec rm -rf {} \;
    COMMAND find ARGS ${CMAKE_SOURCE_DIR} -depth -type f -iname *.a -exec rm -rf {} \;
    COMMAND find ARGS ${CMAKE_SOURCE_DIR} -depth -type d -iname html -exec rm -rf {} \;
    COMMAND find ARGS ${CMAKE_SOURCE_DIR} -depth -type d -iname latex -exec rm -rf {} \;
    COMMAND rm ARGS -rf config.guess config.sub COPYING INSTALL missing depcomp configure
    VERBATIM )
ENDIF(UNIX)

