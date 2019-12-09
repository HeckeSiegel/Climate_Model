project(fpda_rrtm_wrapper)
enable_language(Fortran C)
cmake_minimum_required (VERSION 2.8.10)

if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set (CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install" CACHE PATH "default install path" FORCE )
endif()

if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
    "Choose the type of build, options are: None Debug Release."
    FORCE)
endif (NOT CMAKE_BUILD_TYPE)

set (CMAKE_C_FLAGS               " -cpp --std=c99                                                                           " )
set (CMAKE_C_FLAGS_RELEASE       " -O3 -ffast-math                                                                          " )
set (CMAKE_C_FLAGS_DEBUG         " -O0 -g -pedantic -pg  -W -Wall -Wuninitialized                                           " )

set (CMAKE_Fortran_FLAGS         " -cpp -ffree-line-length-none                                                             " )
set (CMAKE_Fortran_FLAGS_RELEASE " -fno-backtrace -fno-range-check -O3                                                      " )
set (CMAKE_Fortran_FLAGS_DEBUG   " -finit-real=nan -fcheck=all -fbacktrace -O0 -g -pedantic -Wsurprising -pg -fbounds-check -fno-range-check -W -Wall -Wuninitialized -g" )


set(INSTALL_LIB_DIR     lib     CACHE PATH "Installation directory for libraries")
set(INSTALL_BIN_DIR     bin     CACHE PATH "Installation directory for executables")
set(INSTALL_INCLUDE_DIR include CACHE PATH "Installation directory for header files")

add_subdirectory (gnuplot_i)
add_subdirectory (rrtm_lw)
add_subdirectory (rrtm_sw)

include_directories("${PROJECT_SOURCE_DIR}/rrtm_lw")
include_directories("${PROJECT_SOURCE_DIR}/rrtm_sw")
add_subdirectory (src)

