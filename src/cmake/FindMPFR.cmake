# - Try to find MPFR
# Once done this will define
#  MPFR_FOUND - System has MPFR
#  MPFR_INCLUDE_DIRS - The MPFR include directories
#  MPFR_LIBRARIES - The libraries needed to use MPFR
#  MPFR_DEFINITIONS - Compiler switches required for using MPFR
#  MPFR_VERSION - version MPFR


find_package(PkgConfig)
pkg_check_modules(PC_MPFR QUIET mpfr)
# MPFR header includes GMP header
pkg_check_modules(PC_GMP QUIET gmp)

set(MPFR_DEFINITIONS ${PC_MPFR_CFLAGS_OTHER})

find_path(MPFR_INCLUDE_DIR mpfr.h
          HINTS ${PC_MPFR_INCLUDEDIR} ${PC_MPFR_INCLUDE_DIRS}
                ${MPFR_DIR}/include $ENV{MPFR_DIR}/include
                ${MPFR_ROOT}/include $ENV{MPFR_ROOT}/include
          PATH_SUFFIXES mpfr
          DOC "Directory where the MPFR header files are located"
          )
find_path(GMP_INCLUDE_DIR gmp.h
          HINTS ${PC_GMP_INCLUDEDIR} ${PC_GMP_INCLUDE_DIRS}
                ${GMP_DIR}/include $ENV{GMP_DIR}/include
                ${GMP_ROOT}/include $ENV{GMP_ROOT}/include
          PATH_SUFFIXES gmp
          DOC "Directory where the GMP header files are located"
          )

find_library(MPFR_LIBRARY NAMES mpfr
             HINTS ${PC_MPFR_LIBDIR} ${PC_MPFR_LIBRARY_DIRS}
                    ${MPFR_DIR}/lib $ENV{MPFR_DIR}/lib
                    ${MPFR_ROOT}/lib $ENV{MPFR_ROOT}/lib
             DOC "Directory where the MPFR library is located"
)

# Try compiling and running test program
# Set flags for building test program
set(CMAKE_REQUIRED_INCLUDES  ${MPFR_INCLUDE_DIR} ${GMP_INCLUDE_DIR})
set(CMAKE_REQUIRED_LIBRARIES ${MPFR_LIBRARY})

# Check MPFR version
set(MPFR_CONFIG_TEST_VERSION_C
  "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/mpfr_config_test_version.c")
file(WRITE ${MPFR_CONFIG_TEST_VERSION_C} "
#include <stdio.h>
#include \"mpfr.h\"

int main() {
printf(\"%d.%d.%d\",
  MPFR_VERSION_MAJOR,
  MPFR_VERSION_MINOR,
  MPFR_VERSION_PATCHLEVEL);
return 0;
}
")

  try_run(
    MPFR_CONFIG_TEST_VERSION_EXITCODE
    MPFR_CONFIG_TEST_VERSION_COMPILED
    ${CMAKE_CURRENT_BINARY_DIR}
    ${MPFR_CONFIG_TEST_VERSION_C}
    CMAKE_FLAGS
      "-DINCLUDE_DIRECTORIES:STRING=${CMAKE_REQUIRED_INCLUDES}"
  "-DLINK_LIBRARIES:STRING=${CMAKE_REQUIRED_LIBRARIES}"
    COMPILE_OUTPUT_VARIABLE MPFR_CONFIG_TEST_VERSION_COMPILE_OUTPUT
    RUN_OUTPUT_VARIABLE MPFR_CONFIG_TEST_VERSION_OUTPUT
    )

  if (MPFR_CONFIG_TEST_VERSION_EXITCODE EQUAL 0)
    set(MPFR_VERSION ${MPFR_CONFIG_TEST_VERSION_OUTPUT} CACHE STRING "MPFR version")
    mark_as_advanced(MPFR_VERSION)
  endif()

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set MPFR_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(
    MPFR
    REQUIRED_VARS MPFR_LIBRARY MPFR_INCLUDE_DIR GMP_INCLUDE_DIR
    VERSION_VAR MPFR_VERSION
)

mark_as_advanced(MPFR_INCLUDE_DIRS MPFR_LIBRARY )

set(MPFR_LIBRARIES ${MPFR_LIBRARY} )
set(MPFR_INCLUDE_DIRS ${MPFR_INCLUDE_DIR} ${GMP_INCLUDE_DIR})
