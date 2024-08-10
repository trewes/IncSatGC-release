# This module finds the cryptominisat header and libraries.
#
# User can give CRYPTOMINISAT_ROOT_DIR as a hint stored in the cmake cache.
#
# It sets the following variables:
#  CRYPTOMINISAT_FOUND              - Set to false, or undefined, if cryptominisat isn't found.
#  CRYPTOMINISAT_INCLUDE_DIR        - include directory
#  CRYPTOMINISAT_LIBRARY            - library files

# use given hint directory or look in parent/cryptominisat folder
IF(NOT DEFINED CRYPTOMINISAT_ROOT_DIR)
    SET(CRYPTOMINISAT_ROOT_DIR ${CMAKE_SOURCE_DIR}/../cryptominisat)
ENDIF()

FIND_PATH(CRYPTOMINISAT_INCLUDE_DIR cryptominisat.h
        PATH_SUFFIXES src
        PATHS ${CRYPTOMINISAT_ROOT_DIR}
)
MESSAGE(STATUS "CRYPTOMINISAT Include Dir: ${CRYPTOMINISAT_INCLUDE_DIR}")


FIND_LIBRARY(CRYPTOMINISAT_LIBRARY
        NAMES cryptominisat cryptominisat5 libcryptominisat5.so
        PATH_SUFFIXES lib build build/lib
        PATHS ${CRYPTOMINISAT_ROOT_DIR}
)
MESSAGE(STATUS "CRYPTOMINISAT Library: ${CRYPTOMINISAT_LIBRARY}")


# check whether required things have been found and set CRYPTOMINISAT_FOUND accordingly
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(CRYPTOMINISAT DEFAULT_MSG
        CRYPTOMINISAT_INCLUDE_DIR CRYPTOMINISAT_LIBRARY)

