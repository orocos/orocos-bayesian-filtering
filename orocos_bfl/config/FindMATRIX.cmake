# Locate MATRIX LIB install directory

# This module defines
# MATRIX_FOUND, is set to true


# variables
# ---------
IF (NOT __MATRIXWRAPPER_NEWMAT__)
  SET(__MATRIXWRAPPER_NEWMAT__ OFF CACHE BOOL "define for newmat")
  MARK_AS_ADVANCED(__MATRIXWRAPPER_NEWMAT__)
ENDIF (NOT __MATRIXWRAPPER_NEWMAT__)
SET(__MATRIXWRAPPER_NEWMAT__ OFF)

IF (NOT __MATRIXWRAPPER_BOOST__)
  SET(__MATRIXWRAPPER_BOOST__ OFF CACHE BOOL "define for boost")
  MARK_AS_ADVANCED(__MATRIXWRAPPER_BOOST__)
ENDIF (NOT __MATRIXWRAPPER_BOOST__)
SET(__MATRIXWRAPPER_BOOST__ OFF)

IF (NOT __MATRIXWRAPPER_EIGEN__)
  SET(__MATRIXWRAPPER_EIGEN__ OFF CACHE BOOL "define for eigen")
  MARK_AS_ADVANCED(__MATRIXWRAPPER_EIGEN__)
ENDIF (NOT __MATRIXWRAPPER_EIGEN__)
SET(__MATRIXWRAPPER_EIGEN__ OFF)

# install path
# ------------
IF(NOT MATRIX_LIB)
  SET( MATRIX_LIB boost CACHE STRING "Which matrix library to use: newmat, boost or eigen")
ENDIF(NOT MATRIX_LIB)
MESSAGE("Searching for matrix lib ${MATRIX_LIB}")


# find libs
# ---------
IF (MATRIX_LIB STREQUAL "newmat")
  FIND_LIBRARY(NEWMAT newmat )
  FIND_PATH(NEWMAT_FOUND newmat/include.h )
  IF ( NEWMAT AND NEWMAT_FOUND )
    MESSAGE("-- Looking for Newmat - found")
    SET( MATRIX_INCLUDE "${NEWMAT_FOUND}")
    SET( MATRIX_LIBS "${NEWMAT}")
    MESSAGE( "-- Newmat includes ${MATRIX_INCLUDE}")
    MESSAGE( "-- Newmat libs     ${MATRIX_LIBS}")
    SET(__MATRIXWRAPPER_NEWMAT__ ON)
  ELSE ( NEWMAT AND NEWMAT_FOUND )
    MESSAGE( FATAL_ERROR "Looking for Newmat - not found")
  ENDIF ( NEWMAT AND NEWMAT_FOUND )
ELSE (MATRIX_LIB STREQUAL "newmat")


IF (MATRIX_LIB STREQUAL "boost")
  FIND_PATH(BOOST_FOUND boost/numeric/ublas/matrix.hpp )
  IF ( BOOST_FOUND )
    MESSAGE("-- Looking for Boost - found")
    SET( MATRIX_INCLUDE "${BOOST_FOUND}")
    SET( MATRIX_LIBS "")
    MESSAGE( "-- Boost includes ${MATRIX_INCLUDE}")
    MESSAGE( "-- Boost libs     ${MATRIX_LIBS}")
    SET(__MATRIXWRAPPER_BOOST__ ON)
  ELSE ( BOOST_FOUND )
    MESSAGE(FATAL_ERROR "Looking for Boost in - not found")
  ENDIF ( BOOST_FOUND )
ELSE (MATRIX_LIB STREQUAL "boost")


IF (MATRIX_LIB STREQUAL "eigen")
  IF ( NOT EIGEN_INCLUDE )
    FIND_PATH(EIGEN_FOUND Eigen/Core
      PATH_SUFFIXES eigen3 include/eigen3 include)
    IF ( EIGEN_FOUND )
      MESSAGE("-- Looking for Eigen - found")
      SET( EIGEN_INCLUDE "${EIGEN_FOUND}")
      SET( MATRIX_LIBS "")
    ELSE ( EIGEN_FOUND )
      MESSAGE(FATAL_ERROR "Looking for Eigen in - not found")
    ENDIF ( EIGEN_FOUND )
  ENDIF ( NOT EIGEN_INCLUDE )

  SET(__MATRIXWRAPPER_EIGEN__ ON)
  SET( MATRIX_INCLUDE "${EIGEN_INCLUDE}")
  MESSAGE( "-- Eigen includes ${EIGEN_INCLUDE}")
ELSE (MATRIX_LIB STREQUAL "eigen")

MESSAGE( FATAL_ERROR "No valid matrix lib specified. Please choose lti, newmat, boost or eigen")

ENDIF (MATRIX_LIB STREQUAL "eigen")
ENDIF (MATRIX_LIB STREQUAL "boost")
ENDIF (MATRIX_LIB STREQUAL "newmat")
