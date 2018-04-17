# Locate RNG LIB install directory

# This module defines
# RNG_FOUND, is set to true

# variables
# ---------
IF (NOT __RNGWRAPPER_BOOST__)
  SET(__RNGWRAPPER_BOOST__ OFF CACHE BOOL "define for boost")
  MARK_AS_ADVANCED(__RNGWRAPPER_BOOST__)
ENDIF (NOT __RNGWRAPPER_BOOST__)
SET(__RNGWRAPPER_BOOST__ OFF)

# install path
# ------------
IF(NOT RNG_LIB)
  SET( RNG_LIB boost CACHE STRING "Which rng library to use: boost")
ENDIF(NOT RNG_LIB)
MESSAGE("Searching for rng lib ${RNG_LIB}")


IF (RNG_LIB STREQUAL "boost")
  FIND_PATH(BOOST_FOUND boost/random/mersenne_twister.hpp )
  IF ( BOOST_FOUND )
    MESSAGE("-- Looking for Boost - found")
    SET( RNG_INCLUDE "${BOOST_FOUND}")
    SET( RNG_LIBS "") 
    MESSAGE( "-- Boost includes ${RNG_INCLUDE}")
    MESSAGE( "-- Boost libs     ${RNG_LIBS}")
    SET(__RNGWRAPPER_BOOST__ ON)	
  ELSE ( BOOST_FOUND )
    MESSAGE(FATAL_ERROR "Looking for Boost - not found")
  ENDIF ( BOOST_FOUND )
ELSE (RNG_LIB STREQUAL "boost")
 
MESSAGE( FATAL_ERROR "No valid rng lib specified. Please choose boost")
 
ENDIF (RNG_LIB STREQUAL "boost")
