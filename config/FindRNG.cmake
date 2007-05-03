# Locate RNG LIB install directory

# This module defines
# RNG_INSTALL where to find include, lib, bin, etc.
# RNG_FOUND, is set to true

MESSAGE("Searching for rng lib ${RNG_LIB}")

SET(__RNGWRAPPER_LTI__ OFF)
SET(__RNGWRAPPER_BOOST__ OFF)


IF (RNG_LIB STREQUAL "lti")
  SET(LTI_FOUND LTI_FOUND-NOTFOUND)
  MARK_AS_ADVANCED(LTI_FOUND)
  FIND_FILE(LTI_FOUND ltiMatrix.h ${RNG_INSTALL}/include/ltilib/ )
  IF ( LTI_FOUND )
    MESSAGE("-- Looking for Lti in ${RNG_INSTALL}/include/ltilib/ - found")
    SET( RNG_INCLUDE "${RNG_INSTALL}/include/")
    SET( RNG_LIBS "-L${RNG_INSTALL}/lib/ltilib -lltir") 
    MESSAGE( "-- Lti includes ${RNG_INCLUDE}")
    MESSAGE( "-- Lti libs     ${RNG_LIBS}")
    SET(__RNGWRAPPER_LTI__ ON)	
  ELSE ( LTI_FOUND )
    MESSAGE(FATAL_ERROR "Looking for Lti in ${RNG_INSTALL}/include/ltilib/ - not found")
  ENDIF ( LTI_FOUND )
ELSE (RNG_LIB STREQUAL "lti")


IF (RNG_LIB STREQUAL "boost")
  SET(BOOST_FOUND BOOST_FOUND-NOTFOUND)
  MARK_AS_ADVANCED(BOOST_FOUND)
  FIND_FILE(BOOST_FOUND mersenne_twister.hpp ${RNG_INSTALL}/include/boost/random/ )
  IF ( BOOST_FOUND )
    MESSAGE("-- Looking for Boost in ${RNG_INSTALL}/include/boost/ - found")
    SET( RNG_INCLUDE "${RNG_INSTALL}/include/")
    SET( RNG_LIBS "") 
    MESSAGE( "-- Boost includes ${RNG_INCLUDE}")
    MESSAGE( "-- Boost libs     ${RNG_LIBS}")
    SET(__RNGWRAPPER_BOOST__ ON)	
  ELSE ( BOOST_FOUND )
    MESSAGE(FATAL_ERROR "Looking for Boost in ${RNG_INSTALL}/include/boost/ - not found")
  ENDIF ( BOOST_FOUND )
ELSE (RNG_LIB STREQUAL "boost")


MESSAGE( FATAL_ERROR "No valid rng lib specified. Please choose lti or boost")


ENDIF (RNG_LIB STREQUAL "boost")
ENDIF (RNG_LIB STREQUAL "lti")
