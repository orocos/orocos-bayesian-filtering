# Locate CPPUNIT library install directory

# This module defines
# CPPUNIT_INSTALL where to find include, lib, bin, etc.
# CPPUNIT_FOUND, is set to true

MESSAGE("Searching for cppunit lib ${CPPUNIT_LIB}")

FIND_LIBRARY(CPPUNIT cppunit )
FIND_PATH(CPPUNIT_FOUND  cppunit/Asserter.h )
IF ( CPPUNIT AND CPPUNIT_FOUND )
  MESSAGE("-- Looking for Cppunit - found")
  SET( CPPUNIT_INCLUDE "${CPPUNIT_FOUND}")
  SET( CPPUNIT_LIBS "${CPPUNIT}") 
  MESSAGE( "-- Cppunit includes ${CPPUNIT_INCLUDE}")
  MESSAGE( "-- Cppunit libs     ${CPPUNIT_LIBS}")
ELSE ( CPPUNIT AND CPPUNIT_FOUND )
  MESSAGE( ERROR )
  MESSAGE("-- Looking for Cppunit - not found")
  MESSAGE("-- You will not be able to build tests. To build tests, first set the path where cppunit can be found using CMAKE_INCLUDE_PATH and CMAKE_LIBRARY_PATH.")
ENDIF ( CPPUNIT AND CPPUNIT_FOUND )
