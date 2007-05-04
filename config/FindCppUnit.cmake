# Locate CPPUNIT library install directory

# This module defines
# CPPUNIT_INSTALL where to find include, lib, bin, etc.
# CPPUNIT_FOUND, is set to true

MESSAGE("Searching for cppunit lib ${CPPUNIT_LIB}")

SET(CPPUNIT_FOUND CPPUNIT_FOUND-NOTFOUND)
MARK_AS_ADVANCED(CPPUNIT_FOUND)
FIND_FILE(CPPUNIT_FOUND Asserter.h ${CPPUNIT_INSTALL}/include/cppunit/ )
IF ( CPPUNIT_FOUND )
  MESSAGE("-- Looking for Cppunit in ${CPPUNIT_INSTALL}/include/cppunit/ - found")
  SET( CPPUNIT_INCLUDE "${CPPUNIT_INSTALL}/include/")
  SET( CPPUNIT_LIBS "-L${CPPUNIT_INSTALL}/lib -lcppunit") 
  MESSAGE( "-- Cppunit includes ${CPPUNIT_INCLUDE}")
  MESSAGE( "-- Cppunit libs     ${CPPUNIT_LIBS}")
ELSE ( CPPUNIT_FOUND )
  MESSAGE( ERROR "Looking for Cppunit in ${CPPUNIT_INSTALL}/include/cppunit/ - not found")
ENDIF ( CPPUNIT_FOUND )
