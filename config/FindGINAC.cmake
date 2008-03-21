# Locate GINAC install directory

# This module defines
# GINAC_INSTALL where to find include, lib, bin, etc.
# GINAC_FOUND, is set to true

OPTION( GINAC_SUPPORT "Turn me off to disable Ginac support" OFF )
IF(NOT GINAC_INSTALL)
  SET( GINAC_INSTALL /usr CACHE PATH "The Ginac lib installation directory.")
ENDIF(NOT GINAC_INSTALL)

IF (GINAC_SUPPORT)

  MESSAGE("Searching for ginac lib")

  SET(GINAC_FOUND GINAC_FOUND-NOTFOUND)
  MARK_AS_ADVANCED(GINAC_FOUND)
  FIND_FILE(GINAC_FOUND ginac.h ${GINAC_INSTALL}/include/ginac/ )
  IF ( GINAC_FOUND )
    MESSAGE("-- Looking for Ginac in ${GINAC_INSTALL}/include/ginac/ - found")
    SET( GINAC_INCLUDE "-I${GINAC_INSTALL}/include/")
    SET( GINAC_LIBS "-L ${GINAC_INSTALL}/lib/ -lginac") 
    MESSAGE( "-- Ginac includes ${GINAC_INCLUDE}")
    MESSAGE( "-- Ginac libs     ${GINAC_LIBS}")
  ELSE ( GINAC_FOUND )
    MESSAGE(FATAL_ERROR "Looking for Ginac in ${GINAC_INSTALL}/include/ginac/ - not found")
  ENDIF ( GINAC_FOUND )

ENDIF (GINAC_SUPPORT)