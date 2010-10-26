# Locate GINAC install directory

# This module defines
# GINAC location and name of library
# GINAC_FOUND, path to ginac header files

OPTION( GINAC_SUPPORT "Turn me off to disable Ginac support" OFF )

IF (GINAC_SUPPORT)

  MESSAGE("Searching for ginac lib")

  FIND_LIBRARY(GINAC ginac)
  FIND_PATH(GINAC_FOUND ginac/ginac.h )
  IF ( GINAC AND GINAC_FOUND )
    MESSAGE("-- Looking for Ginac - found")
    SET( GINAC_INCLUDE "${GINAC_FOUND}")
    SET( GINAC_LIBS "${GINAC}") 
    MESSAGE( "-- Ginac includes ${GINAC_INCLUDE}")
    MESSAGE( "-- Ginac libs     ${GINAC_LIBS}")
  ELSE ( GINAC AND GINAC_FOUND )
    MESSAGE(FATAL_ERROR "Looking for Ginac - not found")
  ENDIF ( GINAC AND GINAC_FOUND )

ENDIF (GINAC_SUPPORT)