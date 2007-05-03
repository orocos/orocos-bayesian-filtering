# Locate BFL install directory

# This module defines
# BLF_INSTALL where to find include, lib, bin, etc.
# BLF_FOUND, is set to true

INCLUDE (${PROJ_SOURCE_DIR}/config/FindPkgConfig.cmake)

IF ( CMAKE_PKGCONFIG_EXECUTABLE )

    MESSAGE( STATUS "Detecting gthread-2.0" )
    PKGCONFIG( "gthread-2.0" GTHREAD_FOUND GTHREAD_INCLUDE_DIRS GTHREAD_DEFINES GTHREAD_LINK_DIRS GTHREAD_LIBS )

    IF( GTHREAD_FOUND )
        MESSAGE("   Includes in: ${GTHREAD_INCLUDE_DIRS}")
        MESSAGE("   Libraries in: ${GTHREAD_LINK_DIRS}")
        MESSAGE("   Libraries: ${GTHREAD_LIBS}")
        MESSAGE("   Defines: ${GTHREAD_DEFINES}")

	INCLUDE_DIRECTORIES( ${GTHREAD_INCLUDE_DIRS} )
    LINK_DIRECTORIES( ${GTHREAD_LINK_DIRS} )

    ENDIF ( GTHREAD_FOUND )

ELSE  ( CMAKE_PKGCONFIG_EXECUTABLE )

    # Can't find pkg-config -- have to search manually
    MESSAGE( FATAL_ERROR "Can't find GTHREAD-2.0 without pkgconfig !")

ENDIF ( CMAKE_PKGCONFIG_EXECUTABLE )
