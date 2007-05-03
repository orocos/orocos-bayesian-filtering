# Locate Orocos::RTT install directory

# This module defines
# OROCOS_RTT_HOME where to find include, lib, bin, etc.
# OROCOS_RTT_FOUND, is set to true

INCLUDE (${PROJ_SOURCE_DIR}/config/FindPkgConfig.cmake)

IF ( CMAKE_PKGCONFIG_EXECUTABLE AND NOT SKIP_BUILD)

    MESSAGE( STATUS "Detecting RTT" )
    
    SET(ENV{PKG_CONFIG_PATH} "${OROCOS_INSTALL}/lib/pkgconfig:${OROCOS_INSTALL}/packages/install/lib/pkgconfig/")
    #MESSAGE( "Setting environment of PKG_CONFIG_PATH to: $ENV{PKG_CONFIG_PATH}")
    MESSAGE( "Searching RTT in ${OROCOS_INSTALL}:" )

    # if OROCOS_RTT_1.2, enable SHARED libraries (see component_rules).
    PKGCONFIG( "orocos-rtt >= 1.1.0" OROCOS_RTT_1.2 OROCOS_RTT_INCLUDE_DIRS_1.2 OROCOS_RTT_DEFINES_1.2 OROCOS_RTT_LINK_DIRS_1.2 OROCOS_RTT_LIBS_1.2 )
    PKGCONFIG( "orocos-rtt >= 1.0.0" OROCOS_RTT OROCOS_RTT_INCLUDE_DIRS OROCOS_RTT_DEFINES OROCOS_RTT_LINK_DIRS OROCOS_RTT_LIBS )

    IF( OROCOS_RTT )
        MESSAGE("   Includes in: ${OROCOS_RTT_INCLUDE_DIRS}")
        MESSAGE("   Libraries in: ${OROCOS_RTT_LINK_DIRS}")
        MESSAGE("   Libraries: ${OROCOS_RTT_LIBS}")
        MESSAGE("   Defines: ${OROCOS_RTT_DEFINES}")

	INCLUDE_DIRECTORIES( ${OROCOS_RTT_INCLUDE_DIRS} )
        LINK_DIRECTORIES( ${OROCOS_RTT_LINK_DIRS} )

	# Detect OS:
	# FIND_XXX has a special treatment if the variable's value is XXX-NOTFOUND
	# This 'value' will trigger a new search in the next cmake run.
	SET(OS_GNULINUX OS_GNULINUX-NOTFOUND)
	SET(OS_XENOMAI OS_XENOMAI-NOTFOUND)
	SET(OS_LXRT OS_LXRT-NOTFOUND)
	#SET( OS_GNULINUX 1 CACHE INTERNAL "")
	#CHECK_INCLUDE_FILE( rtt/os/gnulinux.h OS_GNULINUX ${OROCOS_RTT_INCLUDE_DIRS})
	FIND_FILE( OS_GNULINUX rtt/os/gnulinux.h ${OROCOS_RTT_INCLUDE_DIRS})
	FIND_FILE( OS_GNULINUX rtt/os/gnulinux/fosi.h ${OROCOS_RTT_INCLUDE_DIRS})
        IF(OS_GNULINUX )
           MESSAGE( "Detected GNU/Linux installation: ${OS_GNULINUX}" )
        ENDIF(OS_GNULINUX)

        SET( CMAKE_REQUIRED_INCLUDES ${OROCOS_RTT_INCLUDE_DIRS})
	FIND_FILE( OS_XENOMAI rtt/os/xenomai.h ${OROCOS_RTT_INCLUDE_DIRS})
	FIND_FILE( OS_XENOMAI rtt/os/xenomai/fosi.h ${OROCOS_RTT_INCLUDE_DIRS})
        IF(OS_XENOMAI)
           MESSAGE( "Detected Xenomai installation: ${OS_XENOMAI}" )
	ENDIF(OS_XENOMAI)

	SET( CMAKE_REQUIRED_INCLUDES ${OROCOS_RTT_INCLUDE_DIRS})
 	FIND_FILE( OS_LXRT rtt/os/lxrt.h ${OROCOS_RTT_INCLUDE_DIRS})
 	FIND_FILE( OS_LXRT rtt/os/lxrt/fosi.h ${OROCOS_RTT_INCLUDE_DIRS})
        IF(OS_LXRT)
           MESSAGE( "Detected LXRT installation: ${OS_LXRT}" )
        ENDIF(OS_LXRT)

	IF ( NOT OS_LXRT AND NOT OS_GNULINUX AND NOT OS_XENOMAI )
	  MESSAGE( FATAL_ERROR "Could not detect target OS!" )
	ENDIF ( NOT OS_LXRT AND NOT OS_GNULINUX AND NOT OS_XENOMAI )

	FIND_FILE( CORBA_ENABLED rtt/corba/ControlTaskProxy.hpp ${OROCOS_RTT_INCLUDE_DIRS} )
	IF (CORBA_ENABLED)
	  MESSAGE( "Detected CORBA build of RTT: ${CORBA_ENABLED}" )
	ELSE(CORBA_ENABLED)
	  MESSAGE( "CORBA not detected." )
	ENDIF (CORBA_ENABLED)

    ELSE  ( OROCOS_RTT )
        MESSAGE( FATAL_ERROR "Can't find Orocos Real-Time Toolkit (orocos-rtt.pc)")
    ENDIF ( OROCOS_RTT )

ELSE  ( CMAKE_PKGCONFIG_EXECUTABLE  AND NOT SKIP_BUILD)

    IF (NOT CMAKE_PKGCONFIG_EXECUTABLE)
    # Can't find pkg-config -- have to search manually
    MESSAGE( FATAL_ERROR "Can't find pkg-config ")
    ENDIF (NOT CMAKE_PKGCONFIG_EXECUTABLE)

ENDIF ( CMAKE_PKGCONFIG_EXECUTABLE  AND NOT SKIP_BUILD)
