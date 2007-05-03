CMAKE_MINIMUM_REQUIRED(VERSION 2.2)
#MARK_AS_ADVANCED( FORCE CMAKE_BACKWARDS_COMPATIBILITY )

###################################################
#                                                 #
#       Start project customization section       #
#                                                 #
###################################################

PROJECT(orocos-bfl)
SET( VERSION 0.4.9 CACHE STRING "Version of Orocos BFL" )
MARK_AS_ADVANCED(VERSION)
MESSAGE( "Orocos BFL version ${VERSION}" )
MARK_AS_ADVANCED(LIBRARY_OUTPUT_PATH)
MARK_AS_ADVANCED(EXECUTABLE_OUTPUT_PATH)
MARK_AS_ADVANCED(CMAKE_BACKWARDS_COMPATIBILITY)

# Make compilation not verbose by default
SET(CMAKE_VERBOSE_MAKEFILE FALSE)


# project dirs
# ------------
IF(NOT CMAKE_INSTALL_PREFIX)
  SET( CMAKE_INSTALL_PREFIX /usr/local/ CACHE PATH "Installation directory" FORCE)
  MESSAGE( "Setting installation directory to ${CMAKE_INSTALL_PREFIX}" )
ENDIF(NOT CMAKE_INSTALL_PREFIX)

SET( PROJ_SOURCE_DIR ${orocos-bfl_SOURCE_DIR} )

SET( PROJ_BINARY_DIR ${orocos-bfl_BINARY_DIR} )


# create defines for matrix and rng
# ---------------------------------
IF (NOT __MATRIXWRAPPER_NEWMAT__)
  SET(__MATRIXWRAPPER_NEWMAT__ OFF CACHE BOOL "define for newmat")
  MARK_AS_ADVANCED(__MATRIXWRAPPER_NEWMAT__)
ENDIF (NOT __MATRIXWRAPPER_NEWMAT__)

IF (NOT __MATRIXWRAPPER_LTI__)
  SET(__MATRIXWRAPPER_LTI__ OFF CACHE BOOL "define for lti")
  MARK_AS_ADVANCED(__MATRIXWRAPPER_LTI__)
ENDIF (NOT __MATRIXWRAPPER_LTI__)

IF (NOT __RNGWRAPPER_BOOST__)
  SET(__RNGWRAPPER_BOOST__ OFF CACHE BOOL "define for boost")
  MARK_AS_ADVANCED(__RNGWRAPPER_BOOST__)
ENDIF (NOT __RNGWRAPPER_BOOST__)

IF (NOT __RNGWRAPPER_LTI__)
  SET(__RNGWRAPPER_LTI__ OFF CACHE BOOL "define for lti")
  MARK_AS_ADVANCED(__RNGWRAPPER_LTI__)
ENDIF (NOT __RNGWRAPPER_LTI__)


# install dir of matrix, rng and ginac
# ------------------------------------
IF(NOT MATRIX_LIB)
  SET( MATRIX_LIB lti CACHE STRING "Which matrix library to use: lti or newmat")
ENDIF(NOT MATRIX_LIB)
IF(NOT MATRIX_INSTALL)
  SET( MATRIX_INSTALL /usr CACHE PATH "The Matrix lib installation directory.")
ENDIF(NOT MATRIX_INSTALL)

IF(NOT RNG_LIB)
  SET( RNG_LIB lti CACHE STRING "Which rng library to use: lti or boost")
ENDIF(NOT RNG_LIB)
IF(NOT RNG_INSTALL)
  SET( RNG_INSTALL /usr CACHE PATH "The rng lib installation directory.")
ENDIF(NOT RNG_INSTALL)

IF (NOT GINAC_SUPPORT)
  SET(GINAC_SUPPORT OFF CACHE BOOL "try to compile ginac dependant classes ")
ENDIF (NOT GINAC_SUPPORT)

IF(NOT GINAC_INSTALL)
  SET( GINAC_INSTALL /usr CACHE PATH "The Ginac lib installation directory.")
ENDIF(NOT GINAC_INSTALL)


# Set the build type (affects debugging symbols and optimization)
# ------------------
IF ( NOT CMAKE_BUILD_TYPE )
  SET( CMAKE_BUILD_TYPE "Debug" CACHE STRING "Build type: None (Use CMAKE_C_FLAGS and CMAKE_CXX_FLAGS), Release, Debug, RelWithDebInfo, MinSizeRel." FORCE )
  MESSAGE( "Setting build type to '${CMAKE_BUILD_TYPE}'" )
ELSE ( NOT CMAKE_BUILD_TYPE )
  MESSAGE( "Build type set to '${CMAKE_BUILD_TYPE}' by user." )
  IF ( NOT CMAKE_BUILD_TYPE STREQUAL "BFL" )
    MARK_AS_ADVANCED(FORCE CMAKE_C_FLAGS_BFL CMAKE_CXX_FLAGS_BFL)
  ELSE ( NOT CMAKE_BUILD_TYPE STREQUAL "BFL" )
    MARK_AS_ADVANCED(CLEAR CMAKE_C_FLAGS_BFL CMAKE_CXX_FLAGS_BFL)
  ENDIF ( NOT CMAKE_BUILD_TYPE STREQUAL "BFL" )
ENDIF ( NOT CMAKE_BUILD_TYPE )

# MESSAGE("Looking for compiler:")
INCLUDE (${CMAKE_ROOT}/Modules/CMakeDetermineCXXCompiler.cmake)



###################################################
#                                                 #
#       Do some global set-up & checking          #
#                                                 #
###################################################

#
# This is the first in the include path.
INCLUDE_DIRECTORIES( ${PROJ_SOURCE_DIR} )

# Put your 'flags' inhere if you want them to apear in the
# orocos-BFL.pc file. Various .cmake scripts below do this.
SET( ENV{OROCOS_BFL_CFLAGS} "")
SET( ENV{OROCOS_BFL_LINKFLAGS} "")

INCLUDE(${PROJ_SOURCE_DIR}/config/global_setup.cmake)

## Find MATRIX LIB
INCLUDE (${PROJ_SOURCE_DIR}/config/FindMATRIX.cmake)

## Find RNG LIB
INCLUDE (${PROJ_SOURCE_DIR}/config/FindRNG.cmake)

## Find GINAC
INCLUDE (${PROJ_SOURCE_DIR}/config/FindGINAC.cmake)

## Uninstall target
CONFIGURE_FILE(
  "${CMAKE_CURRENT_SOURCE_DIR}/config/cmake_uninstall.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
  IMMEDIATE @ONLY)

ADD_CUSTOM_TARGET(uninstall
  "${CMAKE_COMMAND}" -P "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake")


###########################################################
#                                                         #
# Look for dependencies required by individual components #
#                                                         #
###########################################################

INCLUDE(config/check_depend.cmake)
INCLUDE(config/bfl_macros.cmake)

# First unset the LOCAL_LIBRARIES
# -------------------------------
SET(ENV{SELECTED_LIBS} "")
SET(ENV{SELECTED_DIRS} "")
ADD_SUBDIRECTORY( src )
ADD_SUBDIRECTORY( tests )
ADD_SUBDIRECTORY( examples )

# configure files
# ---------------
MESSAGE("Selected Libraries: $ENV{SELECTED_LIBS}")
SET(SELECTED_LIBS $ENV{SELECTED_LIBS})
SET(OROCOS_BFL_LINKFLAGS "$ENV{OROCOS_BFL_LINKFLAGS}")
SET(OROCOS_BFL_CFLAGS "$ENV{OROCOS_BFL_CFLAGS}")
CONFIGURE_FILE( orocos-bfl.pc.in orocos-bfl.pc @ONLY)
INSTALL_FILES( /lib/pkgconfig FILES orocos-bfl.pc)
CONFIGURE_FILE(Doxyfile.in Doxyfile @ONLY)



############################################################
##                                                         #
## Setup PACKAGING   THIS IS EXPERIMENTAL/BETA/NOT USED    #
##                                                         #
############################################################

IF ( "${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" STREQUAL "2.4" )
OPTION( CPACK_PACKAGES "Set to ON to build the packages. Requires cmake >2.4" OFF )
MARK_AS_ADVANCED( CPACK_PACKAGES )
IF (CPACK_PACKAGES)

INCLUDE(InstallRequiredSystemLibraries)

# ATTENTION: There is sometimes a _SOURCE_ version of an
# option as well, set both if necessary !

# Create .tar.gz and .tar.tbz2 files:
SET(CPACK_GENERATOR "TBZ2")
SET(CPACK_SOURCE_GENERATOR "TBZ2")

# The plain 'package' target works correctly.
SET(CPACK_IGNORE_FILES        "/CVS/;/.svn/;.swp$;.#;/#;/build/")
# Since the 'package_source' target does a bold copy, define a list of
# files which should be excluded. Note that 'subpattern' matching is used,
# thus to exclude a directory use /mydir/
SET(CPACK_SOURCE_IGNORE_FILES "/CVS/;/.svn/;.swp$;.#;/#;/build/;~")

SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Orocos BFL")
SET(CPACK_PACKAGE_VENDOR "The Orocos Community")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/COPYING")
SET(CPACK_PACKAGE_MAJOR_VERSION "0")
SET(CPACK_PACKAGE_MINOR_VERSION "1")
SET(CPACK_PACKAGE_VERSION_PATCH "1")

SET(CPACK_PACKAGE_EXECUTABLES "OrocosExec" "Orocos Executable")

SET(CPACK_PACKAGE_INSTALL_DIRECTORY "${CMAKE_INSTALL_PREFIX}")
SET(CPACK_SOURCE_TOPLEVEL_TAG "src")
SET(CPACK_PACKAGE_FILE_NAME "orocos-bfl")

INCLUDE(CPack)
ENDIF (CPACK_PACKAGES)
ELSE ( "${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" STREQUAL "2.4" )

   MESSAGE ( "Disabling packaging for version ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}")

ENDIF ( "${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" STREQUAL "2.4" )



#############################
#			    #
# Building Doxygen documents#
#                           #
#############################

ADD_CUSTOM_TARGET(docs "doxygen" "Doxyfile")

## Debugging
OPTION( DISPLAY_FLAGS "Display cmake variables at end of configure step." OFF)
MARK_AS_ADVANCED(DISPLAY_FLAGS)
IF( DISPLAY_FLAGS )
INCLUDE(config/cmake_logging.cmake)
ENDIF( DISPLAY_FLAGS )
