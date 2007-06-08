CMAKE_MINIMUM_REQUIRED(VERSION 2.2)
#MARK_AS_ADVANCED( FORCE CMAKE_BACKWARDS_COMPATIBILITY )



###################################################
#                                                 #
#       Start project customization section       #
#                                                 #
###################################################

PROJECT(orocos-bfl)
SET(VERSION_MAJOR "0")
SET(VERSION_MINOR "4")
SET(VERSION_PATCH "9")
SET(VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")
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


# get svn revision number
# ------------
CONFIGURE_FILE(revision.in revision @ONLY)
EXECUTE_PROCESS(COMMAND ./revision OUTPUT_VARIABLE SVN_REVISION)
MESSAGE("SVN REVISION ${SVN_REVISION}")

IF (WIN32)
	ADD_DEFINITIONS( /D_SCL_SECURE_NO_DEPRECATE)
ENDIF (WIN32)

# Set the build type (affects debugging symbols and optimization)
# ------------------
IF ( NOT CMAKE_BUILD_TYPE )
  SET( CMAKE_BUILD_TYPE "Debug" CACHE STRING "Build type: None (Use CMAKE_C_FLAGS and CMAKE_CXX_FLAGS), Release, Debug, RelWithDebInfo, MinSizeRel." FORCE )
  MESSAGE( "Build type set to '${CMAKE_BUILD_TYPE}'" )
ELSE ( NOT CMAKE_BUILD_TYPE )
  MESSAGE( "Build type set to '${CMAKE_BUILD_TYPE}'" )
  IF ( NOT CMAKE_BUILD_TYPE STREQUAL "BFL" )
    MARK_AS_ADVANCED(FORCE CMAKE_C_FLAGS_BFL CMAKE_CXX_FLAGS_BFL)
  ELSE ( NOT CMAKE_BUILD_TYPE STREQUAL "BFL" )
    MARK_AS_ADVANCED(CLEAR CMAKE_C_FLAGS_BFL CMAKE_CXX_FLAGS_BFL)
  ENDIF ( NOT CMAKE_BUILD_TYPE STREQUAL "BFL" )
ENDIF ( NOT CMAKE_BUILD_TYPE )

# MESSAGE("Looking for compiler:")
INCLUDE (${CMAKE_ROOT}/Modules/CMakeDetermineCXXCompiler.cmake)
INCLUDE( ${CMAKE_ROOT}/Modules/CheckIncludeFileCXX.cmake )
INCLUDE( ${CMAKE_ROOT}/Modules/CheckIncludeFile.cmake )
INCLUDE( ${PROJ_SOURCE_DIR}/config/DependentOption.cmake )


###################################################
#                                                 #
#       Do some global set-up & checking          #
#                                                 #
###################################################

# Put your 'flags' inhere if you want them to apear in the
# orocos-BFL.pc file. Various .cmake scripts below do this.
SET( ENV{OROCOS_BFL_CFLAGS} "")
SET( ENV{OROCOS_BFL_LINKFLAGS} "")

## Shared or static library
IF (NOT LIBRARY_TYPE)
  SET( LIBRARY_TYPE shared CACHE STRING "What type of library to build: (shared, static or both)")
ENDIF (NOT LIBRARY_TYPE)


## Find MATRIX LIB
# ----------------
IF (NOT __MATRIXWRAPPER_NEWMAT__)
  SET(__MATRIXWRAPPER_NEWMAT__ OFF CACHE BOOL "define for newmat")
  MARK_AS_ADVANCED(__MATRIXWRAPPER_NEWMAT__)
ENDIF (NOT __MATRIXWRAPPER_NEWMAT__)
IF (NOT __MATRIXWRAPPER_LTI__)
  SET(__MATRIXWRAPPER_LTI__ OFF CACHE BOOL "define for lti")
  MARK_AS_ADVANCED(__MATRIXWRAPPER_LTI__)
ENDIF (NOT __MATRIXWRAPPER_LTI__)
IF (NOT __MATRIXWRAPPER_BOOST__)
  SET(__MATRIXWRAPPER_BOOST__ OFF CACHE BOOL "define for boost")
  MARK_AS_ADVANCED(__MATRIXWRAPPER_BOOST__)
ENDIF (NOT __MATRIXWRAPPER_BOOST__)
IF(NOT MATRIX_LIB)
  SET( MATRIX_LIB lti CACHE STRING "Which matrix library to use: boost, lti or newmat")
ENDIF(NOT MATRIX_LIB)
IF(NOT MATRIX_INSTALL)
  SET( MATRIX_INSTALL /usr CACHE PATH "The Matrix lib installation directory.")
ENDIF(NOT MATRIX_INSTALL)
INCLUDE (${PROJ_SOURCE_DIR}/config/FindMATRIX.cmake)

## Find RNG LIB
# ----------------
IF (NOT __RNGWRAPPER_BOOST__)
  SET(__RNGWRAPPER_BOOST__ OFF CACHE BOOL "define for boost")
  MARK_AS_ADVANCED(__RNGWRAPPER_BOOST__)
ENDIF (NOT __RNGWRAPPER_BOOST__)
IF (NOT __RNGWRAPPER_LTI__)
  SET(__RNGWRAPPER_LTI__ OFF CACHE BOOL "define for lti")
  MARK_AS_ADVANCED(__RNGWRAPPER_LTI__)
ENDIF (NOT __RNGWRAPPER_LTI__)
IF(NOT RNG_LIB)
  SET( RNG_LIB lti CACHE STRING "Which rng library to use: lti or boost")
ENDIF(NOT RNG_LIB)
IF(NOT RNG_INSTALL)
  SET( RNG_INSTALL /usr CACHE PATH "The rng lib installation directory.")
ENDIF(NOT RNG_INSTALL)
INCLUDE (${PROJ_SOURCE_DIR}/config/FindRNG.cmake)

## Find GINAC
# ----------------
OPTION( GINAC_SUPPORT "Turn me off to disable Ginac support" OFF )
IF(NOT GINAC_INSTALL)
  SET( GINAC_INSTALL /usr CACHE PATH "The Ginac lib installation directory.")
ENDIF(NOT GINAC_INSTALL)
INCLUDE (${PROJ_SOURCE_DIR}/config/FindGINAC.cmake)

## Find CPP Unit
# ----------------
IF(NOT CPPUNIT_INSTALL)
  SET( CPPUNIT_INSTALL /usr CACHE PATH "The CppUnit lib installation directory.")
ENDIF(NOT CPPUNIT_INSTALL)
INCLUDE (${PROJ_SOURCE_DIR}/config/FindCppUnit.cmake)
DEPENDENT_OPTION( BUILD_TESTS "Turn me off to disable the BFL functionality testing" ON "CPPUNIT_FOUND" OFF)
OPTION( BUILD_EXAMPLES "Turn me off to disable the compilation of the BFL examples" ON )


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




#############################
#			    #
# Building Doxygen documents#
#                           #
#############################
CONFIGURE_FILE(Doxyfile.in Doxyfile @ONLY)
ADD_CUSTOM_TARGET(docs "doxygen" "Doxyfile")
