CMAKE_MINIMUM_REQUIRED(VERSION 2.2)

# project version
# ---------------
PROJECT(orocos-bfl)
SET(VERSION_MAJOR "0")
SET(VERSION_MINOR "6")
SET(VERSION_PATCH "0")
SET(VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")
MARK_AS_ADVANCED(VERSION)
MESSAGE( "Orocos BFL version ${VERSION}" )
MARK_AS_ADVANCED(LIBRARY_OUTPUT_PATH)
MARK_AS_ADVANCED(EXECUTABLE_OUTPUT_PATH)
MARK_AS_ADVANCED(CMAKE_BACKWARDS_COMPATIBILITY)


# Make compilation not verbose by default
# ---------------------------------------
SET(CMAKE_VERBOSE_MAKEFILE FALSE)


# project dirs
# ------------
IF(NOT CMAKE_INSTALL_PREFIX)
  SET( CMAKE_INSTALL_PREFIX /usr/local/ CACHE PATH "Installation directory" FORCE)
  MESSAGE( "Setting installation directory to ${CMAKE_INSTALL_PREFIX}" )
ENDIF(NOT CMAKE_INSTALL_PREFIX)
SET( PROJ_SOURCE_DIR ${orocos-bfl_SOURCE_DIR} )
SET( PROJ_BINARY_DIR ${orocos-bfl_BINARY_DIR} )


# win 32 options
# --------------
IF (WIN32)
  ADD_DEFINITIONS( /D_SCL_SECURE_NO_DEPRECATE)
ENDIF (WIN32)


# Set the build type (affects debugging symbols and optimization)
# ---------------------------------------------------------------
IF ( NOT CMAKE_BUILD_TYPE )
  SET( CMAKE_BUILD_TYPE "Debug" CACHE STRING "Build type: None, Release, Debug, RelWithDebInfo, MinSizeRel." FORCE )
  MARK_AS_ADVANCED(FORCE CMAKE_C_FLAGS_BFL CMAKE_CXX_FLAGS_BFL)
ENDIF ( NOT CMAKE_BUILD_TYPE )
MESSAGE( "Build type set to '${CMAKE_BUILD_TYPE}'" )


# MESSAGE("Looking for compiler:")
# --------------------------------
INCLUDE (${CMAKE_ROOT}/Modules/CMakeDetermineCXXCompiler.cmake)
INCLUDE( ${PROJ_SOURCE_DIR}/config/DependentOption.cmake )


# Shared or static library
# ------------------------
IF (NOT LIBRARY_TYPE)
  SET( LIBRARY_TYPE shared CACHE STRING "What type of library to build: (shared, static or both)")
ENDIF (NOT LIBRARY_TYPE)


# find libs 
# ----------
INCLUDE (${PROJ_SOURCE_DIR}/config/check_depend.cmake)
INCLUDE (${PROJ_SOURCE_DIR}/config/bfl_macros.cmake)
INCLUDE (${PROJ_SOURCE_DIR}/config/FindMATRIX.cmake)
INCLUDE (${PROJ_SOURCE_DIR}/config/FindRNG.cmake)
INCLUDE (${PROJ_SOURCE_DIR}/config/FindGINAC.cmake)
INCLUDE (${PROJ_SOURCE_DIR}/config/FindCppUnit.cmake)


# options to build tests and examples
# -----------------------------------
DEPENDENT_OPTION( BUILD_TESTS "Turn me off to disable the BFL functionality testing" ON "CPPUNIT_FOUND" OFF)
OPTION( BUILD_EXAMPLES "Turn me off to disable the compilation of the BFL examples" ON )


# Uninstall target
# ----------------
CONFIGURE_FILE(
  "${CMAKE_CURRENT_SOURCE_DIR}/config/cmake_uninstall.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
  IMMEDIATE @ONLY)
ADD_CUSTOM_TARGET(uninstall
  "${CMAKE_COMMAND}" -P "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake")


# add subfolders
# --------------
SET(ENV{SELECTED_LIBS} "")
SET(ENV{SELECTED_DIRS} "")
ADD_SUBDIRECTORY( src )
ADD_SUBDIRECTORY( tests )
ADD_SUBDIRECTORY( examples )


# .pc file
# --------
CONFIGURE_FILE( orocos-bfl.pc.in orocos-bfl.pc @ONLY)
INSTALL_FILES( /lib/pkgconfig FILES orocos-bfl.pc)


# Building Doxygen documents
# --------------------------
CONFIGURE_FILE(revision.in revision @ONLY)
EXECUTE_PROCESS(COMMAND ./revision OUTPUT_VARIABLE SVN_REVISION)
MESSAGE("SVN REVISION ${SVN_REVISION}")
CONFIGURE_FILE(Doxyfile.in Doxyfile @ONLY)
ADD_CUSTOM_TARGET(docs "doxygen" "Doxyfile")
