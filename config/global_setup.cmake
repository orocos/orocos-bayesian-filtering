#
# Include cmake modules required to look for dependencies
#
INCLUDE( ${CMAKE_ROOT}/Modules/CheckIncludeFileCXX.cmake )
INCLUDE( ${CMAKE_ROOT}/Modules/CheckIncludeFile.cmake )

INCLUDE( ${PROJ_SOURCE_DIR}/config/DependentOption.cmake )

# An option for tests, to make it easy to turn off all tests
DEPENDENT_OPTION( BUILD_TESTS "Turn me off to disable compilation of all tests" ON "CPPUNIT_FOUND" OFF)
OPTION( BUILD_EXAMPLES "Turn me off to disable compilation of all examples" ON )



