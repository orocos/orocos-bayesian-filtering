# ------------------------- Begin Generic CMake Variable Logging ------------------


# if you are building in-source, this is the same as CMAKE_SOURCE_DIR, otherwise 
# this is the top level directory of your build tree 
MESSAGE( "CMAKE_BINARY_DIR:         " ${CMAKE_BINARY_DIR} )

# this is the directory, from which cmake was started, i.e. the top level source directory 
MESSAGE( "CMAKE_SOURCE_DIR:         " ${CMAKE_SOURCE_DIR} )

# contains the full path to the top level directory of your build tree 
MESSAGE( "PROJECT_BINARY_DIR: " ${PROJECT_BINARY_DIR} )

# contains the full path to the root of your project source directory,
# i.e. to the nearest directory where CMakeLists.txt contains the PROJECT() command 
MESSAGE( "PROJECT_SOURCE_DIR: " ${PROJECT_SOURCE_DIR} )

# set this variable to specify a common place where CMake should put all executable files
# (instead of CMAKE_CURRENT_BINARY_DIR)
MESSAGE( "EXECUTABLE_OUTPUT_PATH: " ${EXECUTABLE_OUTPUT_PATH} )

# set this variable to specify a common place where CMake should put all libraries 
# (instead of CMAKE_CURRENT_BINARY_DIR)
MESSAGE( "LIBRARY_OUTPUT_PATH:     " ${LIBRARY_OUTPUT_PATH} )

# tell CMake to search first in directories listed in CMAKE_MODULE_PATH
# when you use FIND_PACKAGE() or INCLUDE()
MESSAGE( "CMAKE_MODULE_PATH: " ${CMAKE_MODULE_PATH} )

# this is used when searching for include files e.g. using the FIND_PATH() command.
MESSAGE( "CMAKE_INCLUDE_PATH: " ${CMAKE_INCLUDE_PATH} )

# this is used when searching for libraries e.g. using the FIND_LIBRARY() command.
MESSAGE( "CMAKE_LIBRARY_PATH: " ${CMAKE_LIBRARY_PATH} )

# the complete system name, e.g. "Linux-2.4.22", "FreeBSD-5.4-RELEASE" or "Windows 5.1" 
MESSAGE( "CMAKE_SYSTEM: " ${CMAKE_SYSTEM} )

# the short system name, e.g. "Linux", "FreeBSD" or "Windows"
MESSAGE( "CMAKE_SYSTEM_NAME: " ${CMAKE_SYSTEM_NAME} )

# only the version part of CMAKE_SYSTEM 
MESSAGE( "CMAKE_SYSTEM_VERSION: " ${CMAKE_SYSTEM_VERSION} )

# the processor name (e.g. "Intel(R) Pentium(R) M processor 2.00GHz") 
MESSAGE( "CMAKE_SYSTEM_PROCESSOR: " ${CMAKE_SYSTEM_PROCESSOR} )

# is TRUE on all UNIX-like OS's, including Apple OS X and CygWin
MESSAGE( "UNIX: " ${UNIX} )

# is TRUE on Windows, including CygWin 
MESSAGE( "WIN32: " ${WIN32} )

# is TRUE on Apple OS X
MESSAGE( "APPLE: " ${APPLE} )

# is TRUE when using the MinGW compiler in Windows
MESSAGE( "MINGW: " ${MINGW} )

# is TRUE on Windows when using the CygWin version of cmake
MESSAGE( "CYGWIN: " ${CYGWIN} )

# is TRUE on Windows when using a Borland compiler 
MESSAGE( "BORLAND: " ${BORLAND} )

# Microsoft compiler 
MESSAGE( "MSVC: " ${MSVC} )
MESSAGE( "MSVC_IDE: " ${MSVC_IDE} )
MESSAGE( "MSVC60: " ${MSVC60} )
MESSAGE( "MSVC70: " ${MSVC70} )
MESSAGE( "MSVC71: " ${MSVC71} )
MESSAGE( "MSVC80: " ${MSVC80} )
MESSAGE( "CMAKE_COMPILER_2005: " ${CMAKE_COMPILER_2005} )


# set this to true if you don't want to rebuild the object files if the rules have changed, 
# but not the actual source files or headers (e.g. if you changed the some compiler switches) 
MESSAGE( "CMAKE_SKIP_RULE_DEPENDENCY: " ${CMAKE_SKIP_RULE_DEPENDENCY} )

MESSAGE( "CMAKE_C_COMPILER: "${CMAKE_C_COMPILER} )

MESSAGE( "CMAKE_CXX_COMPILER: "${CMAKE_CXX_COMPILER} )

MESSAGE( "CMAKE_COMPILER_IS_GNUCC: "${CMAKE_COMPILER_IS_GNUCC} )

MESSAGE( "CMAKE_COMPILER_IS_GNUCXX: "${CMAKE_COMPILER_IS_GNUCXX} )

# since CMake 2.1 the install rule depends on all, i.e. everything will be built