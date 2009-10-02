IF( OROCOS_PLUGIN )


  #CONFIGURE_FILE( bfl_toolkit.hpp.in src/bindings/rtt/bfl_toolkit.hpp)
  INSTALL_FILES( /include/bfl/bindings/rtt FILES bfl_toolkit.hpp)
  INCLUDE_DIRECTORIES(${PROJ_SOURCE_DIR}/src)

  FILE( GLOB SRCS [^.]*.cpp [^.]*.cxx)
  FILE( GLOB HPPS [^.]*.hpp [^.]*.h)
  
  INCLUDE_DIRECTORIES(${PROJ_SOURCE_DIR} ${PROJ_BINARY_DIR}/src/bindings/rtt
      ${OROCOS_RTT_INCLUDE_DIRS} ${OROCOS_INSTALL}/include/rtt )
  ADD_DEFINITIONS(-D_REENTRANT)
  INCLUDE_DIRECTORIES(/usr/include/rtt)
  ADD_LIBRARY(orocos-bfl_toolkit-${OROCOS_TARGET} SHARED ${SRCS})
  SET_TARGET_PROPERTIES(orocos-bfl_toolkit-${OROCOS_TARGET} PROPERTIES
  SOVERSION "${bfl_VERSION_MAJOR}.${bfl_VERSION_MINOR}" 
  VERSION "${bfl_VERSION}") 
  INSTALL_TARGETS( /lib/rtt/${OROCOS_TARGET}/plugins orocos-bfl_toolkit-${OROCOS_TARGET})
  INSTALL_FILES( /include/bfl/bindings/rtt FILES ${HPPS})
  TARGET_LINK_LIBRARIES(orocos-bfl_toolkit-${OROCOS_TARGET} orocos-bfl orocos-rtt-${OROCOS_TARGET})

  CONFIGURE_FILE( bfl_toolkit.pc.in src/bindings/rtt/orocos-bfl_toolkit-${OROCOS_TARGET}.pc @ONLY)
  INSTALL_FILES( /lib/pkgconfig FILES orocos-bfl_toolkit-${OROCOS_TARGET}.pc)

  INCLUDE_DIRECTORIES( $ENV{GLOBAL_INCLUDE_DIRS} )

ENDIF( OROCOS_PLUGIN )
IF( OROCOS_PLUGIN )


  #CONFIGURE_FILE( bfl_toolkit.hpp.in src/bindings/rtt/bfl_toolkit.hpp)
  INSTALL_FILES( /include/bfl/bindings/rtt FILES bfl_toolkit.hpp)
  INCLUDE_DIRECTORIES(${PROJ_SOURCE_DIR}/src)

  FILE( GLOB SRCS [^.]*.cpp [^.]*.cxx)
  FILE( GLOB HPPS [^.]*.hpp [^.]*.h)
  
  INCLUDE_DIRECTORIES(${PROJ_SOURCE_DIR} ${PROJ_BINARY_DIR}/src/bindings/rtt
      ${OROCOS_RTT_INCLUDE_DIRS} ${OROCOS_INSTALL}/include/rtt )
  ADD_DEFINITIONS(-D_REENTRANT)
  INCLUDE_DIRECTORIES(/usr/include/rtt)
  ADD_LIBRARY(orocos-bfl_toolkit-${OROCOS_TARGET} SHARED ${SRCS})
  SET_TARGET_PROPERTIES(orocos-bfl_toolkit-${OROCOS_TARGET} PROPERTIES
  SOVERSION "${bfl_VERSION_MAJOR}.${bfl_VERSION_MINOR}" 
  VERSION "${bfl_VERSION}") 
  INSTALL_TARGETS( /lib/rtt/${OROCOS_TARGET}/plugins orocos-bfl_toolkit-${OROCOS_TARGET})
  INSTALL_FILES( /include/bfl/bindings/rtt FILES ${HPPS})
  TARGET_LINK_LIBRARIES(orocos-bfl_toolkit-${OROCOS_TARGET} orocos-bfl orocos-rtt-${OROCOS_TARGET})

  CONFIGURE_FILE( bfl_toolkit.pc.in src/bindings/rtt/orocos-bfl_toolkit-${OROCOS_TARGET}.pc @ONLY)
  INSTALL_FILES( /lib/pkgconfig FILES orocos-bfl_toolkit-${OROCOS_TARGET}.pc)

  INCLUDE_DIRECTORIES( $ENV{GLOBAL_INCLUDE_DIRS} )

ENDIF( OROCOS_PLUGIN )
