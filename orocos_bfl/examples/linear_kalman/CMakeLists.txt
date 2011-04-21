IF (BUILD_EXAMPLES)
  ADD_EXECUTABLE( test_linear_kalman  test_linear_kalman.cpp ../mobile_robot.cpp
../nonlinearanalyticconditionalgaussianmobile.cpp) 

  TARGET_LINK_LIBRARIES( test_linear_kalman orocos-bfl ${MATRIX_LIBS} ${RNG_LIBS}  ${GINAC_LIBS})

  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/src/wrappers/matrix ${CMAKE_SOURCE_DIR}/src ${RNG_INCLUDE} ${MATRIX_INCLUDE})

  INSTALL(TARGETS test_linear_kalman
          DESTINATION ${CMAKE_INSTALL_PREFIX}/bin/bfl 
          PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
ENDIF (BUILD_EXAMPLES)

