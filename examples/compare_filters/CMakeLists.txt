IF (BUILD_EXAMPLES)
  ADD_EXECUTABLE( test_compare_filters test_compare_filters.cpp nonlinearanalyticconditionalgaussianmobile.cpp ../mobile_robot.cpp )

  TARGET_LINK_LIBRARIES( test_compare_filters orocos-bfl ${MATRIX_LIBS} ${RNG_LIBS}  ${GINAC_LIBS})

  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/src/wrappers/matrix ${CMAKE_SOURCE_DIR}/src ${RNG_INCLUDE} ${MATRIX_INCLUDE})
ENDIF (BUILD_EXAMPLES)

