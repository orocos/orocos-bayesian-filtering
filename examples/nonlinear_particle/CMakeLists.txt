IF (BUILD_EXAMPLES)
  ADD_EXECUTABLE( test_nonlinear_particle test_nonlinear_particle.cpp nonlinearSystemPdf.cpp nonlinearMeasurementPdf.cpp ../mobile_robot.cpp )

  TARGET_LINK_LIBRARIES( test_nonlinear_particle orocos-bfl ${MATRIX_LIBS} ${RNG_LIBS}  ${GINAC_LIBS})

  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/src/wrappers/matrix ${CMAKE_SOURCE_DIR}/src ${RNG_INCLUDE} ${MATRIX_INCLUDE})
ENDIF (BUILD_EXAMPLES)

