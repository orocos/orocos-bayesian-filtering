IF (BUILD_TESTS)

  ENABLE_TESTING()

  INCLUDE_DIRECTORIES( ${PROJ_SOURCE_DIR}/src ${PROJ_BINARY_DIR}/src/wrappers/matrix ${MATRIX_INCLUDE} ${RNG_INCLUDE} ${CPPUNIT_INCLUDE})
  LINK_DIRECTORIES( ${PROJ_BINARY_DIR}/src )
  ADD_DEFINITIONS("-g")

  ADD_EXECUTABLE( test_ekf ${DO_EXCLUDE} test-runner.cpp approxEqual.cpp ekf_test.cpp  ekf_test.hpp )
  TARGET_LINK_LIBRARIES( test_ekf orocos-bfl ${MATRIX_LIBS} ${RNG_LIBS} ${GINAC_LIBS} ${CPPUNIT_LIBS} )
  ADD_TEST( test_ekf test_ekf )

  ADD_EXECUTABLE( test_pdf ${DO_EXCLUDE} test-runner.cpp approxEqual.cpp pdf_test.cpp  pdf_test.hpp )
  TARGET_LINK_LIBRARIES( test_pdf orocos-bfl ${MATRIX_LIBS} ${RNG_LIBS} ${GINAC_LIBS} ${CPPUNIT_LIBS} )
  ADD_TEST( test_pdf test_pdf )

  ADD_EXECUTABLE( test_sample ${DO_EXCLUDE} test-runner.cpp  approxEqual.cpp sample_test.cpp sample_test.hpp)
  TARGET_LINK_LIBRARIES( test_sample orocos-bfl ${MATRIX_LIBS} ${RNG_LIBS} ${GINAC_LIBS} ${CPPUNIT_LIBS} )
  ADD_TEST( test_sample test_sample )

  ADD_EXECUTABLE( test_matrixwrapper ${DO_EXCLUDE} test-runner.cpp  approxEqual.cpp matrixwrapper_test.cpp matrixwrapper_test.hpp)
  TARGET_LINK_LIBRARIES( test_matrixwrapper orocos-bfl ${MATRIX_LIBS} ${RNG_LIBS} ${GINAC_LIBS} ${CPPUNIT_LIBS} )
  ADD_TEST( test_matrixwrapper test_matrixwrapper )

  ADD_EXECUTABLE( test_model ${DO_EXCLUDE} test-runner.cpp approxEqual.cpp model_test.cpp model_test.hpp)
  TARGET_LINK_LIBRARIES( test_model orocos-bfl ${MATRIX_LIBS} ${RNG_LIBS} ${GINAC_LIBS} ${CPPUNIT_LIBS} )
  ADD_TEST( test_model test_model )

  ADD_EXECUTABLE( test_smoother ${DO_EXCLUDE} test-runner.cpp approxEqual.cpp smoother_test.cpp smoother_test.hpp
                  ${PROJ_SOURCE_DIR}/examples/mobile_robot.cpp ${PROJ_SOURCE_DIR}/examples/compare_filters/nonlinearanalyticconditionalgaussianmobile.cpp)
  TARGET_LINK_LIBRARIES( test_smoother orocos-bfl ${MATRIX_LIBS} ${RNG_LIBS} ${GINAC_LIBS} ${CPPUNIT_LIBS} )
  ADD_TEST( test_smoother test_smoother )

  ADD_EXECUTABLE( test_complete_filter ${DO_EXCLUDE} test-runner.cpp approxEqual.cpp complete_filter_test.cpp complete_filter_test.hpp
                  ${PROJ_SOURCE_DIR}/examples/mobile_robot.cpp ${PROJ_SOURCE_DIR}/examples/compare_filters/nonlinearanalyticconditionalgaussianmobile.cpp
                  ${PROJ_SOURCE_DIR}/examples/discrete_filter/conditionalUniformMeasPdf1d.cpp )
  TARGET_LINK_LIBRARIES( test_complete_filter orocos-bfl ${MATRIX_LIBS} ${RNG_LIBS} ${GINAC_LIBS} ${CPPUNIT_LIBS} )
  ADD_TEST( test_complete_filter test_complete_filter )

  ADD_EXECUTABLE( test_dataAssociationFilter ${DO_EXCLUDE} test-runner.cpp approxEqual.cpp dataAssociationFilterTest.cpp dataAssociationFilterTest.hpp
                  ${PROJ_SOURCE_DIR}/examples/mobile_robot.cpp ${PROJ_SOURCE_DIR}/examples/compare_filters/nonlinearanalyticconditionalgaussianmobile.cpp
                  ${PROJ_SOURCE_DIR}/examples/discrete_filter/conditionalUniformMeasPdf1d.cpp )
  TARGET_LINK_LIBRARIES( test_dataAssociationFilter orocos-bfl ${MATRIX_LIBS}  ${RNG_LIBS} ${GINAC_LIBS} ${CPPUNIT_LIBS})
  ADD_TEST( test_dataAssociationFilter test_dataAssociationFilter )

IF (GINAC_SUPPORT)
  ADD_EXECUTABLE( test_model_ginac ${DO_EXCLUDE} test-runner.cpp approxEqual.cpp model_test_ginac.cpp model_test_ginac.hpp)
  TARGET_LINK_LIBRARIES( test_model_ginac orocos-bfl ${MATRIX_LIBS} ${RNG_LIBS} ${GINAC_LIBS} ${CPPUNIT_LIBS} )
  ADD_TEST( test_model_ginac test_model_ginac )
ENDIF (GINAC_SUPPORT)




  ADD_CUSTOM_TARGET(check ctest -V 
                    DEPENDS test_ekf test_pdf test_sample test_matrixwrapper test_model test_smoother test_complete_filter
                    WORKING_DIRECTORY ${PROJ_BINARY_DIR}/tests)

IF (GINAC_SUPPORT)
  ADD_CUSTOM_TARGET(check ctest -V 
                    DEPENDS test_model_ginac
                    WORKING_DIRECTORY ${PROJ_BINARY_DIR}/tests)
ENDIF (GINAC_SUPPORT)



ENDIF (BUILD_TESTS)

