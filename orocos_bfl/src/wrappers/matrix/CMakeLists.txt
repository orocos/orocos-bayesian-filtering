GLOBAL_ADD_INCLUDE_DIR( ${CMAKE_SOURCE_DIR}/src/wrappers/matrix )
GLOBAL_ADD_INCLUDE_DIR( ${CMAKE_BINARY_DIR}/src/wrappers/matrix )
GLOBAL_ADD_INCLUDE_DIR( ${MATRIX_INCLUDE} )
GLOBAL_ADD_INCLUDE( bfl/wrappers/ ${CMAKE_BINARY_DIR}/src/wrappers/config.h)
GLOBAL_ADD_INCLUDE( bfl/wrappers/matrix matrix_wrapper.h vector_wrapper.h matrix_BOOST.h vector_BOOST.h
                                        matrix_NEWMAT.h vector_NEWMAT.h matrix_LTI.h vector_LTI.h
                                        matrix_EIGEN.h vector_EIGEN.h)
GLOBAL_ADD_SRC ( wrappers/matrix/matrix_BOOST.cpp  wrappers/matrix/vector_BOOST.cpp
                 wrappers/matrix/matrix_NEWMAT.cpp wrappers/matrix/vector_NEWMAT.cpp  
                 wrappers/matrix/matrix_LTI.cpp wrappers/matrix/vector_LTI.cpp  
                 wrappers/matrix/matrix_EIGEN.cpp wrappers/matrix/vector_EIGEN.cpp
                 wrappers/matrix/matrix_wrapper.cpp )

