GLOBAL_ADD_INCLUDE_DIR( ${CMAKE_SOURCE_DIR}/src/smoother )
GLOBAL_ADD_INCLUDE( bfl/smoother backwardfilter.h backwardfilter.cpp rauchtungstriebel.h particlesmoother.h particlesmoother.cpp)
GLOBAL_ADD_SRC ( smoother/rauchtungstriebel.cpp)
