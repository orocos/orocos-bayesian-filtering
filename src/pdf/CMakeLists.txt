GLOBAL_ADD_INCLUDE_DIR( ${CMAKE_SOURCE_DIR}/src/pdf )
GLOBAL_ADD_INCLUDE( bfl/pdf pdf.h conditionalpdf.h discretepdf.h discreteconditionalpdf.h mcpdf.h mcpdf.cpp 
                            gaussian.h conditionalgaussian.h conditionalgaussian_additivenoise.h analyticconditionalgaussian.h 
                            analyticconditionalgaussian_additivenoise.h linearanalyticconditionalgaussian.h filterproposaldensity.h 
                            EKF_proposaldensity.h)
GLOBAL_ADD_SRC ( pdf/analyticconditionalgaussian.cpp pdf/conditionalgaussian.cpp pdf/gaussian.cpp 
                 pdf/conditionalgaussian_additivenoise.cpp pdf/analyticconditionalgaussian_additivenoise.cpp 
                 pdf/linearanalyticconditionalgaussian.cpp pdf/filterproposaldensity.cpp pdf/EKF_proposaldensity.cpp 
                 pdf/discretepdf.cpp pdf/discreteconditionalpdf.cpp)

IF (GINAC_SUPPORT)
 GLOBAL_ADD_INCLUDE( bfl/pdf  nonlinearanalyticconditionalgaussian_ginac.h)
 GLOBAL_ADD_SRC( pdf/nonlinearanalyticconditionalgaussian_ginac.cpp)
ENDIF (GINAC_SUPPORT)
