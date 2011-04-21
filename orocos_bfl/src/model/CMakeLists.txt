GLOBAL_ADD_INCLUDE_DIR( ${CMAKE_SOURCE_DIR}/src/model )
GLOBAL_ADD_INCLUDE( bfl/model systemmodel.h systemmodel.cpp measurementmodel.h measurementmodel.cpp analyticsystemmodel_gaussianuncertainty.h 
                              analyticmeasurementmodel_gaussianuncertainty.h linearanalyticmeasurementmodel_gaussianuncertainty.h 
                              linearanalyticsystemmodel_gaussianuncertainty.h discretesystemmodel.h 
                              linearanalyticmeasurementmodel_gaussianuncertainty_implicit.h)
GLOBAL_ADD_SRC (model/analyticmeasurementmodel_gaussianuncertainty.cpp model/analyticsystemmodel_gaussianuncertainty.cpp
                model/linearanalyticsystemmodel_gaussianuncertainty.cpp model/linearanalyticmeasurementmodel_gaussianuncertainty.cpp
                model/discretesystemmodel.cpp model/linearanalyticmeasurementmodel_gaussianuncertainty_implicit.cpp)

IF (GINAC_SUPPORT)
 GLOBAL_ADD_INCLUDE( bfl/model nonlinearanalyticmeasurementmodel_gaussianuncertainty_ginac.h
                               nonlinearanalyticsystemmodel_gaussianuncertainty_ginac.h)
 GLOBAL_ADD_SRC( model/nonlinearanalyticmeasurementmodel_gaussianuncertainty_ginac.cpp
                 model/nonlinearanalyticsystemmodel_gaussianuncertainty_ginac.cpp)
ENDIF (GINAC_SUPPORT)
