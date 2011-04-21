GLOBAL_ADD_INCLUDE_DIR( ${CMAKE_SOURCE_DIR}/src/dataAssociation )
GLOBAL_ADD_INCLUDE( bfl/dataAssociation dataAssociationFilter.h dataAssociationFilterMCPdf.h dataAssociationFilterGaussian.h)
GLOBAL_ADD_SRC ( dataAssociation/dataAssociationFilterGaussian.cpp )
