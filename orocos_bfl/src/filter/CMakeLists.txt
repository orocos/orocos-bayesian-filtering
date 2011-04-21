GLOBAL_ADD_INCLUDE_DIR( ${CMAKE_SOURCE_DIR}/src/filter )
GLOBAL_ADD_INCLUDE( bfl/filter filter.h filter.cpp  particlefilter.h particlefilter.cpp bootstrapfilter.h bootstrapfilter.cpp
                        asirfilter.h asirfilter.cpp kalmanfilter.h extendedkalmanfilter.h iteratedextendedkalmanfilter.h
                        EKparticlefilter.h SRiteratedextendedkalmanfilter.h innovationCheck.h histogramfilter.h histogramfilter.cpp
                        mixtureParticleFilter.h mixtureParticleFilter.cpp mixtureBootstrapFilter.h mixtureBootstrapFilter.cpp)
GLOBAL_ADD_SRC (filter/iteratedextendedkalmanfilter.cpp filter/kalmanfilter.cpp filter/extendedkalmanfilter.cpp
                filter/EKparticlefilter.cpp filter/SRiteratedextendedkalmanfilter.cpp filter/innovationCheck.cpp )

