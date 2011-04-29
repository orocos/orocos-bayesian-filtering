INSTALL_PREFIX=${PWD}/install
EXTRA_CMAKE_FLAGS=-DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX} \
                         -DBUILD_TESTS=OFF \
						 -DCMAKE_INCLUDE_PATH=$(BOOST_INCLUDE) \
						 -DOROCOS_PLUGIN=OFF \
                       	 -DCMAKE_BUILD_TYPE=Release \
						 -DMATRIX_LIB=boost -DRNG_LIB=boost
                        

default: all
	cd build; make install


include $(shell rospack find mk)/cmake.mk

