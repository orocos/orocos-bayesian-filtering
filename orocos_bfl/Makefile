default: install

SVN_REVISION= -r HEAD
SVN_DIR = bfl-svn
SVN_URL = http://svn.mech.kuleuven.be/repos/orocos/branches/bfl/branch-0.8
INSTALL_DIR = bfl-boost
CMAKE = cmake
OROCOS_INSTALL=`rospack find rtt_1.10`/rtt
#BOOST_INCLUDE =$(shell rosboost-cfg --include_dirs)
EXTRA_CMAKE_FLAGS = -DBUILD_TESTS=OFF \
						 -DCMAKE_INSTALL_PREFIX=`rospack find bfl_0.7`/$(INSTALL_DIR)/ \
						 -DCMAKE_INCLUDE_PATH=$(BOOST_INCLUDE) \
						 -DOROCOS_PLUGIN=OFF \
          				 -DOROCOS_INSTALL=$(OROCOS_INSTALL) \
                       	 -DCMAKE_BUILD_TYPE=Release \
						 -DMATRIX_LIB=boost -DRNG_LIB=boost

include $(shell rospack find mk)/svn_checkout.mk
include $(shell rospack find mk)/cmake.mk

install: $(SVN_DIR) patched
	mkdir -p $(SVN_DIR)/build-boost
	cd $(SVN_DIR)/build-boost && $(CMAKE) $(EXTRA_CMAKE_FLAGS) ..
	cd $(SVN_DIR)/build-boost && make $(ROS_PARALLEL_JOBS) && make install
	touch installed

wiped: Makefile
	make wipe
	touch wiped

clean:
	rm -rf $(INSTALL_DIR)
	rm -rf $(SVN_DIR)/build-boost
	rm -rf installed

wipe: 	clean
	rm -rf $(SVN_DIR)
