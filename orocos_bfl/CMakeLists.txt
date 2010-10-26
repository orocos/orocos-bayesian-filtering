cmake_minimum_required(VERSION 2.6.0)
include($ENV{ROS_ROOT}/core/rosbuild/rosbuild.cmake)

rosbuild_init()

#set(CMAKE_BUILD_TYPE "RTT")
add_subdirectory(bfl-svn)

