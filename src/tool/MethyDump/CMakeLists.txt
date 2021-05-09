cmake_minimum_required(VERSION 3.19)
project(MethyDump)

set(CMAKE_CXX_STANDARD 14)

link_directories("${CMAKE_SOURCE_DIR}/lib")
include_directories("${CMAKE_SOURCE_DIR}/include")

SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)

add_executable(MethyDump main.cpp my_config.h)