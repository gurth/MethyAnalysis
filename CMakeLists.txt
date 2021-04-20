cmake_minimum_required(VERSION 3.16)
project(methyanalysis)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake;${CMAKE_MODULE_PATH}")

include(libzlog)

add_subdirectory("${PROJECT_SOURCE_DIR}/src/tool/MethyProfile")