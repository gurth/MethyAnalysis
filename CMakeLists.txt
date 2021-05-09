cmake_minimum_required(VERSION 3.16)
project(methyanalysis)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake;${CMAKE_MODULE_PATH}")

IF (WIN32 OR (CMAKE_SYSTEM_NAME MATCHES "CYGWIN"))
ELSE()
    include(libzlog)
ENDIF()

add_subdirectory("${PROJECT_SOURCE_DIR}/src/tool/MethyProfile")

add_subdirectory("${PROJECT_SOURCE_DIR}/src/MethyProfile-GUI")

add_subdirectory("${PROJECT_SOURCE_DIR}/src/tool/MethyDump")