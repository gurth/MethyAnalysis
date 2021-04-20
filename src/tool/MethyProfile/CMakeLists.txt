cmake_minimum_required(VERSION 3.16)
project(methyprofile)

set(CMAKE_CXX_STANDARD 14)

link_directories("${CMAKE_SOURCE_DIR}/lib")
include_directories("${CMAKE_SOURCE_DIR}/include")

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake;${CMAKE_MODULE_PATH}")

IF (WIN32)
    MESSAGE(STATUS "Now is windows")
    ADD_DEFINITIONS(-D _WIN32_PLATFORM_)
ELSEIF (APPLE OR UNIX)
    MESSAGE(STATUS "Now is UNIX-like OS's. Including APPLE OSX and CygWin")
    ADD_DEFINITIONS(-D _UNIX_PLATFORM_)
    find_package(Threads)
ENDIF ()

ADD_DEFINITIONS(-D CMAKE_SOURCE_DIR="${CMAKE_SOURCE_DIR}")

SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)

add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/plug-in/libsave.so")

add_executable(methyprofile main.cpp bed.cpp bed.h config.h merror.h profile_node.h)

target_link_libraries(methyprofile dl ${CMAKE_THREAD_LIBS_INIT} libzlog.so)