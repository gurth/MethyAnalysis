cmake_minimum_required(VERSION 3.16)
project(methyprofile)

set(CMAKE_CXX_STANDARD 14)

link_directories("${CMAKE_SOURCE_DIR}/lib")
include_directories("${CMAKE_SOURCE_DIR}/include")

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake;${CMAKE_MODULE_PATH}")

IF (WIN32 OR (CMAKE_SYSTEM_NAME MATCHES "CYGWIN"))
    ADD_DEFINITIONS(-D _WIN32_PLATFORM_)
else()
    ADD_DEFINITIONS(-D _UNIX_PLATFORM_)
    ADD_DEFINITIONS(-D ENABLE_LOG)
    find_package(Threads)
ENDIF ()

ADD_DEFINITIONS(-D CMAKE_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
ADD_DEFINITIONS(-D CMAKE_CURRENT_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

ADD_DEFINITIONS(-D MY_CONFIG_PATH="${CMAKE_CURRENT_SOURCE_DIR}/my_config.h")

SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)

add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/plug-in/libsave")

add_executable(methyprofile main.cpp bed.cpp bed.h profile_node.h my_config.h)

IF (WIN32 OR (CMAKE_SYSTEM_NAME MATCHES "CYGWIN"))
    target_link_libraries(methyprofile dl ${CMAKE_THREAD_LIBS_INIT})
ELSE()
    target_link_libraries(methyprofile dl ${CMAKE_THREAD_LIBS_INIT} libzlog.so)
ENDIF()