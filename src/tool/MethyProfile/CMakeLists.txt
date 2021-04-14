cmake_minimum_required(VERSION 3.16)
project(methyprofile)

set(CMAKE_CXX_STANDARD 14)

link_directories("./lib")
include_directories("./include")

IF (WIN32)
    MESSAGE(STATUS "Now is windows")
    ADD_DEFINITIONS(-D _WIN32_PLATFORM_)
ELSEIF (APPLE OR UNIX)
    MESSAGE(STATUS "Now is UNIX-like OS's. Including APPLE OSX and CygWin")
    ADD_DEFINITIONS(-D _UNIX_PLATFORM_)
    find_package(Threads)
ENDIF ()

add_executable(methyprofile main.cpp bed.cpp bed.h config.h merror.h profile_node.h)

target_link_libraries(methyprofile dl ${CMAKE_THREAD_LIBS_INIT})