cmake_minimum_required(VERSION 3.16)
project(methyprofile)

set(CMAKE_CXX_STANDARD 14)

link_directories("./lib")
include_directories("./include")

find_package(Threads)

add_executable(methyprofile main.cpp bed.cpp bed.h config.h)

target_link_libraries(methyprofile ${CMAKE_THREAD_LIBS_INIT})