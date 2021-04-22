cmake_minimum_required(VERSION 3.16)

project(libsave.so)

add_custom_target(libsave.so ALL
        COMMAND ${CMAKE_MAKE_PROGRAM}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMENT "Building CXX library etc/libsave.so")