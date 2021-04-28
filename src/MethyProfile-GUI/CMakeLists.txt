cmake_minimum_required(VERSION 3.5)

project(MethyProfile-GUI LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# QtCreator supports the following variables for Android, which are identical to qmake Android variables.
# Check http://doc.qt.io/qt-5/deployment-android.html for more information.
# They need to be set before the find_package(Qt5 ...) call.

#if(ANDROID)
#    set(ANDROID_PACKAGE_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/android")
#    if (ANDROID_ABI STREQUAL "armeabi-v7a")
#        set(ANDROID_EXTRA_LIBS
#            ${CMAKE_CURRENT_SOURCE_DIR}/path/to/libcrypto.so
#            ${CMAKE_CURRENT_SOURCE_DIR}/path/to/libssl.so)
#    endif()
#endif()

IF (WIN32)
    # Your QT path
    set(CMAKE_PREFIX_PATH "C:\\Qt\\Qt\\5.15.2\\mingw81_64")

    set(Qt5_DIR "${CMAKE_PREFIX_PATH}/lib/cmake/Qt5")
    set(Qt5Widgets_DIR "${CMAKE_PREFIX_PATH}/lib/cmake/Qt5Widgets")

    ADD_DEFINITIONS(-D _WIN32_PLATFORM_)

ELSEIF (APPLE OR UNIX)
IF (CMAKE_SYSTEM_NAME MATCHES "CYGWIN")
        set(CMAKE_PREFIX_PATH "/cygdrive/c/Qt/Qt/5.15.2/mingw81_64")
        set(Qt${QT_VERSION_MAJOR}_DIR "${CMAKE_PREFIX_PATH}/lib/cmake/Qt5")
        set(Qt${QT_VERSION_MAJOR}Widgets_DIR "${CMAKE_PREFIX_PATH}/lib/cmake/Qt5Widgets")

        ADD_DEFINITIONS(-D _WIN32_PLATFORM_)
    ELSE()
        ADD_DEFINITIONS(-D _UNIX_PLATFORM_)
    ENDIF()

    find_package(Threads)
ENDIF ()

find_package(Qt5 COMPONENTS Widgets REQUIRED)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake;${CMAKE_MODULE_PATH}")

file(GLOB_RECURSE QRC_SOURCE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.qrc)

set(QRC_SOURCE_FILES main.qrc)
qt5_add_resources(QRC_FILES ${QRC_SOURCE_FILES})

set(VCS_ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../..)
ADD_DEFINITIONS(-D CMAKE_SOURCE_DIR="${VCS_ROOT_DIR}")
ADD_DEFINITIONS(-D CMAKE_CURRENT_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

ADD_DEFINITIONS(-D MY_CONFIG_PATH="${CMAKE_CURRENT_SOURCE_DIR}/my_config.h")

link_directories("${VCS_ROOT_DIR}/lib")
include_directories("${VCS_ROOT_DIR}/include")

if(ANDROID)
  add_library(MethyProfile-GUI SHARED
    main.cpp
    mainwindow.cpp
    mainwindow.h
    mainwindow.ui
    ../tool/MethyProfile/bed.h
    ../tool/MethyProfile/bed.cpp
    interface.h
    interface.cpp
    my_config.h
    my_error.hpp
    ${QRC_FILES}
  )
else()
  add_executable(MethyProfile-GUI
    main.cpp
    mainwindow.cpp
    mainwindow.h
    mainwindow.ui
    ../tool/MethyProfile/bed.h
    ../tool/MethyProfile/bed.cpp
    interface.h
    interface.cpp
    my_config.h
    my_error.hpp
    ${QRC_FILES}
  )
endif()

SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)

IF (WIN32 OR (CMAKE_SYSTEM_NAME MATCHES "CYGWIN"))
  target_link_libraries(MethyProfile-GUI PRIVATE Qt5::Widgets)
ELSE()
  target_link_libraries(MethyProfile-GUI PRIVATE Qt5::Widgets ${CMAKE_THREAD_LIBS_INIT} dl libzlog.so)
ENDIF ()
