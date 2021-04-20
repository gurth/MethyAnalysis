include(ExternalProject)
 
set(ZLOG_ROOT ${CMAKE_SOURCE_DIR}/3rdParty/zlog)
set(ZLOG_GIT_TAG  1.2.15)
set(ZLOG_GIT_URL      https://github.com/HardySimpson/zlog.git)
set(ZLOG_CONFIGURE    echo "[  EX] Building libzlog..." )
set(ZLOG_MAKE         cd "${ZLOG_ROOT}" && make PREFIX="${CMAKE_SOURCE_DIR}")
set(ZLOG_INSTALL      cd "${ZLOG_ROOT}" && make PREFIX="${CMAKE_SOURCE_DIR}" install)
 
ExternalProject_Add(zlog
        SOURCE_DIR        ${ZLOG_ROOT}
        PREFIX            ${ZLOG_ROOT}
        CONFIGURE_COMMAND ${ZLOG_CONFIGURE}
        BUILD_COMMAND     ${ZLOG_MAKE}
        INSTALL_COMMAND   ${ZLOG_INSTALL}
)