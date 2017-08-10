INCLUDE(FindPkgConfig)
PKG_CHECK_MODULES(PC_UWSPR uwspr)

FIND_PATH(
    UWSPR_INCLUDE_DIRS
    NAMES uwspr/api.h
    HINTS $ENV{UWSPR_DIR}/include
        ${PC_UWSPR_INCLUDEDIR}
    PATHS ${CMAKE_INSTALL_PREFIX}/include
          /usr/local/include
          /usr/include
)

FIND_LIBRARY(
    UWSPR_LIBRARIES
    NAMES gnuradio-uwspr
    HINTS $ENV{UWSPR_DIR}/lib
        ${PC_UWSPR_LIBDIR}
    PATHS ${CMAKE_INSTALL_PREFIX}/lib
          ${CMAKE_INSTALL_PREFIX}/lib64
          /usr/local/lib
          /usr/local/lib64
          /usr/lib
          /usr/lib64
)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(UWSPR DEFAULT_MSG UWSPR_LIBRARIES UWSPR_INCLUDE_DIRS)
MARK_AS_ADVANCED(UWSPR_LIBRARIES UWSPR_INCLUDE_DIRS)

