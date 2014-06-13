# Try to find the STAPL library and headers
# Usage of this module is as follows
#
# == Using any header-only components of STAPL: ==
#
#     find_package( STAPL )
#     if(STAPL_FOUND)
#         include_directories(${STAPL_INCLUDE_DIRS})
#         add_executable(foo foo.cc)
#     endif()
#
# == Using the binary STAPL library ==
#
#     find_package( STAPL )
#     if(STAPL_FOUND)
#         include_directories(${STAPL_INCLUDE_DIRS})
#         add_executable(foo foo.cc)
#         target_link_libraries(foo ${STAPL_LIBRARIES})
#     endif()
#
# You can provide a minimum version number that should be used.
# If you provide this version number and specify the REQUIRED attribute,
# this module will fail if it can't find a STAPL of the specified version
# or higher. If you further specify the EXACT attribute, then this module 
# will fail if it can't find a STAPL with a version eaxctly as specified.
#
# ===========================================================================
# Variables used by this module which can be used to change the default
# behaviour, and hence need to be set before calling find_package:
#
#  STAPL_ROOT_DIR        The preferred installation prefix for searching for 
#                        STAPL. Set this if the module has problems finding 
#                        the proper STAPL installation.
#
# If you don't supply STAPL_ROOT_DIR, the module will search on the standard
# system paths. On UNIX, the module will also try to find the STAPL-config
# program in the PATH, and if found will use the prefix supplied by this
# program as a HINT on where to find the STAPL headers and libraries.
#
# You can re-run CMake with a different version of STAPL_ROOT_DIR to 
# force a new search for STAPL using the new version of STAPL_ROOT_DIR.
#
# ============================================================================
# Variables set by this module:
#
#  STAPL_FOUND           System has STAPL.
#
#  STAPL_INCLUDE_DIRS    STAPL include directories: not cached.
#
#  STAPL_LIBRARIES       Link to these to use the STAPL library: not cached.
#
# ===========================================================================
#
#        BELOW: NOT IMPLEMENTED YET...
#
# ===========================================================================
# If STAPL is installed in a non-standard way, e.g. a non GNU-style install
# of <prefix>/{lib,include}, then this module may fail to locate the headers
# and libraries as needed. In this case, the following cached variables can
# be editted to point to the correct locations.
#
#  STAPL_INCLUDE_DIR    The path to the STAPL include directory: cached
#
#  STAPL_LIBRARY        The path to the STAPL library: cached
# 
# You should not need to set these in the vast majority of cases
#
  

#include(ResolveCompilerPaths)

#------------------------------------------------------------------------------
add_feature(STAPL_ROOT "Path to the root folder of STAPL")
set(STAPL_ROOT_DIR )
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# If the user specifies STAPL_ROOT via environment variable or with -DSTAPL_ROOT on cmake line
set(_stapl_root_hints ${STAPL_ROOT} $ENV{STAPL_ROOT})
if(NOT "${STAPL_ROOT}" STREQUAL "")
    list(APPEND _stapl_root_hints ${${CMAKE_PROJECT_NAME}_BINARY_DIR}/${STAPL_ROOT})
endif()
if("${_stapl_root_hints}" STREQUAL "")
    list(APPEND _stapl_root_hints ${${CMAKE_PROJECT_NAME}_SOURCE_DIR})
    list(APPEND _stapl_root_hints ${INCLUDE_BIN_PATHS})
    foreach(_syspath ${SYSTEM_PATHS})
        foreach(_module ${MODULE_PATHS})
            list(APPEND _stapl_root_hints ${_syspath}/${_module})
        endforeach()
    endforeach()
endif()
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
#message(STATUS "STAPL HINTS : ${_stapl_root_hints}")
find_path(STAPL_ROOT_DIR include/utility/stapl_location.h
    HINTS ${_stapl_root_hints}
    PATH_SUFFIXES stapl
    DOC "Path to the STAPL headers"
)
find_library(STAPL_LIB_STAPL stapl
    HINTS ${_stapl_root_hints}
    PATHS ${STAPL_ROOT_DIR}
    PATH_SUFFIXES lib lib64
    DOC "Path to the STAPL headers"
)
find_library(STAPL_LIB_STAPL_RT stapl_rt
    HINTS ${_stapl_root_hints}
    PATHS ${STAPL_ROOT_DIR}
    PATH_SUFFIXES lib lib64
    DOC "Path to the STAPL headers"
)
#------------------------------------------------------------------------------

#message(STATUS "\n\nSTAPL_LIB_STAPL : ${STAPL_LIB_STAPL}\n\n")
#------------------------------------------------------------------------------
if(DEFINED STAPL_ROOT_DIR AND NOT "${STAPL_ROOT_DIR}" STREQUAL "STAPL_ROOT_DIR-NOTFOUND")

    set(STAPL_FOUND TRUE)
    set(STAPL_INCLUDE_DIR ${STAPL_ROOT_DIR}
                  ${STAPL_ROOT_DIR}/include
                          ${STAPL_ROOT_DIR}/stapl
                          ${STAPL_ROOT_DIR}/tools
                          ${STAPL_ROOT_DIR}/tools/libstdc++/${_stl_version}
                )
                

    get_filename_component(_stapl_lib_path ${STAPL_LIB_STAPL} PATH ABSOLUTE)
    set(STAPL_LIBRARY ${_stapl_lib_path})
    set(STAPL_INCLUDE_DIRS ${STAPL_INCLUDE_DIR})
    set(STAPL_LIBRARIES ${STAPL_LIB_STAPL} ${STAPL_LIB_STAPL_RT})

    #message(STATUS "STAPL_ROOT_DIR : ${STAPL_ROOT_DIR}")
    #message(STATUS "STAPL_INCLUDE_DIR : ${STAPL_INCLUDE_DIR}")
    #message(STATUS "STAPL_LIBRARY : ${STAPL_LIBRARY}")
    #message(STATUS "STAPL_LIBRARIES : ${STAPL_LIBRARIES}")

else()
    set(STAPL_FOUND FALSE)
endif()
#------------------------------------------------------------------------------




