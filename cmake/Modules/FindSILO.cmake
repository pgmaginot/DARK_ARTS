# Try to find the SILO library and headers
# Usage of this module is as follows
#
# == Using any header-only components of SILO: ==
#
#     find_package( SILO )
#     if(SILO_FOUND)
#         include_directories(${SILO_INCLUDE_DIRS})
#         add_executable(foo foo.cc)
#     endif()
#
# == Using the binary SILO library ==
#
#     find_package( SILO )
#     if(SILO_FOUND)
#         include_directories(${SILO_INCLUDE_DIRS})
#         add_executable(foo foo.cc)
#         target_link_libraries(foo ${SILO_LIBRARIES})
#     endif()
#
# You can provide a minimum version number that should be used.
# If you provide this version number and specify the REQUIRED attribute,
# this module will fail if it can't find a SILO of the specified version
# or higher. If you further specify the EXACT attribute, then this module 
# will fail if it can't find a SILO with a version eaxctly as specified.
#
# ===========================================================================
# Variables used by this module which can be used to change the default
# behaviour, and hence need to be set before calling find_package:
#
#  SILO_ROOT_DIR        The preferred installation prefix for searching for 
#                        SILO. Set this if the module has problems finding 
#                        the proper SILO installation.
#
# If you don't supply SILO_ROOT_DIR, the module will search on the standard
# system paths. On UNIX, the module will also try to find the SILO-config
# program in the PATH, and if found will use the prefix supplied by this
# program as a HINT on where to find the SILO headers and libraries.
#
# You can re-run CMake with a different version of SILO_ROOT_DIR to 
# force a new search for SILO using the new version of SILO_ROOT_DIR.
#
# ============================================================================
# Variables set by this module:
#
#  SILO_FOUND           System has SILO.
#
#  SILO_INCLUDE_DIRS    SILO include directories: not cached.
#
#  SILO_LIBRARIES       Link to these to use the SILO library: not cached.
#
# ===========================================================================
#
#        BELOW: NOT IMPLEMENTED YET...
#
# ===========================================================================
# If SILO is installed in a non-standard way, e.g. a non GNU-style install
# of <prefix>/{lib,include}, then this module may fail to locate the headers
# and libraries as needed. In this case, the following cached variables can
# be editted to point to the correct locations.
#
#  SILO_INCLUDE_DIR    The path to the SILO include directory: cached
#
#  SILO_LIBRARY        The path to the SILO library: cached
# 
# You should not need to set these in the vast majority of cases
#

#include(ResolveCompilerPaths)

#------------------------------------------------------------------------------
add_feature(SILO_ROOT "Path to the root folder of SILO")
set(SILO_ROOT_DIR SILO_ROOT_DIR-NOTFOUND)
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# If the user specifies SILO_ROOT via environment variable or with -DSILO_ROOT on cmake line
set(_silo_root_hints ${SILO_ROOT} $ENV{SILO_ROOT})
if(NOT "${SILO_ROOT}" STREQUAL "")
    list(APPEND _silo_root_hints ${${CMAKE_PROJECT_NAME}_BINARY_DIR}/${SILO_ROOT})
endif()
if("${_silo_root_hints}" STREQUAL "")
    list(APPEND _silo_root_hints ${${CMAKE_PROJECT_NAME}_SOURCE_DIR})
    list(APPEND _silo_root_hints ${INCLUDE_BIN_PATHS})
    foreach(_syspath ${SYSTEM_PATHS})
        foreach(_module ${MODULE_PATHS})
            list(APPEND _silo_root_hints ${_syspath}/${_module})
        endforeach()
    endforeach()
endif()
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
#message(STATUS "SILO HINTS : ${_silo_root_hints}")
find_path(SILO_INCLUDE_DIR silo.h
    HINTS ${_silo_root_hints}
    PATH_SUFFIXES include
    DOC "Path to the SILO headers"
)

find_library(SILO_LIBRARY_SILOH5 siloh5
    HINTS ${_silo_root_hints}
    PATH_SUFFIXES lib lib64
    DOC "Path to the SILO headers"
)
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
if(DEFINED SILO_INCLUDE_DIR AND NOT "${SILO_INCLUDE_DIR}" STREQUAL "SILO_INCLUDE_DIR-NOTFOUND")

    set(SILO_FOUND TRUE)
    get_filename_component(_silo_root_dir ${SILO_INCLUDE_DIR} PATH ABSOLUTE)
    get_filename_component(_silo_lib_path ${SILO_LIBRARY_SILOH5} PATH ABSOLUTE)
    set(SILO_ROOT_DIR ${_silo_root_dir})
    set(SILO_LIBRARY ${_silo_lib_path})
    set(SILO_LIBRARIES ${SILO_LIBRARY_SILOH5})

    #get_filename_component(_silo_lib_name ${SILO_LIBRARY_SILOH5} NAME_WE)
    #string(REPLACE "lib" "" _silo_lib_name ${_silo_lib_name})
    #set(SILO_LIBRARIES ${_silo_lib_name})


    set(SILO_INCLUDE_DIRS ${SILO_INCLUDE_DIR})

    #message(STATUS "SILO_ROOT_DIR : ${SILO_ROOT_DIR}")
    #message(STATUS "SILO_INCLUDE_DIR : ${SILO_INCLUDE_DIR}")
    #message(STATUS "SILO_LIBRARY : ${SILO_LIBRARY}")
    #message(STATUS "SILO_LIBRARIES : ${SILO_LIBRARIES}")
    
    #resolve_libraries(SILO_LIBRARIES ${SILO_LIBRARY})
    
else()
    set(SILO_FOUND FALSE)
endif()
#------------------------------------------------------------------------------




