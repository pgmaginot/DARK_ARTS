# Try to find the HDF5 library and headers
# Usage of this module is as follows
#
# == Using any header-only components of HDF5: ==
#
#     find_package( HDF5 )
#     if(HDF5_FOUND)
#         include_directories(${HDF5_INCLUDE_DIRS})
#         add_executable(foo foo.cc)
#     endif()
#
# == Using the binary HDF5 library ==
#
#     find_package( HDF5 )
#     if(HDF5_FOUND)
#         include_directories(${HDF5_INCLUDE_DIRS})
#         add_executable(foo foo.cc)
#         target_link_libraries(foo ${HDF5_LIBRARIES})
#     endif()
#
# You can provide a minimum version number that should be used.
# If you provide this version number and specify the REQUIRED attribute,
# this module will fail if it can't find a HDF5 of the specified version
# or higher. If you further specify the EXACT attribute, then this module 
# will fail if it can't find a HDF5 with a version eaxctly as specified.
#
# ===========================================================================
# Variables used by this module which can be used to change the default
# behaviour, and hence need to be set before calling find_package:
#
#  HDF5_ROOT_DIR        The preferred installation prefix for searching for 
#                        HDF5. Set this if the module has problems finding 
#                        the proper HDF5 installation.
#
# If you don't supply HDF5_ROOT_DIR, the module will search on the standard
# system paths. On UNIX, the module will also try to find the HDF5-config
# program in the PATH, and if found will use the prefix supplied by this
# program as a HINT on where to find the HDF5 headers and libraries.
#
# You can re-run CMake with a different version of HDF5_ROOT_DIR to 
# force a new search for HDF5 using the new version of HDF5_ROOT_DIR.
#
# ============================================================================
# Variables set by this module:
#
#  HDF5_FOUND           System has HDF5.
#
#  HDF5_INCLUDE_DIRS    HDF5 include directories: not cached.
#
#  HDF5_LIBRARIES       Link to these to use the HDF5 library: not cached.
#
# ===========================================================================
#
#        BELOW: NOT IMPLEMENTED YET...
#
# ===========================================================================
# If HDF5 is installed in a non-standard way, e.g. a non GNU-style install
# of <prefix>/{lib,include}, then this module may fail to locate the headers
# and libraries as needed. In this case, the following cached variables can
# be editted to point to the correct locations.
#
#  HDF5_INCLUDE_DIR    The path to the HDF5 include directory: cached
#
#  HDF5_LIBRARY        The path to the HDF5 library: cached
# 
# You should not need to set these in the vast majority of cases
#
  


#------------------------------------------------------------------------------
add_feature(HDF5_ROOT "Path to the root folder of HDF5")
#set(HDF5_ROOT_DIR )
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# If the user specifies HDF5_ROOT via environment variable or with -DHDF5_ROOT on cmake line
set(_hdf5_root_hints ${HDF5_ROOT} $ENV{HDF5_ROOT})
if(NOT "${HDF5_ROOT}" STREQUAL "")
    list(APPEND _hdf5_root_hints ${${CMAKE_PROJECT_NAME}_BINARY_DIR}/${HDF5_ROOT})
endif()
if("${_hdf5_root_hints}" STREQUAL "")
    list(APPEND _hdf5_root_hints ${${CMAKE_PROJECT_NAME}_SOURCE_DIR})
    list(APPEND _hdf5_root_hints ${INCLUDE_BIN_PATHS})
    foreach(_syspath ${SYSTEM_PATHS})
        foreach(_module ${MODULE_PATHS})
            list(APPEND _hdf5_root_hints ${_syspath}/${_module})
        endforeach()
    endforeach()
endif()
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
#message(STATUS "HDF5 HINTS : ${_hdf5_root_hints}")
find_path(HDF5_INCLUDE_DIRS hdf5.h
    HINTS ${_hdf5_root_hints}
    PATH_SUFFIXES include
    DOC "Path to the HDF5 headers"
)

find_library(HDF5_LIB_HDF5 hdf5
    HINTS ${_hdf5_root_hints}
    PATH_SUFFIXES lib lib64
    DOC "Path to the HDF5 library"
)

find_library(HDF5_LIB_HDF5_HL hdf5_hl
    HINTS ${_hdf5_root_hints}
    PATH_SUFFIXES lib lib64
    DOC "Path to the HDF5 hl library"
)
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
if(DEFINED HDF5_INCLUDE_DIRS AND NOT "${HDF5_INCLUDE_DIRS}" STREQUAL "HDF5_INCLUDE_DIRS-NOTFOUND")

    set(HDF5_FOUND TRUE)
    
    get_filename_component(_hdf5_root_dir ${HDF5_INCLUDE_DIRS} PATH ABSOLUTE)
    get_filename_component(_hdf5_lib_path ${HDF5_LIB_HDF5} PATH ABSOLUTE)

    set(HDF5_ROOT_DIR ${_hdf5_root_dir})
    set(HDF5_LIBRARY ${_hdf5_lib_path})
    set(HDF5_LIBRARIES ${HDF5_LIB_HDF5} ${HDF5_LIB_HDF5_HL})

    #get_filename_component(_hdf5_lib_name ${HDF5_LIB_HDF5} NAME_WE)
    #get_filename_component(_hdf5_lib_name_hl ${HDF5_LIB_HDF5_HL} NAME_WE)
    #string(REPLACE "lib" "" _hdf5_lib_name ${_hdf5_lib_name})
    #string(REPLACE "lib" "" _hdf5_lib_name_hl ${_hdf5_lib_name_hl})
    #set(HDF5_LIBRARIES ${_hdf5_lib_name} ${_hdf5_lib_name_hl})

    set(HDF5_INCLUDE_DIR ${HDF5_INCLUDE_DIRS})

    #message(STATUS "HDF5_ROOT_DIR : ${HDF5_ROOT_DIR}")
    #message(STATUS "HDF5_INCLUDE_DIRS : ${HDF5_INCLUDE_DIRS}")
    #message(STATUS "HDF5_LIBRARY : ${HDF5_LIBRARY}")
    #message(STATUS "HDF5_LIBRARIES : ${HDF5_LIBRARIES}")
    
else()
    set(HDF5_FOUND FALSE)
endif()
#------------------------------------------------------------------------------



