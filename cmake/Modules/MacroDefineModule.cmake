#
# Modified/Copied from example in Geant4 version 9.6.2:
#
# http://geant4.web.cern.ch/geant4/
#
# - Macros for organizing and specializing code in  modules
#
# This file defines the following macros for  developers needing to
# define the sources, headers and library dependencies for a standard
#  granular library module, plus utlities for specializing source
# file properties
#
# DEFINE_MODULE      - define a standard  Granular Library
#                             Module
#
# ======================
# A  Module is defined as a directory containing subdirectories
#
#   include - holds all header files for the module
#   src     - holds all source files for the module
#
# DEFINE_MODULE will take the name of the module, a list of header
# files, a list of source files and dependencies:
#
# DEFINE_MODULE(NAME name
#                      HEADERS hdr1 hdr2 ...
#                      SOURCES src1 src2 ...
#                      GRANULAR_DEPENDENCIES dep1 dep2 ...
#                      GLOBAL_DEPENDENCIES dep1 dep2 ...
#                      LINK_LIBRARIES lib1 lib2 ...)
#
# It assumes that it will be called from a CMake located in the module
# root directory (i.e. the directory containing the include and src
# subdirectories). It uses this location to define absolute paths to the
# headers and sources, and will ignore any absolute paths passed in HEADERS
# or SOURCES so that generated files can be passed in.
#
# The macro defines the variables
#
# MODULENAME = name
# ${MODULENAME}_HEADERS = List of absolute paths to files given in HEADERS
# ${MODULENAME}_SOURCES = List of absolute paths to files given in SOURCES
# ${MODULENAME}_GRANULAR_DEPS  = List of granular libraries on which this
#                                  module depends
# ${MODULENAME}_GLOBAL_DEPS    = List of global libraries on which this module
#                                  depends
# ${MODULENAME}_LINK_LIBRARIES = List of external libraries to be linked to
#                                  this module
#
# It will also add the module include directory to the list of directories
# using include_directories
#
#
# ADD_COMPILE_DEFINITIONS - add compile defintions to a list of files
# ================================
# ADD_COMPILE_DEFINITIONS(SOURCES src1 src2 ...
#                                COMPILE_DEFINITIONS def1 def 2)
#
# Here, SOURCES is the list of source files to which compile definitions
# will be added. COMPILE_DEFINITIONS gives the list of definitions that
# should be added. These definitions will be appended to any existing
# definitions given to the sources.
#

# - Include guard
if(__macrodefinemodule_isloaded)
  return()
endif()
set(__macrodefinemodule_isloaded YES)

include(CMakeMacroParseArguments)

#-----------------------------------------------------------------------
# macro define_module(NAME <name>
#                            HEADERS <header1> <header2> ... <headerN>
#                            SOURCES <source1> <source2> ... <sourceN>
#                            GRANULAR_DEPENDENCIES <dep1> ... <depN>
#                            GLOBAL_DEPENDENCIES <dep1> ... <depN>
#                            LINK_LIBRARIES <lib1> ... <lib2>)
#       Define a  Module's sources and what internal and external
#       libraries it links to.
#
macro(define_module)
  cmake_parse_arguments(DEFMOD
    ""
    "NAME"
    "HEADERS;SOURCES;GRANULAR_DEPENDENCIES;GLOBAL_DEPENDENCIES;LINK_LIBRARIES"
    ${ARGN}
    )

  set(MODULENAME ${DEFMOD_NAME})

  get_filename_component(${MODULENAME}_BASEDIR ${CMAKE_CURRENT_LIST_FILE} PATH)
  set(${MODULENAME}_SRCDIR ${${MODULENAME}_BASEDIR}/src)
  set(${MODULENAME}_INCDIR ${${MODULENAME}_BASEDIR}/include)

  # We now create absolute paths to the headers and sources,
  # ignoring any already absolute paths
  foreach(_HDR ${DEFMOD_HEADERS})
    if(IS_ABSOLUTE ${_HDR})
      list(APPEND ${MODULENAME}_HEADERS ${_HDR})
    else()
      list(APPEND ${MODULENAME}_HEADERS ${${MODULENAME}_INCDIR}/${_HDR})
    endif()
  endforeach()

  foreach(_SRC ${DEFMOD_SOURCES})
    if(IS_ABSOLUTE ${_SRC})
      list(APPEND ${MODULENAME}_SOURCES ${_SRC})
    else()
      list(APPEND ${MODULENAME}_SOURCES ${${MODULENAME}_SRCDIR}/${_SRC})
    endif()
  endforeach()

  foreach(_LIB ${DEFMOD_GRANULAR_DEPENDENCIES})
    list(APPEND ${MODULENAME}_GRANULAR_DEPENDENCIES ${_LIB})
  endforeach()

  foreach(_LIB ${DEFMOD_GLOBAL_DEPENDENCIES})
    list(APPEND ${MODULENAME}_GLOBAL_DEPENDENCIES ${_LIB})
  endforeach()

  foreach(_LIB ${DEFMOD_LINK_LIBRARIES})
    list(APPEND ${MODULENAME}_LINK_LIBRARIES ${_LIB})
  endforeach()

  include_directories(${${MODULENAME}_INCDIR})
endmacro()

#-----------------------------------------------------------------------
# macro add_compile_definitions(SOURCES <source1> ... <sourceN>
#                                      COMPILE_DEFINITIONS <def1> ... <defN>
#                                      )
#       Add extra compile definitions to a specific list of sources.
#       Macroized to handle the need to specify absolute paths.
#
macro(add_compile_definitions)
  cmake_parse_arguments(ADDDEF
    ""
    ""
    "SOURCES;COMPILE_DEFINITIONS"
    ${ARGN}
    )

  # We assume that the sources have been added at the level of a
  # a sources.cmake, so are inside the src subdir of the sources.cmake
  get_filename_component(_ACD_BASE_PATH ${CMAKE_CURRENT_LIST_FILE} PATH)

  # Now for each file, add the definitions
  foreach(_acd_source ${ADDDEF_SOURCES})
    # Extract any existing compile definitions
    get_source_file_property(
      _acd_existing_properties
      ${_ACD_BASE_PATH}/src/${_acd_source}
      COMPILE_DEFINITIONS)

    if(_acd_existing_properties)
      set(_acd_new_defs ${_acd_existing_properties}
        ${ADDDEF_COMPILE_DEFINITIONS})
    else()
      set(_acd_new_defs ${ADDDEF_COMPILE_DEFINITIONS})
    endif()

    # quote compile defs because this must epand to space separated list
    set_source_files_properties(${_ACD_BASE_PATH}/src/${_acd_source}
      PROPERTIES COMPILE_DEFINITIONS "${_acd_new_defs}")
  endforeach()
endmacro()

