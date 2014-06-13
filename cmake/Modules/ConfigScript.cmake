#
# Modified from example in Geant4 version 9.6.2:
#
# http://geant4.web.cern.ch/geant4/
#
# - Script for configuring and installing ${PROJECT_NAME}-config script
#
# The ${PROJECT_NAME}-config script provides an sh based interface to provide
# information on the ${PROJECT_NAME} installation, including installation prefix,
# version number, compiler and linker flags.
#
# The script is generated from a template file and then installed to the
# known bindir as an executable.
#
# Paths are always hardcoded in the build tree version as this is never
# intended to be relocatable.
# The Install Tree script uses self-location based on that in
# absolute paths are encoded.
#
#

#-----------------------------------------------------------------------
# function get_system_include_dirs
#          return list of directories our C++ compiler searches
#          by default.
#
#          The idea comes from CMake's inbuilt technique to do this
#          for the Eclipse and CodeBlocks generators, but we implement
#          our own function because the CMake functionality is internal
#          so we can't rely on it.
function(get_system_include_dirs _dirs)
  # Only for GCC, Clang and Intel
  if("${CMAKE_CXX_COMPILER_ID}" MATCHES GNU OR "${CMAKE_CXX_COMPILER_ID}" MATCHES Clang OR "${CMAKE_CXX_COMPILER_ID}" MATCHES Intel)
    # Proceed
    file(WRITE "${CMAKE_BINARY_DIR}/CMakeFiles/${PROJECT_NAME}dummy" "\n")

    # Save locale, them to "C" english locale so we can parse in English
    set(_orig_lc_all      $ENV{LC_ALL})
    set(_orig_lc_messages $ENV{LC_MESSAGES})
    set(_orig_lang        $ENV{LANG})

    set(ENV{LC_ALL}      C)
    set(ENV{LC_MESSAGES} C)
    set(ENV{LANG}        C)

    execute_process(COMMAND ${CMAKE_CXX_COMPILER} ${CMAKE_CXX_COMPILER_ARG1} -v -E -x c++ -dD ${PROJECT_NAME}dummy
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/CMakeFiles
      ERROR_VARIABLE _cxxOutput
      OUTPUT_VARIABLE _cxxStdout
      )

    file(REMOVE "${CMAKE_BINARY_DIR}/CMakeFiles/${PROJECT_NAME}dummy")

    # Parse and extract search dirs
    set(_resultIncludeDirs )
    if( "${_cxxOutput}" MATCHES "> search starts here[^\n]+\n *(.+ *\n) *End of (search) list" )
      string(REGEX MATCHALL "[^\n]+\n" _includeLines "${CMAKE_MATCH_1}")
      foreach(nextLine ${_includeLines})
        string(REGEX REPLACE "\\(framework directory\\)" "" nextLineNoFramework "${nextLine}")
        string(STRIP "${nextLineNoFramework}" _includePath)
        list(APPEND _resultIncludeDirs "${_includePath}")
      endforeach()
    endif()

    # Restore original locale
    set(ENV{LC_ALL}      ${_orig_lc_all})
    set(ENV{LC_MESSAGES} ${_orig_lc_messages})
    set(ENV{LANG}        ${_orig_lang})

    set(${_dirs} ${_resultIncludeDirs} PARENT_SCOPE)
  else()
    set(${_dirs} "" PARENT_SCOPE)
  endif()
endfunction()


#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
# ${PROJECT_NAME}_VERSION, ${PROJECT_NAME}_DESCRIPTION should be
# set in main (top-level) CMakeLists.txt
set(PROJECT_DESCRIPTION ${${PROJECT_NAME}_DESCRIPTION})
set(PROJECT_VERSION ${${PROJECT_NAME}_VERSION})
set(PROJECT_LIBRARIES ${${PROJECT_NAME}_INSTALL_LIBRARIES})

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------


#-----------------------------------------------------------------------
# Only create script if we have a global library build...
#
if(UNIX)
  # Get implicit search paths
  get_system_include_dirs(_cxx_compiler_dirs)


  # Configure the script
  # - BUILD TREE
  # Ouch, the include path will be LONG, but at least we always have
  # absolute paths...
  set(CONFIG_SELF_LOCATION "# BUILD TREE IS NON-RELOCATABLE")
  set(CONFIG_INSTALL_PREFIX "${PROJECT_BINARY_DIR}")
  set(CONFIG_INSTALL_EXECPREFIX \"\")
  set(CONFIG_LIBDIR ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

  get_property(_buildtree_include_dirs GLOBAL PROPERTY
    BUILDTREE_INCLUDE_DIRS)

    if(NOT "${_buildtree_include_dirs}" STREQUAL "")
        list(REMOVE_DUPLICATES _buildtree_include_dirs)
          foreach(_dir ${_buildtree_include_dirs})
            set(CONFIG_INCLUDE_DIRS "$CONFIG_INCLUDE_DIRS} \\
            ${_dir}")
          endforeach()
    endif()

    if(NOT "${${PROJECT_NAME}_INCLUDE_DIRS}" STREQUAL "")
        list(REMOVE_DUPLICATES ${PROJECT_NAME}_INCLUDE_DIRS)
        foreach(_dir ${${PROJECT_NAME}_INCLUDE_DIRS})
        set(CONFIG_INCLUDE_DIRS "${CONFIG_INCLUDE_DIRS} \\
            ${_dir}")
        endforeach()
    endif()


    #list(REMOVE_DUPLICATES CONFIG_INCLUDE_DIRS)
    #message(STATUS "CONFIG_INCLUDE_DIRS \t:\t ${CONFIG_INCLUDE_DIRS}")

    foreach(_lib ${${PROJECT_NAME}_LINKED_LIBRARIES})
    set(CONFIG_LINKED_LIBS "${CONFIG_LINKED_LIBS} \\
        ${_lib}")
    endforeach()

  # Configure the build tree script
  # If we're on CMake 2.8 and above, we try to use file(COPY) to create an
  # executable script
  # Not sure if version check is o.k., but I'll be shocked if we ever see
  # a CMake 2.7 in the wild...
  if(${CMAKE_VERSION} VERSION_GREATER 2.7)
    configure_file(
      ${CMAKE_SOURCE_DIR}/cmake/Templates/config.in
      ${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${LC_PROJECT_NAME}-config
      @ONLY
      )

    file(COPY
      ${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${LC_PROJECT_NAME}-config
      DESTINATION ${PROJECT_BINARY_DIR}
      FILE_PERMISSIONS
        OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_EXECUTE
        WORLD_READ WORLD_EXECUTE
      )

  else()
    # Changing permissions is awkward, so just configure and document
    # that you have to do 'sh ${PROJECT_NAME}-config' in this case.
    configure_file(
      ${CMAKE_SOURCE_DIR}/cmake/Templates/config.in
      ${PROJECT_BINARY_DIR}/${LC_PROJECT_NAME}-config
      @ONLY
      )

  endif()

  # - Install Tree
  # Much easier :-)
  # Non-Relocatable case...
  if(CMAKE_INSTALL_IS_NONRELOCATABLE)
    # Hardcoded paths
    set(CONFIG_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")
    set(CONFIG_INSTALL_EXECPREFIX \"\")
    set(CONFIG_LIBDIR "${CMAKE_INSTALL_FULL_LIBDIR}")
    set(CONFIG_INCLUDE_DIRS "${CMAKE_INSTALL_FULL_INCLUDEDIR}")
  else()
    # Calculate base of self contained install based on relative path from
    # CMAKE_INSTALL_FULL_BINDIR to CMAKE_INSTALL_PREFIX.
    file(RELATIVE_PATH _bin_to_prefix ${CMAKE_INSTALL_FULL_BINDIR} ${CMAKE_INSTALL_PREFIX})
    # Strip any trailing path separators just for neatness.
    string(REGEX REPLACE "[/\\]$" "" _bin_to_prefix "${_bin_to_prefix}")

    set(CONFIG_INSTALL_PREFIX "$scriptloc/${_bin_to_prefix}")
    set(CONFIG_INSTALL_EXECPREFIX \"\")
    set(CONFIG_LIBDIR "\${prefix}/${CMAKE_INSTALL_LIBDIR}")
    set(CONFIG_INCLUDE_DIRS "\${prefix}/${CMAKE_INSTALL_INCLUDEDIR}")
  endif()

  # Configure the install tree script
  configure_file(
    ${CMAKE_SOURCE_DIR}/cmake/Templates/config.in
    ${PROJECT_BINARY_DIR}/InstallTreeFiles/${LC_PROJECT_NAME}-config
    @ONLY
    )

  # Install it
  install(FILES ${PROJECT_BINARY_DIR}/InstallTreeFiles/${LC_PROJECT_NAME}-config
    DESTINATION ${CMAKE_INSTALL_BINDIR}
    PERMISSIONS
      OWNER_READ OWNER_WRITE OWNER_EXECUTE
      GROUP_READ GROUP_EXECUTE
      WORLD_READ WORLD_EXECUTE
    COMPONENT Development
    )
endif()

