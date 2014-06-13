#
# Modified/Copied from example in Geant4 version 9.6.2:
#
# http://geant4.web.cern.ch/geant4/
#
# - Define useful macros for building and installing  library targets
#
# This file defines the following macros for  developers needing to
# add shared and static library targets.
#
# LIBRARY_TARGET        - define standard  library targets
#
# The macro will take the name of the library and its sources, defining
# static and shared targets depending on the value of BUILD_SHARED_LIBS and
# BUILD_STATIC_LIBS. Install targets are also created.
#
# A custom compile definition "DEVELOPER_<CONFIG>" is set on
# each target using the target property COMPILE_DEFINITIONS_<CONFIG>
# target property

if(__MACROLIBRARYTARGETS_ISLOADED)
  return()
endif()
set(__MACROLIBRARYTARGETS_ISLOADED TRUE)

include(CMakeMacroParseArguments)

#-----------------------------------------------------------------------
# function compile_definitions_config(<target>)
#          Set a custom compile definition for a  target on a
#          per configuration basis:
#            For mode <CONFIG>, define DEVELOPER_<CONFIG>
#
function(compile_definitions_config _target)
  if(NOT TARGET ${_target})
    message(FATAL_ERROR "compile_definitions_config passed target '${_target}' which is not a valid CMake target")
  endif()

  if(CMAKE_CONFIGURATION_TYPES)
    # - Multimode tools
    foreach(_mode ${CMAKE_CONFIGURATION_TYPES})
      string(TOUPPER ${_mode} _mode_upper)
      set_property(TARGET ${_target}
        APPEND PROPERTY COMPILE_DEFINITIONS_${_mode_upper} DEVELOPER_${_mode_upper}
        )
    endforeach()
  elseif(CMAKE_BUILD_TYPE)
    # - Single mode tools, only if set
    string(TOUPPER ${CMAKE_BUILD_TYPE} _mode_upper)
    set_property(TARGET ${_target}
      APPEND PROPERTY COMPILE_DEFINITIONS_${_mode_upper} DEVELOPER_${_mode_upper}
      )
  endif()
endfunction()


#-----------------------------------------------------------------------
# - LIBRARY_TARGET
# General build and install of a  library target
#
MACRO(LIBRARY_TARGET)
  CMAKE_PARSE_ARGUMENTS(LIBTARGET
    ""
    "NAME" "SOURCES;LINK_LIBRARIES;LINK_LIBRARIES"
    ${ARGN}
    )

  if(BUILD_SHARED_LIBS)
    # Add the shared library target and link its dependencies
    # WIN32 first
    if(WIN32)
      # We have to generate the def export file from an archive library.
      # If we're building Static libraries already, use that existing
      # target, otherwise, build a temporary uninstalled archive...
      if(BUILD_STATIC_LIBS)
        set(_archive ${LIBTARGET_NAME}-static)
      else()
        add_library(_${LIBTARGET_NAME}-archive STATIC EXCLUDE_FROM_ALL ${LIBTARGET_SOURCES})
        set(_archive _${LIBTARGET_NAME}-archive)
      endif()

      # - Add the config specific compile definitions
      compile_definitions_config(${_archive})

      # - Create the .def file for this library
      # Note that we have to pass the actual full path to the library
      # to the command. CMake unfortunately won't generate this for us.
      # Note also that because we're likely to be on a platform with
      # multiconfig build tools. we use the CFG_INTDIR to locate the
      # archive we need...
      add_custom_command(OUTPUT _${LIBTARGET_NAME}.def
        COMMAND genwindef -o _${LIBTARGET_NAME}.def -l ${LIBTARGET_NAME} ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}/${CMAKE_CFG_INTDIR}/${_archive}.lib
        DEPENDS ${_archive} genwindef)

      # - Now we can build the DLL
      # We create it from a dummy empty C++ file plus the def file.
      file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/_${LIBTARGET_NAME}.cpp
        "// empty _${LIBTARGET_NAME}.cpp\n")

      add_library(${LIBTARGET_NAME} SHARED _${LIBTARGET_NAME}.cpp
        _${LIBTARGET_NAME}.def)

      # - Link the DLL.
      # We link it to the archive, and the supplied libraries,
      # but then remove the archive from the LINK_INTERFACE.
      target_link_libraries(${LIBTARGET_NAME}
        ${_archive}
        ${LIBTARGET_LINK_LIBRARIES}
        ${LIBTARGET_LINK_LIBRARIES})

      set_target_properties(${LIBTARGET_NAME}
        PROPERTIES LINK_INTERFACE_LIBRARIES "${LIBTARGET_LINK_LIBRARIES};${LIBTARGET_LINK_LIBRARIES}")

    else()
      # - We build a Shared library in the usual fashion...
      add_library(${LIBTARGET_NAME} SHARED ${LIBTARGET_SOURCES})
      compile_definitions_config(${LIBTARGET_NAME})
      target_link_libraries(${LIBTARGET_NAME}
        ${LIBTARGET_LINK_LIBRARIES}
        ${LIBTARGET_LINK_LIBRARIES})
    endif()

    # This property is set to prevent concurrent builds of static and
    # shared libs removing each others files.
    set_target_properties(${LIBTARGET_NAME}
      PROPERTIES CLEAN_DIRECT_OUTPUT 1)

    # Set the INSTALL_NAME_DIR of the library to its final installation
    # location (Only affects Mac OS X). This will only affect the library
    # when installed, BUT it does hard code this in. One should still be
    # able to bundle up the libraries later as CMake should build the
    # library with headerpad_max_install_names
    set_target_properties(${LIBTARGET_NAME}
      PROPERTIES INSTALL_NAME_DIR ${CMAKE_INSTALL_FULL_LIBDIR})

    # Install the library - note the use of RUNTIME, LIBRARY and ARCHIVE
    # this helps with later DLL builds.
    # Export to standard depends file for later install
    # NEEDS WORK TO REMOVE HARDCODED LIB/BIN DIR
    install(TARGETS ${LIBTARGET_NAME}
      EXPORT LibraryDepends
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
      ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development)

    # Append the library target to a global property so that build tree
    # export of library dependencies can pick up all targets
    set_property(GLOBAL APPEND
      PROPERTY EXPORTED_TARGETS ${LIBTARGET_NAME})
  endif()

  #
  # As above, but for static rather than shared library
  if(BUILD_STATIC_LIBS)
    # We have to distinguish the static from shared lib, so use -static in
    # name. Link its dependencies, and ensure we actually link to the
    # -static targets (We should strictly do this for the external
    # libraries as well if we want a pure static build).
    add_library(${LIBTARGET_NAME}-static STATIC ${LIBTARGET_SOURCES})
    compile_definitions_config(${LIBTARGET_NAME}-static)

    set(LIBTARGET_LINK_LIBRARIES_STATIC )
    foreach(_tgt ${LIBTARGET_LINK_LIBRARIES})
      list(APPEND LIBTARGET_LINK_LIBRARIES_STATIC ${_tgt}-static)
    endforeach()

    target_link_libraries(${LIBTARGET_NAME}-static
      ${LIBTARGET_LINK_LIBRARIES_STATIC}
      ${LIBTARGET_LINK_LIBRARIES})

    # But we can rename the output library to the correct name
    # On WIN32 we *retain* the -static postfix because otherwise
    # we'll conflict with the .lib from the DLL build...
    # We could also install differently...
    if(NOT WIN32)
      set_target_properties(${LIBTARGET_NAME}-static
        PROPERTIES OUTPUT_NAME ${LIBTARGET_NAME})
    endif()

    set_target_properties(${LIBTARGET_NAME}-static
      PROPERTIES CLEAN_DIRECT_OUTPUT 1)

    install(TARGETS ${LIBTARGET_NAME}-static
      EXPORT LibraryDepends
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
      ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development)

    set_property(GLOBAL APPEND
      PROPERTY EXPORTED_TARGETS ${LIBTARGET_NAME}-static)
  endif()
ENDMACRO()

#-----------------------------------------------------------------------
# - HEADER_MODULE_TARGET
# Build and install for a header only  module.
#
MACRO(HEADER_MODULE_TARGET)
  CMAKE_PARSE_ARGUMENTS(HEADERMOD
    ""
    "COMPONENT"
    ""
    ${ARGN}
    )

  # Only has one component, and we just have to pick out the headers
  include(${HEADERMOD_COMPONENT})

  # Header install?
  install(FILES ${${MODULENAME}_HEADERS}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
    COMPONENT Development)

  # Store the include path of the component so that the build tree
  # config file can pick up all needed header paths
  set_property(GLOBAL APPEND
    PROPERTY BUILDTREE_INCLUDE_DIRS ${${MODULENAME}_INCDIR})
ENDMACRO()

#-----------------------------------------------------------------------
# - GRANULAR_LIBRARY_TARGET
# Build and install for a  module (granular) library
#
MACRO(GRANULAR_LIBRARY_TARGET)
  CMAKE_PARSE_ARGUMENTS(GRANLIB
    ""
    "COMPONENT"
    ""
    ${ARGN}
    )

  # Granular lib only has one component, but we must pick out
  # the granular dependencies
  include(${GRANLIB_COMPONENT})

  # Add the library target, using variables set by the inclusion of
  # the component file
  LIBRARY_TARGET(NAME ${MODULENAME}
    SOURCES ${${MODULENAME}_SOURCES} ${${MODULENAME}_HEADERS}
    LINK_LIBRARIES ${${MODULENAME}_GRANULAR_DEPENDENCIES}
    LINK_LIBRARIES ${${MODULENAME}_LINK_LIBRARIES})

  # Header install?
  install(FILES ${${MODULENAME}_HEADERS}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
    COMPONENT Development)

  # Store the include path of the component so that the build tree
  # config file can pick up all needed header paths
  set_property(GLOBAL APPEND
    PROPERTY BUILDTREE_INCLUDE_DIRS ${${MODULENAME}_INCDIR})
ENDMACRO()

#-----------------------------------------------------------------------
# - GLOBAL_LIBRARY_TARGET
# Build and install of a  category (global) library
#
MACRO(GLOBAL_LIBRARY_TARGET)
  CMAKE_PARSE_ARGUMENTS(GLOBLIB
    ""
    "NAME"
    "COMPONENTS"
    ${ARGN}
    )

  # We loop over the component sources one at a time,
  # appending properties as we go.
  foreach(_comp ${GLOBLIB_COMPONENTS})
    #include(${_comp})
    #message(STATUS "GLOBLIB_COMPONENT : ${_comp} : ${MODULENAME} : ${GLOBLIB_NAME}")
    # In case we have a global lib with one component, ensure name gets set
    if(NOT GLOBLIB_NAME)
      set(GLOBLIB_NAME ${MODULENAME})
    endif()

    list(APPEND ${GLOBLIB_NAME}_GLOBAL_SOURCES ${${MODULENAME}_SOURCES})
    list(APPEND ${GLOBLIB_NAME}_GLOBAL_HEADERS ${${MODULENAME}_HEADERS})

    list(APPEND ${GLOBLIB_NAME}_GLOBAL_DEPENDENCIES
      ${${MODULENAME}_GLOBAL_DEPENDENCIES})

    list(APPEND ${GLOBLIB_NAME}_LINK_LIBRARIES
      ${${MODULENAME}_LINK_LIBRARIES})

    list(APPEND ${GLOBLIB_NAME}_BUILDTREE_INCLUDES ${${MODULENAME}_INCDIR})
  endforeach()

  # Filter out duplicates in GLOBAL_DEPENDENCIES and LINK_LIBRARIES
  if(${GLOBLIB_NAME}_GLOBAL_DEPENDENCIES)
    list(REMOVE_DUPLICATES ${GLOBLIB_NAME}_GLOBAL_DEPENDENCIES)
  endif()
  if(${GLOBLIB_NAME}_LINK_LIBRARIES)
    list(REMOVE_DUPLICATES ${GLOBLIB_NAME}_LINK_LIBRARIES)
  endif()

  # Now add the library target
  LIBRARY_TARGET(NAME ${GLOBLIB_NAME}
    SOURCES
    ${${GLOBLIB_NAME}_GLOBAL_SOURCES}
    ${${GLOBLIB_NAME}_GLOBAL_HEADERS}
    LINK_LIBRARIES
    ${${GLOBLIB_NAME}_GLOBAL_DEPENDENCIES}
    LINK_LIBRARIES
    ${${GLOBLIB_NAME}_LINK_LIBRARIES})

  # Header install?
  install(FILES ${${GLOBLIB_NAME}_GLOBAL_HEADERS}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
    COMPONENT Development)

  # Store the include path of the component so that the build tree
  # config file can pick up all needed header paths
  set_property(GLOBAL APPEND
    PROPERTY BUILDTREE_INCLUDE_DIRS ${${GLOBLIB_NAME}_BUILDTREE_INCLUDES})

ENDMACRO()



