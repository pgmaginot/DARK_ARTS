#==============================================================================
# Top Level CMakeLists.txt for DARK ARTS arts Build
# Dfem And Runge Kutta for Advanced Radiative Transfer in Slabs
# 3rd January 2014 Jonathan Madsen
#
#





#==============================================================================
# - Enforce an out-of-source builds before anything else
#
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(STATUS "")
    message(STATUS "")
    message(STATUS "\t DARK_ARTS requires an out-of-source build.")
    message(STATUS "\t Please remove these files from ${CMAKE_BINARY_DIR} first:")
    message(STATUS "\t\t CMakeCache.txt")
    message(STATUS "\t\t CMakeFiles")
    message(STATUS "\t Once these files are removed, create a separate directory")
    message(STATUS "\t inside or outside of the source code folder ${CMAKE_PROJECT_DIR}")
    message(STATUS "\t and run CMake from there")
    message(STATUS "")
    message(STATUS "")
    message(FATAL_ERROR "in-source build detected")
endif()





#==============================================================================
# - Define CMake requirements and override make rules as needed
#
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
set(CMAKE_USER_MAKE_RULES_OVERRIDE_CXX
    ${CMAKE_SOURCE_DIR}/cmake/Modules/MakeRules_cxx.cmake)





#==============================================================================
# - Project definition
#
project(DARK_ARTS)
#-------------------------------------------------------
# Lowercase project name (used for EXE name and scripts)
set(LC_PROJECT_NAME ${PROJECT_NAME})
string(TOLOWER ${LC_PROJECT_NAME} LC_PROJECT_NAME)





#==============================================================================
# - Prepend our own CMake Modules to the search path
# This allows us to customize, add, and factor out functionality
# NB: if our custom modules include others that we don't supply, those in the
# base path will be used, so watch for incompatibilities!!
#
#-------------------------------------------------------
set(CMAKE_MODULE_PATH
    ${PROJECT_SOURCE_DIR}/cmake/Modules
    ${CMAKE_MODULE_PATH})





#==============================================================================
include(GenericCMakeOptions)
include(GenericCMakeFunctions)
include(MacroUtilities)
#-------------------------------------------
SET_PROJECT_VERSION(0    # Major
                    0    # Minor
                    0    # Patch
                    "DARK_ARTS")
#-------------------------------------------
undefset(${PROJECT_NAME}_LINKED_LIBRARIES)
undefset(${PROJECT_NAME}_LIBRARY_DIR)
undefset(${PROJECT_NAME}_INCLUDE_DIRS)
undefset(${PROJECT_NAME}_SOURCE_FILES)
undefset(${PROJECT_NAME}_HEADER_FILES)
undefset(${PROJECT_NAME}_EXCLUDE_FILES)
undefset(PROJECT_FOLDERS)
set(${PROJECT_NAME}_LIB_NAME ${PROJECT_NAME}-${${PROJECT_NAME}_VERSION})
# ${PROJECT_NAME}_EXE_NAME and ${PROJECT_NAME}_EXE_PATH
# is set in source/CMakeLists.txt and propagated to this scope





############################################
#   OPTIONS (ON/OFF)
############################################
option(SORTED "Sort the folder groups into another level of include/ and src/ folders" ON)
option(VERBOSE "Turn on verbose Makefile" OFF)
option(USE_OPENGL "Turn of visualization with OpenGL" OFF)
if(VERBOSE)
    set(CMAKE_VERBOSE_MAKEFILE ON)
else(VERBOSE)
    set(CMAKE_VERBOSE_MAKEFILE OFF)
endif()
OPTION_SET_BUILD_RPATH(ON)
OPTION_RETURN_TYPE_COMPILER_WARNINGS(ON)
OPTION_ACCELERATE_FRAMEWORK(ON)
OPTION_FAST_MATH(OFF)
ENABLE_OPENMP_SUPPORT(OFF)




#==============================================================================
# - Provide an 'uninstall' target.
include(CMakeUninstallTarget)
#-------------------------------------------------------
# - Provide options to control how libraries are built
include(LibraryBuildOptions)
#-------------------------------------------------------
# - Provide standard install directories permitting customization.
include(InstallDirs)





#==============================================================================
#    Configure required components
include(RequiredComponents)
#-------------------------------------------------------
#    Configure optional components (specifically the testing components)
# - currently configures nothing
include(OptionalComponents)





#==============================================================================
# Add the source and environments subdirectories
add_subdirectory(src)
#-------------------------------------------------------
include_directories(${${PROJECT_NAME}_INCLUDE_DIRS})
remove_duplicates(${PROJECT_NAME}_LINKED_LIBRARIES)
remove_duplicates(PROJECT_FOLDERS)




#==============================================================================
# Generate the ${PROJECT_NAME}.sh/.csh files for adding project to
# PATH and {DY}LD_LIBARY_PATH environment variables
include(SourceScript)
# Generate the ${LC_PROJECT_NAME}-config executable script for getting the
# compile flags, link flags, etc.
include(ConfigScript)





#message(STATUS "Linked libraries : ${${PROJECT_NAME}_LINKED_LIBRARIES}")
#message(STATUS "${PROJECT_NAME}_COMPILE_DEFINITIONS : ${${PROJECT_NAME}_COMPILE_DEFINITIONS}")

foreach(_definition ${${PROJECT_NAME}_COMPILE_DEFINITIONS})
	add_definitions(-D${_definition})
endforeach()

#file(GLOB resources ${PROJECT_SOURCE_DIR}/examples/*.xml)

#############################################
#    LIBRARIES and EXECUTABLES
#############################################
add_library(${${PROJECT_NAME}_LIB_NAME} ${${PROJECT_NAME}_HEADER_FILES} ${${PROJECT_NAME}_SOURCE_FILES})
target_link_libraries(${${PROJECT_NAME}_LIB_NAME} ${${PROJECT_NAME}_LINKED_LIBRARIES} ${EXTRA_LIBS})

set(${PROJECT_NAME}_INSTALL_EXECUTABLES ${${PROJECT_NAME}_EXE_NAME})
set(${PROJECT_NAME}_INSTALL_LIBRARIES ${${PROJECT_NAME}_LIB_NAME})

add_executable(${${PROJECT_NAME}_EXE_NAME} ${${PROJECT_NAME}_EXE_PATH}/${${PROJECT_NAME}_EXE_NAME}.cc ${resources})
target_link_libraries(${${PROJECT_NAME}_EXE_NAME} ${${PROJECT_NAME}_LIB_NAME} ${${PROJECT_NAME}_LINKED_LIBRARIES} ${EXTRA_LIBS})





#==============================================================================
# Sort the files into source groups for when using an IDE
foreach(_folder ${PROJECT_FOLDERS})
    set(sg_name ${_folder})
    string(REPLACE "." "\\\\" _folder ${_folder})
    string(TOUPPER ${sg_name} sg_name)

    #message(STATUS "Project folders : ${_folder} for ${sg_name}")

    if(SORTED)
        if(NOT "${SG_${sg_name}_I}" STREQUAL "")
            set(inc_folder "${_folder}\\include")
        else()
            set(inc_folder "${_folder}")
        endif()
        set(src_folder "${_folder}")

        source_group(${inc_folder} FILES ${SG_${sg_name}_H})
        source_group(${src_folder} FILES ${SG_${sg_name}_I})

    else()
        set(sg_folder "${_folder}")
        source_group(${sg_folder} FILES ${SG_${sg_name}_H})
        source_group(${sg_folder} FILES ${SG_${sg_name}_I})

    endif()

endforeach()




#==============================================================================
# Put example XML files in Resources folder for Xcode IDE
if(XCODE)
	message(STATUS "Adding example XML files to Resources folder")
    set_source_files_properties(
        ${resources}
        PROPERTIES
        MACOSX_PACKAGE_LOCATION Resources
    )
endif()




#==============================================================================
# TESTING
#-------------------------------------------
if(BUILD_TESTS)
	enable_testing()
	add_subdirectory(testing)
endif()




#==============================================================================
# Final output - show what's been enabled so that user knows what's
# happening - also useful for later problem solving!
#
#-------------------------------------------
PRINT_ENABLED_FEATURES()
PRINT_DISABLED_FEATURES()
message(STATUS "${PROJECT_NAME} Build Type : ${CMAKE_BUILD_TYPE}")
#-------------------------------------------
set(${PROJECT_NAME}_INSTALL_HEADERS )
foreach(_header ${${PROJECT_NAME}_HEADER_FILES})
    get_filename_component(_install_header ${_header} NAME)
    list(APPEND ${PROJECT_NAME}_INSTALL_HEADERS ${_install_header})
endforeach()



#==============================================================================
# Install under CMAKE_INSTALL_PREFIX
#==============================================================================
install(TARGETS ${${PROJECT_NAME}_INSTALL_EXECUTABLES} DESTINATION ${CMAKE_INSTALL_BINDIR})
install(TARGETS ${${PROJECT_NAME}_INSTALL_LIBRARIES} DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(FILES ${${PROJECT_NAME}_HEADER_FILES} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
#install(FILES ${${PROJECT_NAME}_DATA_FILES} DESTINATION ${CMAKE_INSTALL_DATADIR})

