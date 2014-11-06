# - Setup core required and optional components of ${CMAKE_PROJECT_NAME}
#
# Here we provide options to enable and configure required and optional
# components, which may require third party libraries.
#
# We don't configure User Interface options here because these require
# a higher degree of configuration so to keep things neat these have their
# own Module.
#
# Options configured here:
#
#

########################################################################################
# /usr/local for Linux, /opt/local for Mac
set(SYSTEM_PATHS /usr/local /opt/local /usr)
# Module paths for Linux machines using module load utility
set(MODULE_PATHS $ENV{LOADEDMODULES})
set(HOME_FOLDER $ENV{HOME})
set(INCLUDE_BIN_PATHS $ENV{PATH})
if(NOT "${MODULE_PATHS}" STREQUAL "")
    string(REPLACE ":" ";" MODULE_PATHS ${MODULE_PATHS})
endif()

# Additional possible search paths to look for some root directories
# need to remove bin so that include can be added on end, need to replace colons with semicolons
# for CMake to recognize the variables as separate

string(REPLACE ":" ";" TMP_INCLUDE_BIN_PATHS ${INCLUDE_BIN_PATHS})
set(INCLUDE_BIN_PATHS )
foreach(_path ${TMP_INCLUDE_BIN_PATHS})
    get_filename_component(_path_no_bin ${_path} PATH ABSOLUTE)
    list(APPEND INCLUDE_BIN_PATHS ${_path_no_bin})
endforeach()
list(REMOVE_DUPLICATES INCLUDE_BIN_PATHS)

########################################################################################

include(GenericCMakeFunctions)





########################################################################################
#
#        SILO
#
########################################################################################
option(USE_SILO "Enable the usage of SILO" OFF)

if(USE_SILO)
    set(SILO_FOUND FALSE)

    find_package(SILO)
    if(SILO_FOUND)
        list(APPEND ${PROJECT_NAME}_INCLUDE_DIRS ${SILO_INCLUDE_DIR})
        list(APPEND ${PROJECT_NAME}_LINKED_LIBRARIES ${SILO_LIBRARIES})
        include_directories(${SILO_INCLUDE_DIR})
        #link_directories(${SILO_LIBRARY})

        message(STATUS "SILO FOUND : ${SILO_ROOT_DIR}")
        #message(STATUS "")
        #message(STATUS "SILO_ROOT_DIR : ${SILO_ROOT_DIR}")
        #message(STATUS "SILO_INCLUDE_DIR : ${SILO_INCLUDE_DIR}")
        #message(STATUS "SILO_LIBRARY : ${SILO_LIBRARY}")
        #message(STATUS "SILO_LIBRARIES : ${SILO_LIBRARIES}")
        #message(STATUS "")

        set(${PROJECT_NAME}_SILO_INCLUDE "${SILO_INCLUDE_DIR}")
        string(REPLACE ";" ":" ${PROJECT_NAME}_SILO_INCLUDE ${${PROJECT_NAME}_SILO_INCLUDE})

        #add_definitions(-D_SILO_ -DENABLED_SILO)
        add_package_definitions(SILO)

    else()
        message(WARNING "NO SILO_ROOT FOUND -- ${SILO_ROOT} -- Please set SILO_ROOT")
        remove_package_definitions(SILO)
    endif()
else()
	remove_package_definitions(SILO)
endif()



########################################################################################
#
#        HDF5
#
########################################################################################
option(USE_HDF5 "Enable the usage of HDF5" ${USE_SILO})

if(USE_HDF5)
    set(HDF5_FOUND FALSE)

    find_package(HDF5)
    if(HDF5_FOUND)
        list(APPEND ${PROJECT_NAME}_INCLUDE_DIRS ${HDF5_INCLUDE_DIR})
        list(APPEND ${PROJECT_NAME}_LINKED_LIBRARIES ${HDF5_LIBRARIES})
        include_directories(${HDF5_INCLUDE_DIR})
        #link_directories(${HDF5_LIBRARY})

        message(STATUS "HDF5 FOUND : ${HDF5_ROOT_DIR}")
        #message(STATUS "")
        #message(STATUS "HDF5_ROOT_DIR : ${HDF5_ROOT_DIR}")
        #message(STATUS "HDF5_INCLUDE_DIR : ${HDF5_INCLUDE_DIR}")
        #message(STATUS "HDF5_LIBRARY : ${HDF5_LIBRARY}")
        #message(STATUS "HDF5_LIBRARIES : ${HDF5_LIBRARIES}")
        #message(STATUS "")

        set(${PROJECT_NAME}_HDF5_INCLUDE "${HDF5_INCLUDE_DIR}")
        string(REPLACE ";" ":" ${PROJECT_NAME}_HDF5_INCLUDE ${${PROJECT_NAME}_HDF5_INCLUDE})
		add_package_definitions(HDF5)
		
    else()
        message(WARNING "NO HDF5_ROOT FOUND -- ${HDF5_ROOT} -- Please set HDF5_ROOT")
        remove_package_definitions(HDF5)
    endif()
else()
	remove_package_definitions(HDF5)
endif()


########################################################################################
#
#        OPENGL using X11
#
########################################################################################
OPTION(USE_OPENGL_X11 "Enable the use of OpenGL real-time visualization with X11" OFF)

if(USE_OPENGL_X11)

    IF(APPLE)
        OPTION_NODEPREC_COMPILER_WARNINGS(ON)
    ELSE()
        OPTION_NODEPREC_COMPILER_WARNINGS(OFF)
    ENDIF()

    set(DEFAULT_CMAKE_FIND_FRAMEWORK ${CMAKE_FIND_FRAMEWORK})
    set(CMAKE_FIND_FRAMEWORK NEVER)

    find_package(OpenGL REQUIRED)
    if(OPENGL_FOUND)
        #set(OPENGL_INCLUDE_DIR /opt/local/include)
        #set(OPENGL_LIBRARIES /opt/local/lib/libGL.dylib /opt/local/lib/libGLU.dylib)
        list(APPEND ${PROJECT_NAME}_INCLUDE_DIRS ${OPENGL_INCLUDE_DIR})
        list(APPEND ${PROJECT_NAME}_LINKED_LIBRARIES ${OPENGL_LIBRARIES})
        include_directories(${OPENGL_INCLUDE_DIR})

        message(STATUS "OpenGL FOUND")
        message(STATUS "OpenGL INCLUDE : ${OPENGL_INCLUDE_DIR}")
        message(STATUS "OpenGL LIBRARY : ${OPENGL_LIBRARIES}")

    endif()

    find_package(GLUT REQUIRED)
    if(GLUT_FOUND)
        #set(GLUT_INCLUDE_DIR /opt/local/include)
        #set(GLUT_LIBRARIES /opt/local/lib/libglut.dylib)
        list(APPEND ${PROJECT_NAME}_INCLUDE_DIRS ${GLUT_INCLUDE_DIR})
        list(APPEND ${PROJECT_NAME}_LINKED_LIBRARIES ${GLUT_LIBRARIES})
        include_directories(${GLUT_INCLUDE_DIR})

        message(STATUS "GLUT FOUND")
        message(STATUS "GLUT INCLUDE : ${GLUT_INCLUDE_DIR}")
        message(STATUS "GLUT LIBRARY : ${GLUT_LIBRARIES}")

    endif()

    #find_package(GLEW REQUIRED)
    #if(GLEW_FOUND)
    #    list(APPEND ${PROJECT_NAME}_INCLUDE_DIRS ${GLEW_INCLUDE_DIRS})
    #    list(APPEND ${PROJECT_NAME}_LINKED_LIBRARIES ${GLEW_LIBRARIES})
    #    include_directories(${GLEW_INCLUDE_DIR})

    #    message(STATUS "GLEW FOUND")
    #    message(STATUS "GLEW INCLUDE : ${GLEW_INCLUDE_DIR}")
    #    message(STATUS "GLEW LIBRARY : ${GLEW_LIBRARIES}")

    #endif()


    find_package(X11 REQUIRED)
    if(X11_FOUND)
        list(APPEND ${PROJECT_NAME}_INCLUDE_DIRS ${X11_INCLUDE_DIR})
        list(APPEND ${PROJECT_NAME}_LINKED_LIBRARIES ${X11_LIBRARIES})
        include_directories(${X11_INCLUDE_DIR})

        message(STATUS "X11 FOUND")
        message(STATUS "X11 INCLUDE : ${X11_INCLUDE_DIR}")
        message(STATUS "X11 LIBRARY : ${X11_LIBRARIES}")

    endif()

	add_package_definitions(OPENGL_X11)
	
    set(CMAKE_FIND_FRAMEWORK ${DEFAULT_CMAKE_FIND_FRAMEWORK})

else()
	remove_package_definitions(OPENGL_X11)
endif()


########################################################################################
#
#        TBB - Intel Thread Building Blocks
#
########################################################################################
option(USE_TBB "Enable Intel Thread Building Blocks (TBB)" OFF)

if(USE_TBB)
    set(TBB_FOUND FALSE)

    find_package(TBB)
    if(TBB_FOUND)
        list(APPEND ${PROJECT_NAME}_INCLUDE_DIRS ${TBB_INCLUDE_DIRS})
        list(APPEND ${PROJECT_NAME}_LINKED_LIBRARIES ${TBB_LIBRARIES})
        include_directories(${TBB_INCLUDE_DIRS})
        #link_directories(${HDF5_LIBRARY})

        message(STATUS "TBB FOUND : ${TBB_ROOT_DIR}")
        #message(STATUS "")
        #message(STATUS "TBB_ROOT_DIR : ${TBB_ROOT_DIR}")
        #message(STATUS "TBB_INCLUDE_DIR : ${TBB_INCLUDE_DIR}")
        #message(STATUS "TBB_LIBRARY : ${TBB_LIBRARY}")
        #message(STATUS "TBB_LIBRARIES : ${TBB_LIBRARIES}")
        #message(STATUS "")

        set(${PROJECT_NAME}_TBB_INCLUDE "${TBB_INCLUDE_DIR}")
        string(REPLACE ";" ":" ${PROJECT_NAME}_TBB_INCLUDE ${${PROJECT_NAME}_TBB_INCLUDE})

		add_package_definitions(TBB)
    else()
        message(WARNING "\n\tNO TBB_ROOT FOUND -- ${TBB_ROOT} -- Please set TBB_ROOT\n")
        remove_package_definitions(TBB)
    endif()
else()
	remove_package_definitions(TBB)
endif()


