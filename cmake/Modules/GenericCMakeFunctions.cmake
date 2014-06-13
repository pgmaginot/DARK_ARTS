

include(CMakeMacroParseArguments)


######################################################################################
#   OPTIONS TAKING STRING
######################################################################################
function(ParseCommandLineOptions)

    get_cmake_property(CACHE_VARS CACHE_VARIABLES)
    foreach(CACHE_VAR ${CACHE_VARS})

        get_property(CACHE_VAR_HELPSTRING CACHE ${CACHE_VAR} PROPERTY HELPSTRING)
        if(CACHE_VAR_HELPSTRING STREQUAL "No help, variable specified on the command line.")
            get_property(CACHE_VAR_TYPE CACHE ${CACHE_VAR} PROPERTY TYPE)
            if(CACHE_VAR_TYPE STREQUAL "UNINITIALIZED")
                set(CACHE_VAR_TYPE)
            else()
                set(CACHE_VAR_TYPE :${CACHE_VAR_TYPE})
            endif()
            set(CMAKE_ARGS "${CMAKE_ARGS} -D${CACHE_VAR}${CACHE_VAR_TYPE}=\"${${CACHE_VAR}}\"")
        endif()
    endforeach()
    message("CMAKE_ARGS: ${CMAKE_ARGS}")

endfunction()

######################################################################################
#   Propagate to parent scope
######################################################################################
macro(Propagate SET_VARIABLE)
    set(${SET_VARIABLE} ${${SET_VARIABLE}} PARENT_SCOPE)
endmacro()


######################################################################################
#   Propagate to parent scope
######################################################################################
macro(PropagateSG SET_VARIABLE)
    set(SG_${SET_VARIABLE}_H ${SG_${SET_VARIABLE}_H} PARENT_SCOPE)
    set(SG_${SET_VARIABLE}_I ${SG_${SET_VARIABLE}_I} PARENT_SCOPE)
    set(PROJECT_FOLDERS ${PROJECT_FOLDERS} ${SET_VARIABLE})
    set(PROJECT_FOLDERS ${PROJECT_FOLDERS} PARENT_SCOPE)
    #message(STATUS "Adding ${SET_VARIABLE} to PROJECT_FOLDERS : ${PROJECT_FOLDERS}")
endmacro()


######################################################################################
#   Set project version
######################################################################################
macro(SET_PROJECT_VERSION _MAJOR _MINOR _PATCH _DESCRIPTION)
    set(${PROJECT_NAME}_VERSION_MAJOR "${_MAJOR}")
    set(${PROJECT_NAME}_VERSION_MINOR "${_MINOR}")
    set(${PROJECT_NAME}_VERSION_PATCH "${_PATCH}")
    set(${PROJECT_NAME}_VERSION "${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}.${${PROJECT_NAME}_VERSION_PATCH}")
    set(${PROJECT_NAME}_DESCRIPTION "${_DESCRIPTION}")
endmacro()


######################################################################################
#    Undefined set
######################################################################################
macro(undefset _ARG)
    if(NOT DEFINED ${_ARG})
        set(${_ARG} )
    endif()
endmacro()


######################################################################################
#    Add Package definitions 
######################################################################################
macro(ADD_PACKAGE_DEFINITIONS PKG)
	if(NOT "${PKG}" STREQUAL "")
		#add_definitions(-D_${PKG}_ -DUSING_${PKG})
		list(APPEND ${PROJECT_NAME}_COMPILE_DEFINITIONS _${PKG}_ USING_${PKG})
	endif()
endmacro()


######################################################################################
#    Remove Package definitions 
######################################################################################
macro(REMOVE_PACKAGE_DEFINITIONS PKG)
	if(NOT "${PKG}" STREQUAL "")
		#remove_definitions(-D_${PKG}_ -DUSING_${PKG})
		if(NOT "${${PROJECT_NAME}_COMPILE_DEFINITIONS}" STREQUAL "")
			list(REMOVE_ITEM ${PROJECT_NAME}_COMPILE_DEFINITIONS _${PKG}_ USING_${PKG})
		endif()
	endif()
endmacro()



######################################################################################
#    Remove duplicates if string exists
######################################################################################
macro(REMOVE_DUPLICATES _ARG)
    if(NOT "${${_ARG}}" STREQUAL "")
        #message(STATUS "Removing duplicates from ${_ARG} -- (${${_ARG}})")
        list(REMOVE_DUPLICATES ${_ARG})
    endif()
endmacro()

######################################################################################
#   Load source files
######################################################################################
macro(AddSubdirectoryFiles   NAME
                             NAME_ABBREV
                             SUB_DIR
                             CURRENT_DIR
                             GLOB_HEADER_DIR
                             HEADER_EXT
                             SOURCE_EXT
                             GLOB_HEADERS GLOB_SOURCES
                        )
	set(TOP_DIR ${${PROJECT_NAME}_SOURCE_DIR})
	set(DIR_PATH ${SUB_DIR})
	get_filename_component(TOP_DIR_NAME ${TOP_DIR} NAME)
    get_filename_component(PARENT_DIR ${CURRENT_DIR} PATH)
    #get_filename_component(PARENT_DIR ${PARENT_DIR} NAME)
    string(REPLACE "${TOP_DIR}/" "" PARENT_DIR ${PARENT_DIR})		
    message(STATUS "\tAdding ${PARENT_DIR}/${SUB_DIR}")

    # Include the directory
    include_directories(${PROJECT_SOURCE_DIR}/${SUB_DIR} ${${GLOB_HEADER_DIR}})
    # Link the directory for object files
    link_directories(${PROJECT_BINARY_DIR}/${SUB_DIR})

    # Get files
    file(GLOB ${PROJECT_NAME}_${NAME}_HEADERS ${CURRENT_DIR}/*.${HEADER_EXT})
    file(GLOB ${PROJECT_NAME}_${NAME}_SOURCES ${CURRENT_DIR}/*.${SOURCE_EXT})

    #message("\n\t${PROJECT_NAME}_${NAME}_HEADERS are : ${${PROJECT_NAME}_${NAME}_HEADERS}\n")
    #message("\n\t${PROJECT_NAME}_${NAME}_SOURCES are : ${${PROJECT_NAME}_${NAME}_SOURCES}\n")

    # Append files to global set
    set(${GLOB_HEADERS} ${${GLOB_HEADERS}} ${${PROJECT_NAME}_${NAME}_HEADERS} PARENT_SCOPE)
    set(${GLOB_SOURCES} ${${GLOB_SOURCES}} ${${PROJECT_NAME}_${NAME}_SOURCES} PARENT_SCOPE)

    #message("\n\tGLOBAL Headers are : ${${GLOB_HEADERS}}\n")
    #message("\n\tGLOBAL Sources are : ${${GLOB_SOURCES}}\n")
    
    # set the source group
    set(SG_${NAME_ABBREV}_H ${${PROJECT_NAME}_${NAME}_HEADERS} PARENT_SCOPE)
    set(SG_${NAME_ABBREV}_I ${${PROJECT_NAME}_${NAME}_SOURCES} PARENT_SCOPE)

    #set(${SG_HEADER_SETS} ${${SG_HEADER_SETS}} ${SG_${NAME_ABBREV}_H} PARENT_SCOPE)
    #set(${SG_SOURCE_SETS} ${${SG_SOURCE_SETS}} ${SG_${NAME_ABBREV}_I} PARENT_SCOPE)

    #list(APPEND ${SG_HEADER_NAMES} "${SUB_DIR}/include")
    #list(APPEND ${SG_SOURCE_NAMES} "${SUB_DIR}/src")

    set(${GLOB_HEADER_DIR} ${${GLOB_HEADER_DIR}} ${CURRENT_DIR} PARENT_SCOPE)

endmacro()


######################################################################################
#   Load source files
######################################################################################
macro(AddSubdirectorySourceSorted  NAME
                             NAME_ABBREV
                             SUB_DIR
                             CURRENT_DIR
                             GLOB_HEADER_DIR
                             HEADER_FOLDER
                             SOURCE_FOLDER
                             HEADER_EXT
                             SOURCE_EXT
                             GLOB_HEADERS GLOB_SOURCES
                             #SG_HEADER_SETS SG_SOURCE_SETS
                             #SG_HEADER_NAMES SG_SOURCE_NAMES
                        )

    message(STATUS "\tAdding ${NAME}")

    # Include the directory
    include_directories(${PROJECT_SOURCE_DIR}/${SUB_DIR}/${HEADER_FOLDER} ${${GLOB_HEADER_DIR}})
    # Link the directory for object files
    link_directories(${PROJECT_BINARY_DIR}/${SUB_DIR}/${SOURCE_FOLDER})

    # Get files
    file(GLOB ${PROJECT_NAME}_${NAME}_HEADERS ${CURRENT_DIR}/${HEADER_FOLDER}/*.${HEADER_EXT})
    file(GLOB ${PROJECT_NAME}_${NAME}_SOURCES ${CURRENT_DIR}/${SOURCE_FOLDER}/*.${SOURCE_EXT})

    message("\n\tHeaders are : ${${PROJECT_NAME}_${NAME}_HEADERS}\n")
    message("\n\tSources are : ${${PROJECT_NAME}_${NAME}_SOURCES}\n")

    # Append files to global set
    set(${GLOB_HEADERS} ${${GLOB_HEADERS}} ${${PROJECT_NAME}_${NAME}_HEADERS} PARENT_SCOPE)
    set(${GLOB_SOURCES} ${${GLOB_SOURCES}} ${${PROJECT_NAME}_${NAME}_SOURCES} PARENT_SCOPE)

    #message("\n\tGLOBAL Headers are : ${${GLOB_HEADERS}}\n")
    #message("\n\tGLOBAL Sources are : ${${GLOB_SOURCES}}\n")

    # set the source group
    set(SG_${NAME_ABBREV}_H ${${PROJECT_NAME}_${NAME}_HEADERS} PARENT_SCOPE)
    set(SG_${NAME_ABBREV}_I ${${PROJECT_NAME}_${NAME}_SOURCES} PARENT_SCOPE)

endmacro()













