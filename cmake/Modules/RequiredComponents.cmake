# - Setup core required and required components of 
#
#
# Components configured here:
#
#    ${PROJECT_NAME}_LINKED_LIBRARIES
#    ${PROJECT_NAME}_INCLUDE_DIRS
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




########################################################################################
#
#        BOOST
#
########################################################################################
option(USE_BOOST "Enable the use of BOOST" ON)
add_feature(BOOST_ROOT "The root location of Boost" $ENV{BOOST_ROOT})
add_feature(BOOST_ROOT_DIR "The root location of Boost" $ENV{BOOST_ROOT})

if(USE_BOOST)
    ####### POSSIBLE BOOST LIBRARIES FOR  #######
    set(USING_BOOST_LIBRARIES
                                                # atomic
                                                # chrono
                                                # context
                                                # date_time
                                                # exception
                                                # filesystem
                                                # graph
                                                # graph_parallel
                                                # iostreams
                                                # locale
                                                # math
                                                # math_c99f
                                                # math_c99l
                                                # math_c99
                                                # math_tr1f
                                                # math_tr1l
                                                # math_tr1
                                                # prg_exec_monitor
                                                # mpi
                                                # python
                                                # regex
                                                # program_options
                                                random
                                                # serialization
                                                # signals
                                                # system
                                                # test
                                                # thread
                                                # timer
                                                # unit_test_framework
                                                # wave
                                                # wserialization
                        )

    if(NOT "$ENV{BOOST_ROOT}" STREQUAL "" OR NOT "${BOOST_ROOT}" STREQUAL "")
        if(NOT "$ENV{BOOST_ROOT}" STREQUAL "")
            set(Boost_DIR $ENV{BOOST_ROOT})
        elseif(NOT "${BOOST_ROOT}" STREQUAL "")
            set(Boost_DIR ${BOOST_ROOT})
        endif()
    endif()

    set(Boost_NO_BOOST_CMAKE         OFF)
    find_package(Boost 1.53 COMPONENTS ${USING_BOOST_LIBRARIES})
    if(Boost_FOUND)
        include_directories(${Boost_INCLUDE_DIRS})
        #link_directories(${Boost_LIBRARY_DIRS})
        list(APPEND ${PROJECT_NAME}_LINKED_LIBRARIES ${Boost_LIBRARIES})
        list(APPEND ${PROJECT_NAME}_INCLUDE_DIRS ${Boost_INCLUDE_DIRS})

        #message(STATUS "")
        #message(STATUS "BOOST FOUND")
        #set(Boost_DIR ${BOOST_ROOT})
        #message(STATUS "Boost FOUND : ${Boost_DIR}")
        #message(STATUS "Boost_INCLUDE_DIRS : ${Boost_INCLUDE_DIRS}")
        #message(STATUS "Boost_LIBRARY_DIRS : ${Boost_LIBRARY_DIRS}")
        #message(STATUS "Boost_LIBRARIES : ${Boost_LIBRARIES}")
        #message(STATUS "")
        add_package_definitions(BOOST)

    else()
        message(WARNING "BOOST NOT FOUND -- ${BOOST_ROOT} -- Please set BOOST_ROOT")
        #remove_package_definitions(BOOST)        
    endif(Boost_FOUND)
else()
    remove_package_definitions(BOOST)
endif()

