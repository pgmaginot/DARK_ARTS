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
#        MPI
#
########################################################################################
option(USE_MPI "Turn on MPI" ON)
if(USE_MPI)
    ENABLE_MPI_SUPPORT(ON)
    if(MPI_FOUND)
        message(STATUS "MPI FOUND : ${MPI_ROOT}")

        #if(XCODE OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
            # set the C and C++ compilers to MPI compilers
        #    message(STATUS "")
        #    message(STATUS "\tXCODE support for MPI")
        #    message(STATUS "")
        #    set(CMAKE_C_COMPILER ${MPI_CC_COMPILER})
        #    set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
        #    set($ENV{OMPI_CXX} clang++)
        #    set($ENV{OMPI_CC} clang)
        #endif()

        include_directories(${MPI_INCLUDE_PATH})
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_CC_COMPILE_FLAGS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_LINK_FLAGS}")


        list(APPEND ${PROJECT_NAME}_LINKED_LIBRARIES ${MPI_CC_LIBRARIES} ${MPI_CXX_LIBRARIES})
        list(APPEND ${PROJECT_NAME}_INCLUDE_DIRS ${MPI_INCLUDE_PATH})

        message(STATUS "MPI_INCLUDE_PATH : ${MPI_INCLUDE_PATH}")
        message(STATUS "MPI  CC Libraries : ${MPI_CC_LIBRARIES}")
        message(STATUS "MPI CXX Libraries : ${MPI_CXX_LIBRARIES}")
		message(STATUS "Adding MPI package definitions...")
        add_package_definitions(MPI)
        
    else()
        message(WARNING "MPI NOT FOUND - ${PROJECT_NAME} is best run with MPI")
        remove_package_definitions(MPI)
    endif()
else()
	remove_package_definitions(MPI)
endif()

########################################################################################
#
#        PETSc
#
########################################################################################

option(USE_PETSC "Use the PETSc library" ON)
if(USE_PETSC)
  message(STATUS "Trying to tie into the PETSc library")
#  SET( PETSC_DIR "$ENV{PETSC_DIR}" )
  
 # IF( NOT PETSC_DIR )
 #   MESSAGE( FATAL_ERROR "Please point the environment variable PETSC_DIR to the include directory of your PETSc installation.")
 # ENDIF()
  
  find_package(PETSc REQUIRED)
  
  if( NOT PETSC_FOUND)
    MESSAGE( FATAL_ERROR "Could Not find PETSc")
  else()
    MESSAGE(STATUS "Found PETSc")
  ENDIF()
  
  # grove-01.ne.tamu.edu {~/Research/dark_arts_local/build }202 :which mpicc
  # /usr/local/openmpi-1.6.3-gcc-4.7.2/bin/mpicc

  
  message( STATUS "Here are the PETSc includes: ${PETSC_INCLUDES}")
  message( STATUS "PETSc was compiled with this compiler: ${PETSC_COMPILER}" )
  message( STATUS "PETSc version: ${PETSC_VERSION}" )
  
#  PETSC_INCLUDES     - the PETSc include directories
#  PETSC_LIBRARIES    - Link these to use PETSc
#  PETSC_COMPILER     - Compiler used by PETSc, helpful to find a compatible MPI
#  PETSC_DEFINITIONS  - Compiler switches for using PETSc
#  PETSC_MPIEXEC      - Executable for running MPI programs
#  PETSC_VERSION 

  list(APPEND ${PROJECT_NAME}_LINKED_LIBRARIES ${PETSC_LIBRARIES})
  list(APPEND ${PROJECT_NAME}_INCLUDE_DIRS ${PETSC_INCLUDES})
  
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${PETSC_DEFINITIONS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${PETSC_DEFINITIONS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${PETSC_DEFINITIONS}")

endif()


########################################################################################
#
#        BOOST
#
########################################################################################
option(USE_BOOST "Enable the use of BOOST" OFF)
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
                                                # random
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

