

# - Include guard
#if(__genericcmakeoptions_isloaded)
#  return()
#endif()
#set(__genericcmakeoptions_isloaded YES)



######################################################################################
# Library Build Type Options
######################################################################################

MACRO(OPTION_LIBRARY_BUILD_STATIC DEFAULT)
    OPTION(BUILD_STATIC "Build static libraries" ${DEFAULT})

    IF(BUILD_STATIC)
        LIST(APPEND NP_LIB_TYPES STATIC)
        MESSAGE(STATUS "Building Static Libraries for ${CMAKE_PROJECT_NAME}")
    ELSE(BUILD_STATIC)
        MESSAGE(STATUS "NOT Building Static Libraries for ${CMAKE_PROJECT_NAME}")
    ENDIF(BUILD_STATIC)
ENDMACRO(OPTION_LIBRARY_BUILD_STATIC)

MACRO(OPTION_LIBRARY_BUILD_SHARED DEFAULT)
    OPTION(BUILD_SHARED "Build shared libraries" ${DEFAULT})

    IF(BUILD_SHARED)
        LIST(APPEND NP_LIB_TYPES SHARED)
        MESSAGE(STATUS "Building Shared Libraries for ${CMAKE_PROJECT_NAME}")
    ELSE(BUILD_SHARED)
        MESSAGE(STATUS "NOT Building Shared Libraries for ${CMAKE_PROJECT_NAME}")
    ENDIF(BUILD_SHARED)
ENDMACRO(OPTION_LIBRARY_BUILD_SHARED)

######################################################################################
# RPATH Relink Options
######################################################################################

MACRO(OPTION_SET_BUILD_RPATH DEFAULT)
    OPTION(SET_BUILD_RPATH "Set the build RPATH to local directories, relink to install directories at install time" ${DEFAULT})

    IF(SET_BUILD_RPATH)
        if(APPLE)
            set(CMAKE_MACOSX_RPATH 1)
        endif(APPLE)
        
        # the RPATH to be used when installing
        if(DEFINED LIBRARY_INSTALL_DIR)
            SET(CMAKE_INSTALL_RPATH "${LIBRARY_INSTALL_DIR}")
        else()
            SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
        endif()

        MESSAGE(STATUS "Setting build RPATH for ${CMAKE_PROJECT_NAME} to ${CMAKE_INSTALL_RPATH}")
        # use, i.e. don't skip the full RPATH for the build tree
        SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

        # when building, don't use the install RPATH already
        # (but later on when installing)
        SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
        
        # add the automatically determined parts of the RPATH
        # which point to directories outside the build tree to the install RPATH
        SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
        
        # the RPATH to be used when installing, but only if it's not a system directory
        LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
        IF("${isSystemDir}" STREQUAL "-1")
           SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
        ENDIF("${isSystemDir}" STREQUAL "-1")
        
    ENDIF(SET_BUILD_RPATH)
ENDMACRO(OPTION_SET_BUILD_RPATH)

######################################################################################
# Create software version code file
######################################################################################

MACRO(OPTION_CREATE_VERSION_FILE DEFAULT OUTPUT_FILES)
    OPTION(CREATE_VERSION_FILE "Creates a version.cc file using the setlocalversion script" ${DEFAULT})
    IF(CREATE_VERSION_FILE)
        MESSAGE(STATUS "Generating git information for ${CMAKE_PROJECT_NAME}")
        FOREACH(VERSION_FILE ${OUTPUT_FILES})
            EXECUTE_PROCESS(COMMAND "${FIVETEN_SCRIPT_DIR}/setlocalversion.sh" OUTPUT_FILE ${VERSION_FILE})
            MESSAGE(STATUS "- Generating to ${VERSION_FILE}")
        ENDFOREACH(VERSION_FILE ${OUTPUT_FILES})
    ELSE(CREATE_VERSION_FILE)
        MESSAGE(STATUS "NOT generating git information for ${CMAKE_PROJECT_NAME}")
    ENDIF(CREATE_VERSION_FILE)
ENDMACRO(OPTION_CREATE_VERSION_FILE)

######################################################################################
# Fast Math Option
######################################################################################

MACRO(OPTION_FAST_MATH FALSE)
    IF(CMAKE_COMPILER_IS_GNUCXX)
        OPTION(FAST_MATH "Use -ffast-math for GCC 4.0" DEFAULT)

        IF(FAST_MATH)
            SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffast-math")
            MESSAGE(STATUS "Turning on -ffast-math for ${CMAKE_PROJECT_NAME}")
        ELSE(FAST_MATH)
            MESSAGE(STATUS "NOT Turning on -ffast-math for ${CMAKE_PROJECT_NAME}")
        ENDIF(FAST_MATH)
    ELSE(CMAKE_COMPILER_IS_GNUCXX)
        MESSAGE(STATUS "Fast math NOT AVAILABLE - Not a GNU compiler")
    ENDIF(CMAKE_COMPILER_IS_GNUCXX)
ENDMACRO(OPTION_FAST_MATH)

######################################################################################
# Turn on GProf based profiling
######################################################################################

MACRO(OPTION_GPROF DEFAULT)
    IF(CMAKE_COMPILER_IS_GNUCXX)
        OPTION(ENABLE_GPROF "Compile using -g -pg for gprof output" ${DEFAULT})
        IF(ENABLE_GPROF)
            MESSAGE(STATUS "Using gprof output for ${CMAKE_PROJECT_NAME}")
            SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -pg")
            SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -g -pg")
            SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -g -pg")
        ELSE(ENABLE_GPROF)
            MESSAGE(STATUS "NOT using gprof output for ${CMAKE_PROJECT_NAME}")
        ENDIF(ENABLE_GPROF)
    ELSE(CMAKE_COMPILER_IS_GNUCXX)
        MESSAGE(STATUS "gprof generation NOT AVAILABLE - Not a GNU compiler")
    ENDIF(CMAKE_COMPILER_IS_GNUCXX)
ENDMACRO(OPTION_GPROF)

######################################################################################
# Turn on "extra" compiler warnings (SPAMMY WITH BOOST)
######################################################################################

MACRO(OPTION_EXTRA_COMPILER_WARNINGS DEFAULT)
    IF(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        OPTION(EXTRA_COMPILER_WARNINGS "Turn on -Wextra for ${CMAKE_CXX_COMPILER_ID}" ${DEFAULT})
        IF(EXTRA_COMPILER_WARNINGS)
            MESSAGE(STATUS "Turning on extra c/c++ warnings for ${CMAKE_PROJECT_NAME}")
            SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra")
        ELSE(EXTRA_COMPILER_WARNINGS)
            MESSAGE(STATUS "NOT turning on extra c/c++ warnings for ${CMAKE_PROJECT_NAME}")
        ENDIF(EXTRA_COMPILER_WARNINGS)
    ELSE(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        MESSAGE(STATUS "Extra compiler warnings NOT AVAILABLE - Not a GNU/Clang compiler")
    ENDIF(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
ENDMACRO(OPTION_EXTRA_COMPILER_WARNINGS )

######################################################################################
# Turn on effective C++ compiler warnings
######################################################################################

MACRO(OPTION_EFFCXX_COMPILER_WARNINGS DEFAULT)
    IF(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        OPTION(EFFCXX_COMPILER_WARNINGS "Turn on -Weffc++ (effective c++ warnings) for ${CMAKE_CXX_COMPILER_ID}" ${DEFAULT})
        IF(EFFCXX_COMPILER_WARNINGS)
            MESSAGE(STATUS "Turning on Effective c++ warnings for ${CMAKE_PROJECT_NAME}")
            SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Weffc++")
        ELSE(EFFCXX_COMPILER_WARNINGS)
            MESSAGE(STATUS "NOT turning on Effective c++ warnings for ${CMAKE_PROJECT_NAME}")
        ENDIF(EFFCXX_COMPILER_WARNINGS)
    ELSE(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        MESSAGE(STATUS "Effective C++ compiler warnings NOT AVAILABLE - Not a GNU/Clang compiler")
    ENDIF(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
ENDMACRO(OPTION_EFFCXX_COMPILER_WARNINGS)

######################################################################################
# Turn on no depreciated
######################################################################################

MACRO(OPTION_NODEPREC_COMPILER_WARNINGS DEFAULT)
    IF(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        OPTION(NODEPREC_COMPILER_WARNINGS "Turn on -Wno-deprecated -Wno-deprecated-declarations for ${CMAKE_CXX_COMPILER_ID}" ${DEFAULT})
        IF(NODEPREC_COMPILER_WARNINGS)
            MESSAGE(STATUS "Turning on no depreciated function c++ warnings for ${CMAKE_PROJECT_NAME}")
            SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated -Wno-deprecated-declarations")
        ELSE(NODEPREC_COMPILER_WARNINGS)
            MESSAGE(STATUS "NOT Turning on no depreciated function c++ warnings for ${CMAKE_PROJECT_NAME}")
        ENDIF(NODEPREC_COMPILER_WARNINGS)
    ELSE(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        MESSAGE(STATUS "Effective C++ compiler warnings NOT AVAILABLE - Not a GNU/Clang compiler")
    ENDIF(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
ENDMACRO(OPTION_NODEPREC_COMPILER_WARNINGS)

######################################################################################
# Return type compiler warnings
######################################################################################

MACRO(OPTION_RETURN_TYPE_COMPILER_WARNINGS DEFAULT)
    IF(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        OPTION(RETURN_TYPE_COMPILER_WARNINGS "Turn on -Wreturn-type for ${CMAKE_CXX_COMPILER_ID}" ${DEFAULT})
        IF(RETURN_TYPE_COMPILER_WARNINGS)
            SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wreturn-type")
            MESSAGE(STATUS "Turning on return type warnings for ${CMAKE_PROJECT_NAME}")
        ELSE(RETURN_TYPE_COMPILER_WARNINGS)
            MESSAGE(STATUS "NOT turning on return type warnings for ${CMAKE_PROJECT_NAME}")
        ENDIF(RETURN_TYPE_COMPILER_WARNINGS)
    ELSE(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        MESSAGE(STATUS "Return type warnings NOT AVAILABLE - Not a GNU/Clang compiler")
    ENDIF(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
ENDMACRO(OPTION_RETURN_TYPE_COMPILER_WARNINGS)

######################################################################################
# Look for accelerate, if found, add proper includes
######################################################################################

MACRO(OPTION_ACCELERATE_FRAMEWORK DEFAULT)
    IF(APPLE)
        OPTION(ACCELERATE_FRAMEWORK "Use Accelerate Framework for Math (Adds -D_ACCELERATE_ for compiling and Accelerate Framework linking)" ${DEFAULT})
        IF(ACCELERATE_FRAMEWORK)
            FIND_LIBRARY(ACCELERATE_LIBRARY Accelerate REQUIRED PATHS "/System/Library/Frameworks/Accelerate.framework/" )
            SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_ACCELERATE_")
            MESSAGE(STATUS "Turning on Accelerate Framework for ${CMAKE_PROJECT_NAME}")
        ELSE(ACCELERATE_FRAMEWORK)
            MESSAGE(STATUS "NOT turning on Accelerate Framework for ${CMAKE_PROJECT_NAME}")
        ENDIF(ACCELERATE_FRAMEWORK)
    ELSE(APPLE)
        MESSAGE(STATUS "Accelerate Framework NOT AVAILABLE - Not compiling for OS X")
    ENDIF(APPLE)
ENDMACRO(OPTION_ACCELERATE_FRAMEWORK)

######################################################################################
# Enable C++11 for GNU/Clang Compilers
######################################################################################

MACRO(OPTION_ENABLE_CXX11 DEFAULT)
    OPTION(ENABLE_CXX11 "Enable C++11 for GNU/Clang compilers" ${DEFAULT})
    if(ENABLE_CXX11)
        message(STATUS "Enabling C++11 for GNU/Clang compilers")
        if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -lpetsc")
            if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
                set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
            endif()
        endif()
    else()
        message(STATUS "NOT enabling C++11 for GNU/Clang compilers")
    endif()
ENDMACRO()

######################################################################################
# Force 32-bit, regardless of the platform we're on
######################################################################################

MACRO(OPTION_FORCE_32_BIT DEFAULT)
    IF(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
        IF(CMAKE_COMPILER_IS_GNUCXX)
            OPTION(FORCE_32_BIT "Force compiler to use -m32 when compiling" ${DEFAULT})
            IF(FORCE_32_BIT)
                MESSAGE(STATUS "Forcing 32-bit on 64-bit platform (using -m32)")
                SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m32")
            ELSE(FORCE_32_BIT)
                MESSAGE(STATUS "Not forcing 32-bit on 64-bit platform")
            ENDIF(FORCE_32_BIT)
        ELSE(CMAKE_COMPILER_IS_GNUCXX)
            MESSAGE(STATUS "Force 32 bit NOT AVAILABLE - Not using gnu compiler")
        ENDIF(CMAKE_COMPILER_IS_GNUCXX)
    ELSE({CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
        MESSAGE(STATUS "Force 32 bit NOT AVAILABLE - Already on a 32 bit platform")
    ENDIF(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
ENDMACRO(OPTION_FORCE_32_BIT)


####################################################
## ---------------------------------------------- ##
## -                                            - ##
## -            Enable MPI Support              - ##
## -                                            - ##
## ---------------------------------------------- ##
####################################################

# Begin configuring MPI options
macro( ENABLE_MPI_SUPPORT )
    IF(USE_MPI)
        find_package(MPI REQUIRED)

        message(STATUS "Enabling MPI support")
        # Add the MPI-specific compiler and linker flags
        # Also, search for #includes in MPI's paths

        set(CMAKE_C_COMPILE_FLAGS "${CMAKE_C_COMPILE_FLAGS} ${MPI_C_COMPILE_FLAGS}")
        set(CMAKE_C_LINK_FLAGS "${CMAKE_C_LINK_FLAGS} ${MPI_C_LINK_FLAGS}")
        include_directories(${MPI_C_INCLUDE_PATH})

        set(CMAKE_CXX_COMPILE_FLAGS "${CMAKE_CXX_COMPILE_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
        set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} ${MPI_CXX_LINK_FLAGS}")
        include_directories(${MPI_CXX_INCLUDE_PATH})
    else()
      message(STATUS "Not enabling MPI support")
    ENDIF()
endmacro(ENABLE_MPI_SUPPORT)

# Done configuring MPI Options


####################################################
## ---------------------------------------------- ##
## -                                            - ##
## -            Enable OpenMP Support           - ##
## -                                            - ##
## ---------------------------------------------- ##
####################################################

# Begin configuring OpenMP options
macro(enable_openmp_support DEFAULT)

    OPTION(USE_OPENMP "Enable OpenMP multithreading" ${DEFAULT})
    IF(USE_OPENMP)
    
    	if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    		message(WARNING "\n\tClang on Darwin (as of OS X Mavericks) does not have OpenMP Support\n")
    	endif()
    	#else()
	        find_package(OpenMP)
        
        	if(OpenMP_FOUND)
	    	    # Add the OpenMP-specific compiler and linker flags
    	    	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        		set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        		set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
	        endif()
		#endif()
		
    ENDIF(USE_OPENMP)
endmacro(enable_openmp_support)
# Done configuring OpenMP Options

####################################################
## ---------------------------------------------- ##
## -                                            - ##
## -            Configure BOOST                    - ##
## -                                            - ##
## ---------------------------------------------- ##
####################################################

macro(configure_boost VERSION)

    #set(Boost_USE_STATIC_LIBS        ON)
    #set(Boost_USE_MULTITHREADED      ON)
    #set(Boost_USE_STATIC_RUNTIME    OFF)

    # - BoostConfig
    # looks for local installation using env variable WCCPP_EXTERNALS as hints
    # if no custom install in $BOOST_ROOT, looks for installation in standard directtory /usr
    # This module defines
    #      Boost_INCLUDE_DIR: where to find include files
    #      Boost_LIBRARY NAMES : which library to look for
    #      Boost_LIBRARY_DIRS: where to find the libraries
    #      Boost_LIBRARIES: where to find the libraries
    # This module defines also implicitly
    #      Boost_FOUND and set value to 1
    #      Boost_DIR and set value to path where this file is located
    #      Boost_CONFIG and set value to filename

    # define a list of library names to look for (when one is found found, it is a success!!!)
    set(Boost_NAMES ${Boost_NAMES}
                    boost_atomic
                    boost_chrono
                    boost_context
                    boost_date_time
                    boost_exception
                    boost_filesystem
                    boost_graph
                    boost_graph_parallel
                    boost_iostreams
                    boost_locale
                    boost_math
                    boost_math_c99f
                    boost_math_c99l
                    boost_math_c99
                    boost_math_tr1f
                    boost_math_tr1l
                    boost_math_tr1
                    boost_prg_exec_monitor
                    boost_mpi
                    boost_python
                    boost_regex
                    boost_program_options
                    boost_random
                    boost_serialization
                    boost_signals
                    boost_system
                    boost_test
                    boost_thread
                    boost_timer
                    boost_unit_test_framework
                    boost_wave
                    #boost_wserialization
            )

    #set(Boost_NAMES ####### ALL BOOST LIBRARIES #######
                #atomic chrono context date_time exception filesystem graph graph_parallel iostreams locale math math_c99f
                #math_c99l math_c99 math_tr1f math_tr1l math_tr1 prg_exec_monitor mpi python regex program_options random
                #serialization signals system test thread timer unit_test_framework wave wserialization

                ####### BOOST LIBRARIES ON CLUSTER.NE.TAMU.EDU #######
                #atomic chrono context date_time filesystem
                #graph locale math_c99f math_c99l math_c99
                #math_tr1f math_tr1l math_tr1 prg_exec_monitor
                #program_options random regex serialization
                #signals system thread timer unit_test_framework
                #wave wserialization
    #)

    # check if there is a local installation
    #   message(STATUS "inside BoostConfig:  Boost_FOUND is ${Boost_FOUND}, no system installation found")
    #message(STATUS "looking for local installation in $ENV{BOOST_ROOT}...")
    find_path(Boost_INCLUDE_DIR boost/graph/adjacency_list_io.hpp PATHS $ENV{BOOST_ROOT}/include )
    if(IS_DIRECTORY ${Boost_INCLUDE_DIR})
        #message(STATUS "inside BoostConfig: local installation in ${Boost_INCLUDE_DIR}")
          set(Boost_FOUND TRUE)
    else(IS_DIRECTORY ${Boost_INCLUDE_DIR})
        message(STATUS "inside BoostConfig:  no local installation in $ENV{BOOST_ROOT}")
          set(Boost_FOUND FALSE)
    endif(IS_DIRECTORY ${Boost_INCLUDE_DIR})

    # find the path where the library is
    foreach(_boost_lib ${Boost_NAMES})
        find_library(_boost_library NAMES ${_boost_lib} PATHS $ENV{BOOST_ROOT}/lib)
        if(EXISTS ${_boost_library})
            #get_filename_component(Boost_LIBRARY_DIRS ${Boost_LIBRARIES} PATH)
            #message(STATUS "inside BoostConfig: ${_boost_lib} library found in ${Boost_LIBRARY_DIRS}")
            list(APPEND Boost_LIBRARIES ${_boost_library})
        else()
            message(STATUS "inside BoostConfig: no ${_boost_lib} library found in local installation in $ENV{BOOST_ROOT}")
            set(Boost_FOUND FALSE)
        endif()
        unset(_boost_library CACHE)
    endforeach()


    if(Boost_FOUND)
          message(STATUS "found local installation in $ENV{BOOST_ROOT} for Boost")
    endif(Boost_FOUND)


    if(NOT Boost_FOUND)
      message(STATUS "no local installation found for Boost, looking for system installation using /usr/share/cmake/FindBoost.cmake...")
      # check in system installation /usr/share/cmake/FindBoost.cmake
      FIND_PACKAGE(Boost ${VERSION})
      if( (NOT EXISTS ${Boost_LIBRARY_DIRS}) OR (NOT EXISTS ${Boost_INCLUDE_DIR}) )
        message(STATUS "inside BoostConfig: Boost system installation not found")
        set(Boost_FOUND FALSE)
      endif( (NOT EXISTS ${Boost_LIBRARY_DIRS}) OR (NOT EXISTS ${Boost_INCLUDE_DIR}) )
    endif(NOT Boost_FOUND)

    #message(STATUS "inside BoostConfig:  Boost_FOUND is ${Boost_FOUND}")
    #message(STATUS "inside BoostConfig:  Boost_INCLUDE_DIR is ${Boost_INCLUDE_DIR}")
    #message(STATUS "inside BoostConfig:  Boost_LIBRARY_DIRS is ${Boost_LIBRARY_DIRS}")
    #message(STATUS "inside BoostConfig:  Boost_LIBRARIES is ${Boost_LIBRARIES}")
    #message(STATUS "inside BoostConfig:  Boost_NAMES is ${Boost_NAMES}")
endmacro()
