#
# Modified/Copied from example in Geant4 version 9.6.2:
#
# http://geant4.web.cern.ch/geant4/
#
# - Set up backward compatible ${PROJECT_NAME} GNU make toolchain
#
# The GNU make based buildsystem for ${PROJECT_NAME} provides a toolchain for
# users building simple ${PROJECT_NAME} applications. The old style build and
# install of ${PROJECT_NAME} provides a customized set of non-standard install
# paths with use of the toolchain dependent on environment variables
# pointing to the install paths.
#
# This script processes information on the CMake install paths, system
# and compiler to determine the following variables for backward
# compatibility:
#
#  ${PROJECT_NAME}_SYSTEM      Old style system name, e.g. 'Linux', 'Darwin'
#                     or 'WIN32'
#
#  ${PROJECT_NAME}_COMPILER    Old system compiler id, e.g. 'g++', 'VC'.
#
#  ${PROJECT_NAME}INSTALL          Location of 'config' subdirectory which contains
#                     all the GNU make toolchain fragments
#
#  ${PROJECT_NAME}INCLUDE          Old style path to location of ${PROJECT_NAME} headers
#
#  ${PROJECT_NAME}LIB              Old style library directory path. Rather than
#                     containing the actual libraries, it is expected to
#                     contain subdirectories named
#                     ${PROJECT_NAME}_SYSTEM-${PROJECT_NAME}_COMPILER
#
# These variables are used in a CMake configuration file which is used
# to generate shell scripts (C and Bourne flavour) the user can source
# to set up their environment for use of the old toolchain.
# These replace the old 'env.(c)sh' scripts to allow users to work with
# the new CMake built libraries transparently if their application
# relies on the old style toolchain.
#
# The scripts are generated for both the build and install trees so that
# developers wishing to write test applications do not have to install
# their fresh build of ${PROJECT_NAME}.
#
# Compatibility with the library path style:
#
#  <prefix>/lib/${PROJECT_NAME}SYSTEM-${PROJECT_NAME}COMPILER
#
# is provided by installing a directory '${PROJECT_NAME}-<version>' in the
# <prefix>/lib directory and creating a symbolic link inside here
# pointing up one directory level.
# This will not work on Windows however, and here users are recommended
# to use Visual Studio directly, or to use CMake for application
# configuration.
#

#-----------------------------------------------------------------------
# - Functions and Macros to help configuration of shell scripts.
#-----------------------------------------------------------------------
# macro _tc_shell_setup(<shell>)
#       Set shell parameters such as program, family and common builtins
#       for supplied shell (e.g. 'bourne' or 'cshell'
#
macro(_tc_shell_setup SHELL_FAMILY)
  if(${SHELL_FAMILY} STREQUAL "bourne")
    set(${PROJECT_NAME}_TC_SHELL_PROGRAM "/bin/sh")
    set(${PROJECT_NAME}_TC_SHELL_FAMILY "Bourne shell")
    set(${PROJECT_NAME}_TC_UNSET_COMMAND "unset")
    set(${PROJECT_NAME}_TC_SHELL_EXTENSION ".sh")
  elseif(${SHELL_FAMILY} STREQUAL "cshell")
    set(${PROJECT_NAME}_TC_SHELL_PROGRAM "/bin/csh")
    set(${PROJECT_NAME}_TC_SHELL_FAMILY "C shell")
    set(${PROJECT_NAME}_TC_UNSET_COMMAND "unsetenv")
    set(${PROJECT_NAME}_TC_SHELL_EXTENSION ".csh")
  else()
    message(FATAL_ERROR "Unsupported shell '${SHELL_FAMILY}'")
  endif()
endmacro()

#-----------------------------------------------------------------------
# function _tc_selflocate(<output> <shell> <script> <variable name>)
#          Set output to string containing shell commands needed to
#          locate the directory in which script is located if the
#          script is sourced. This derived location is set as the
#          value of the shell variable name.
#
function(_tc_selflocate TEMPLATE_NAME SHELL_FAMILY SCRIPT_NAME LOCATION_VARIABLE)
  if(${SHELL_FAMILY} STREQUAL "bourne")
    set(${TEMPLATE_NAME}
      "# Self locate script when sourced
if [ -z \"\$BASH_VERSION\" ]; then
  # Not bash, so rely on sourcing from correct location
  if [ ! -f ${SCRIPT_NAME}${${PROJECT_NAME}_TC_SHELL_EXTENSION} ]; then
    echo 'ERROR: ${SCRIPT_NAME}${${PROJECT_NAME}_TC_SHELL_EXTENSION} could NOT self-locate ${PROJECT_NAME} installation'
    echo 'This is most likely because you are using ksh, zsh or similar'
    echo 'To fix this issue, cd to the directory containing this script'
    echo 'and source it in that directory.'
    return 1
  fi
  ${LOCATION_VARIABLE}=\$\(pwd\)
else
  sls_sourced_dir=\$\(dirname \${BASH_ARGV[0]}\)
  ${LOCATION_VARIABLE}=$\(cd \$sls_sourced_dir > /dev/null ; pwd\)
fi
      "
      PARENT_SCOPE
      )
    # For bourne shell, set the values of the guard variables
    set(TC_IF_SELFLOCATED "" PARENT_SCOPE)
    set(TC_ENDIF_SELFLOCATED "" PARENT_SCOPE)


  elseif(${SHELL_FAMILY} STREQUAL "cshell")
    set(${TEMPLATE_NAME}
      "# Self locate script when sourced
# If sourced interactively, we can use $_ as this should be
#
#   source path_to_script_dir/${SCRIPT_NAME}${${PROJECT_NAME}_TC_SHELL_EXTENSION}
#
unset sls_sourced_dir
unset ${LOCATION_VARIABLE}

set ARGS=($_)
if (\"$ARGS\" != \"\") then
  if (\"$ARGS[2]\" =~ */${SCRIPT_NAME}${${PROJECT_NAME}_TC_SHELL_EXTENSION}) then
    set sls_sourced_dir=\"`dirname \${ARGS[2]}`\"
  endif
endif

if (! \$?sls_sourced_dir) then
  # Oh great, we were sourced non-interactively. This means that $_
  # won't be set, so we need an external source of information on
  # where the script is located.
  # We obtain this in one of two ways:
  #   1) Current directory:
  #     cd script_dir ; source ${SCRIPT_NAME}${${PROJECT_NAME}_TC_SHELL_EXTENSION}
  #
  #   2) Supply the directory as an argument to the script:
  #     source script_dir/${SCRIPT_NAME}${${PROJECT_NAME}_TC_SHELL_EXTENSION} script_dir
  #
  if ( -e ${SCRIPT_NAME}${${PROJECT_NAME}_TC_SHELL_EXTENSION} ) then
    set sls_sourced_dir=\"`pwd`\"
  else if ( \"\$1\" != \"\" )  then
    if ( -e \${1}/${SCRIPT_NAME}${${PROJECT_NAME}_TC_SHELL_EXTENSION} ) then
      set sls_sourced_dir=\${1}
    else
      echo \"ERROR \${1} does not contain a ${PROJECT_NAME} installation\"
    endif
  endif
endif

if (! \$?sls_sourced_dir) then
  echo \"ERROR: ${SCRIPT_NAME}${${PROJECT_NAME}_TC_SHELL_EXTENSION} could NOT self-locate ${PROJECT_NAME} installation\"
  echo \"because it was sourced (i.e. embedded) in another script.\"
  echo \"This is due to limitations of (t)csh but can be worked around by providing\"
  echo \"the directory where ${SCRIPT_NAME}${${PROJECT_NAME}_TC_SHELL_EXTENSION} is located\"
  echo \"to it, either via cd-ing to the directory before sourcing:\"
  echo \"  cd where_script_is ; source ${SCRIPT_NAME}${${PROJECT_NAME}_TC_SHELL_EXTENSION}\"
  echo \"or by supplying the directory as an argument to the script:\"
  echo \"  source where_script_is/${SCRIPT_NAME}${${PROJECT_NAME}_TC_SHELL_EXTENSION} where_script_is\"
  echo \" \"
  exit 1
endif

set ${LOCATION_VARIABLE}=\"`cd \${sls_sourced_dir} > /dev/null ; pwd`\"
"
      PARENT_SCOPE
      )

    # For C-shell, set the values of the guard variables
    set(TC_IF_SELFLOCATED "" PARENT_SCOPE)
   set(TC_ENDIF_SELFLOCATED "" PARENT_SCOPE)
  endif()
endfunction()

#-----------------------------------------------------------------------
# function _tc_setenv_command(<output> <shell> <name> <value>)
#          Set output to a string whose value is the shell command to
#          set an environment variable with name and value
#
function(_tc_setenv_command TEMPLATE_NAME SHELL_FAMILY VARIABLE_NAME VARIABLE_VALUE)
  if(${SHELL_FAMILY} STREQUAL "bourne")
    set(${TEMPLATE_NAME}
      "export ${VARIABLE_NAME}=${VARIABLE_VALUE}"
      PARENT_SCOPE
      )
  elseif(${SHELL_FAMILY} STREQUAL "cshell")
    set(${TEMPLATE_NAME}
      "setenv ${VARIABLE_NAME} ${VARIABLE_VALUE}"
      PARENT_SCOPE
      )
  endif()
endfunction()

#-----------------------------------------------------------------------
# function _tc_setenv_ifnotset_command(<output> <shell> <name> <value>)
#          Set output to a string whose value is the shell command to
#          set an environment variable with name and value if the
#          variable is not already set
#
function(_tc_setenv_ifnotset_command TEMPLATE_NAME SHELL_FAMILY VARIABLE_NAME VARIABLE_VALUE)
  # -- bourne
  if(${SHELL_FAMILY} STREQUAL "bourne")
    # Have to make this section verbatim to get correct formatting
    set(${TEMPLATE_NAME}
      "
if test \"x\$${VARIABLE_NAME}\" = \"x\" ; then
  export ${VARIABLE_NAME}=${VARIABLE_VALUE}
fi
"
      PARENT_SCOPE
      )
  # -- cshell
  elseif(${SHELL_FAMILY} STREQUAL "cshell")
    # Again, verbatim to get correct formatting...
    set(${TEMPLATE_NAME}
      "
if ( ! \${?${VARIABLE_NAME}} ) then
  setenv ${VARIABLE_NAME} ${VARIABLE_VALUE}
endif
"
       PARENT_SCOPE
       )
  endif()
endfunction()

#-----------------------------------------------------------------------
# function _tc_prepend_path(<output> <shell> <name> <value>)
#          Set output to a string whose value is the shell command to
#          prepend supplied value to the path style environment variable
#          name (e.g. 'PATH')
#
    #if [[ ! \":$${PATH_VARIABLE}:\" == *\":${APPEND_VARIABLE}:\"* ]]; then
    #    export ${PATH_VARIABLE}=${APPEND_VARIABLE}:\${${PATH_VARIABLE}}
    #fi
#
function(_tc_prepend_path TEMPLATE_NAME SHELL_FAMILY PATH_VARIABLE
  APPEND_VARIABLE)
  # -- bourne block
  if(${SHELL_FAMILY} STREQUAL "bourne")
    # We have to make this section verbatim
    set(${TEMPLATE_NAME}
    "
if test \"x\$${PATH_VARIABLE}\" = \"x\" ; then
  export ${PATH_VARIABLE}=${APPEND_VARIABLE}
else
  export ${PATH_VARIABLE}=${APPEND_VARIABLE}:\${${PATH_VARIABLE}}
fi
        "
    PARENT_SCOPE
    )
  # -- cshell block
  elseif(${SHELL_FAMILY} STREQUAL "cshell")
    # Again, this is verbatim so final output is formatted correctly
    set(${TEMPLATE_NAME}
      "
if ( ! \${?${PATH_VARIABLE}} ) then
  setenv ${PATH_VARIABLE} ${APPEND_VARIABLE}
else
  setenv ${PATH_VARIABLE} ${APPEND_VARIABLE}:\${${PATH_VARIABLE}}
endif
      "
      PARENT_SCOPE
      )
  endif()
endfunction()

#-----------------------------------------------------------------------
# MACRO(_tc_configure_tc_variables)
# Macro to perform the actual setting of the low level toolchain variables
# which need to be set in the final shell files.
# We do this in a separate macro so that we can wrap it in different ways for
# the install and build trees.
#
macro(_tc_configure_tc_variables SHELL_FAMILY SCRIPT_NAME)
  # - Set up the requested shell
  _tc_shell_setup(${SHELL_FAMILY})

  # - Locate self
  _tc_selflocate(${PROJECT_NAME}_TC_LOCATE_SELF_COMMAND ${SHELL_FAMILY} ${SCRIPT_NAME} ${PROJECT_NAME}make_root)


  # - Standard Setup and Paths
  _tc_setenv_command(${PROJECT_NAME}_TC_${PROJECT_NAME}SYSTEM ${SHELL_FAMILY} ${PROJECT_SYSTEM ${${PROJECT_NAME}SYSTEM})
  _tc_setenv_command(${PROJECT_NAME}_TC_${PROJECT_NAME}INSTALL ${SHELL_FAMILY} ${PROJECT_NAME}INSTALL ${${PROJECT_NAME}INSTALL})
  _tc_setenv_command(${PROJECT_NAME}_TC_${PROJECT_NAME}INCLUDE ${SHELL_FAMILY} ${PROJECT_NAME}INCLUDE ${${PROJECT_NAME}INCLUDE})
  _tc_setenv_command(${PROJECT_NAME}_TC_${PROJECT_NAME}LIB ${SHELL_FAMILY} ${PROJECT_NAME}LIB ${${PROJECT_NAME}LIB})

  if(${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
    _tc_prepend_path(${PROJECT_NAME}_TC_${PROJECT_NAME}LIB_PATH_SETUP ${SHELL_FAMILY} DYLD_LIBRARY_PATH ${${PROJECT_NAME}LIB_DIR})
  else()
    _tc_prepend_path(${PROJECT_NAME}_TC_${PROJECT_NAME}LIB_PATH_SETUP ${SHELL_FAMILY} LD_LIBRARY_PATH ${${PROJECT_NAME}LIB_DIR})
  endif()

  _tc_setenv_ifnotset_command(${PROJECT_NAME}_TC_${PROJECT_NAME}WORKDIR_SETUP ${SHELL_FAMILY} ${PROJECT_NAME}WORKDIR ${${PROJECT_NAME}WORKDIR_DEFAULT})
  _tc_prepend_path(${PROJECT_NAME}_TC_${PROJECT_NAME}WORKDIR_PATH_SETUP ${SHELL_FAMILY} PATH
  \${${PROJECT_NAME}WORKDIR}/bin/\${${PROJECT_NAME}SYSTEM})

  # - ${PROJECT_NAME} Library build setup
  # We prefer shared libs if these are built, otherwise fall back to static
  # On Win32, we also want DLLs?
  if(BUILD_SHARED_LIBS)
    _tc_setenv_command(${PROJECT_NAME}_TC_${PROJECT_NAME}LIB_BUILD_SHARED ${SHELL_FAMILY} ${PROJECT_NAME}LIB_BUILD_SHARED 1)
    if(WIN32)
      _tc_setenv_command(${PROJECT_NAME}_TC_${PROJECT_NAME}LIB_USE_DLL ${SHELL_FAMILY} ${PROJECT_NAME}LIB_USE_DLL 1)
    endif()
  else()
    _tc_setenv_command(${PROJECT_NAME}_TC_${PROJECT_NAME}LIB_BUILD_STATIC ${SHELL_FAMILY} ${PROJECT_NAME}LIB_BUILD_STATIC 1)
  endif()


  # - Resource file paths
  set(${PROJECT_NAME}_TC_DATASETS )
  foreach(_ds ${${PROJECT_NAME}_EXPORTED_DATASETS})
    _tc_setenv_command(_dssetenvcmd ${SHELL_FAMILY} ${${_ds}_ENVVAR} ${${_ds}_PATH})
    set(${PROJECT_NAME}_TC_DATASETS "${${PROJECT_NAME}_TC_DATASETS}${_dssetenvcmd}\n")
  endforeach()



endmacro()



#-----------------------------------------------------------------------
# Implementation section
#-----------------------------------------------------------------------
# Configure shell scripts for BUILD TREE
# This means we have to point to libraries in the build tree, but
# includes and resource files will be in the source tree
# This script never needs to be relocatable, so we don't need to use the
# self location functionality.
# N.B. IT WILL NOT WORK when building with VS/Xcode or any multiconfig
# buildtool because we cannot reconcile the output paths these use with
# those expected by the old toolchain...
#
set(${PROJECT_NAME}SYSTEM  "${${PROJECT_NAME}_SYSTEM}-${${PROJECT_NAME}_COMPILER}")
set(${PROJECT_NAME}INSTALL ${PROJECT_SOURCE_DIR})
set(${PROJECT_NAME}INCLUDE ${PROJECT_SOURCE_DIR}/this_is_a_deliberate_dummy_path)
set(${PROJECT_NAME}LIB ${PROJECT_BINARY_DIR})
set(${PROJECT_NAME}LIB_DIR ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
set(${PROJECT_NAME}WORKDIR_DEFAULT "\$HOME/${PROJECT_NAME}_workdir")



#-----------------------------------------------------------------------
# Configure shell scripts for INSTALL TREE
# This means we have to point things to their final location when
# installed. These paths are all determined by the CMAKE_INSTALL_FULL
# directories and others.
# If we are relocatable, then the structure we will have is
# +- CMAKE_INSTALL_PREFIX
#    +- LIBDIR/${PROJECT_NAME}-VERSION (${PROJECT_NAME}LIB)
#    +- INCLUDEDIR/${PROJECT_NAME}     (${PROJECT_NAME}INCLUDE)
#    +- DATAROOTDIR/${PROJECT_NAME}-VERSION/
#       +- ${PROJECT_NAME}make              (${PROJECT_NAME}INSTALL!)
#          +- ${PROJECT_NAME}make.(c)sh
#          +- config/

# - Construct universal backward compatible INSTALL TREE PATHS.
set(${PROJECT_NAME}SYSTEM  "${${PROJECT_NAME}_SYSTEM}-${${PROJECT_NAME}_COMPILER}")
set(${PROJECT_NAME}INSTALL "\"\$${PROJECT_NAME}make_root\"")

# - Now need relative paths between '${PROJECT_NAME}INSTALL' and include/lib dirs
# - Include dir
file(RELATIVE_PATH
  ${PROJECT_NAME}MAKE_TO_INCLUDEDIR
  ${CMAKE_INSTALL_FULL_DATAROOTDIR}/${PROJECT_NAME}-${${PROJECT_NAME}_VERSION}/${PROJECT_NAME}make
  ${CMAKE_INSTALL_FULL_INCLUDEDIR}/${PROJECT_NAME}
  )
set(${PROJECT_NAME}INCLUDE "\"`cd \$${PROJECT_NAME}make_root/${${PROJECT_NAME}MAKE_TO_INCLUDEDIR} > /dev/null \; pwd`\"")

# - Lib dir
file(RELATIVE_PATH
  ${PROJECT_NAME}MAKE_TO_LIBDIR
  ${CMAKE_INSTALL_FULL_DATAROOTDIR}/${PROJECT_NAME}-${${PROJECT_NAME}_VERSION}/${PROJECT_NAME}make
  ${CMAKE_INSTALL_FULL_LIBDIR}
  )
set(${PROJECT_NAME}LIB "\"`cd \$${PROJECT_NAME}make_root/${${PROJECT_NAME}MAKE_TO_LIBDIR}/${PROJECT_NAME}-${${PROJECT_NAME}_VERSION} > /dev/null \; pwd`\"")
set(${PROJECT_NAME}LIB_DIR "\"`cd \$${PROJECT_NAME}make_root/${${PROJECT_NAME}MAKE_TO_LIBDIR} > /dev/null \; pwd`\"")

set(${PROJECT_NAME}WORKDIR_DEFAULT "\$HOME/${PROJECT_NAME}_workdir")




#-----------------------------------------------------------------------
# TEMPORARY
# Configure environment setup script for install of ${PROJECT_NAME}
# Temporarily here to keep all shell setup in one place.
# Later, should be refactored into its own module, with module containing
# all the shell tools above.
#
# - Script base name (without extension
set(_scriptbasename ${PROJECT_NAME})

# - Relative path between bindir (where script is) and library directory
file(RELATIVE_PATH
  ${PROJECT_NAME}ENV_BINDIR_TO_LIBDIR
  ${CMAKE_INSTALL_FULL_BINDIR}
  ${CMAKE_INSTALL_FULL_LIBDIR}
  )


# - Configure for each shell
foreach(_shell bourne;cshell)
  # Setup the shell
  _tc_shell_setup(${_shell})

  # Set script full name
  set(_scriptfullname ${_scriptbasename}${${PROJECT_NAME}_TC_SHELL_EXTENSION})

  # Set locate self command
  _tc_selflocate(ENV_SELFLOCATE_COMMAND
    ${_shell}
    ${_scriptbasename}
    ${PROJECT_NAME}_envbindir
    )

  # Set path, which should be where the script itself is installed
  _tc_prepend_path(ENV_BINPATH_SETUP
    ${_shell}
    PATH
    "\"\$${PROJECT_NAME}_envbindir\""
    )

  # Set library path, based on relative paths between bindir and libdir
  if(${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
    set(_libpathname DYLD_LIBRARY_PATH)
  else()
    set(_libpathname LD_LIBRARY_PATH)
  endif()

  _tc_prepend_path(ENV_LIBPATH_SETUP
    ${_shell}
    ${_libpathname}
    #"\"$${PROJECT_NAME}_envbindir/${${PROJECT_NAME}ENV_BINDIR_TO_LIBDIR}\""
    "\"`cd $${PROJECT_NAME}_envbindir/${${PROJECT_NAME}ENV_BINDIR_TO_LIBDIR} > /dev/null ; pwd`\""

    )

  # - Set data paths
  set(${PROJECT_NAME}_ENV_DATASETS )
  foreach(_ds ${${PROJECT_NAME}_EXPORTED_DATASETS})
    _tc_setenv_command(_dssetenvcmd ${_shell} ${${_ds}_ENVVAR} ${${_ds}_PATH})
    set(${PROJECT_NAME}_ENV_DATASETS "${${PROJECT_NAME}_ENV_DATASETS}${_dssetenvcmd}\n")
  endforeach()

  # Configure the file
  configure_file(
    ${CMAKE_SOURCE_DIR}/cmake/Templates/env-skeleton.in
    ${PROJECT_BINARY_DIR}/InstallTreeFiles/${_scriptfullname}
    @ONLY
    )

  # Install it to the required location
  install(FILES
    ${PROJECT_BINARY_DIR}/InstallTreeFiles/${_scriptfullname}
    DESTINATION ${CMAKE_INSTALL_BINDIR}
    PERMISSIONS
    OWNER_READ OWNER_WRITE OWNER_EXECUTE
    GROUP_READ GROUP_EXECUTE
    WORLD_READ WORLD_EXECUTE
    COMPONENT Runtime
    )
endforeach()

