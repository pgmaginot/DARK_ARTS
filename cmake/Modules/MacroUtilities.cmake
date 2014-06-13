#
# Modified/Copied from example in Geant4 version 9.6.2:
#
# http://geant4.web.cern.ch/geant4/
#
# MacroUtilities - useful macros and functions for generic tasks
#
# CMake Extensions
# ----------------
# macro set_ifnot(<var> <value>)
#       If variable var is not set, set its value to that provided
#
# function enum_option(<option>
#                      VALUES <value1> ... <valueN>
#                      TYPE   <valuetype>
#                      DOC    <docstring>
#                      [DEFAULT <elem>]
#                      [CASE_INSENSITIVE])
#          Declare a cache variable <option> that can only take values
#          listed in VALUES. TYPE may be FILEPATH, PATH or STRING.
#          <docstring> should describe that option, and will appear in
#          the interactive CMake interfaces. If DEFAULT is provided,
#          <elem> will be taken as the zero-indexed element in VALUES
#          to which the value of <option> should default to if not
#          provided. Otherwise, the default is taken as the first
#          entry in VALUES. If CASE_INSENSITIVE is present, then
#          checks of the value of <option> against the allowed values
#          will ignore the case when performing string comparison.
#
#
# General
# --------------
# function add_feature(<NAME> <DOCSTRING>)
#          Add a  feature, whose activation is specified by the
#          existence of the variable <NAME>, to the list of enabled/disabled
#          features, plus a docstring describing the feature
#
# function print_enabled_features()
#          Print enabled  features plus their docstrings.
#
#

# - Include guard
if(__macroutilities_isloaded)
  return()
endif()
set(__macroutilities_isloaded YES)

include(MacroDefineModule)
include(MacroLibraryTargets)

#-----------------------------------------------------------------------
# CMAKE EXTENSIONS
#-----------------------------------------------------------------------
# macro set_ifnot(<var> <value>)
#       If variable var is not set, set its value to that provided
#
macro(set_ifnot _var _value)
  if(NOT ${_var})
    set(${_var} ${_value})
  endif()
endmacro()

#-----------------------------------------------------------------------
# function enum_option(<option>
#                      VALUES <value1> ... <valueN>
#                      TYPE   <valuetype>
#                      DOC    <docstring>
#                      [DEFAULT <elem>]
#                      [CASE_INSENSITIVE])
#          Declare a cache variable <option> that can only take values
#          listed in VALUES. TYPE may be FILEPATH, PATH or STRING.
#          <docstring> should describe that option, and will appear in
#          the interactive CMake interfaces. If DEFAULT is provided,
#          <elem> will be taken as the zero-indexed element in VALUES
#          to which the value of <option> should default to if not
#          provided. Otherwise, the default is taken as the first
#          entry in VALUES. If CASE_INSENSITIVE is present, then
#          checks of the value of <option> against the allowed values
#          will ignore the case when performing string comparison.
#
function(enum_option _var)
  set(options CASE_INSENSITIVE)
  set(oneValueArgs DOC TYPE DEFAULT)
  set(multiValueArgs VALUES)
  cmake_parse_arguments(_ENUMOP "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  # - Validation as needed arguments
  if(NOT _ENUMOP_VALUES)
    message(FATAL_ERROR "enum_option must be called with non-empty VALUES\n(Called for enum_option '${_var}')")
  endif()

  # - Set argument defaults as needed
  if(_ENUMOP_CASE_INSENSITIVE)
    set(_ci_values )
    foreach(_elem ${_ENUMOP_VALUES})
      string(TOLOWER "${_elem}" _ci_elem)
      list(APPEND _ci_values "${_ci_elem}")
    endforeach()
    set(_ENUMOP_VALUES ${_ci_values})
  endif()

  set_ifnot(_ENUMOP_TYPE STRING)
  set_ifnot(_ENUMOP_DEFAULT 0)
  list(GET _ENUMOP_VALUES ${_ENUMOP_DEFAULT} _default)

  if(NOT DEFINED ${_var})
    set(${_var} ${_default} CACHE ${_ENUMOP_TYPE} "${_ENUMOP_DOC} (${_ENUMOP_VALUES})")
  else()
    set(_var_tmp ${${_var}})
    if(_ENUMOP_CASE_INSENSITIVE)
      string(TOLOWER ${_var_tmp} _var_tmp)
    endif()

    list(FIND _ENUMOP_VALUES ${_var_tmp} _elem)
    if(_elem LESS 0)
      message(FATAL_ERROR "Value '${${_var}}' for variable ${_var} is not allowed\nIt must be selected from the set: ${_ENUMOP_VALUES} (DEFAULT: ${_default})\n")
    else()
      # - convert to lowercase
      if(_ENUMOP_CASE_INSENSITIVE)
        set(${_var} ${_var_tmp} CACHE ${_ENUMOP_TYPE} "${_ENUMOP_DOC} (${_ENUMOP_VALUES})" FORCE)
      endif()
    endif()
  endif()
endfunction()

#-----------------------------------------------------------------------
# GENERAL
#-----------------------------------------------------------------------
# function add_feature(<NAME> <DOCSTRING>)
#          Add a  feature, whose activation is specified by the
#          existence of the variable <NAME>, to the list of enabled/disabled
#          features, plus a docstring describing the feature
#
function(ADD_FEATURE _var _description)
  if(${_var})
    set_property(GLOBAL APPEND PROPERTY ENABLED_FEATURES ${_var})
  else()
    set_property(GLOBAL APPEND PROPERTY DISABLED_FEATURES ${_var})
  endif()

  set_property(GLOBAL PROPERTY ${_var}_DESCRIPTION "${_description}")
endfunction()

#-----------------------------------------------------------------------
# function print_enabled_features()
#          Print enabled  features plus their docstrings.
#
function(print_enabled_features)
  set(_currentFeatureText "The following  features are enabled:")
  get_property(_enabledFeatures GLOBAL PROPERTY ENABLED_FEATURES)

  foreach(_feature ${_enabledFeatures})
    set(_currentFeatureText "${_currentFeatureText}\n\t${_feature}")

    get_property(_desc GLOBAL PROPERTY ${_feature}_DESCRIPTION)

    if(_desc)
      set(_currentFeatureText "${_currentFeatureText}: ${_desc}")
      set(_desc NOTFOUND)
    endif(_desc)
  endforeach(_feature)

  message(STATUS "${_currentFeatureText}\n")
endfunction()


#-----------------------------------------------------------------------
# function print_disabled_features()
#          Print disabled features plus their docstrings.
#
function(print_disabled_features)
  set(_currentFeatureText "The following ${PROJECT_NAME} features are NOT enabled (may be implicitly defined/derived):")
  get_property(_disabledFeatures GLOBAL PROPERTY DISABLED_FEATURES)

  foreach(_feature ${_disabledFeatures})
    set(_currentFeatureText "${_currentFeatureText}\n\t${_feature}")

    get_property(_desc GLOBAL PROPERTY ${_feature}_DESCRIPTION)

    if(_desc)
      set(_currentFeatureText "${_currentFeatureText}: ${_desc}")
      set(_desc NOTFOUND)
    endif(_desc)
  endforeach(_feature)

  message(STATUS "${_currentFeatureText}\n")
endfunction()

