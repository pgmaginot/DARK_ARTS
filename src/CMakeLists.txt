#
# NOTE:
#   All subdirectories use this format to add the files inside the
#   directory, the macro "AddSubdirectoryFiles" is defined in
#   cmake/Modules/GenericCMakeFunctions.cmake
#
#
# AddSubdirectoryFiles(DEFINITIONS   # this is used internally by function
#     DEFINITIONS # this is what you want to propagate
#           # if the folder is embedded more than one folder
#           # inside this folder (i.e. this folder is source/
#           # and one lower is source/definitions
#           # then this should be DEFINITIONS.{SUBFOLDER}
#           # where {SUBFOLDER} is the name of the first
#						# argument. Do this for each subfolder (e.g.)
#						# DEFINITIONS.{SUBFOLDER1}.{SUBFOLDER2}.{etc...}
#     definitions		# the exact name of the current folder
#     ${CMAKE_CURRENT_SOURCE_DIR}   # this is always the same
#     ${PROJECT_NAME}_INCLUDE_DIRS  # this is always the same
#     "hh"		# this is the extension of the header files
#     "cc"      # this is the extension of the implementation files
#     ${PROJECT_NAME}_HEADER_FILES  # this is always the same
#     ${PROJECT_NAME}_SOURCE_FILES) # this is always the same
#
# NOTE ON PROPAGATION:
#   There are two functions "PropagateSG" and "Propagate" that move variables
#   to the PARENT_SCOPE. "Propagate" should be used for standard variables
#   "PropagateSG" is a specialization of "Propagate" that comes from how
#   "AddSubdirectoryFiles" groups headers and implementation files. The second
#   argument of "AddSubdirectoryFiles" creates two variables SG_{ARGUMENT}_I and
#   SG_{ARGUMENT}_H for implementation and header files, respectively.
#   "PropagateSG" should only be used with the second argument of
#   in "AddSubdirectoryFiles"
#
#   SG stands for source group. Creating separate source groups are used
#   in IDEs (i.e. Xcode, Visual Studio, KDevelop, etc.) to sort the files
#   into subfolders instead of putting them all in one folder
#
#   Propagation means putting the variable in the scope of the file/function
#   that called it, i.e. variables are local and you must propagate the
#   the variable to the PARENT_SCOPE for the calling file/function to see
#   the variable. When you have nested folders, you must propagate the variable
#   to each level of the parent scope until the top-level CMakeLists.txt
#   file sees the variable. E.g. a variable in
#   definitions/subfolder1/subfolder2/CMakeLists.txt
#   must have a propagate statement in:
#     - source/definitions/subfolder1/subfolder2/CMakeLists.txt (handled by AddSubdirectoryFiles so it can be omitted)
#     - source/definitions/subfolder1/CMakeLists.txt
#     - source/definitions/CMakeLists.txt
#     - source/CMakeLists.txt
#   This final propagation in source/CMakeLists.txt makes it so the top-level
#   CMakeLists.txt sees the variable
#
#
# FINAL NOTE ON NESTED FOLDERS:
#   The macro "AddSubdirectorySource" modified the variables
#     - ${PROJECT_NAME}_HEADER_FILES
#	  - ${PROJECT_NAME}_SOURCE_FILES
#     - ${PROJECT_NAME}_INCLUDE_DIRS
#     - PROJECT_FOLDERS
#   These must be propagated also because the modification happens locally
#   source/rng/CMakeLists.txt shows how to do this when a subdirectory is added
#   and "AddSubdirectorySource" is also called in the same CMakeLists.txt
#


set(${PROJECT_NAME}_EXE_PATH ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)
set(${PROJECT_NAME}_EXE_NAME dark_arts PARENT_SCOPE)

#----------------------------------------
# Exceptions
#----------------------------------------
add_subdirectory(exceptions)
PropagateSG(EXCEPTIONS)


#----------------------------------------
# Input
#----------------------------------------
add_subdirectory(input)
PropagateSG(INPUT)

#----------------------------------------
# Quadrature rules, fem, and math things
#----------------------------------------
add_subdirectory(quadrature_fem_math)
PropagateSG(QUADRATURE_FEM_MATH)

#----------------------------------------
# Cell / Geometry Data
#----------------------------------------

add_subdirectory(cells)
PropagateSG(CELLS)

#----------------------------------------
# SDIRK data
#----------------------------------------

add_subdirectory(time_data)
PropagateSG(TIME_DATA)

#----------------------------------------
# Intensity and Temperature Storage Arrays
#----------------------------------------

add_subdirectory(variable_storage)
PropagateSG(VARIABLE_STORAGE)

#----------------------------------------
# Opacities and Heat Capacities
#----------------------------------------

add_subdirectory(materials)
PropagateSG(MATERIALS)

#----------------------------------------
# Time Marcher
#----------------------------------------

add_subdirectory(time_marcher)
PropagateSG(TIME_MARCHER)

#----------------------------------------
# Temperature and Radiation Solvers
#----------------------------------------

add_subdirectory(solvers)
PropagateSG(SOLVERS)

#----------------------------------------
# Stuff unique to Performing Transport Sweeps
#----------------------------------------

add_subdirectory(transport_sweep)
PropagateSG(transport_sweep)

#----------------------------------------
# Propagate variables to parents scope
#----------------------------------------
Propagate(${PROJECT_NAME}_HEADER_FILES)
Propagate(${PROJECT_NAME}_SOURCE_FILES)
Propagate(${PROJECT_NAME}_INCLUDE_DIRS)
Propagate(${PROJECT_NAME}_LINKED_LIBRARIES)
Propagate(PROJECT_FOLDERS)





###################################################
# --------------------------------------------
# Everything below here is a relic of PGM copying this from the meh 
# 2014 Code-a-thon repository
# --------------------------------------------
###################################################

#----------------------------------------
# Data
#----------------------------------------
# add_subdirectory(data)
# PropagateSG(DATA)

#----------------------------------------
# Definitions
#----------------------------------------
# add_subdirectory(definitions)
# PropagateSG(DEFINITIONS)

#----------------------------------------
# Geometry
#----------------------------------------
# add_subdirectory(geometry)
# PropagateSG(GEOMETRY)
# PropagateSG(GEOMETRY.CONSTRUCTS)
# PropagateSG(GEOMETRY.SURFACE)
# PropagateSG(GEOMETRY.VOLUME)
# PropagateSG(GEOMETRY.MESH)

#----------------------------------------
# Input/Output
#----------------------------------------
# add_subdirectory(io)
# PropagateSG(IO)


#----------------------------------------
# Integrators
#----------------------------------------
# add_subdirectory(integrators)
#PropagateSG(INTEGRATORS)


#----------------------------------------
# Management
#----------------------------------------
# add_subdirectory(management)
# PropagateSG(MANAGEMENT)


#----------------------------------------
# Material
#----------------------------------------
#add_subdirectory(material)
#PropagateSG(MATERIAL)


#----------------------------------------
# Parallel
#----------------------------------------
#add_subdirectory(parallel)
#PropagateSG(PARALLEL)


#----------------------------------------
# Tools
#----------------------------------------
# add_subdirectory(tools)
# PropagateSG(TOOLS)
# PropagateSG(TOOLS.MATRIX)



#----------------------------------------
# Visualization
#----------------------------------------
# add_subdirectory(visualization)
# PropagateSG(VISUALIZATION)




