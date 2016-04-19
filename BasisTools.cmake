# ============================================================================
# Copyright (c) 2011-2012 University of Pennsylvania
# Copyright (c) 2013-2016 Andreas Schuh
# All rights reserved.
#
# See COPYING file for license information or visit
# https://cmake-basis.github.io/download.html#license
# ============================================================================

##############################################################################
# @file  BasisTools.cmake
# @brief Definition of functions and macros used by BASIS project.
#
# This is the main module that is included by BASIS projects. Most of the other
# BASIS CMake modules are included by this main module and hence do not need
# to be included separately. In particular, all CMake modules which are part
# of BASIS and whose name does not include the prefix "Basis" are not
# supposed to be included directly by a project that makes use of BASIS.
# Only the modules with the prefix "Basis" should be included directly.
#
# @ingroup CMakeAPI
##############################################################################

# ----------------------------------------------------------------------------
# include guard
if (__BASIS_TOOLS_INCLUDED)
  return ()
else ()
  set (__BASIS_TOOLS_INCLUDED TRUE)
endif ()

# ----------------------------------------------------------------------------
# externally developed modules

# ExternalData.cmake module - yet only part of ITK, not CMake
include ("${CMAKE_CURRENT_LIST_DIR}/ExternalData.cmake")

# the module for the topological sort of modules according to their
# inter-dependencies was copied from the ITK v4 project
include ("${CMAKE_CURRENT_LIST_DIR}/TopologicalSort.cmake")

# ----------------------------------------------------------------------------
# BASIS modules
include ("${CMAKE_CURRENT_LIST_DIR}/CommonTools.cmake")
include ("${CMAKE_CURRENT_LIST_DIR}/DocTools.cmake")
include ("${CMAKE_CURRENT_LIST_DIR}/InterpTools.cmake")
include ("${CMAKE_CURRENT_LIST_DIR}/InstallationTools.cmake")
include ("${CMAKE_CURRENT_LIST_DIR}/MatlabTools.cmake")
include ("${CMAKE_CURRENT_LIST_DIR}/ProjectTools.cmake")
include ("${CMAKE_CURRENT_LIST_DIR}/RevisionTools.cmake")
include ("${CMAKE_CURRENT_LIST_DIR}/SlicerTools.cmake")
include ("${CMAKE_CURRENT_LIST_DIR}/TargetTools.cmake")
include ("${CMAKE_CURRENT_LIST_DIR}/ExportTools.cmake")
include ("${CMAKE_CURRENT_LIST_DIR}/UtilitiesTools.cmake")

# Issues/limitations of CMake's "function" command complicate the definition
# of a custom set_target_properties function which can be used to collect
# "global" information about targets imported from the CMake package
# configuration of project dependencies. The workaround in the custom
# set_target_properties function defined in the ImportTools.cmake is
# extremely inefficient and slows down the configuration step a lot
# (cf. https://github.com/cmake-basis/BASIS/issues/494).
# 
# The only need for collecting this information for all (executable)
# targets imported from dependencies is for generating the executable
# target info table for the BASIS Utilities (cf. UtilitiesTools.cmake).
# Hence, when these are not used, the ImportTools.cmake are not needed.
# Further, when a project does not consist of modules, the imported
# targets are available in the scope of the project.
#
# A project has to set BASIS_IMPORT_TARGETS to TRUE in its root CMakeLists.txt
# file before basis_project_begin() or basis_project_impl(), respectively.
