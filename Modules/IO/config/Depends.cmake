# ==============================================================================
# Medical Image Registration ToolKit (MIRTK)
#
# Copyright 2013-2015 Imperial College London
# Copyright 2013-2015 Andreas Schuh
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# ==============================================================================

##############################################################################
# @file  Depends.cmake
# @brief Find additional dependencies.
#
# This file is included by basis_project_impl() after it included the
# BasisProject.cmake file of the project and collected information about its
# modules. Furthermore, it is included before it looks for the packages
# specified as arguments to the basis_project() command. At this point, the
# CMake project is not initialized yet and most BASIS variables are not set.
#
# Use this file to find additional dependencies or dependencies which are
# only required by a project if certain modules are enabled.
#
# Use case: If this project has a module which is a Slicer Extension
# and the project itself shall be build as Slicer Extension only if the
# module is enabled, the Slicer package configuration file has to be
# included here such that the Slicer settings are imported before any
# module is configured. This is done by using the command
#
# Another use case would be that you want to specify advanced options to
# the basis_find_package() function which you cannot specify as part of
# the dependencies arguments of the basis_project() function.
#
# Example:
# @code
# if (SlicerExtension_ENABLED)
#   # requires Slicer if the SlicerExtension module is enabled
#   basis_find_package (Slicer REQUIRED)
#   basis_use_package (Slicer)
# endif ()
# @endcode
#
# @ingroup BasisSettings
##############################################################################

option(WITH_NiftiCLib "Require system installation of NIfTI C library" OFF)
mark_as_advanced(WITH_NiftiCLib)

if (WITH_NiftiCLib)
  basis_find_package(NiftiCLib REQUIRED)
else ()
  add_subdirectory(
    "${TOPLEVEL_PROJECT_SOURCE_DIR}/ThirdParty/NiftiCLib"
    "${TOPLEVEL_PROJECT_BINARY_DIR}/ThirdParty/NiftiCLib"
  )
  set(NiftiCLib_INCLUDE_DIRS
    "${TOPLEVEL_PROJECT_SOURCE_DIR}/ThirdParty/NiftiCLib/znzlib"
    "${TOPLEVEL_PROJECT_SOURCE_DIR}/ThirdParty/NiftiCLib/niftilib"
  )
  set(NiftiCLib_LIBRARIES nifticlib)
  set(NiftiCLib_FOUND TRUE)
endif ()
include_directories(${NiftiCLib_INCLUDE_DIRS})
