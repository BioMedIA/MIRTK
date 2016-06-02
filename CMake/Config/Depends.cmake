# ==============================================================================
# Medical Image Registration ToolKit (MIRTK)
#
# Copyright 2013-2016 Imperial College London
# Copyright 2013-2016 Andreas Schuh
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

# Directory with source files of third-party libraries
set(MIRTK_THIRDPARTY_DIR "${MIRTK_SOURCE_DIR}/ThirdParty")


# ------------------------------------------------------------------------------
# libLBFGS (optional)
set(LBFGS_SOURCE_DIR "${MIRTK_THIRDPARTY_DIR}/LBFGS")

# Set default for WITH_LibLBFGS option
if (NOT DEFINED WITH_LBFGS_DEFAULT AND EXISTS "${LBFGS_SOURCE_DIR}/CMakeLists.txt")
  set(WITH_LBFGS_DEFAULT ON)
endif ()

# When third-party submodule is available and library is needed
if ((DEFINED WITH_LBFGS AND WITH_LBFGS) OR (NOT DEFINED WITH_LBFGS AND WITH_LBFGS_DEFAULT))
  basis_check_if_package_is_needed_by_modules(is_needed LBFGS)
  if (is_needed AND EXISTS "${LBFGS_SOURCE_DIR}/CMakeLists.txt")
    message(STATUS "Configuring module ThirdParty/LBFGS...")
    set(LBFGS_INSTALL_HEADERS     FALSE)
    set(LBFGS_INSTALL_STATIC_LIBS TRUE) # needed because of CMake install(EXPORT)
    set(LBFGS_INSTALL_LIBRARY_DIR ${INSTALL_LIBRARY_DIR})
    set(LBFGS_lib_TARGET_NAME     LibLBFGS)
    basis_add_export_target(export_args ${LBFGS_lib_TARGET_NAME} FALSE ${LBFGS_INSTALL_STATIC_LIBS})
    if (export_args)
      list(GET export_args 1 LBFGS_EXPORT_NAME)
    else ()
      set(LBFGS_EXPORT_NAME)
    endif ()
    add_subdirectory("${LBFGS_SOURCE_DIR}" "${MIRTK_BINARY_DIR}/ThirdParty/LBFGS")
    set_target_properties(${LBFGS_lib_TARGET_NAME} PROPERTIES OUTPUT_NAME MIRTKlbfgs)
    set(LBFGS_FOUND       TRUE)
    set(LBFGS_LIBRARIES   ${LBFGS_lib_TARGET_NAME})
    set(LBFGS_INCLUDE_DIR "${LBFGS_SOURCE_DIR}/include")
    set(LBFGS_DIR         "${MIRTK_BINARY_DIR}/ThirdParty/LBFGS")
    message(STATUS "Configuring module ThirdParty/LBFGS... - done")
  endif ()
  unset(is_needed)
endif ()
