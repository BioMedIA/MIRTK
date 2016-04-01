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

################################################################################
# @file  Settings.cmake
# @brief Non-default project settings.
#
# This file is included by basis_project_impl(), after it looked for the
# required and optional dependencies and the CMake variables related to the
# project directory structure were defined (see Directories.cmake file in
# @c BINARY_CONFIG_DIR). It is also included before the BasisSettings.cmake
# file.
#
# In particular build options should be added in this file using CMake's
# <a href="http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:option">
# option()</a> command. Further, any common settings related to using a found
# dependency can be set here if the basis_use_package() command was enable
# to import the required configuration of a particular external package.
#
# @ingroup BasisSettings
################################################################################

# By default, build applications
set(BUILD_APPLICATIONS_DEFAULT ON)

# By default, do not regenerate .rst files of command help pages
option(BUILD_DOCUMENTATION_SOURCES "Regenerate help of commands in reStructuredText format when missing" OFF)

# By default, most modules are enabled already...
mark_as_advanced(BUILD_ALL_MODULES)

# Choose between shared or static linkage (shared is recommended)
# Note: BUILD_SHARED_LIBS option added by mirtk_project_begin
set(MIRTK_BUILD_SHARED_LIBS_DEFAULT ON)

# Enable profiling of program execution (cf. Common/include/mirtkProfiling.h)
# Note: WITH_PROFILING option and MIRTK_WITH_PROFILING compile definition added
#       by mirtk_project_begin when WITH_PROFILING is ON
set(MIRTK_WITH_PROFILING_DEFAULT ON)

# Testing is yet very limited, hence mark this option as advanced for now
mark_as_advanced(BUILD_TESTING)

# Always use/find VTK when anyway required by one of the enabled modules
if (MODULE_PointSet)
  if (NOT WITH_VTK)
    message("VTK required by PointSet module, setting WITH_VTK to ON for all modules")
    basis_update_value(WITH_VTK ON)
  endif ()
endif ()

# Always use IO module when Image or PointSet module and Applications enabled
if (BUILD_APPLICATIONS AND (MODULE_Image OR MODULE_PointSet))
  if (NOT MODULE_IO)
    message("I/O module required by Applications using the enabled MODULE_Image or MODULE_PointSet,"
            " setting also MODULE_IO to ON")
    basis_update_value(MODULE_IO ON)
  endif ()
endif ()

# Installation directories
if (WIN32)
  set(INSTALL_CMAKE_MODULES_DIR "${INSTALL_SHARE_DIR}/CMake")
  set(MIRTK_TOOLS_DIR           "${INSTALL_LIBEXEC_DIR}/Tools")
else ()
  set(INSTALL_CMAKE_MODULES_DIR "${INSTALL_SHARE_DIR}/cmake")
  set(MIRTK_TOOLS_DIR           "${INSTALL_LIBEXEC_DIR}/tools")
endif ()
