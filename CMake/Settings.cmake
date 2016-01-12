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
option(BUILD_APPLICATIONS "Request build and installation of command-line tools" ON)

# By default, do not regenerate .rst files of command help pages
option(BUILD_DOCUMENTATION_SOURCES "Regenerate help of commands in reStructuredText format when missing" OFF)

# By default, most modules are enabled already...
mark_as_advanced(BUILD_ALL_MODULES)

# Choose build configuration from CMake list of common configuration types
if (CMAKE_CONFIGURATION_TYPES)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS ${CMAKE_CONFIGURATION_TYPES})
  set_property(CACHE CMAKE_CONFIGURATION_TYPES PROPERTY TYPE INTERNAL)
endif ()

# Choose between shared or static linkage (shared is recommended)
option(BUILD_SHARED_LIBS "Request build of shared libraries" ON)
mark_as_advanced(BUILD_SHARED_LIBS)

# Enable profiling of program execution (cf. Common/include/mirtkProfiling.h)
option(WITH_PROFILING "Enable profiling of program execution" ON)
if (WITH_PROFILING)
  add_definitions(-DMIRTK_WITH_PROFILING)
endif ()

# Testing is yet very limited, hence mark this option as advanced for now
mark_as_advanced(BUILD_TESTING)

# CodeBlocks generators add this cache entry, hide it when set
if (CMAKE_CODEBLOCKS_EXECUTABLE)
  mark_as_advanced(FORCE CMAKE_CODEBLOCKS_EXECUTABLE)
else ()
  mark_as_advanced(CLEAR CMAKE_CODEBLOCKS_EXECUTABLE)
endif ()

# Always use/find VTK when anyway required by one of the enabled modules
if (MODULE_PointSet)
  if (NOT WITH_VTK)
    message("VTK required by PointSet module, setting WITH_VTK to ON for all modules")
    set(WITH_VTK ON CACHE BOOL "Request build with VTK library" FORCE)
  endif ()
endif ()
