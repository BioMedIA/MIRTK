# ============================================================================
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
# ============================================================================

##############################################################################
# @file  CMakeLists.txt
# @brief Root build configuration file.
#
# @sa <a href="https://cmake-basis.github.io">CMake BASIS</a>
##############################################################################

# ----------------------------------------------------------------------------
# Minimum required CMake version
cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)

# ----------------------------------------------------------------------------
# Include BASIS policies, settings, macros, and functions
set(BASIS_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMake/Basis")
set(MIRTK_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMake/Modules")
set(CMAKE_MODULE_PATH "${MIRTK_MODULE_PATH}" "${BASIS_MODULE_PATH}")

include(mirtkPolicies NO_POLICY_SCOPE)

include(BasisSettings)
include(BasisTools)
include(mirtkTools)

# ----------------------------------------------------------------------------
# Default settings
#
# These variables have to be set before (mirtk|basis)_project_begin.
# See CMake/Config/Settings.cmake file for other settings.

# By default, require these optional dependencies such that the library
# features which depend on these external libraries are available
if (UNIX)
  set(WITH_ZLIB_DEFAULT ON)
endif ()

# ThirdParty/Eigen
if (NOT Eigen3_DIR)
  if (Eigen3_INCLUDE_DIR)
    set(Eigen3_DIR "${Eigen3_INCLUDE_DIR}")
  elseif (EIGEN3_INCLUDE_DIR)
    set(Eigen3_DIR "${EIGEN3_INCLUDE_DIR}")
  elseif (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/Eigen/signature_of_eigen3_matrix_library")
    set(Eigen3_DIR "${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/Eigen")
  endif ()
endif ()

# ThirdParty/Boost
if (NOT Boost_DIR)
  if (BOOST_INCLUDEDIR)
    set(Boost_DIR "${BOOST_INCLUDEDIR}")
  elseif (BOOST_ROOT)
    set(Boost_DIR "${BOOST_ROOT}/include")
  elseif (BOOSTROOT)
    set(Boost_DIR "${BOOSTROOT}/include")
  elseif (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/Boost/boost/version.hpp")
    set(Boost_DIR "${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/Boost")
  endif ()
endif ()

# ----------------------------------------------------------------------------
# Support use of ccache to speed up builds (-D WITH_CCACHE=ON)
option(WITH_CCACHE "Use ccache if available" OFF)
mark_as_advanced(WITH_CCACHE)
if (WITH_CCACHE)
  if (CMAKE_CXX_COMPILER MATCHES "ccache")
    message(STATUS "Detected C++ compiler is in ccache path, ignoring WITH_CCACHE option")
  else ()
    find_program(CCACHE_COMMAND "ccache")
    if (NOT CCACHE_COMMAND)
      message(FATAL_ERROR "Use of ccache requested, but CCACHE_COMMAND not found!")
    endif ()
    message(STATUS "Using ccache for compilation of C++ code")
    mark_as_advanced(CCACHE_COMMAND)
    # For Unix Makefiles and Ninja generators
    if (CMAKE_VERSION VERSION_LESS 3.4)
      set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_COMMAND}")
    else ()
      set(CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_COMMAND}")
    endif ()
    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
      add_compile_options("-Qunused-arguments")
    endif ()
  endif ()
else ()
  message(STATUS "Not using ccache for compilation of C++ code (WITH_CCACHE=${WITH_CCACHE})")
endif ()

# ----------------------------------------------------------------------------
# Configure build system

# Initialize project
mirtk_project_begin()

# Reset global properties used to collect information about MIRTK commands
# and the LD_LIBRARY_PATH of external libraries required by these
set_property(GLOBAL PROPERTY MIRTK_COMMANDS     "")
set_property(GLOBAL PROPERTY MIRTK_LIBRARY_PATH "")

# Process modules
foreach (MODULE IN LISTS PROJECT_MODULES_ENABLED)
  basis_add_module(${MODULE})
  basis_use_module(${MODULE})
endforeach ()

# Process subdirectories
foreach (SUBDIR IN LISTS PROJECT_SUBDIRS)
  basis_add_subdirectory(${SUBDIR})
endforeach ()

# Finalize project
mirtk_project_end()
