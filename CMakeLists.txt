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
# @sa <a href="http://opensource.andreasschuh.com/cmake-basis">CMake BASIS</a>
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

# Reset global properties used to collect information about MIRTK commands
# and the LD_LIBRARY_PATH of external libraries required by these
set_property(GLOBAL PROPERTY MIRTK_COMMANDS     "")
set_property(GLOBAL PROPERTY MIRTK_LIBRARY_PATH "")

# ----------------------------------------------------------------------------
# Configure build system

# Initialize project
mirtk_project_begin()

# Main "mirtk" script
#
# This wrapper script sets up the environment and executes the specified
# subcommand. Do not use mirtk_add_executable which is only used for
# commands that are executed by this wrapper script.
if (BUILD_APPLICATIONS)
  basis_add_executable("${PROJECT_TOOLS_DIR}/mirtk.py")
endif ()

# Add auxiliary command used to (re-)format help output
#
# Note: Added as first application because mirtk_add_executable adds a
#       dependency on this target for each command to generate the documentation
#       when BUILD_DOCUMENTATION is ON.
if (BUILD_APPLICATIONS AND TARGET mirtk)
  basis_add_executable(help-rst "${PROJECT_TOOLS_DIR}/help-rst.py" LIBEXEC)
  basis_set_target_properties(help-rst PROPERTIES
    OUTPUT_DIRECTORY "${BINARY_LIBEXEC_DIR}/mirtk"
  )
  if (NOT INSTALL_LIBEXEC_DIR MATCHES "/mirtk/*$")
    basis_set_target_properties(help-rst PROPERTIES
      INSTALL_DIRECTORY "${INSTALL_LIBEXEC_DIR}/mirtk"
    )
  endif ()
  basis_add_dependencies(help-rst mirtk)
endif ()

# Process modules
foreach (MODULE IN LISTS PROJECT_MODULES_ENABLED)
  basis_add_module(${MODULE})
  basis_use_module(${MODULE})
endforeach ()

# Process subdirectories
foreach (SUBDIR IN LISTS PROJECT_SUBDIRS)
  basis_add_subdirectory(${SUBDIR})
endforeach ()

# Set MIRTK_LIBRARY_PATH for "mirtk" script
get_property(ldpath GLOBAL PROPERTY MIRTK_LIBRARY_PATH)
list(REMOVE_DUPLICATES ldpath)
if (UNIX)
  list(REMOVE_ITEM ldpath "/lib")
  list(REMOVE_ITEM ldpath "/lib64")
  list(REMOVE_ITEM ldpath "/usr/lib")
  list(REMOVE_ITEM ldpath "/usr/lib64")
  list(REMOVE_ITEM ldpath "/usr/local/lib")
  list(REMOVE_ITEM ldpath "/opt/local/lib")
endif ()
basis_set_target_properties(mirtk PROPERTIES
  SCRIPT_DEFINITIONS
    "basis_list_to_delimited_string(MIRTK_LIBRARY_PATH_CONFIG \"', '\" NOAUTOQUOTE \"${ldpath}\")"
)
unset(ldpath)

# Install CMake modules
basis_add_subdirectory(CMake)

# Finalize project
mirtk_project_end()

# ------------------------------------------------------------------------------
# Generate auto-completion scripts
if (BUILD_APPLICATIONS AND UNIX)
  include(mirtkGenerateBashCompletions)
  mirtk_generate_bash_completions()
endif ()
