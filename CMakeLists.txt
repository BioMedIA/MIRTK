# ============================================================================
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
# ============================================================================

##############################################################################
# @file  CMakeLists.txt
# @brief Root build configuration file.
#
# @note This package utilizes <a href="http://opensource.andreasschuh.com/cmake-basis">CMake BASIS</a>.
#
##############################################################################

# ----------------------------------------------------------------------------
# Minimum required CMake version
cmake_minimum_required(VERSION 2.8.12)

# ----------------------------------------------------------------------------
# Include BASIS policies, settings, macros, and functions
set(BASIS_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMake/Basis")
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMake" ${BASIS_MODULE_PATH})

include("${BASIS_MODULE_PATH}/BasisSettings.cmake" NO_POLICY_SCOPE)
include("${BASIS_MODULE_PATH}/BasisTools.cmake")

include(mirtk_add_library)
if (BUILD_APPLICATIONS)
  include(mirtk_add_executable)
endif ()

# ----------------------------------------------------------------------------
# Override default BASIS settings and hide unused cache entries
set(BUILD_MODULES_BY_DEFAULT    OFF       CACHE INTERNAL "" FORCE)
set(BASIS_BUILD_ONLY            OFF       CACHE INTERNAL "" FORCE)
set(BASIS_ALL_DOC               OFF       CACHE INTERNAL "" FORCE)
set(BASIS_COMPILE_MATLAB        ON        CACHE INTERNAL "" FORCE)
set(BASIS_COMPILE_SCRIPTS       OFF       CACHE INTERNAL "" FORCE)
set(BASIS_CREATE_EXPORTS_FILE   ON        CACHE INTERNAL "" FORCE)
set(BASIS_INSTALL_RPATH         ON        CACHE INTERNAL "" FORCE)
set(BASIS_INSTALL_SCHEME        "default" CACHE INTERNAL "" FORCE)
set(BASIS_INSTALL_SITE_PACKAGES OFF       CACHE INTERNAL "" FORCE)
set(BASIS_REGISTER              ON        CACHE INTERNAL "" FORCE)
set(BASIS_REVISION_INFO         ON        CACHE INTERNAL "" FORCE)
set(BASIS_SUPERBUILD_MODULES    OFF       CACHE INTERNAL "" FORCE)

foreach (lang IN ITEMS BASH CXX PYTHON PERL)
  set(BUILD_BASIS_UTILITIES_FOR_${lang} OFF CACHE INTERNAL "" FORCE)
endforeach ()
unset(lang)

# ----------------------------------------------------------------------------
# Configure build system

# Initialize project
basis_project_begin()

# Main "mirtk" script
#
# This wrapper script sets up the environment and executes the specified
# subcommand. Do not use mirtk_add_executable which is only used for
# commands that are executed by this wrapper script.
if (BUILD_APPLICATIONS)
  if (WIN32)
    # TODO: Either implement it in Python or as Windows Batch command (.cmd)
    message(WARNING "MIRTK wrapper script not implemented for Windows")
  else ()
    basis_add_executable("${PROJECT_TOOLS_DIR}/mirtk.sh")
  endif ()
endif ()

# Add auxiliary command used to (re-)format help output
#
# Note: Added as first application because mirtk_add_executable adds a
#       dependency on this target for each command to generate the documentation
#       when BUILD_DOCUMENTATION is ON.
if (BUILD_APPLICATIONS)
  basis_add_executable(help-rst "${PROJECT_TOOLS_DIR}/help-rst.py" LIBEXEC)
  basis_set_target_properties(help-rst PROPERTIES OUTPUT_NAME "mirtk-help-rst")
  basis_set_project_property(APPEND PROPERTY MIRTK_COMMANDS help-rst)
  basis_add_dependencies(help-rst mirtk)
endif ()

# Process modules
if (NOT PROJECT_IS_MODULE)
  foreach (MODULE IN LISTS PROJECT_MODULES_ENABLED)
    basis_add_module(${MODULE})
    basis_use_module(${MODULE})
  endforeach ()
endif ()

# Process subdirectories
foreach (SUBDIR IN LISTS PROJECT_SUBDIRS)
  basis_add_subdirectory(${SUBDIR})
endforeach ()

# Finalize project
basis_project_end()

# Hide unused BASIS installation directory settings
foreach (lang IN ITEMS BASH JYTHON PYTHON MATLAB PERL PYTHON)
  set_property(CACHE INSTALL_${lang}_SITE_DIR PROPERTY TYPE INTERNAL)
endforeach ()

# ------------------------------------------------------------------------------
# Generate Bash completions file
if (UNIX)
  get_property(COMMANDS GLOBAL PROPERTY MIRTK_COMMANDS)
  list(REMOVE_ITEM COMMANDS help-rst)
  if (COMMANDS)
    message(STATUS "Writing Bash completions script...")
    basis_list_to_delimited_string(COMMANDS " " ${COMMANDS})
    file(WRITE "${PROJECT_BINARY_DIR}/mirtk_bash_completions.sh"
"
# Auto-generated by MIRTK root CMakeLists.txt during CMake configure.
#
# To enable Bash completion of \"mirtk\" subcommands, add the following to your ~/.bashrc file:
#
# if [ -f \"\$MIRTK_DIR/share/mirtk_bash_completions.sh\" ]; then
#   . \"\$MIRTK_DIR/share/mirtk_bash_completions.sh\"
# fi
#
mirtk_bash_completions()
{
  local commands='${COMMANDS}'
  local noreply='true'
  if [ \$COMP_CWORD -eq 2 ]; then
    if [[ \${COMP_WORDS[1]} == help ]] || [[ \${COMP_WORDS[1]} == help-rst ]]; then
      COMPREPLY=( \$(compgen -W \"\$commands\" -- \${COMP_WORDS[COMP_CWORD]}) )
      noreply='false'
    fi
  elif [ \$COMP_CWORD -eq 1 ]; then
    COMPREPLY=( \$(compgen -W \"help help-rst \$commands\" -- \${COMP_WORDS[COMP_CWORD]}) )
    noreply='false'
  fi
  if [[ \$noreply == true ]]; then
    local IFS=$'\n'
    COMPREPLY=($(compgen -o plusdirs -f -- \"\${COMP_WORDS[COMP_CWORD]}\"))
    compopt -o filenames +o nospace 2>/dev/null || compgen -f /non-existing-dir/ > /dev/null
  fi
}
complete -F mirtk_bash_completions mirtk
"
    )
    install(FILES "${PROJECT_BINARY_DIR}/mirtk_bash_completions.sh" DESTINATION "${INSTALL_SHARE_DIR}")
    message(STATUS "Writing Bash completions script... done")
  endif ()
endif ()
