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

if (COMMAND mirtk_add_executable)
  return()
endif ()

# Name of subdirectory of MIRTK command executables.
# This variable is also used by the Applications module and therefore global.
if (WIN32)
  set(MIRTK_TOOLS_SUBDIR "Tools")
else ()
  set(MIRTK_TOOLS_SUBDIR "tools")
endif ()

# ------------------------------------------------------------------------------
## Add build target for executable MIRTK command
function(mirtk_add_executable target_name)
  # Parse arguments
  if (NOT PROJECT_NAME)
    message(FATAL_ERROR "mirtk_add_executable called outside project scope!")
  endif ()
  if (NOT COMMAND cmake_parse_arguments)
    include("${CMAKE_ROOT}/Modules/CMakeParseArguments.cmake")
  endif ()
  cmake_parse_arguments(TARGET "" "" "SOURCES;DEPENDS;OPTIONAL" ${ARGN})
  if (TARGET_UNPARSED_ARGUMENTS)
    message(FATAL_ERROR "mirtk_add_executable called with unrecognized arguments: ${TARGET_UNPARSED_ARGUMENTS}!")
  endif ()
  if (target_name MATCHES "^(commands|commands-help)$")
    message(FATAL_ERROR "Internally used target name cannot be used for command: ${target_name}")
  endif ()
  # Add executable target
  set(LANGUAGE CXX)
  set(BINARY_LIBEXEC_DIR  "${BINARY_LIBEXEC_DIR}/${MIRTK_TOOLS_SUBDIR}")
  if (MIRTK_TOOLS_DIR)
    set(INSTALL_LIBEXEC_DIR "${MIRTK_TOOLS_DIR}")
  else ()
    set(INSTALL_LIBEXEC_DIR "${INSTALL_LIBEXEC_DIR}/${MIRTK_TOOLS_SUBDIR}")
  endif ()
  if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${target_name}.cc")
    basis_add_executable(${target_name}.cc ${TARGET_SOURCES} LIBEXEC)
  elseif (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${target_name}.cxx")
    basis_add_executable(${target_name}.cxx ${TARGET_SOURCES} LIBEXEC)
  elseif (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${target_name}.cpp")
    basis_add_executable(${target_name}.cpp ${TARGET_SOURCES} LIBEXEC)
  elseif (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${target_name}.py" OR
          EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${target_name}.py.in")
    basis_add_executable(${target_name}.py ${TARGET_SOURCES} LIBEXEC FINAL)
    set(LANGUAGE PYTHON)
  elseif (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${target_name}.sh" OR
          EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${target_name}.sh.in")
    basis_add_executable(${target_name}.sh ${TARGET_SOURCES} LIBEXEC FINAL)
    set(LANGUAGE BASH)
  else ()
    message(FATAL_ERROR "Source file for MIRTK command ${target_name} not found!"
                        " Must have extension '.cc', '.cxx', '.cpp', '.sh', or '.py'."
                        " Exclude extension from mirtk_add_executable argument.")
  endif ()
  set_property(GLOBAL APPEND PROPERTY MIRTK_COMMANDS ${target_name})
  if (LANGUAGE STREQUAL CXX)
    basis_get_target_uid(target_uid "${target_name}")
    # Add target dependencies
    if (TARGET_DEPENDS OR TARGET_OPTIONAL)
      if (NOT COMMAND mirtk_target_dependencies)
        include("${MIRTK_MODULE_PATH}/mirtkTargetDependencies.cmake")
      endif ()
      mirtk_target_dependencies(${target_name} ${TARGET_DEPENDS})
      mirtk_target_dependencies(${target_name} ${TARGET_OPTIONAL} OPTIONAL)
    endif ()
    # Set RPATH of installed binary
    #
    # Would be set automatically by basis_project_end/mirtk_project_end when
    # BASIS_INSTALL_RPATH is set to TRUE/ON in mirtk_project_begin instead.
    # This requires, however, the use of basis_target_link_libraries instead of
    # the CMake target_link_libraries command for basis_get_target_link_libraries
    # used by basis_set_target_install_rpath to work. This function is very
    # costly. Hence, as we know what RPATH we want without inspecting the list
    # of all link dependencies, we rather set the INSTALL_RPATH ourselves.
    if (CMAKE_HOST_APPLE)
      set(ORIGIN "@loader_path")
    else ()
      set(ORIGIN "\$ORIGIN")
    endif ()
    basis_get_relative_path(rpath
      "${CMAKE_INSTALL_PREFIX}/${INSTALL_LIBEXEC_DIR}"
      "${CMAKE_INSTALL_PREFIX}/${INSTALL_LIBRARY_DIR}"
    )
    string(REGEX REPLACE "/+$" "" rpath "${rpath}")
    set_target_properties(${target_uid} PROPERTIES INSTALL_RPATH "${ORIGIN}/${rpath}")
  endif ()
endfunction()
