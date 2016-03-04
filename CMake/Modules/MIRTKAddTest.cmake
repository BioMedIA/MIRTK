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

if (COMMAND mirtk_add_test)
  return()
endif ()

include(CMakeParseArguments)
include("${CMAKE_CURRENT_LIST_DIR}/MIRTKGetTargetName.cmake")

# ------------------------------------------------------------------------------
## Add build target for executable MIRTK command
function(mirtk_add_test target_name)
  # Parse arguments
  if (NOT PROJECT_NAME)
    message(FATAL_ERROR "mirtk_add_test called outside project scope!")
  endif ()
  cmake_parse_arguments(TARGET "" "" "SOURCES;DEPENDS;OPTIONAL" ${ARGN})
  if (TARGET_UNPARSED_ARGUMENTS)
    message(FATAL_ERROR "mirtk_add_test called with unrecognized arguments: ${TARGET_UNPARSED_ARGUMENTS}!")
  endif ()
  # Add optional dependencies if targets exist
  foreach (dep IN LISTS TARGET_OPTIONAL)
    mirtk_get_target_name(dep_name ${dep})
    if (TARGET ${dep_name})
      list(APPEND TARGET_DEPENDS ${dep_name})
    endif ()
  endforeach ()
  # Add test executable build target and CTest target
  if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/test${target_name}.cc")
    basis_add_test(test${target_name}
      SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/test${target_name}.cc
        ${TARGET_SOURCES}
      LINK_DEPENDS
        ${GTEST_LIBRARY}
        ${GTEST_MAIN_LIBRARY}
        ${TARGET_DEPENDS}
    )
  else ()
    message(FATAL_ERROR "Source file for MIRTK unit test \"test${target_name}.cc\" not found!")
  endif ()
endfunction()
