# ==============================================================================
# Medical Image Registration ToolKit (MIRTK)
#
# Copyright 2016 Imperial College London
# Copyright 2016 Andreas Schuh
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

if (COMMAND mirtk_configure_module)
  return()
endif ()

include("${CMAKE_CURRENT_LIST_DIR}/MIRTKProjectBegin.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/MIRTKProjectEnd.cmake")

# ------------------------------------------------------------------------------
## Configure MIRTK module, to be called in top-level CMakeLists.txt of module
#
# This macro initializes the MIRTK module configuration by calling
# mirtk_project_begin. It then processes the subdirectories to add the
# build targets. Finally, the module configuration is finalized by calling
# mirtk_project_end. When a module needs to call additional code between
# the mandatory mirtk_project_begin and mirtk_project_end commands, it should
# replace this macro call by custom CMake code enclosed by the mirtk_project_begin
# and mirtk_project_end calls.
macro(mirtk_configure_module)
  mirtk_project_begin()
  foreach (SUBDIR IN LISTS PROJECT_SUBDIRS)
    basis_add_subdirectory("${SUBDIR}")
  endforeach ()
  mirtk_project_end()
endmacro ()
