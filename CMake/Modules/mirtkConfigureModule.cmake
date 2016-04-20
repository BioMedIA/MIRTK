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
  if (NOT COMMAND mirtk_project_begin)
    include("${MIRTK_MODULE_PATH}/mirtkProjectBegin.cmake")
  endif ()
  if (NOT COMMAND mirtk_project_end)
    include("${MIRTK_MODULE_PATH}/mirtkProjectEnd.cmake")
  endif ()
  mirtk_project_begin()
  if (GTest_FOUND AND UNIX)
    option(GTest_NO_PTHREADS "Whether gtest library does not require pthreads" OFF)
    mark_as_advanced(GTest_NO_PTHREADS)
    if (NOT GTest_NO_PTHREADS)
      if (DEFINED CMAKE_THREAD_PREFER_PTHREAD)
        set(_CMAKE_THREAD_PREFER_PTHREAD ${CMAKE_THREAD_PREFER_PTHREAD})
      endif ()
      set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
      find_package(Threads QUIET)
      if (DEFINED _CMAKE_THREAD_PREFER_PTHREAD)
        set(CMAKE_THREAD_PREFER_PTHREAD ${_CMAKE_THREAD_PREFER_PTHREAD})
        unset(_CMAKE_THREAD_PREFER_PTHREAD)
      endif ()
    endif ()
  endif ()
  foreach (SUBDIR IN LISTS PROJECT_SUBDIRS)
    basis_add_subdirectory("${SUBDIR}")
  endforeach ()
  mirtk_project_end()
endmacro ()
