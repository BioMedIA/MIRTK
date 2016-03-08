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

if (COMMAND mirtk_get_target_name)
  return()
endif ()

# ------------------------------------------------------------------------------
## Either prefix target name with "mirtk::" namespace identifier or remove it
function (mirtk_get_target_name target_name dep)
  if (TARGET ${dep})
    set(${target_name} ${dep} PARENT_SCOPE)
  elseif (TARGET mirtk::${dep})
    set(${target_name} mirtk::${dep} PARENT_SCOPE)
  elseif (dep MATCHES "^mirtk::(.*)$" AND TARGET "${CMAKE_MATCH_1}")
    set(${target_name} ${CMAKE_MATCH_1} PARENT_SCOPE)
  else ()
    set(${target_name} ${dep} PARENT_SCOPE)
  endif ()
endfunction ()
