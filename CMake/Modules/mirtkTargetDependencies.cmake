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

if (COMMAND mirtk_target_dependencies)
  return()
endif ()

# ------------------------------------------------------------------------------
## @brief Add target dependency
function (mirtk_target_dependencies target_name)
  # Check arguments
  cmake_parse_arguments(ARGN "OPTIONAL" "" "" ${ARGN})
  if (NOT ARGN_UNPARSED_ARGUMENTS)
    return()
  endif ()
  basis_get_target_uid(target_uid ${target_name})
  if (NOT COMMAND mirtk_get_target_name)
    include("${MIRTK_MODULE_PATH}/mirtkGetTargetName.cmake")
  endif ()
  # Get current MIRTK_LIBRARY_PATH used by "mirtk" command launcher
  get_property(ldpath GLOBAL PROPERTY MIRTK_LIBRARY_PATH)
  # Add dependency on link library file, imported library target, or other target
  unset(dep_attr)
  foreach (dep IN LISTS ARGN_UNPARSED_ARGUMENTS)
    if (dep MATCHES "^(general|optimized|debug)$")
      set(dep_attr "${dep}")
    else ()
      mirtk_get_target_name(dep_name "${dep}")
      if (TARGET "${dep_name}")
        get_target_property(dep_type ${dep_name} TYPE)
        # Add dependency on (imported|other) library target
        if (dep_type MATCHES "LIBRARY")
          # Need to use basis_target_link_libraries such that BASIS can add
          # link dependencies of non-imported targets added via add_library
          # instead of basis_add_library to the export set if necessary
          basis_target_link_libraries(${target_uid} ${dep_attr} ${dep_name})
          get_property(dep_imported TARGET ${dep_name} PROPERTY IMPORTED)
          if (dep_imported)
            basis_get_target_location(dep_path ${dep_name} PATH)
            if (dep_path)
              list(APPEND ldpath "${dep_path}")
            endif ()
          endif ()
        # Add dependency on non-library target
        else ()
          add_dependencies(${target_uid} ${dep_name})
        endif ()
      # Add dependency on external link library file
      elseif (NOT ARGN_OPTIONAL OR (ARGN_OPTIONAL AND dep_name MATCHES "\\.(so|a|dylib|lib)(\\.[0-9].*)?$"))
        target_link_libraries(${target_uid} ${dep_attr} ${dep_name})
        if (IS_ABSOLUTE "${dep_name}")
          get_filename_component(dep_path "${dep_name}" PATH)
          if (CMAKE_HOST_WIN32)
            get_filename_component(dll_name "${dep_name}" NAME_WE)
            string(REPLACE "/lib/" "/bin/" dll_path "${dep_path}/")
            string(REPLACE "/Lib/" "/Bin/" dll_path "${dll_path}")
            if (EXISTS "${dll_path}${dll_name}.dll")
              string(REGEX REPLACE "/+$" "" dep_path "${dll_path}")
            endif ()
          endif ()
          list(APPEND ldpath "${dep_path}")
        endif ()
      endif ()
      unset(dep_attr)
    endif ()
  endforeach ()
  if (dep_attr)
    message(FATAL_ERROR "Expected library file path or imported target after ${dep_attr}\nArguments = [${ARGN_UNPARSED_ARGUMENTS}]")
  endif ()
  # Update MIRTK_LIBRARY_PATH
  if (ldpath)
    list(REMOVE_DUPLICATES ldpath)
  endif ()
  set_property(GLOBAL PROPERTY MIRTK_LIBRARY_PATH "${ldpath}")
endfunction ()
