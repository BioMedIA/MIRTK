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

include(CMakeParseArguments)

# ------------------------------------------------------------------------------
## Add MIRTK module library
function(mirtk_add_library)
  if (NOT PROJECT_NAME)
    message(FATAL_ERROR "mirtk_add_library called outside project scope!")
  endif ()
  cmake_parse_arguments(TARGET "" "" "HEADERS;SOURCES;DEPENDS" ${ARGN})
  if (TARGET_UNPARSED_ARGUMENTS)
    message(FATAL_ERROR "mirtk_add_library called with unrecognized arguments: ${TARGET_UNPARSED_ARGUMENTS}!")
  endif ()
  foreach (arg IN ITEMS HEADERS SOURCES DEPENDS)
    if (NOT TARGET_${arg} AND ${arg})
      set(TARGET_${arg} "${${arg}}")
    endif ()
  endforeach ()
  set(target_name "Lib${PROJECT_NAME}")
  set(headers)
  foreach (hdr IN LISTS TARGET_HEADERS)
    if (NOT IS_ABSOLUTE "${hdr}")
      set(hdr "${PROJECT_INCLUDE_DIR}/${hdr}")
    endif ()
    list(APPEND headers "${hdr}")
  endforeach ()
  basis_add_library(${target_name} ${TARGET_SOURCES} ${headers})
  basis_set_target_properties(
    ${target_name} PROPERTIES
      VERSION             "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}"
      SOVERSION           "${PROJECT_SOVERSION}"
      LIBRARY_OUTPUT_NAME "${PROJECT_PACKAGE_NAME}${PROJECT_NAME}"
  )
  foreach (dep IN LISTS TARGET_DEPENDS)
    if (TARGET ${dep})
      get_property(type TARGET ${dep} PROPERTY TYPE)
      if (type MATCHES LIBRARY)
        basis_target_link_libraries(${target_name} ${dep})
      else ()
        basis_add_dependencies(${target_name} ${dep})
      endif ()
    else ()
      basis_target_link_libraries(${target_name} ${dep})
    endif ()
  endforeach ()
endfunction()
