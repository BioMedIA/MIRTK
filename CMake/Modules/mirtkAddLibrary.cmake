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

if (COMMAND mirtk_add_library)
  return()
endif ()

include(CMakeParseArguments)
include("${CMAKE_CURRENT_LIST_DIR}/mirtkGetTargetName.cmake")

# ------------------------------------------------------------------------------
## Add MIRTK module library
function(mirtk_add_library)
  if (NOT PROJECT_NAME)
    message(FATAL_ERROR "mirtk_add_library called outside project scope!")
  endif ()
  cmake_parse_arguments(TARGET "AUTO_REGISTER" "" "HEADERS;SOURCES;DEPENDS" ${ARGN})
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
  set(OUTPUT_NAME "${PROJECT_PACKAGE_NAME}${PROJECT_NAME}")
  basis_add_library(${target_name} ${TARGET_SOURCES} ${headers})
  set_target_properties(
    ${target_name} PROPERTIES
      VERSION             "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}"
      SOVERSION           "${PROJECT_SOVERSION}"
      ARCHIVE_OUTPUT_NAME "${OUTPUT_NAME}"
      LIBRARY_OUTPUT_NAME "${OUTPUT_NAME}"
  )
  basis_get_target_uid(target_uid ${target_name})
  target_include_directories(${target_uid}
    PUBLIC $<BUILD_INTERFACE:${PROJECT_INCLUDE_DIR}>
           $<INSTALL_INTERFACE:${INSTALL_INCLUDE_DIR}>
    PRIVATE ${PROJECT_CODE_DIR}
  )
  if (TARGET_AUTO_REGISTER AND BUILD_SHARED_LIBS)
    target_compile_definitions(${target_uid} PRIVATE MIRTK_AUTO_REGISTER)
  endif ()
  set(keyword)
  foreach (dep IN LISTS TARGET_DEPENDS)
    if (dep MATCHES "^(debug|optimized|general)$")
      set(keyword ${dep})
    else ()
      mirtk_get_target_name(dep_name ${dep})
      if (TARGET ${dep_name})
        get_property(type TARGET ${dep_name} PROPERTY TYPE)
        if (type MATCHES LIBRARY)
          target_link_libraries(${target_name} ${keyword} ${dep_name})
        else ()
          add_dependencies(${target_name} ${dep_name})
        endif ()
      else ()
        target_link_libraries(${target_name} ${keyword} ${dep_name})
      endif ()
      set(keyword)
    endif ()
  endforeach ()
endfunction()
