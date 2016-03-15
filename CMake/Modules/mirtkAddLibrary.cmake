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

# ------------------------------------------------------------------------------
## Add MIRTK module library
function(mirtk_add_library)
  if (NOT PROJECT_NAME)
    message(FATAL_ERROR "mirtk_add_library called outside project scope!")
  endif ()
  if (NOT COMMAND cmake_parse_arguments)
    include("${CMAKE_ROOT}/Modules/CMakeParseArguments.cmake")
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
  # Add library target
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
  basis_get_target_uid(target_uid ${target_name})
  set_target_properties(
    ${target_uid} PROPERTIES
      VERSION             "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}"
      SOVERSION           "${PROJECT_SOVERSION}"
      ARCHIVE_OUTPUT_NAME "${OUTPUT_NAME}"
      LIBRARY_OUTPUT_NAME "${OUTPUT_NAME}"
      DEFINE_SYMBOL       "MIRTK_${PROJECT_NAME}_EXPORTS"
  )
  target_include_directories(${target_uid}
    PUBLIC $<BUILD_INTERFACE:${PROJECT_INCLUDE_DIR}>
           $<INSTALL_INTERFACE:${INSTALL_INCLUDE_DIR}>
    PRIVATE ${PROJECT_CODE_DIR}
  )
  if (TARGET_AUTO_REGISTER AND BUILD_SHARED_LIBS)
    target_compile_definitions(${target_uid} PRIVATE MIRTK_AUTO_REGISTER)
  endif ()
  # Generate export header file and set WINDOWS_EXPORT_ALL_SYMBOLS target property
  include(GenerateExportHeader)
  generate_export_header(${target_uid}
    PREFIX_NAME              "MIRTK_"
    BASE_NAME                "${PROJECT_NAME}"
    EXPORT_MACRO_NAME        "${PROJECT_NAME}_EXPORT"
    NO_EXPORT_MACRO_NAME     "${PROJECT_NAME}_NO_EXPORT"
    DEPRECATED_MACRO_NAME    "${PROJECT_NAME}_DEPRECATED"
    NO_DEPRECATED_MACRO_NAME "${PROJECT_NAME}_NO_DEPRECATED"
    STATIC_DEFINE            "${PROJECT_NAME}_STATIC_DEFINE"
    EXPORT_FILE_NAME         "${BINARY_INCLUDE_DIR}/mirtk${PROJECT_NAME}Export.h"
    DEFINE_NO_DEPRECATED
  )
  if (WIN32 AND BUILD_SHARED_LIBS)
    if (CMAKE_VERSION VERSION_LESS 3.4)
      message(FATAL_ERROR "Build of DLLs on Windows (BUILD_SHARED_LIBS=ON) requires CMake version 3.4 or greater!")
    endif ()
    set_target_properties(${target_uid} PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
  endif ()
  # Add link dependencies and record all directories containing
  # external/imported library files for [DY]LD_LIBRARY_PATH
  # setting in "mirtk" command execution script
  if (NOT COMMAND mirtk_get_target_name)
    include("${MIRTK_MODULE_PATH}/mirtkGetTargetName.cmake")
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
          get_property(is_imported TARGET ${dep_name} PROPERTY IMPORTED)
          if (is_imported)
            basis_get_target_location(dep_path ${dep_name} PATH)
            if (dep_path)
              set_property(GLOBAL APPEND PROPERTY MIRTK_LIBRARY_PATH ${dep_path})
            endif ()
          endif ()
        else ()
          add_dependencies(${target_name} ${dep_name})
        endif ()
      else ()
        target_link_libraries(${target_name} ${keyword} ${dep_name})
        if (IS_ABSOLUTE "${dep_name}")
          get_filename_component(dep_path "${dep_name}" PATH)
          set_property(GLOBAL APPEND PROPERTY MIRTK_LIBRARY_PATH ${dep_path})
        endif ()
      endif ()
      set(keyword)
    endif ()
  endforeach ()
  # Remove duplicates from MIRTK_LIBRARY_PATH
  get_property(ldpath GLOBAL PROPERTY MIRTK_LIBRARY_PATH)
  list(REMOVE_DUPLICATES ldpath)
  set_property(GLOBAL PROPERTY MIRTK_LIBRARY_PATH "${ldpath}")
  include(GenerateExportHeader)
  generate_export_header(${target_uid}
    PREFIX_NAME              "MIRTK_"
    BASE_NAME                "${PROJECT_NAME}"
    EXPORT_MACRO_NAME        "${PROJECT_NAME}_EXPORT"
    NO_EXPORT_MACRO_NAME     "${PROJECT_NAME}_NO_EXPORT"
    DEPRECATED_MACRO_NAME    "${PROJECT_NAME}_DEPRECATED"
    NO_DEPRECATED_MACRO_NAME "${PROJECT_NAME}_NO_DEPRECATED"
    STATIC_DEFINE            "${PROJECT_NAME}_STATIC_DEFINE"
    EXPORT_FILE_NAME         "${BINARY_INCLUDE_DIR}/mirtk${PROJECT_NAME}Export.h"
    DEFINE_NO_DEPRECATED
  )
endfunction()
