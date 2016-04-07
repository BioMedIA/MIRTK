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
  # Parse arguments
  if (NOT PROJECT_NAME)
    message(FATAL_ERROR "mirtk_add_library called outside project scope!")
  endif ()
  if (NOT COMMAND cmake_parse_arguments)
    include("${CMAKE_ROOT}/Modules/CMakeParseArguments.cmake")
  endif ()
  cmake_parse_arguments(TARGET "AUTO_REGISTER" "NAME_VAR;UID_VAR" "HEADERS;SOURCES;DEPENDS" ${ARGN})
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
      set(hdr "${PROJECT_INCLUDE_DIR}/${PROJECT_PACKAGE_NAME_L}/${hdr}")
    endif ()
    list(APPEND headers "${hdr}")
  endforeach ()
  basis_add_library(${target_name} ${TARGET_SOURCES} ${headers})
  basis_get_target_uid(target_uid ${target_name})
  if (TARGET_NAME_VAR)
    set(${TARGET_NAME_VAR} ${target_name} PARENT_SCOPE)
  endif ()
  if (TARGET_UID_VAR)
    set(${TARGET_UID_VAR} ${target_uid} PARENT_SCOPE)
  endif ()
  set_target_properties(
    ${target_uid} PROPERTIES
      VERSION       "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}"
      SOVERSION     "${PROJECT_SOVERSION}"
      OUTPUT_NAME   "${PROJECT_PACKAGE_NAME}${PROJECT_NAME}"
      DEFINE_SYMBOL "MIRTK_${PROJECT_NAME}_EXPORTS"
      DEBUG_POSTFIX D
  )
  target_include_directories(${target_uid}
    PUBLIC $<BUILD_INTERFACE:${PROJECT_INCLUDE_DIR}>
           $<INSTALL_INTERFACE:${INSTALL_INCLUDE_DIR}>
    PRIVATE ${PROJECT_CODE_DIR}
  )
  if (BUILD_SHARED_LIBS)
    if (TARGET_AUTO_REGISTER)
      target_compile_definitions(${target_uid} PRIVATE MIRTK_AUTO_REGISTER)
    endif ()
  else ()
    if (WIN32)
      set_target_properties(${target_uid} PROPERTIES SUFFIX "_s.lib")
    endif ()
    set_target_properties(${target_uid} PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
  endif ()
  # VTK 7.0.0 adds its CMake modules directory to CMAKE_MODULE_PATH
  # and this directory contains a GenerateExportHeader.cmake module
  include("${BASIS_MODULE_PATH}/GenerateExportHeader.cmake")
  # Generate MIRTK<Module>Export.h file
  generate_export_header(${target_uid}
    PREFIX_NAME              "MIRTK_"
    BASE_NAME                "${PROJECT_NAME}"
    EXPORT_MACRO_NAME        "${PROJECT_NAME}_EXPORT"
    NO_EXPORT_MACRO_NAME     "${PROJECT_NAME}_NO_EXPORT"
    DEPRECATED_MACRO_NAME    "${PROJECT_NAME}_DEPRECATED"
    NO_DEPRECATED_MACRO_NAME "${PROJECT_NAME}_NO_DEPRECATED"
    STATIC_DEFINE            "${PROJECT_NAME}_STATIC_DEFINE"
    EXPORT_FILE_NAME         "${BINARY_INCLUDE_DIR}/${PROJECT_PACKAGE_NAME_L}/${PROJECT_NAME}Export.h"
  )
  if (BUILD_SHARED_LIBS)
    if (WIN32)
      if (CMAKE_VERSION VERSION_LESS 3.4)
        message(FATAL_ERROR "Build of DLLs on Windows (BUILD_SHARED_LIBS=ON) requires CMake version 3.4 or greater!")
      endif ()
      set_target_properties(${target_uid} PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
    endif ()
    if (MSVC)
      target_compile_options(${target_uid} PRIVATE /bigobj)
    endif ()
  else ()
    target_compile_definitions(${target_uid} PRIVATE MIRTK_${PROJECT_NAME}_STATIC_DEFINE)
  endif ()
  # Add link dependencies and record all directories containing
  # external/imported library files for [DY]LD_LIBRARY_PATH
  # setting in "mirtk" command execution script
  if (TARGET_DEPENDS)
    if (NOT COMMAND mirtk_target_dependencies)
      include("${MIRTK_MODULE_PATH}/mirtkTargetDependencies.cmake")
    endif ()
    mirtk_target_dependencies(${target_name} ${TARGET_DEPENDS})
  endif ()
endfunction()
