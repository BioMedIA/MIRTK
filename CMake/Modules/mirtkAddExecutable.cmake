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
  set(BINARY_LIBEXEC_DIR "${TOPLEVEL_BINARY_LIBEXEC_DIR}/mirtk")
  if (TOPLEVEL_INSTALL_LIBEXEC_DIR MATCHES "/mirtk/*$")
    set(INSTALL_LIBEXEC_DIR "${TOPLEVEL_INSTALL_LIBEXEC_DIR}")
  else ()
    set(INSTALL_LIBEXEC_DIR "${TOPLEVEL_INSTALL_LIBEXEC_DIR}/mirtk")
  endif ()
  if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${target_name}.cc")
    basis_add_executable(${target_name}.cc ${TARGET_SOURCES} LIBEXEC)
  elseif (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${target_name}.cxx")
    basis_add_executable(${target_name}.cxx ${TARGET_SOURCES} LIBEXEC)
  elseif (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${target_name}.cpp")
    basis_add_executable(${target_name}.cpp ${TARGET_SOURCES} LIBEXEC)
  elseif (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${target_name}.py" OR
          EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${target_name}.py.in")
    basis_add_executable(${target_name}.py ${TARGET_SOURCES} LIBEXEC)
    set(LANGUAGE PYTHON)
  elseif (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${target_name}.sh" OR
          EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${target_name}.sh.in")
    basis_add_executable(${target_name}.sh ${TARGET_SOURCES} LIBEXEC)
    set(LANGUAGE BASH)
  else ()
    message(FATAL_ERROR "Source file for MIRTK command ${target_name} not found!"
                        " Must have extension '.cc', '.cxx', '.cpp', '.sh', or '.py'."
                        " Exclude extension from mirtk_add_executable argument.")
  endif ()
  set_property(GLOBAL APPEND PROPERTY MIRTK_COMMANDS ${target_name})
  basis_get_target_uid(target_uid "${target_name}")
  if (LANGUAGE STREQUAL CXX)
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
    # Add custom target to generate command documentation page
    if (BUILD_DOCUMENTATION AND BUILD_DOCUMENTATION_SOURCES AND NOT "^${target_name}$" STREQUAL "^help-rst$")
      basis_get_target_uid(help_rst_target help-rst)
      if (TARGET ${help_rst_target})
        basis_get_target_uid(mirtk_target mirtk)
        basis_get_target_location(mirtk_path mirtk ABSOLUTE)
        add_custom_target(${target_uid}-description
          COMMAND "${mirtk_path}" help-rst "${target_name}" -generated -description -noheaders
                  -output "${TOPLEVEL_PROJECT_DOC_DIR}/commands/_descriptions/${target_name}.rst"
                  -output-brief-description "${TOPLEVEL_PROJECT_DOC_DIR}/commands/_summaries/${target_name}.rst"
          DEPENDS ${mirtk_target} ${help_rst_target} ${target_uid}
          COMMENT "Extracting command description from mirtk ${target_name} -help"
        )
        add_custom_target(${target_uid}-help
          COMMAND "${mirtk_path}" help-rst "${target_name}" -generated -orphan
                  "-include-description" "_descriptions/${target_name}.rst"
                  -output "${TOPLEVEL_PROJECT_DOC_DIR}/commands/${target_name}.rst"
          DEPENDS ${mirtk_target} ${help_rst_target} ${target_uid} ${target_uid}-description
          COMMENT "Generating documentation page for command ${target_name}"
        )
        basis_get_target_uid(commands_help_target commands-help)
        if (NOT TARGET ${commands_help_target})
          add_custom_target(${commands_help_target})
        endif ()
        add_dependencies(${commands_help_target} ${target_uid}-help)
      else ()
        message(WARNING "Target help-rst missing! Skipping auto-generation of command documentation.")
      endif ()
    endif ()
    # Add link dependencies and record all directories containing
    # external/imported library files for [DY]LD_LIBRARY_PATH
    # setting in "mirtk" command execution script
    if (TARGET_DEPENDS OR TARGET_OPTIONAL)
      if (NOT COMMAND mirtk_get_target_name)
        include("${MIRTK_MODULE_PATH}/mirtkGetTargetName.cmake")
      endif ()
      # Add required dependencies
      set(keyword)
      foreach (dep IN LISTS TARGET_DEPENDS)
        if (dep MATCHES "^(debug|optimized|general)$")
          set(keyword ${dep})
        else ()
          mirtk_get_target_name(dep_name ${dep})
          if (TARGET ${dep_name})
            get_property(type TARGET ${dep_name} PROPERTY TYPE)
            if (type MATCHES LIBRARY)
              target_link_libraries(${target_uid} ${keyword} ${dep_name})
              get_property(is_imported TARGET ${dep_name} PROPERTY IMPORTED)
              if (is_imported)
                basis_get_target_location(dep_path ${dep_name} PATH)
                if (dep_path)
                  set_property(GLOBAL APPEND PROPERTY MIRTK_LIBRARY_PATH ${dep_path})
                endif ()
              endif ()
            else ()
              add_dependencies(${target_uid} ${dep_name})
            endif ()
          else ()
            target_link_libraries(${target_uid} ${keyword} ${dep_name})
            if (IS_ABSOLUTE "${dep_name}")
              get_filename_component(dep_path "${dep_name}" PATH)
              set_property(GLOBAL APPEND PROPERTY MIRTK_LIBRARY_PATH ${dep_path})
            endif ()
          endif ()
          set(keyword)
        endif ()
      endforeach ()
      # Add optional dependencies if targets exist
      set(keyword)
      foreach (dep IN LISTS TARGET_OPTIONAL)
        if (dep MATCHES "^(debug|optimized|general)$")
          set(keyword ${dep})
        else ()
          mirtk_get_target_name(dep_name ${dep})
          if (TARGET ${dep_name})
            get_property(type TARGET ${dep_name} PROPERTY TYPE)
            if (type MATCHES LIBRARY)
              target_link_libraries(${target_uid} ${keyword} ${dep_name})
              get_property(is_imported TARGET ${dep_name} PROPERTY IMPORTED)
              if (is_imported)
                basis_get_target_location(dep_path ${dep_name} PATH)
                if (dep_path)
                  set_property(GLOBAL APPEND PROPERTY MIRTK_LIBRARY_PATH ${dep_path})
                endif ()
              endif ()
            else ()
              add_dependencies(${target_uid} ${dep_name})
            endif ()
          elseif (dep_name MATCHES "\\.(so|a|dylib|lib)(\\.[0-9].*)?$")
            target_link_libraries(${target_uid} ${keyword} ${dep_name})
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
    endif ()
  endif ()
endfunction()
