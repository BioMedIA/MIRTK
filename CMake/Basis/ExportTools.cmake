# ============================================================================
# Copyright (c) 2011-2012 University of Pennsylvania
# Copyright (c) 2013-2016 Andreas Schuh
# All rights reserved.
#
# See COPYING file for license information or visit
# https://cmake-basis.github.io/download.html#license
# ============================================================================

##############################################################################
# @file  ExportTools.cmake
# @brief Functions and macros for the export of targets.
#
# @ingroup CMakeTools
##############################################################################

if (__BASIS_EXPORTIMPORTTOOLS_INCLUDED)
  return ()
else ()
  set (__BASIS_EXPORTIMPORTTOOLS_INCLUDED TRUE)
endif ()


## @addtogroup CMakeUtilities
#  @{


# ----------------------------------------------------------------------------
## @brief Add target to export set.
#
# Targets are added to the export set named after the top-level project.
# This is necessary because CMake does not allow us to split the exports
# into different sets when there are inter-dependencies between library
# targets of the modules.
#
# @param[out] EXPORT_OPTION Export option for install() command including
#                           the EXPORT option name. Set to an empty string
#                           if target is not installed (see @p ARGN).
# @param[in]  TARGET_UID    UID of target to be exported.
# @param[in]  IS_TEST       Whether given target is a test executable or library.
# @param[in]  ARGN          Optional installation destinations. The actual
#                           values are ignored and argument @c TRUE should
#                           be used if no specific destination paths are given,
#                           but the target is also to be included in the export
#                           set used for installation. If ARGN is empty, the
#                           target is only added to the build tree export set.
#                           When @p ARGN is @c 0, @c FALSE or @c OFF, the target
#                           is not added to list of installation exports.
function (basis_add_export_target EXPORT_OPTION TARGET_UID IS_TEST)
  set (EXPORT_SET "${TOPLEVEL_PROJECT_NAME}")
  if (IS_TEST)
    basis_set_project_property (PROJECT "${EXPORT_SET}" APPEND PROPERTY TEST_EXPORT_TARGETS "${TARGET_UID}")
    set (${EXPORT_OPTION} "" PARENT_SCOPE)
  else ()
    basis_set_project_property (PROJECT "${EXPORT_SET}" APPEND PROPERTY EXPORT_TARGETS "${TARGET_UID}")
    if (ARGN AND NOT "^${ARGN}$" STREQUAL "^(0|false|FALSE|off|OFF)$")
      basis_set_project_property (PROJECT "${EXPORT_SET}" APPEND PROPERTY INSTALL_EXPORT_TARGETS "${TARGET_UID}")
      set (${EXPORT_OPTION} "EXPORT;${EXPORT_SET}" PARENT_SCOPE)
    else ()
      set (${EXPORT_OPTION} "" PARENT_SCOPE)
    endif ()
  endif ()
endfunction ()

# ----------------------------------------------------------------------------
## @brief Add target to custom export set.
#
# Targets are added to the export set named after the top-level project.
# This is necessary because CMake does not allow us to split the exports
# into different sets when there are inter-dependencies between library
# targets of the modules.
#
# @param[in]  TARGET_UID UID of target to add to the export set.
# @param[in]  IS_TEST    Whether given target is a test executable or library.
function (basis_add_custom_export_target TARGET_UID IS_TEST)
  set (EXPORT_SET "${TOPLEVEL_PROJECT_NAME}")
  if (IS_TEST)
    basis_set_project_property (PROJECT "${EXPORT_SET}" APPEND PROPERTY TEST_EXPORT_TARGETS "${TARGET_UID}")
  else ()
    basis_set_project_property (PROJECT "${EXPORT_SET}" APPEND PROPERTY CUSTOM_EXPORT_TARGETS "${TARGET_UID}")
  endif ()
endfunction ()

# ----------------------------------------------------------------------------
## @brief Get soname of object file.
#
# This function extracts the soname from object files in the ELF format on
# systems where the objdump command is available. On all other systems,
# an empty string is returned.
#
# @param [out] SONAME  The soname of the object file.
# @param [in]  OBJFILE Object file in ELF format.
function (basis_get_soname SONAME OBJFILE)
  # get absolute path of object file
  basis_get_target_uid (TARGET_UID ${OBJFILE})
  if (TARGET TARGET_UID)
    basis_get_target_location (OBJFILE ${TARGET_UID} ABSOLUTE)
    string (REPLACE "$<${BASIS_GE_CONFIG}>" "${CMAKE_BUILD_TYPE}" OBJFILE "${OBJFILE}")
  else ()
    get_filename_component (OBJFILE "${OBJFILE}" ABSOLUTE)
  endif ()
  # usually CMake did this already
  find_program (CMAKE_OBJDUMP NAMES objdump DOC "The objdump command")
  # run objdump and extract soname
  execute_process (
    COMMAND ${CMAKE_OBJDUMP} -p "${OBJFILE}"
    COMMAND sed -n "-e's/^[[:space:]]*SONAME[[:space:]]*//p'"
    RESULT_VARIABLE STATUS
    OUTPUT_VARIABLE SONAME_OUT
    ERROR_QUIET
  )
  # return
  if (STATUS EQUAL 0)
    set (${SONAME} "${SONAME_OUT}" PARENT_SCOPE)
  else ()
    set (${SONAME} "" PARENT_SCOPE)
  endif ()
endfunction ()

# ----------------------------------------------------------------------------
## @brief Generate header of exports file.
function (basis_export_header CODE)
  set (C "# Generated by BASIS\n\n")
  set (C "${C}if (\"\${CMAKE_MAJOR_VERSION}.\${CMAKE_MINOR_VERSION}\" LESS 2.8)\n")
  set (C "${C}  message (FATAL_ERROR \"CMake >= 2.8.4 required\")\n")
  set (C "${C}endif ()\n")
  set (C "${C}cmake_policy (PUSH)\n")
  set (C "${C}cmake_policy (VERSION 2.8.4)\n")
  set (C "${C}#----------------------------------------------------------------\n")
  set (C "${C}# Generated CMake target import file.\n")
  set (C "${C}#----------------------------------------------------------------\n")
  set (C "${C}\n# Commands may need to know the format version.\n")
  set (C "${C}set (CMAKE_IMPORT_FILE_VERSION 1)\n")
  set (${CODE} "${C}" PARENT_SCOPE)
endfunction ()

# ----------------------------------------------------------------------------
## @brief Add code to compute prefix relative to @c INSTALL_CONFIG_DIR.
function (basis_export_prefix CODE)
  set (C "\n# Compute the installation prefix relative to this file.\n")
  set (C "${C}get_filename_component (_IMPORT_PREFIX \"\${CMAKE_CURRENT_LIST_FILE}\" PATH)\n")
  string (REGEX REPLACE "[/\\]" ";" DIRS "${INSTALL_CONFIG_DIR}")
  foreach (D IN LISTS DIRS)
    set (C "${C}get_filename_component (_IMPORT_PREFIX \"\${_IMPORT_PREFIX}\" PATH)\n")
  endforeach ()
  set (${CODE} "${${CODE}}${C}" PARENT_SCOPE)
endfunction ()

# ----------------------------------------------------------------------------
## @brief Add code to add import targets.
function (basis_export_import_targets CODE)
  set (C)
  foreach (T IN LISTS ARGN)
    basis_get_fully_qualified_target_uid (UID "${T}")
    set (C "${C}\n# Create import target \"${UID}\"\n")
    get_target_property (BASIS_TYPE ${T} "BASIS_TYPE")
    if (BASIS_TYPE MATCHES "EXECUTABLE")
      set (C "${C}add_executable (${UID} IMPORTED)\n")
    elseif (BASIS_TYPE MATCHES "LIBRARY|MODULE|MEX")
      if (BASIS_TYPE MATCHES "MEX|MCC")
        set (TYPE SHARED)
      elseif (BASIS_TYPE MATCHES "SCRIPT")
        set (TYPE UNKNOWN)
      else ()
        string (REGEX REPLACE "_LIBRARY" "" TYPE "${BASIS_TYPE}")
      endif ()
      set (C "${C}add_library (${UID} ${TYPE} IMPORTED)\n")
    else ()
      message (FATAL_ERROR "Cannot export target ${T} of type ${BASIS_TYPE}! Use NOEXPORT option.")
    endif ()
    set (C "${C}set_target_properties (${UID} PROPERTIES BASIS_TYPE \"${BASIS_TYPE}\")\n")
  endforeach ()
  set (${CODE} "${${CODE}}${C}" PARENT_SCOPE)
endfunction ()

# ----------------------------------------------------------------------------
## @brief Add code to set properties of imported targets for build tree.
function (basis_export_build_properties CODE)
  set (C)
  if (CMAKE_BUILD_TYPE)
    set (CONFIG "${CMAKE_BUILD_TYPE}")
  else ()
    set (CONFIG "noconfig")
  endif ()
  string (TOUPPER "${CONFIG}" CONFIG_U)
  foreach (T IN LISTS ARGN)
    basis_get_fully_qualified_target_uid (UID "${T}")
    get_target_property (BASIS_TYPE ${T} BASIS_TYPE)
    set (C "${C}\n# Import target \"${UID}\" for configuration \"${CONFIG}\"\n")
    set (C "${C}set_property (TARGET ${UID} APPEND PROPERTY IMPORTED_CONFIGURATIONS ${CONFIG_U})\n")
    set (C "${C}set_target_properties (${UID} PROPERTIES\n")
    if (BASIS_TYPE MATCHES "SCRIPT|LIBRARY|MEX")
      get_target_property (LANGUAGE ${T} LANGUAGE)
      if (LANGUAGE)
        set (C "${C}  IMPORTED_LINK_INTERFACE_LANGUAGES_${CONFIG_U} \"${LANGUAGE}\"\n")
      endif ()
      get_target_property (LINK_DEPENDS ${T} LINK_DEPENDS)
      set (LINK_UIDS)
      foreach (LINK_DEPEND IN LISTS LINK_DEPENDS)
        basis_get_fully_qualified_target_uid (LINK_UID ${LINK_DEPEND})
        list (APPEND LINK_UIDS "${LINK_UID}")
      endforeach ()
      set (C "${C}  IMPORTED_LINK_INTERFACE_LIBRARIES_${CONFIG_U} \"${LINK_UIDS}\"\n")
    endif ()
    basis_get_target_location (LOCATION ${T} ABSOLUTE)
    string (REPLACE "$<${BASIS_GE_CONFIG}>" "${CONFIG}" LOCATION "${LOCATION}")
    set (C "${C}  IMPORTED_LOCATION_${CONFIG_U} \"${LOCATION}\"\n")
    set (C "${C}  )\n")
  endforeach ()
  set (${CODE} "${${CODE}}${C}" PARENT_SCOPE)
endfunction ()

# ----------------------------------------------------------------------------
## @brief Add code to set properties of imported targets for installation.
function (basis_export_install_properties CODE)
  set (C)
  if (CMAKE_BUILD_TYPE)
    set (CONFIG "${CMAKE_BUILD_TYPE}")
  else ()
    set (CONFIG "noconfig")
  endif ()
  string (TOUPPER "${CONFIG}" CONFIG_U)
  foreach (T IN LISTS ARGN)
    basis_get_fully_qualified_target_uid (UID "${T}")
    get_target_property (BASIS_TYPE ${T} BASIS_TYPE)
    set (C "${C}\n# Import target \"${UID}\" for configuration \"${CONFIG}\"\n")
    set (C "${C}set_property (TARGET ${UID} APPEND PROPERTY IMPORTED_CONFIGURATIONS ${CONFIG_U})\n")
    set (C "${C}set_target_properties (${UID} PROPERTIES\n")
    basis_get_target_location (LOCATION ${T} POST_INSTALL_RELATIVE)
    set (C "${C}  IMPORTED_LOCATION_${CONFIG_U} \"\${_IMPORT_PREFIX}/${LOCATION}\"\n")
    if (BASIS_TYPE MATCHES "SCRIPT|LIBRARY|MEX")
      get_target_property (LANGUAGE ${T} LANGUAGE)
      if (LANGUAGE)
        set (C "${C}  IMPORTED_LINK_INTERFACE_LANGUAGES_${CONFIG_U} \"${LANGUAGE}\"\n")
      endif ()
      get_target_property (LINK_DEPENDS ${T} LINK_DEPENDS)
      set (LINK_UIDS)
      foreach (LINK_DEPEND IN LISTS LINK_DEPENDS)
        basis_get_fully_qualified_target_uid (LINK_UID ${LINK_DEPEND})
        list (APPEND LINK_UIDS "${LINK_UID}")
      endforeach ()
      set (C "${C}  IMPORTED_LINK_INTERFACE_LIBRARIES_${CONFIG_U} \"${LINK_UIDS}\"\n")
    endif ()
    set (C "${C}  )\n")
  endforeach ()
  set (${CODE} "${${CODE}}${C}" PARENT_SCOPE)
endfunction ()

# ----------------------------------------------------------------------------
## @brief Add footer of exports file.
function (basis_export_footer CODE)
  set (C "\n# Cleanup temporary variables.\n")
  set (C "${C}set (_IMPORT_PREFIX)\n")
  set (C "${C}\n# Commands beyond this point should not need to know the version.\n")
  set (C "${C}set (CMAKE_IMPORT_FILE_VERSION)\n")
  set (C "${C}cmake_policy (POP)\n")
  set (${CODE} "${${CODE}}${C}" PARENT_SCOPE)
endfunction ()

# ----------------------------------------------------------------------------
## @brief Export all targets added by basis_add_* commands.
function (basis_export_targets)
  # parse arguments
  CMAKE_PARSE_ARGUMENTS (ARGN "" "FILE;CUSTOM_FILE" "" ${ARGN})

  if (NOT ARGN_FILE)
    message (FATAL_ERROR "basis_export_targets(): FILE option is required!")
  endif ()
  if (NOT ARGN_CUSTOM_FILE)
    message (FATAL_ERROR "basis_export_targets(): CUSTOM_FILE option is required!")
  endif ()

  if (IS_ABSOLUTE ${ARGN_FILE})
    message (FATAL_ERROR "basis_export_targets(): FILE option argument must be a relative path!")
  endif ()
  if (IS_ABSOLUTE ${ARGN_CUSTOM_FILE})
    message (FATAL_ERROR "basis_export_targets(): CUSTOM_FILE option argument must be a relative path!")
  endif ()

  # --------------------------------------------------------------------------
  # export non-custom targets
  basis_get_project_property (EXPORT_TARGETS)
  if (EXPORT_TARGETS)
    # add link dependencies to build tree export set
    foreach (EXPORT_TARGET IN LISTS EXPORT_TARGETS)
      get_target_property (LINK_DEPENDS ${EXPORT_TARGET} BASIS_LINK_DEPENDS)
      foreach (LINK_DEPEND IN LISTS LINK_DEPENDS)
        if (TARGET ${LINK_DEPEND})
          list (FIND  EXPORT_TARGETS ${LINK_DEPEND} IDX)
          if (IDX EQUAL -1)
            get_target_property (IMPORTED ${LINK_DEPEND} IMPORTED)
            if (NOT IMPORTED)
              message (AUTHOR_WARNING
                "Adding missing link dependency ${LINK_DEPEND} of ${EXPORT_TARGET} to list of export() targets."
                " Add the target ${LINK_DEPEND} to the export set using basis_add_export_target() after the"
                " add_library(${LINK_DEPEND}) command and before basis_project_end(). If the target is added"
                " using basis_add_library(${LINK_DEPEND}), remove the NOEXPORT option or add EXPORT instead when"
                " the global variable BASIS_EXPORT_DEFAULT is set to FALSE (actual value is ${BASIS_EXPORT_DEFAULT})."
              )
              list (APPEND EXPORT_TARGETS ${LINK_DEPEND})
            endif ()
          endif ()
        endif ()
      endforeach ()
    endforeach ()
    list (REMOVE_DUPLICATES EXPORT_TARGETS)
    # set namespace of exported import targets
    if (BASIS_USE_TARGET_UIDS AND BASIS_USE_FULLY_QUALIFIED_UIDS)
      set (NAMESPACE_OPT)
    elseif (TOPLEVEL_PROJECT_NAMESPACE_CMAKE)
      set (NAMESPACE_OPT NAMESPACE "${TOPLEVEL_PROJECT_NAMESPACE_CMAKE}${BASIS_NAMESPACE_DELIMITER_CMAKE}")
    endif ()
    # export build tree targets
    export (
      TARGETS   ${EXPORT_TARGETS}
      FILE      "${BINARY_LIBCONF_DIR}/${ARGN_FILE}"
      ${NAMESPACE_OPT}
    )
    # install export set
    basis_get_project_property (INSTALL_EXPORT_TARGETS)
    if (INSTALL_EXPORT_TARGETS AND NOT BASIS_BUILD_ONLY)
      foreach (COMPONENT "${BASIS_RUNTIME_COMPONENT}" "${BASIS_LIBRARY_COMPONENT}")
        install (
          EXPORT      "${PROJECT_NAME}"
          DESTINATION "${INSTALL_CONFIG_DIR}"
          FILE        "${ARGN_FILE}"
          COMPONENT   "${COMPONENT}"
          ${NAMESPACE_OPT}
        )
      endforeach ()
    endif ()
  endif ()

  # --------------------------------------------------------------------------
  # export custom targets and/or test targets
  basis_get_project_property (CUSTOM_EXPORT_TARGETS)
  basis_get_project_property (TEST_EXPORT_TARGETS)

  if (CUSTOM_EXPORT_TARGETS OR TEST_EXPORT_TARGETS)

    # write exports for build tree
    basis_export_header           (CONTENT)
    basis_export_import_targets   (CONTENT ${CUSTOM_EXPORT_TARGETS} ${TEST_EXPORT_TARGETS})
    basis_export_build_properties (CONTENT ${CUSTOM_EXPORT_TARGETS} ${TEST_EXPORT_TARGETS})
    basis_export_footer           (CONTENT)

    file (WRITE "${BINARY_LIBCONF_DIR}/${ARGN_CUSTOM_FILE}" "${CONTENT}")
    unset (CONTENT)

    # write exports for installation - excluding test targets
    if (CUSTOM_EXPORT_TARGETS AND NOT BASIS_BUILD_ONLY)
      set (INSTALL_EXPORT_FILE "${PROJECT_BINARY_DIR}/CMakeFiles/Export/${INSTALL_CONFIG_DIR}/${ARGN_CUSTOM_FILE}")

      basis_export_header (CONTENT)
      basis_export_prefix (CONTENT)
      basis_export_import_targets (CONTENT ${CUSTOM_EXPORT_TARGETS})
      basis_export_install_properties (CONTENT ${CUSTOM_EXPORT_TARGETS})
      basis_export_footer (CONTENT)

      file (WRITE "${INSTALL_EXPORT_FILE}" "${CONTENT}")
      unset (CONTENT)

      install (
        FILES       "${INSTALL_EXPORT_FILE}"
        DESTINATION "${INSTALL_CONFIG_DIR}"
      )
    endif ()

  endif ()
endfunction ()


## @}
# end of Doxygen group
