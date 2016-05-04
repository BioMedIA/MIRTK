# ============================================================================
# Copyright (c) 2011-2012 University of Pennsylvania
# Copyright (c) 2013-2016 Andreas Schuh
# All rights reserved.
#
# See COPYING file for license information or visit
# https://cmake-basis.github.io/download.html#license
# ============================================================================

##############################################################################
# @file  InstallationTools.cmake
# @brief CMake functions used for installation.
#
# @ingroup CMakeTools
##############################################################################

# ----------------------------------------------------------------------------
# include guard
if (__BASIS_INSTALLATIONTOOLS_INCLUDED)
  return ()
else ()
  set (__BASIS_INSTALLATIONTOOLS_INCLUDED TRUE)
endif ()


## @addtogroup CMakeUtilities
# @{


# ============================================================================
# Installation
# ============================================================================

# ----------------------------------------------------------------------------
## @brief Specify rules to run at install time.
#
# This function replaces CMake's
# <a href="https://cmake.org/cmake/help/v3.5/command/install.html">install()</a> command.
#
# @sa https://cmake.org/cmake/help/v3.5/command/install.html
#
# @ingroup CMakeAPI
function (basis_install)
  include(CMakeParseArguments)
  cmake_parse_arguments (ARGN "" "DIRECTORY" "TARGETS;FILES;PROGRAMS" ${ARGN})
  set (_errmsg "Options TARGETS, FILES, PROGRAMS, and DIRECTORY are mutually exclusive!")
  if (ARGN_DIRECTORY)
    if (ARGN_FILES OR ARGN_PROGRAMS OR ARGN_TARGETS)
      message (FATAL_ERROR "${_errmsg}")
    endif ()
    install (DIRECTORY ${ARGN_DIRECTORY} ${ARGN_UNPARSED_ARGUMENTS})
  elseif (ARGN_FILES)
    if (ARGN_PROGRAMS OR ARGN_TARGETS OR  ARGN_DIRECTORY)
      message (FATAL_ERROR "${_errmsg}")
    endif ()
    install (FILES ${ARGN_FILES} ${ARGN_UNPARSED_ARGUMENTS})
  elseif (ARGN_PROGRAMS)
    if (ARGN_FILES OR ARGN_TARGETS OR  ARGN_DIRECTORY)
      message (FATAL_ERROR "${_errmsg}")
    endif ()
    install (PROGRAMS ${ARGN_PROGRAMS} ${ARGN_UNPARSED_ARGUMENTS})
  elseif (ARGN_TARGETS)
    if (ARGN_FILES OR ARGN_PROGRAMS OR ARGN_DIRECTORY)
      message (FATAL_ERROR "${_errmsg}")
    endif ()
    install (TARGETS ${ARGN_TARGETS} ${ARGN_UNPARSED_ARGUMENTS})
  else ()
    install (${ARGN_UNPARSED_ARGUMENTS})
  endif ()
endfunction ()

# ----------------------------------------------------------------------------
## @brief Install content of source directory excluding typical files.
#
# Files which are excluded are typical backup files, system files, files
# from revision control systems, and CMakeLists.txt files.
#
# Example:
# @code
# basis_install_directory("${INSTALL_DATA_DIR}")
# basis_install_directory(. "${INSTALL_DATA_DIR}")
# basis_install_directory("${CMAKE_CURRENT_SOURCE_DIR}" "${INSTALL_DATA_DIR}")
# basis_install_directory(images "${INSTALL_DATA_DIR}/images")
# @endcode
#
# @param [in] ARGN The first two arguments are extracted from the beginning
#                  of this list in the named order (without option name),
#                  and the remaining arguments are passed on to CMake's
#                  <a href="http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:install">
#                  <tt>install(DIRECTORY)</tt></a> command.
# @par
# <table border="0">
#   <tr>
#     @tp @b SOURCE @endtp
#     <td>Source directory. Defaults to current source directory
#         if only one argument, the @p DESTINATION, is given./td>
#   </tr>
#   <tr>
#     @tp @b DESTINATION @endtp
#     <td>Destination directory.</td>
#   </tr>
# </table>
#
# @sa http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:install
#
# @ingroup CMakeAPI
function (basis_install_directory)
  if (ARGC EQUAL 1)
    set (SOURCE      "${CMAKE_CURRENT_SOURCE_DIR}")
    set (DESTINATION "${ARGV0}")
    set (OPTIONS     "${ARGV}")
    list (REMOVE_AT OPTIONS 0)
  elseif (ARGC GREATER 1)
    set (SOURCE      "${ARGV0}")
    set (DESTINATION "${ARGV1}")
    set (OPTIONS     "${ARGV}")
    list (REMOVE_AT OPTIONS 0 1)
  else ()
    message (FATAL_ERROR "Too few arguments given!")
  endif ()
  # check arguments
  if (NOT IS_ABSOLUTE "${SOURCE}")
    get_filename_component (SOURCE "${SOURCE}" ABSOLUTE)
  endif ()
  string (REGEX REPLACE "/+$" "" SOURCE "${SOURCE}")
  if (NOT IS_ABSOLUTE "${DESTINATION}")
    set (DESTINATION "${CMAKE_INSTALL_PREFIX}/${DESTINATION}")
  endif ()
  basis_sanitize_for_regex (PROJECT_SOURCE_DIR_RE "${PROJECT_SOURCE_DIR}")
  if ("${DESTINATION}" MATCHES "^${PROJECT_SOURCE_DIR_RE}(/|$)")
    message (FATAL_ERROR "Installation directory ${DESTINATION} is inside the project source tree!")
  endif ()
  # parse options
  set (MATCH)
  set (EXCLUDE)
  set (COMPONENT)
  set (CONFIGURATIONS)
  set (RE)
  set (OPT)
  set (NUM)
  foreach (O IN LISTS OPTIONS)
    if (O MATCHES "^(FILES_MATCHING|.*PERMISSIONS|OPTIONAL|DESTINATION|DIRECTORY)$")
      message (FATAL_ERROR "Option ${O} not supported by basis_install_directory()!")
    endif ()
    if (OPT MATCHES "^(PATTERN|REGEX)$")
      if (NUM EQUAL 1)
        if (O MATCHES "^EXCLUDE$")
          list (APPEND EXCLUDE "${RE}")
        else ()
          list (APPEND MATCH   "${RE}")
        endif ()
        set (OPT)
        set (RE)
      else ()
        if (OPT MATCHES "PATTERN")
          string (REGEX REPLACE "(^|[^\\])\\*" "&#9734" O "${O}")
          basis_sanitize_for_regex (O "${O}")
          string (REPLACE "&#9734" ".*" O "${O}")
        endif ()
        set (RE  "${O}")
        set (NUM 1)
      endif ()
    elseif (OPT MATCHES "^COMPONENT$")
      set (COMPONENT "${O}")
      set (OPT)
    elseif (OPT MATCHES "^CONFIGURATIONS$")
      if (O MATCHES "^(PATTERN|REGEX|COMPONENT|CONFIGURATIONS)$")
        set (OPT)
      else ()
        list (APPEND CONFIGURATIONS "${O}")
        math (EXPR NUM "${NUM} + 1")
      endif ()
    endif ()
    if (O MATCHES "^(PATTERN|REGEX|COMPONENT|CONFIGURATIONS)$")
      if (OPT)
        message (FATAL_ERROR "basis_install_directory(): Option ${OPT} is missing an argument!")
      endif ()
      set (OPT ${O})
      set (NUM 0)
    endif ()
  endforeach ()
  if (OPT MATCHES "CONFIGURATIONS" AND NOT NUM GREATER 0)
    message (FATAL_ERROR "basis_install_directory(): Missing argument(s) for ${OPT} option!")
  elseif (OPT MATCHES "PATTERN|REGEX")
    if (NUM EQUAL 1)
      list (APPEND MATCH "${RE}")
    else ()
      message (FATAL_ERROR "basis_install_directory(): Missing argument(s) for ${OPT} option!")
    endif ()
  elseif (OPT MATCHES "COMPONENT")
    message (FATAL_ERROR "basis_install_directory(): Missing argument for ${OPT} option!")
  endif ()
  list (APPEND EXCLUDE "CMakeLists.txt$" "/.svn/" "/.git/" ".DS_Store$" ".*~$") # default excludes
  basis_list_to_delimited_string (MATCH   "|" NOAUTOQUOTE ${MATCH})
  basis_list_to_delimited_string (EXCLUDE "|" NOAUTOQUOTE ${EXCLUDE})
  string (REPLACE "\\" "\\\\" MATCH   "${MATCH}")
  string (REPLACE "\"" "\\\"" MATCH   "${MATCH}")
  string (REPLACE "\\" "\\\\" EXCLUDE "${EXCLUDE}")
  string (REPLACE "\"" "\\\"" EXCLUDE "${EXCLUDE}")
  # Add installation instructions. Note that install(DIRECTORY) is here not
  # used to avoid the generation of empty directories
  set (OPTIONS)
  if (CONFIGURATIONS)
    list (APPEND OPTIONS CONFIGURATIONS ${CONFIGURATIONS})
  endif ()
  if (COMPONENT)
    list (APPEND OPTIONS COMPONENT ${COMPONENT})
  endif ()
  install (CODE
 "# -----------------------------------------------------------------------
  # basis_install_directory(): ${SOURCE}
  set (BASIS_INSTALL_DIRECTORY_FILES)
  set (BASIS_INSTALL_DIRECTORY_SOURCE      \"${SOURCE}\")
  set (BASIS_INSTALL_DIRECTORY_DESTINATION \"\$ENV{DESTDIR}${DESTINATION}\")
  set (BASIS_INSTALL_DIRECTORY_MATCH       \"${MATCH}\")
  set (BASIS_INSTALL_DIRECTORY_EXCLUDE     \"${EXCLUDE}\")
  file (GLOB_RECURSE BASIS_INSTALL_DIRECTORY_ALL_FILES \"\${BASIS_INSTALL_DIRECTORY_SOURCE}/*\")
  foreach (BASIS_INSTALL_DIRECTORY_FILE IN LISTS BASIS_INSTALL_DIRECTORY_ALL_FILES)
    if (NOT BASIS_INSTALL_DIRECTORY_MATCH                                            OR
            BASIS_INSTALL_DIRECTORY_FILE MATCHES \"\${BASIS_INSTALL_DIRECTORY_MATCH}\" AND
        NOT BASIS_INSTALL_DIRECTORY_FILE MATCHES \"\${BASIS_INSTALL_DIRECTORY_EXCLUDE}\")
      list (APPEND BASIS_INSTALL_DIRECTORY_FILES \"\${BASIS_INSTALL_DIRECTORY_FILE}\")
   endif ()
  endforeach ()
  foreach (BASIS_INSTALL_DIRECTORY_FILE IN LISTS BASIS_INSTALL_DIRECTORY_FILES)
    file (RELATIVE_PATH BASIS_INSTALL_DIRECTORY_FILE \"\${BASIS_INSTALL_DIRECTORY_SOURCE}\" \"\${BASIS_INSTALL_DIRECTORY_FILE}\")
    execute_process (
      COMMAND \"${CMAKE_COMMAND}\" -E compare_files
                  \"\${BASIS_INSTALL_DIRECTORY_SOURCE}/\${BASIS_INSTALL_DIRECTORY_FILE}\"
                  \"\${BASIS_INSTALL_DIRECTORY_DESTINATION}/\${BASIS_INSTALL_DIRECTORY_FILE}\"
      RESULT_VARIABLE RC
      OUTPUT_QUIET
      ERROR_QUIET
    )
    if (RC EQUAL 0)
      message (STATUS \"Up-to-date: \${BASIS_INSTALL_DIRECTORY_DESTINATION}/\${BASIS_INSTALL_DIRECTORY_FILE}\")
    else ()
      message (STATUS \"Installing: \${BASIS_INSTALL_DIRECTORY_DESTINATION}/\${BASIS_INSTALL_DIRECTORY_FILE}\")
      execute_process (
        COMMAND \"${CMAKE_COMMAND}\" -E copy_if_different
            \"\${BASIS_INSTALL_DIRECTORY_SOURCE}/\${BASIS_INSTALL_DIRECTORY_FILE}\"
            \"\${BASIS_INSTALL_DIRECTORY_DESTINATION}/\${BASIS_INSTALL_DIRECTORY_FILE}\"
        RESULT_VARIABLE RC
        OUTPUT_QUIET
        ERROR_QUIET
      )
      if (RC EQUAL 0)
        list (APPEND CMAKE_INSTALL_MANIFEST_FILES \"\${BASIS_INSTALL_DIRECTORY_DESTINATION}/\${BASIS_INSTALL_DIRECTORY_FILE}\")
      else ()
        message (STATUS \"Failed to install \${BASIS_INSTALL_DIRECTORY_DESTINATION}/\${BASIS_INSTALL_DIRECTORY_FILE}\")
      endif ()
    endif ()
  endforeach ()
  # -----------------------------------------------------------------------"
    ${OPTIONS}
  )
endfunction ()

# ----------------------------------------------------------------------------
## @brief Add installation rule for project template.
function (basis_install_template TEMPLATE DESTINATION)
  if (NOT IS_ABSOLUTE "${TEMPLATE}")
    set (TEMPLATE "${CMAKE_CURRENT_SOURCE_DIR}/${TEMPLATE}")
  endif ()
  install (
    DIRECTORY   "${TEMPLATE}/"
    DESTINATION "${DESTINATION}"
    PATTERN     *~             EXCLUDE
    PATTERN     .svn           EXCLUDE
    PATTERN     .git           EXCLUDE
    PATTERN     .hg            EXCLUDE
    PATTERN     .DS_Store      EXCLUDE
  )
endfunction ()

# ----------------------------------------------------------------------------
## @brief Add installation rule to create a symbolic link.
#
# Note that the installation rule will only be effective on a Unix-like
# system, i.e., one which supports the creation of a symbolic link.
#
# @param [in] OLD  The value of the symbolic link.
# @param [in] NEW  The name of the symbolic link.
#
# @returns Adds installation rule to create the symbolic link @p NEW.
#
# @ingroup CMakeAPI
function (basis_install_link OLD NEW)
  set (CMD_IN
    "
    set (OLD \"@OLD@\")
    set (NEW \"@NEW@\")

    if (NOT IS_ABSOLUTE \"\${OLD}\")
      set (OLD \"\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/\${OLD}\")
    endif ()
    if (NOT IS_ABSOLUTE \"\${NEW}\")
      set (NEW \"\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/\${NEW}\")
    endif ()

    if (IS_SYMLINK \"\${NEW}\")
      file (REMOVE \"\${NEW}\")
    endif ()

    if (EXISTS \"\${NEW}\")
      message (STATUS \"Skipping: \${NEW} -> \${OLD}\")
    else ()
      message (STATUS \"Installing: \${NEW} -> \${OLD}\")

      get_filename_component (SYMDIR \"\${NEW}\" PATH)

      file (RELATIVE_PATH OLD \"\${SYMDIR}\" \"\${OLD}\")

      if (NOT EXISTS \${SYMDIR})
        file (MAKE_DIRECTORY \"\${SYMDIR}\")
      endif ()

      execute_process (
        COMMAND \"${CMAKE_COMMAND}\" -E create_symlink \"\${OLD}\" \"\${NEW}\"
        RESULT_VARIABLE RETVAL
      )

      if (NOT RETVAL EQUAL 0)
        message (ERROR \"Failed to create (symbolic) link \${NEW} -> \${OLD}\")
      else ()
        list (APPEND CMAKE_INSTALL_MANIFEST_FILES \"\${NEW}\")
      endif ()
    endif ()
    "
  )

  string (CONFIGURE "${CMD_IN}" CMD @ONLY)
  install (CODE "${CMD}")
endfunction ()

# ----------------------------------------------------------------------------
## @brief Adds installation rules to create default symbolic links.
#
# This function creates for each main executable a symbolic link directly
# in the directory @c INSTALL_PREFIX/bin if @c INSTALL_SINFIX is TRUE and the
# software is installed on a Unix-like system, i.e., one which
# supports the creation of symbolic links.
#
# @returns Adds installation command for creation of symbolic links in the
#          installation tree.
function (basis_install_links)
  if (NOT UNIX)
    return ()
  endif ()

  # main executables
  basis_get_project_property (TARGETS PROPERTY TARGETS)
  foreach (TARGET_UID ${TARGETS})
    get_target_property (IMPORTED ${TARGET_UID} IMPORTED)

    if (NOT IMPORTED)
      get_target_property (BASIS_TYPE ${TARGET_UID} BASIS_TYPE)
      get_target_property (LIBEXEC    ${TARGET_UID} LIBEXEC)
      get_target_property (IS_TEST    ${TARGET_UID} TEST)

      if (BASIS_TYPE MATCHES "EXECUTABLE" AND NOT LIBEXEC AND NOT IS_TEST)
        get_target_property (SYMLINK_NAME ${TARGET_UID} SYMLINK_NAME)
        if (NOT "${SYMLINK_NAME}" MATCHES "^none$|^None$|^NONE$")
          get_target_property (SYMLINK_PREFIX ${TARGET_UID} SYMLINK_PREFIX)
          get_target_property (SYMLINK_SUFFIX ${TARGET_UID} SYMLINK_SUFFIX)
          get_target_property (INSTALL_DIR    ${TARGET_UID} RUNTIME_INSTALL_DIRECTORY)
          if (NOT INSTALL_DIR)
            get_target_property (INSTALL_DIR  ${TARGET_UID} INSTALL_DIRECTORY)
          endif ()

          basis_get_target_location (OUTPUT_NAME ${TARGET_UID} NAME)

          if (NOT SYMLINK_NAME)
            set (SYMLINK_NAME "${OUTPUT_NAME}")
          endif ()
          if (SYMLINK_PREFIX)
            set (SYMLINK_NAME "${SYMLINK_PREFIX}${SYMLINK_NAME}")
          endif ()
          if (SYMLINK_SUFFIX)
            set (SYMLINK_NAME "${SYMLINK_NAME}${SYMLINK_SUFFIX}")
          endif ()

          # avoid creation of symbolic link if there would be a conflict with
          # the subdirectory in bin/ where the actual executables are installed
          if (INSTALL_SINFIX AND "${SYMLINK_NAME}" STREQUAL "${BASIS_INSALL_SINFIX}")
            message (STATUS \"Skipping: ${INSTALL_DIR}/${OUTPUT_NAME} -> ${INSTALL_PREFIX}/bin/${SYMLINK_NAME}\")
          else ()
            basis_install_link (
              "${INSTALL_DIR}/${OUTPUT_NAME}"
              "bin/${SYMLINK_NAME}"
            )
          endif ()
        endif ()
      endif ()
    endif ()
  endforeach ()
endfunction ()

# ============================================================================
# Package registration
# ============================================================================

# ----------------------------------------------------------------------------
## @brief Register installed package with CMake.
#
# This function adds an entry to the CMake registry for packages with the
# path of the directory where the package configuration file is located in
# order to help CMake find the package.
#
# The uninstaller whose template can be found in cmake_uninstaller.cmake.in
# is responsible for removing the registry entry again.
function (basis_register_package)
  set (PKGDIR "${CMAKE_INSTALL_PREFIX}/${INSTALL_CONFIG_DIR}")
  set (PKGUID "${TOPLEVEL_PROJECT_PACKAGE_UID}")
  if (WIN32)
    install (CODE
      "execute_process (
         COMMAND reg add \"HKCU\\\\Software\\\\Kitware\\\\CMake\\\\Packages\\\\${PROJECT_PACKAGE_CONFIG_PREFIX}\" /v \"${PKGUID}\" /d \"${PKGDIR}\" /t REG_SZ /f
         RESULT_VARIABLE RT
         ERROR_VARIABLE  ERR
         OUTPUT_QUIET
       )
       if (RT EQUAL 0)
         message (STATUS \"Register:   Added HKEY_CURRENT_USER\\\\Software\\\\Kitware\\\\CMake\\\\Packages\\\\${PROJECT_PACKAGE_CONFIG_PREFIX}\\\\${PKGUID}\")
       else ()
         string (STRIP \"\${ERR}\" ERR)
         message (STATUS \"Register:   Failed to add registry entry: \${ERR}\")
       endif ()"
    )
  elseif (IS_DIRECTORY "$ENV{HOME}")
    file (WRITE "${BINARY_CONFIG_DIR}/${PROJECT_PACKAGE_CONFIG_PREFIX}RegistryFile" "${PKGDIR}")
    install (
      FILES       "${BINARY_CONFIG_DIR}/${PROJECT_PACKAGE_CONFIG_PREFIX}RegistryFile"
      DESTINATION "$ENV{HOME}/.cmake/packages/${PROJECT_PACKAGE_CONFIG_PREFIX}"
      RENAME      "${PKGUID}"
    )
  endif ()
endfunction ()

# ============================================================================
# Deinstallation
# ============================================================================

# ----------------------------------------------------------------------------
## @brief Add uninstall target.
#
# @returns Adds the custom target @c uninstall and code to
#          <tt>cmake_install.cmake</tt> to install an uninstaller.
function (basis_add_uninstall)
  # add uninstall target
  configure_file (
    ${BASIS_MODULE_PATH}/cmake_uninstall.cmake.in
    ${PROJECT_BINARY_DIR}/cmake_uninstall.cmake
    @ONLY
  )
  add_custom_target (
    uninstall
    COMMAND ${CMAKE_COMMAND} -P "${PROJECT_BINARY_DIR}/cmake_uninstall.cmake"
    COMMENT "Uninstalling..."
  )
endfunction ()


## @}
# end of Doxygen group
