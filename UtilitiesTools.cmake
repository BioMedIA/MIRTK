# ============================================================================
# Copyright (c) 2011-2012 University of Pennsylvania
# Copyright (c) 2013-2016 Andreas Schuh
# All rights reserved.
#
# See COPYING file for license information or visit
# https://cmake-basis.github.io/download.html#license
# ============================================================================

##############################################################################
# @file  UtilitiesTools.cmake
# @brief CMake functions related to configuration of BASIS utilities.
#
# @ingroup CMakeTools
##############################################################################

## @addtogroup CMakeUtilities
#  @{


# ============================================================================
# auto-detect which utilities are used
# ============================================================================

# ----------------------------------------------------------------------------
## @brief Check whether the BASIS utilities are used within a given source file.
#
# This function matches the source code against specific import patterns which
# are all valid imports of the BASIS utilities for the respective programming
# language of the specified file. If the BASIS utilities are used within the
# specified source file, the variable named @p VAR is set to @c TRUE.
# Otherwise, it is set to @c FALSE.
#
# @param [out] VAR         Whether the BASIS utilites are used.
# @param [in]  SOURCE_FILE Path of source file to check.
# @param [in]  ARGN        Source code language. If not specified, the
#                          programming language is determined automatically.
function (basis_utilities_check VAR SOURCE_FILE)
  set (UTILITIES_USED FALSE)
  if (ARGC EQUAL 2)
    basis_get_source_language (LANGUAGE "${SOURCE_FILE}")
  elseif (ARGC EQUAL 3)
    set (LANGUAGE "${ARGN}")
  else ()
    message (FATAL_ERROR "Too many arguments given!")
  endif ()
  # --------------------------------------------------------------------------
  # make file path absolute and append .in suffix if necessary
  get_filename_component (SOURCE_FILE "${SOURCE_FILE}" ABSOLUTE)
  if (NOT EXISTS "${SOURCE_FILE}" AND NOT SOURCE_FILE MATCHES "\\.in$" AND EXISTS "${SOURCE_FILE}.in")
    set (SOURCE_FILE "${SOURCE_FILE}.in")
  endif ()
  # --------------------------------------------------------------------------
  # C++
  if (LANGUAGE MATCHES "CXX")
    # read script file
    file (READ "${SOURCE_FILE}" SOURCE)
    # match use/require statements
    basis_library_prefix (PREFIX ${LANGUAGE})
    basis_sanitize_for_regex (RE "[ \\t]*#[ \\t]*include[ \\t]+[<\"](${PREFIX})?basis.h[\">]") # e.g., #include "basis.h", #include <pkg/basis.h>
    if (SCRIPT MATCHES "(^|\n)[ \t]*${RE}([ \t]*//.*|[ \t]*)(\n|$)")
      set (UTILITIES_USED TRUE)
    endif ()
  # --------------------------------------------------------------------------
  # Python/Jython
  elseif (LANGUAGE MATCHES "[JP]YTHON")
    # read script file
    file (READ "${SOURCE_FILE}" SCRIPT)
    # deprecated BASIS_PYTHON_UTILITIES macro
    if (SCRIPT MATCHES "(^|\n|;)[ \t]*\@BASIS_PYTHON_UTILITIES\@")
      message (FATAL_ERROR "Script ${SOURCE_FILE} uses the deprecated BASIS macro \@BASIS_PYTHON_UTILITIES\@!")
    endif ()
    basis_sanitize_for_regex (PYTHON_PACKAGE "${PROJECT_NAMESPACE_PYTHON}")
    # match use of package-specific utilities
    if (SCRIPT MATCHES "[^a-zA-Z._]${PYTHON_PACKAGE}.basis([.; \t\n]|$)") # e.g., basis = <package>.basis, <package>.basis.exedir()
      set (UTILITIES_USED TRUE)
    else ()
      # match import statements
      foreach (RE IN ITEMS
        "import[ \\t]+basis"                                      # e.g., import basis
        "import[ \\t]+${PYTHON_PACKAGE}\\.basis"                  # e.g., import <package>.basis
        "import[ \\t]+\\.\\.?basis"                               # e.g., import .basis, import ..basis
        "from[ \\t]+${PYTHON_PACKAGE}[ \\t]+import[ \\t]+basis"   # e.g., from <package> import basis
        "from[ \\t]+${PYTHON_PACKAGE}.basis[ \\t]+import[ \\t].*" # e.g., from <package>.basis import which
        "from[ \\t]+\\.\\.?[ \\t]+import[ \\t]+basis"             # e.g., from . import basis", "from .. import basis
        "from[ \\t]+\\.\\.?basis[ \\t]+import[ \\t].*"            # e.g., from .basis import which, WhichError, from ..basis import which
      ) # foreach RE
        if (SCRIPT MATCHES "(^|\n|;)[ \t]*${RE}([ \t]*as[ \t]+.*)?([ \t]*#.*|[ \t]*)(;|\n|$)")
          set (UTILITIES_USED TRUE)
          break ()
        endif ()
      endforeach ()
    endif ()
  # --------------------------------------------------------------------------
  # Perl
  elseif (LANGUAGE MATCHES "PERL")
    # read script file
    file (READ "${SOURCE_FILE}" SCRIPT)
    # deprecated BASIS_PERL_UTILITIES macro
    if (SCRIPT MATCHES "(^|\n|;)[ \t]*\@BASIS_PERL_UTILITIES\@")
      message (FATAL_ERROR "Script ${SOURCE_FILE} uses the deprecated BASIS macro \@BASIS_PERL_UTILITIES\@!")
    endif ()
    # match use/require statements
    basis_sanitize_for_regex (PERL_PACKAGE "${PROJECT_NAMESPACE_PERL}")
    set (RE "(use|require)[ \\t]+${PERL_PACKAGE}::Basis([ \\t]+.*)?") # e.g., use <Package>::Basis qw(:everything);
    if (SCRIPT MATCHES "(^|\n|;)[ \t]*${RE}([ \t]*#.*|[ \t]*)(;|\n|$)")
      set (UTILITIES_USED TRUE)
    endif ()
  # --------------------------------------------------------------------------
  # Bash
  elseif (LANGUAGE MATCHES "BASH")
    # read script file
    file (READ "${SOURCE_FILE}" SCRIPT)
    # deprecated BASIS_BASH_UTILITIES macro
    if (SCRIPT MATCHES "(^|\n|;)[ \t]*\@BASIS_BASH_UTILITIES\@")
      message (FATAL_ERROR "Script ${SOURCE_FILE} uses the deprecated BASIS macro \@BASIS_BASH_UTILITIES\@!")
    endif ()
    # match source/. built-ins
    set (RE "(source|\\.)[ \\t]+\\\"?\\\${BASIS_BASH_UTILITIES}?\\\"?[ \\t]*(\\|\\|.*|&&.*)?(#.*)?") # e.g., . ${BASIS_BASH_UTILITIES} || exit 1
    if (SCRIPT MATCHES "(^|\n|;)[ \t]*(${RE})[ \t]*(;|\n|$)")
      set (UTILITIES_USED TRUE)
    endif ()
  endif ()
  # return
  set (${VAR} "${UTILITIES_USED}" PARENT_SCOPE)
endfunction ()

# ============================================================================
# C++ utilities
# ============================================================================

# ----------------------------------------------------------------------------
## @brief Add build target for BASIS C++ utilities library.
#
# This function is called by basis_add_executable_target() and basis_add_library_target()
# in order to add the "basis" build target for the static project-specific
# BASIS utilities library for C++. If the target was added before, it is
# only used to get the target UID of this build target so the newly added
# executable or library can be linked to it.
#
# The CMake function add_library() checks if the specified source code files
# exist. If a source file is not found, an error is raised by CMake. The BASIS
# utilities can, however, only be configured at the end of the configuration
# step. Therefore, this function simply writes dummy C++ source files in order
# to pass the existence check. The actual source files are configured by the
# function basis_configure_utilities() which is called by basis_project_end().
#
# After writing these dummy source files, a library build target for the
# project-specific BASIS C++ utilities is added. This build target is not
# being build as part of the ALL target in case it is never used by any of
# the build targets of the project. Only if build target links to this
# library, it will be build and installed.
#
# @param [out] UID UID of added build target.
function (basis_add_cxx_utilities_library UID)
  # target UID of "basis" library target
  _basis_make_target_uid (TARGET_UID basis)
  if (NOT TARGET ${TARGET_UID})
    if (PROJECT_IS_SUBPROJECT)
      # a subproject has it's own version of the project-specific BASIS utilities
      # as the targets and functions live in a separate namespace
      set (CODE_DIR    "${BINARY_CODE_DIR}")
      set (INCLUDE_DIR "${BINARY_INCLUDE_DIR}")
      set (OUTPUT_DIR  "${BINARY_ARCHIVE_DIR}")
      set (INSTALL_DIR "${INSTALL_ARCHIVE_DIR}")
    else ()
      # modules, on the other side, share the library with the top-level project
      # the addition of the utilities target is in this case only required because
      # of the install(TARGETS) and install(EXPORT) commands.
      set (CODE_DIR    "${TOPLEVEL_BINARY_CODE_DIR}")
      set (INCLUDE_DIR "${TOPLEVEL_BINARY_INCLUDE_DIR}")
      set (OUTPUT_DIR  "${TOPLEVEL_BINARY_ARCHIVE_DIR}")
      set (INSTALL_DIR "${TOPLEVEL_INSTALL_ARCHIVE_DIR}")
    endif ()
    # write dummy source files
    basis_library_prefix (PREFIX CXX)
    foreach (S IN ITEMS basis.h basis.cxx)
      if (S MATCHES "\\.h$")
        set (S "${INCLUDE_DIR}/${PREFIX}${S}")
      else ()
        set (S "${CODE_DIR}/${S}")
      endif ()
      if (NOT EXISTS "${S}")
        file (WRITE "${S}"
          "#error This dummy source file should have been replaced by the"
          " BASIS CMake function basis_configure_utilities()"
        )
      endif ()
    endforeach ()
    # add library target if not present yet - only build if required
    add_library (${TARGET_UID} STATIC "${CODE_DIR}/basis.cxx")
    # define dependency on non-project specific utilities as the order in
    # which static libraries are listed on the command-line for the linker
    # matters; this will help CMake to get the order right
    target_link_libraries (${TARGET_UID} ${BASIS_CXX_UTILITIES_LIBRARY})
    # set target properties
    set_target_properties (
      ${TARGET_UID}
      PROPERTIES
        BASIS_TYPE                STATIC_LIBRARY
        OUTPUT_NAME               basis
        ARCHIVE_OUTPUT_DIRECTORY  "${OUTPUT_DIR}"
        ARCHIVE_INSTALL_DIRECTORY "${INSTALL_DIR}"
    )
    # export
    basis_add_export_target (EXPORT_OPT ${TARGET_UID} FALSE ${INSTALL_DIR})
    # installation
    install (
      TARGETS ${TARGET_UID} ${EXPORT_OPT}
      ARCHIVE
        DESTINATION "${INSTALL_DIR}"
        COMPONENT   "${BASIS_LIBRARY_COMPONENT}"
    )
    # debug message
    if (BASIS_DEBUG)
      message ("** Added BASIS utilities library ${TARGET_UID}")
    endif ()
  endif ()
  # done
  basis_set_project_property (PROPERTY PROJECT_USES_CXX_UTILITIES TRUE)
  set (${UID} "${TARGET_UID}" PARENT_SCOPE)
endfunction ()

# ============================================================================
# Python utilities
# ============================================================================

# ----------------------------------------------------------------------------
## @brief Add build target for BASIS Python utilities library.
#
# This function is called by basis_target_link_libraries() in order to add the
# "basis" build target for the project-specific BASIS Python utilities.
# If the target was added before, it is only used to get the target UID of
# this build target so the executable or library can be linked to it.
#
# @note The basis_target_link_libraries() function in fact calls
#       basis_add_utilities_library() which calls this function.
#
# @param [out] UID UID of added build target.
function (basis_add_python_utilities_library UID)
  basis_make_target_uid (TARGET_UID basis_py)
  if (NOT TARGET ${TARGET_UID})
    basis_library_prefix (PREFIX PYTHON)
    basis_add_library (.${TARGET_UID} "${BASIS_PYTHON_TEMPLATES_DIR}/basis.py")
    basis_set_target_properties (.${TARGET_UID} PROPERTIES PREFIX ${PREFIX})
    basis_target_link_libraries (.${TARGET_UID} ${BASIS_PYTHON_UTILITIES_LIBRARY})
  endif ()
  basis_set_project_property (PROPERTY PROJECT_USES_PYTHON_UTILITIES TRUE)
  set (${UID} ${TARGET_UID} PARENT_SCOPE)
endfunction ()

# ============================================================================
# Perl utilities
# ============================================================================

# ----------------------------------------------------------------------------
## @brief Add build target for BASIS Perl utilities library.
#
# This function is called by basis_target_link_libraries() in order to add the
# "basis" build target for the project-specific BASIS Perl utilities.
# If the target was added before, it is only used to get the target UID of
# this build target so the executable or library can be linked to it.
#
# @note The basis_target_link_libraries() function in fact calls
#       basis_add_utilities_library() which calls this function.
#
# @param [out] UID UID of added build target.
function (basis_add_perl_utilities_library UID)
  basis_make_target_uid (TARGET_UID Basis_pm)
  if (NOT TARGET ${TARGET_UID})
    basis_library_prefix (PREFIX PERL)
    basis_add_library (.${TARGET_UID} "${BASIS_PERL_TEMPLATES_DIR}/Basis.pm")
    basis_set_target_properties (.${TARGET_UID} PROPERTIES PREFIX ${PREFIX})
    basis_target_link_libraries (.${TARGET_UID} ${BASIS_PERL_UTILITIES_LIBRARY})
  endif ()
  basis_set_project_property (PROPERTY PROJECT_USES_PERL_UTILITIES TRUE)
  set (${UID} ${TARGET_UID} PARENT_SCOPE)
endfunction ()

# ============================================================================
# Bash utilities
# ============================================================================

# ----------------------------------------------------------------------------
## @brief Absolute path of current BASH file.
#
# @note Does not resolve symbolic links.
#
# Example:
# @code
# readonly __MYMODULE=@BASIS_BASH___FILE__@
# @endcode
#
# @ingroup BasisBashUtilities
set (BASIS_BASH___FILE__ "$(cd -- \"$(dirname -- \"\${BASH_SOURCE}\")\" && pwd -P)/$(basename -- \"$BASH_SOURCE\")")

# ----------------------------------------------------------------------------
## @brief Absolute path to directory of current BASH file.
#
# @note Does not resolve symbolic links.
#
# Example:
# @code
# readonly __MYMODULE_dir=@BASIS_BASH___DIR__@
# @endcode
#
# @ingroup BasisBashUtilities
set (BASIS_BASH___DIR__ "$(cd -- \"$(dirname -- \"\${BASH_SOURCE}\")\" && pwd -P)")

# ----------------------------------------------------------------------------
## @brief Add build target for BASIS Bash utilities library.
#
# This function is called by basis_target_link_libraries() in order to add the
# "basis" build target for the project-specific BASIS Bash utilities.
# If the target was added before, it is only used to get the target UID of
# this build target so the executable or library can be linked to it.
#
# @note The basis_target_link_libraries() function in fact calls
#       basis_add_utilities_library() which calls this function.
#
# @param [out] UID UID of added build target.
function (basis_add_bash_utilities_library UID)
  basis_make_target_uid (TARGET_UID basis_sh)
  if (NOT TARGET ${TARGET_UID})
    basis_library_prefix (PREFIX BASH)
    basis_add_library (.${TARGET_UID} "${BASIS_BASH_TEMPLATES_DIR}/basis.sh")
    basis_set_target_properties (.${TARGET_UID} PROPERTIES PREFIX ${PREFIX})
    basis_target_link_libraries (.${TARGET_UID} ${BASIS_BASH_UTILITIES_LIBRARY})
  endif ()
  basis_set_project_property (PROPERTY PROJECT_USES_BASH_UTILITIES TRUE)
  set (${UID} ${TARGET_UID} PARENT_SCOPE)
endfunction ()

# ============================================================================
# configuration
# ============================================================================

# ----------------------------------------------------------------------------
## @brief Add build target for BASIS utilities library.
#
# This function is called by basis_target_link_libraries() in order to add the
# "basis" build target for the project-specific BASIS utilities for a given
# source code language. If the target was added before, it is only used to get
# the target UID of this build target so the executable or library can be
# linked to it.
#
# @param [out] UID      UID of added build target.
# @param [in]  LANGUAGE Programming language of utilities.
macro (basis_add_utilities_library UID LANGUAGE)
  if ("^${LANGUAGE}$" STREQUAL "^CXX$")
    if (NOT TARGET ${BASIS_CXX_UTILITIES_LIBRARY})
      message (FATAL_ERROR "This project makes use of the BASIS C++ utilities"
                           " but BASIS was built without these utilities enabled.")
    endif ()
    basis_add_cxx_utilities_library (${UID})
  elseif ("^${LANGUAGE}$" STREQUAL "^PYTHON$")
    if (NOT TARGET ${BASIS_PYTHON_UTILITIES_LIBRARY})
      message (FATAL_ERROR "This project makes use of the BASIS Python utilities"
                           " but BASIS was built without these utilities enabled.")
    endif ()
    basis_add_python_utilities_library (${UID})
  elseif ("^${LANGUAGE}$" STREQUAL "^PERL$")
    if (NOT TARGET ${BASIS_PERL_UTILITIES_LIBRARY})
      message (FATAL_ERROR "This project makes use of the BASIS Perl utilities"
                           " but BASIS was built without these utilities enabled.")
    endif ()
    basis_add_perl_utilities_library (${UID})
  elseif ("^${LANGUAGE}$" STREQUAL "^BASH$")
    if (NOT TARGET ${BASIS_BASH_UTILITIES_LIBRARY})
      message (FATAL_ERROR "This project makes use of the BASIS Bash utilities"
                           " but BASIS was built without these utilities enabled.")
    endif ()
    basis_add_bash_utilities_library (${UID})
  else ()
    message (FATAL_ERROR "Unsupported language: ${LANGUAGE}")
  endif ()
endmacro ()

# ----------------------------------------------------------------------------
## @brief Determine whether this project uses any of the BASIS Utilities.
function (basis_get_project_uses_utilities RETVAL)
  basis_get_project_property (CXX    PROPERTY PROJECT_USES_CXX_UTILITIES)
  basis_get_project_property (PYTHON PROPERTY PROJECT_USES_PYTHON_UTILITIES)
  basis_get_project_property (PERL   PROPERTY PROJECT_USES_PERL_UTILITIES)
  basis_get_project_property (BASH   PROPERTY PROJECT_USES_BASH_UTILITIES)
  if (CXX OR PYTHON OR PERL OR BASH)
    set (RETVAL TRUE PARENT_SCOPE)
  else ()
    set (RETVAL FALSE PARENT_SCOPE)
  endif ()
endfunction ()

# ----------------------------------------------------------------------------
## @brief Configure BASIS utilities.
#
# This function configures the following source files which can be used
# within the source code of the project. If the BASIS utilities for a specific
# language are not used by any of the project's build targets, no target for
# the build of these utilities is added, unless the
# @c BUILD_BASIS_UTILITIES_FOR_<LANGUAGE> option is set to @c ON. A reason
# for forcing the build of the BASIS utilities is that the libraries should
# be used by other projects which may want to make use of the BASIS Utilities
# to get access to the project attributes.
#
# <table border="0">
#   <tr>
#     @tp @b basis.h @endtp
#     <td>Header file declaring the BASIS utilities for C++.</td>
#   </tr>
#   <tr>
#     @tp @b basis.cxx @endtp
#     <td>Definitions of the constants and functions declared in basis.h.</td>
#   </tr>
#   <tr>
#     @tp @b basis.py @endtp
#     <td>Module defining the BASIS utilities for Python.</td>
#   </tr>
#   <tr>
#     @tp @b Basis.pm @endtp
#     <td>Module defining the BASIS utilities for Perl.</td>
#   </tr>
#   <tr>
#     @tp @b basis.sh @endtp
#     <td>Module defining the BASIS utilities for Bash.</td>
#   </tr>
# </table>
#
# @note Dummy versions of the C++ source files have been written by the
#       function basis_add_utilities_library() beforehand. This is
#       necessary because CMake's add_executable() and add_library() commands
#       raise an error if any of the specified source files does not exist.
function (basis_configure_utilities)
  set (SKIP TRUE)
  foreach (L IN ITEMS CXX PYTHON PERL BASH)
    if (BUILD_BASIS_UTILITIES_FOR_${L})
      set (${L} ON)
    else ()
      basis_get_project_property (${L} PROPERTY PROJECT_USES_${L}_UTILITIES)
    endif ()
    if (${L})
      set (SKIP FALSE)
    endif ()
  endforeach ()
  if (SKIP)
    return ()
  endif ()
  message (STATUS "Configuring BASIS utilities...")
  # --------------------------------------------------------------------------
  # executable target information
  _basis_generate_executable_target_info (${CXX} ${PYTHON} ${PERL} ${BASH})
  # --------------------------------------------------------------------------
  # project ID -- used by print_version() in particular
  set (PROJECT_ID "${PROJECT_PACKAGE_NAME}")
  if (NOT PROJECT_NAME MATCHES "${PROJECT_PACKAGE_NAME_RE}")
    set (PROJECT_ID "${PROJECT_ID}, ${PROJECT_NAME}")
  endif ()
  # --------------------------------------------------------------------------
  # C++
  if (CXX)
    # make sure that library target is added which is not the case yet
    # if the BASIS C++ utilities are not used by any project target, but
    # their build is forced via the BUILD_BASIS_UTILITIES_FOR_CXX option
    basis_add_cxx_utilities_library (TARGET_UID)
    # paths - build tree
    set (BUILD_ROOT_PATH_CONFIG    "${CMAKE_BINARY_DIR}")
    set (RUNTIME_BUILD_PATH_CONFIG "${BINARY_RUNTIME_DIR}")
    set (LIBEXEC_BUILD_PATH_CONFIG "${BINARY_LIBEXEC_DIR}")
    set (LIBRARY_BUILD_PATH_CONFIG "${BINARY_LIBRARY_DIR}")
    set (DATA_BUILD_PATH_CONFIG    "${PROJECT_DATA_DIR}")
    # paths - installation
    file (RELATIVE_PATH RUNTIME_PATH_PREFIX_CONFIG "${CMAKE_INSTALL_PREFIX}/${INSTALL_RUNTIME_DIR}" "${CMAKE_INSTALL_PREFIX}")
    string (REGEX REPLACE "/$|\\$" "" RUNTIME_PATH_PREFIX_CONFIG "${RUNTIME_PATH_PREFIX_CONFIG}")
    file (RELATIVE_PATH LIBEXEC_PATH_PREFIX_CONFIG "${CMAKE_INSTALL_PREFIX}/${INSTALL_LIBEXEC_DIR}" "${CMAKE_INSTALL_PREFIX}")
    string (REGEX REPLACE "/$|\\$" "" LIBEXEC_PATH_PREFIX_CONFIG "${LIBEXEC_PATH_PREFIX_CONFIG}")
    set (RUNTIME_PATH_CONFIG "${INSTALL_RUNTIME_DIR}")
    set (LIBEXEC_PATH_CONFIG "${INSTALL_LIBEXEC_DIR}")
    set (LIBRARY_PATH_CONFIG "${INSTALL_LIBRARY_DIR}")
    set (DATA_PATH_CONFIG    "${INSTALL_DATA_DIR}")
    # namespace
    set (PROJECT_NAMESPACE_CXX_BEGIN "namespace ${PROJECT_PACKAGE_NAME_L} {")
    set (PROJECT_NAMESPACE_CXX_END   "}")
    if (PROJECT_IS_SUBPROJECT)
      set (PROJECT_NAMESPACE_CXX_BEGIN "${PROJECT_NAMESPACE_CXX_BEGIN} namespace ${PROJECT_NAME_L} {")
      set (PROJECT_NAMESPACE_CXX_END   "${PROJECT_NAMESPACE_CXX_END} }")
    endif ()
    # executable target information
    set (EXECUTABLE_TARGET_INFO "${EXECUTABLE_TARGET_INFO_CXX}")
    # configure source files
    basis_library_prefix (PREFIX CXX)
    configure_file ("${BASIS_CXX_TEMPLATES_DIR}/basis.h.in"   "${BINARY_INCLUDE_DIR}/${PREFIX}basis.h" @ONLY)
    configure_file ("${BASIS_CXX_TEMPLATES_DIR}/basis.cxx.in" "${BINARY_CODE_DIR}/basis.cxx"            @ONLY)
    source_group (BASIS FILES "${BINARY_INCLUDE_DIR}/${PREFIX}basis.h" "${BINARY_CODE_DIR}/basis.cxx")
  endif ()
  # --------------------------------------------------------------------------
  # Python
  if (PYTHON)
    # utilities available?
    if (NOT BASIS_UTILITIES_ENABLED MATCHES "PYTHON")
      message (FATAL_ERROR "BASIS Python utilities required by this package"
                           " but BASIS was built without Python utilities."
                           " Rebuild BASIS with Python utilities enabled.")
    endif ()
    # make sure that library target is added which is not the case yet
    # if the BASIS Python utilities are not used by any project target, but
    # their build is forced via the BUILD_BASIS_UTILITIES_FOR_PYTHON option
    basis_add_python_utilities_library (TARGET_UID)
    # set target properties
    set (SCRIPT_DEFINITIONS
      "set (PROJECT_ID \"${PROJECT_ID}\")
       if (BUILD_INSTALL_SCRIPT)
         set (EXECUTABLE_TARGET_INFO \"${EXECUTABLE_TARGET_INFO_PYTHON_I}\")
       else ()
         set (EXECUTABLE_TARGET_INFO \"${EXECUTABLE_TARGET_INFO_PYTHON_B}\")
       endif ()"
    )
    if ("^${PROJECT_NAME}$" STREQUAL "^BASIS$")
      set (SCRIPT_DEFINITIONS "${SCRIPT_DEFINITIONS}\nbasis_set_script_path (_BASIS_PYTHONPATH \"${BINARY_PYTHON_LIBRARY_DIR}\" \"${INSTALL_PYTHON_LIBRARY_DIR}\")")
    elseif (BUNDLE_PROJECTS MATCHES "(^|;)BASIS(;|$)")
      set (SCRIPT_DEFINITIONS "${SCRIPT_DEFINITIONS}\nbasis_set_script_path (_BASIS_PYTHONPATH \"${BASIS_PYTHONPATH}\")")
    else ()
      set (SCRIPT_DEFINITIONS "${SCRIPT_DEFINITIONS}\nset (_BASIS_PYTHONPATH \"${BASIS_PYTHONPATH}\")")
    endif ()
    basis_library_prefix (PREFIX PYTHON)
    basis_set_target_properties (
      .${TARGET_UID}
      PROPERTIES
        SOURCE_DIRECTORY          "${BASIS_PYTHON_TEMPLATES_DIR}"
        BINARY_DIRECTORY          "${BINARY_CODE_DIR}"
        LIBRARY_OUTPUT_DIRECTORY  "${BINARY_PYTHON_LIBRARY_DIR}"
        LIBRARY_INSTALL_DIRECTORY "${INSTALL_PYTHON_LIBRARY_DIR}"
        PREFIX                    "${PREFIX}"
        SCRIPT_DEFINITIONS        "${SCRIPT_DEFINITIONS}"
    )
    basis_finalize_targets (.${TARGET_UID})
  endif ()
  # --------------------------------------------------------------------------
  # Perl
  if (PERL)
    # utilities available?
    if (NOT BASIS_UTILITIES_ENABLED MATCHES "PERL")
      message (FATAL_ERROR "BASIS Perl utilities required by this package"
                           " but BASIS was built without Perl utilities."
                           " Rebuild BASIS with Perl utilities enabled.")
    endif ()
    # make sure that library target is added which is not the case yet
    # if the BASIS Perl utilities are not used by any project target, but
    # their build is forced via the BUILD_BASIS_UTILITIES_FOR_PERL option
    basis_add_perl_utilities_library (TARGET_UID)
    # set target properties
    basis_library_prefix (PREFIX PERL)
    basis_set_target_properties (
      .${TARGET_UID}
      PROPERTIES
        SOURCE_DIRECTORY          "${BASIS_PERL_TEMPLATES_DIR}"
        BINARY_DIRECTORY          "${BINARY_CODE_DIR}"
        LIBRARY_OUTPUT_DIRECTORY  "${BINARY_PERL_LIBRARY_DIR}"
        LIBRARY_INSTALL_DIRECTORY "${INSTALL_PERL_LIBRARY_DIR}"
        PREFIX                    "${PREFIX}"
        SCRIPT_DEFINITIONS
          "set (PROJECT_ID \"${PROJECT_ID}\")
           if (BUILD_INSTALL_SCRIPT)
             set (EXECUTABLE_TARGET_INFO \"${EXECUTABLE_TARGET_INFO_PERL_I}\")
           else ()
             set (EXECUTABLE_TARGET_INFO \"${EXECUTABLE_TARGET_INFO_PERL_B}\")
           endif ()"
    )
    basis_finalize_targets (.${TARGET_UID})
  endif ()
  # --------------------------------------------------------------------------
  # Bash
  if (BASH)
    # utilities available?
    if (NOT UNIX)
      message (WARNING "Package uses BASIS Bash utilities but is build"
                       " on a non-Unix system.")
    endif ()
    if (NOT BASIS_UTILITIES_ENABLED MATCHES "BASH")
      message (FATAL_ERROR "BASIS Bash utilities required by this package"
                           " but BASIS was built without Bash utilities."
                           " Rebuild BASIS with Bash utilities enabled.")
    endif ()
    # make sure that library target is added which is not the case yet
    # if the BASIS Python utilities are not used by any project target, but
    # their build is forced via the BUILD_BASIS_UTILITIES_FOR_BASH option
    basis_add_bash_utilities_library (TARGET_UID)
    # set target properties
    set (SCRIPT_DEFINITIONS
      "set (PROJECT_ID \"${PROJECT_ID}\")
       if (BUILD_INSTALL_SCRIPT)
         set (EXECUTABLE_TARGET_INFO \"${EXECUTABLE_TARGET_INFO_BASH_I}\")
       else ()
         set (EXECUTABLE_TARGET_INFO \"${EXECUTABLE_TARGET_INFO_BASH_B}\")
       endif ()
       set (EXECUTABLE_ALIASES \"${EXECUTABLE_TARGET_INFO_BASH_A}\n\n    # define short aliases for this project's targets\n    ${EXECUTABLE_TARGET_INFO_BASH_S}\")"
    )
    if ("^${PROJECT_NAME}$" STREQUAL "^BASIS$")
      set (SCRIPT_DEFINITIONS "${SCRIPT_DEFINITIONS}\nbasis_set_script_path (_BASIS_BASH_LIBRARY_DIR \"${BINARY_BASH_LIBRARY_DIR}\" \"${INSTALL_BASH_LIBRARY_DIR}\")")
    elseif (BUNDLE_PROJECTS MATCHES "(^|;)BASIS(;|$)")
      set (SCRIPT_DEFINITIONS "${SCRIPT_DEFINITIONS}\nbasis_set_script_path (_BASIS_BASH_LIBRARY_DIR \"${BASIS_BASHPATH}\")")
    else ()
      set (SCRIPT_DEFINITIONS "${SCRIPT_DEFINITIONS}\nset (_BASIS_BASH_LIBRARY_DIR \"${BASIS_BASHPATH}\")")
    endif ()
    basis_library_prefix (PREFIX BASH)
    basis_set_target_properties (
      .${TARGET_UID}
      PROPERTIES
        SOURCE_DIRECTORY          "${BASIS_BASH_TEMPLATES_DIR}"
        BINARY_DIRECTORY          "${BINARY_CODE_DIR}"
        LIBRARY_OUTPUT_DIRECTORY  "${BINARY_BASH_LIBRARY_DIR}"
        LIBRARY_INSTALL_DIRECTORY "${INSTALL_BASH_LIBRARY_DIR}"
        PREFIX                    "${PREFIX}"
        SCRIPT_DEFINITIONS        "${SCRIPT_DEFINITIONS}"
    )
    basis_finalize_targets (.${TARGET_UID})
  endif ()
  message (STATUS "Configuring BASIS utilities... - done")
endfunction ()

# ----------------------------------------------------------------------------
## @brief Generate code for initialization of executable target information.
#
# This macro generates the initialization code of the executable target
# information dictionaries for different supported programming languages.
# In case of C++, the source file has been configured and copied to the binary
# tree in a first configuration pass such that it could be used in basis_add_*()
# commands which check the existence of the arguments immediately.
# As the generation of the initialization code requires a complete list of
# build targets (cached in @c BASIS_TARGETS), this function has to be called
# after all targets have been added and finalized (in case of custom targets).
#
# @param [in] CXX    Request code for C++.
# @param [in] PYTHON Request code for Python.
# @param [in] PERL   Request code for Perl.
# @param [in] BASH   Request code for Bash.
#
# @returns Sets the following variables for each requested language.
#
# @retval EXECUTABLE_TARGET_INFO_CXX      C++ code for both build tree and installation.
# @retval EXECUTABLE_TARGET_INFO_PYTHON_B Python code for build tree.
# @retval EXECUTABLE_TARGET_INFO_PYTHON_I Python code for installation.
# @retval EXECUTABLE_TARGET_INFO_PERL_B   Perl code for build tree.
# @retval EXECUTABLE_TARGET_INFO_PERL_I   Perl code for installation.
# @retval EXECUTABLE_TARGET_INFO_BASH_B   Bash code for build tree.
# @retval EXECUTABLE_TARGET_INFO_BASH_I   Bash code for installation.
# @retval EXECUTABLE_TARGET_INFO_BASH_A   Bash code to set aliases.
# @retval EXECUTABLE_TARGET_INFO_BASH_S   Bash code to set short aliases.
function (_basis_generate_executable_target_info CXX PYTHON PERL BASH)
  # --------------------------------------------------------------------------
  if (NOT CXX AND NOT PYTHON AND NOT PERL AND NOT BASH)
    return ()
  endif ()
  # --------------------------------------------------------------------------
  # local constants
  basis_sanitize_for_regex (PROJECT_NAMESPACE_CMAKE_RE "${PROJECT_NAMESPACE_CMAKE}")
  # --------------------------------------------------------------------------
  # lists of executable targets and their location
  set (EXECUTABLE_TARGETS)
  set (EXECUTABLE_IMPORTED)
  set (BUILD_LOCATIONS)
  set (INSTALL_LOCATIONS)
  # project targets
  basis_get_project_property (TARGETS)
  foreach (TARGET_UID IN LISTS TARGETS)
    basis_get_target_type (TYPE ${TARGET_UID})
    if (TYPE MATCHES "EXECUTABLE")
      get_target_property (IMPORTED ${TARGET_UID} IMPORTED)
      basis_get_target_location (BUILD_LOCATION   ${TARGET_UID} ABSOLUTE)
      basis_get_target_location (INSTALL_LOCATION ${TARGET_UID} POST_INSTALL)
      if (BUILD_LOCATION)
        list (APPEND EXECUTABLE_TARGETS  "${TARGET_UID}")
        list (APPEND EXECUTABLE_IMPORTED "${IMPORTED}")
        list (APPEND BUILD_LOCATIONS     "${BUILD_LOCATION}")
        list (APPEND INSTALL_LOCATIONS   "${INSTALL_LOCATION}")
      else ()
        message (FATAL_ERROR "Failed to determine build location of target ${TARGET_UID}!")
      endif ()
    endif ()
  endforeach ()
  # imported targets
  basis_get_project_property (IMPORTED_TARGETS)
  basis_get_project_property (IMPORTED_TYPES)
  basis_get_project_property (IMPORTED_LOCATIONS)
  set (I 0)
  list (LENGTH IMPORTED_TARGETS N)
  while (I LESS N)
    list (GET IMPORTED_TARGETS   ${I} TARGET)
    list (GET IMPORTED_TYPES     ${I} TYPE)
    list (GET IMPORTED_LOCATIONS ${I} LOCATION)
    if (TYPE MATCHES "EXECUTABLE")
      # get corresponding UID (target may be imported from other module)
      basis_get_target_uid (TARGET_UID ${TARGET})
      # skip already considered executables
      list (FIND EXECUTABLE_TARGETS ${TARGET_UID} IDX)
      if (IDX EQUAL -1)
        if (LOCATION MATCHES "^NOTFOUND$")
          message (WARNING "Imported target ${TARGET} has no location property!")
        else ()
          list (APPEND EXECUTABLE_TARGETS  "${TARGET_UID}")
          list (APPEND EXECUTABLE_IMPORTED TRUE)
          list (APPEND BUILD_LOCATIONS     "${LOCATION}")
          list (APPEND INSTALL_LOCATIONS   "${LOCATION}")
        endif ()
      endif ()
    endif ()
    math (EXPR I "${I} + 1")
  endwhile ()
  # --------------------------------------------------------------------------
  # determine maximum length of target alias for prettier output
  set (MAX_ALIAS_LENGTH 0)
  foreach (TARGET_UID IN LISTS EXECUTABLE_TARGETS)
    basis_get_fully_qualified_target_uid (ALIAS "${TARGET_UID}")
    string (LENGTH "${ALIAS}" LENGTH)
    if (LENGTH GREATER MAX_ALIAS_LENGTH)
      set (MAX_ALIAS_LENGTH ${LENGTH})
    endif ()
  endforeach ()
  # --------------------------------------------------------------------------
  # generate source code
  set (CC)   # C++    - build tree and install tree version, constructor block
  set (PY_B) # Python - build tree version
  set (PY_I) # Python - install tree version
  set (PL_B) # Perl   - build tree version, hash entries
  set (PL_I) # Perl   - install tree version, hash entries
  set (SH_B) # Bash   - build tree version
  set (SH_I) # Bash   - install tree version
  set (SH_A) # Bash   - aliases
  set (SH_S) # Bash   - short aliases

  if (CXX)
    set (CC            "// the following code was automatically generated by the BASIS")
    set (CC "${CC}\n    // CMake function basis_configure_ExecutableTargetInfo()")
  endif ()

  set (I 0)
  list (LENGTH EXECUTABLE_TARGETS N)
  while (I LESS N)
    # ------------------------------------------------------------------------
    # get executable information
    list (GET EXECUTABLE_TARGETS ${I} TARGET_UID)
    list (GET BUILD_LOCATIONS    ${I} BUILD_LOCATION)
    list (GET INSTALL_LOCATIONS  ${I} INSTALL_LOCATION)
    # installation path (relative) to different library paths
    if (INSTALL_LOCATION)
      list (GET EXECUTABLE_IMPORTED ${I} IMPORTED)
      get_target_property (BUNDLED ${TARGET_UID} BUNDLED)
    endif ()
    foreach (L LIBRARY PYTHON_LIBRARY PERL_LIBRARY BASH_LIBRARY)
      if (INSTALL_LOCATION AND (NOT IMPORTED OR BUNDLED))
        file (
          RELATIVE_PATH INSTALL_LOCATION_REL2${L}
            "${CMAKE_INSTALL_PREFIX}/${INSTALL_${L}_DIR}/<package>"
            "${INSTALL_LOCATION}"
        )
      else ()
        set (INSTALL_LOCATION_REL2${L} "${INSTALL_LOCATION}")
      endif ()
    endforeach ()
    # target UID including project namespace
    basis_get_fully_qualified_target_uid (ALIAS "${TARGET_UID}")
    # indentation after dictionary key, i.e., alias
    string (LENGTH "${ALIAS}" ALIAS_LENGTH)
    math (EXPR NUM "${MAX_ALIAS_LENGTH} - ${ALIAS_LENGTH} + 1")
    if (NUM GREATER 1)
      string (RANDOM LENGTH ${NUM} ALPHABET " " S)
    else ()
      set (S " ")
    endif ()
    # ------------------------------------------------------------------------
    # C++
    if (CXX)
      get_filename_component (EXEC_NAME "${BUILD_LOCATION}" NAME)
      get_filename_component (BUILD_DIR "${BUILD_LOCATION}" PATH)
      if (INSTALL_LOCATION)
        get_filename_component (INSTALL_DIR "${INSTALL_LOCATION}" PATH)
      endif ()

      set (CC "${CC}\n")
      set (CC "${CC}\n    // ${TARGET_UID}")
      set (CC "${CC}\n    _exec_names  [\"${ALIAS}\"]${S}= \"${EXEC_NAME}\";")
      set (CC "${CC}\n    _build_dirs  [\"${ALIAS}\"]${S}= \"${BUILD_DIR}\";")
      set (CC "${CC}\n    _install_dirs[\"${ALIAS}\"]${S}= \"${INSTALL_DIR}\";")
    endif ()
    # ------------------------------------------------------------------------
    # Python
    if (PYTHON)
      set (PY_B "${PY_B}    '${ALIAS}':${S}'${BUILD_LOCATION}',\n")
      if (INSTALL_LOCATION)
        set (PY_I "${PY_I}    '${ALIAS}':${S}'${INSTALL_LOCATION_REL2PYTHON_LIBRARY}',\n")
      else ()
        set (PY_I "${PY_I}    '${ALIAS}':${S}'',\n")
      endif ()
    endif ()
    # ------------------------------------------------------------------------
    # Perl
    if (PERL)
      if (PL_B)
        set (PL_B "${PL_B},\n")
      endif ()
      set (PL_B "${PL_B}        '${ALIAS}'${S}=> '${BUILD_LOCATION}'")
      if (PL_I)
        set (PL_I "${PL_I},\n")
      endif ()
      if (INSTALL_LOCATION)
        set (PL_I "${PL_I}        '${ALIAS}'${S}=> '${INSTALL_LOCATION_REL2PERL_LIBRARY}'")
      else ()
        set (PL_I "${PL_I}        '${ALIAS}'${S}=> ''")
      endif ()
    endif ()
    # ------------------------------------------------------------------------
    # Bash
    if (BASH)
      # hash entry
      set (SH_B "${SH_B}\n    _basis_executabletargetinfo_add '${ALIAS}'${S}LOCATION '${BUILD_LOCATION}'")
      if (INSTALL_LOCATION)
        set (SH_I "${SH_I}\n    _basis_executabletargetinfo_add '${ALIAS}'${S}LOCATION '${INSTALL_LOCATION_REL2BASH_LIBRARY}'")
      else ()
        set (SH_I "${SH_I}\n    _basis_executabletargetinfo_add '${ALIAS}'${S}LOCATION ''")
      endif ()
      # alias
      set (SH_A "${SH_A}\n    alias '${ALIAS}'=`get_executable_path '${ALIAS}'`")
      # short alias (if target belongs to this project)
      if (ALIAS MATCHES "^${PROJECT_NAMESPACE_CMAKE_RE}\\.")
        basis_get_target_name (TARGET_NAME "${ALIAS}")
        set (SH_S "${SH_S}\n    alias '${TARGET_NAME}'='${ALIAS}'")
      endif ()
    endif ()
    # ------------------------------------------------------------------------
    # next executable target
    math (EXPR I "${I} + 1")
  endwhile ()
  # --------------------------------------------------------------------------
  # remove unnecessary leading newlines
  string (STRIP "${CC}"   CC)
  string (STRIP "${PY_B}" PY_B)
  string (STRIP "${PY_I}" PY_I)
  string (STRIP "${PL_B}" PL_B)
  string (STRIP "${PL_I}" PL_I)
  string (STRIP "${SH_B}" SH_B)
  string (STRIP "${SH_I}" SH_I)
  string (STRIP "${SH_A}" SH_A)
  string (STRIP "${SH_S}" SH_S)
  # --------------------------------------------------------------------------
  # return
  set (EXECUTABLE_TARGET_INFO_CXX      "${CC}"   PARENT_SCOPE)
  set (EXECUTABLE_TARGET_INFO_PYTHON_B "${PY_B}" PARENT_SCOPE)
  set (EXECUTABLE_TARGET_INFO_PYTHON_I "${PY_I}" PARENT_SCOPE)
  set (EXECUTABLE_TARGET_INFO_PERL_B   "${PL_B}" PARENT_SCOPE)
  set (EXECUTABLE_TARGET_INFO_PERL_I   "${PL_I}" PARENT_SCOPE)
  set (EXECUTABLE_TARGET_INFO_BASH_B   "${SH_B}" PARENT_SCOPE)
  set (EXECUTABLE_TARGET_INFO_BASH_I   "${SH_I}" PARENT_SCOPE)
  set (EXECUTABLE_TARGET_INFO_BASH_A   "${SH_A}" PARENT_SCOPE)
  set (EXECUTABLE_TARGET_INFO_BASH_S   "${SH_S}" PARENT_SCOPE)
endfunction ()


## @}
# end of Doxygen group
