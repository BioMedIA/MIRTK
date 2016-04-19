# ============================================================================
# Copyright (c) 2011-2012 University of Pennsylvania
# Copyright (c) 2013-2016 Andreas Schuh
# All rights reserved.
#
# See COPYING file for license information or visit
# https://cmake-basis.github.io/download.html#license
# ============================================================================

##############################################################################
# @file  BasisConfigSettings.cmake
# @brief Sets basic variables used in CMake package configuration.
#
# It is suggested to use @c _CONFIG as suffix for variable names that are to be
# substituted in the Config.cmake.in template file in order to distinguish
# these variables from the build configuration.
#
# @ingroup BasisSettings
##############################################################################

# ============================================================================
# common configuration settings
# ============================================================================

## @brief Include directories of dependencies.
set (INCLUDE_DIRS_CONFIG)
## @brief Directories of libraries this package depends on.
set (LIBRARY_DIRS_CONFIG)

## @brief Code to set cached &lt;Pkg&gt;_DIR variables in package configuration.
set (DEPENDS_CONFIG)

set (PKGS)
foreach (DEP IN LISTS PROJECT_DEPENDS PROJECT_OPTIONAL_DEPENDS)
  basis_tokenize_dependency ("${DEP}" PKG VER CMPS)
  if (DEPENDS_${PKG}_DIR)
    list (APPEND PKGS ${PKG})
  endif ()
endforeach ()

if (PKGS)
  list (REMOVE_DUPLICATES PKGS)
endif ()

if (PKGS)
  set (DEPENDS_CONFIG "set (_depwarn \" set to different value than during the configuration of ${PROJECT_NAME}.\""
                      "              \" Using different versions of a dependency may cause inconsistencies!\")")
endif ()
foreach (PKG IN LISTS PKGS)
  list (APPEND DEPENDS_CONFIG
    "# ${PKG}"
    "if (DEPENDS_${PKG}_DIR)"
    "  if (NOT DEPENDS_${PKG}_DIR STREQUAL \"${DEPENDS_${PKG}_DIR}\")"
    "    message (WARNING DEPENDS_${PKG}_DIR \${_depwarn})"
    "  endif ()"
    "else ()"
    "  basis_set_or_update_value (DEPENDS_${PKG}_DIR \"${DEPENDS_${PKG}_DIR}\")"
    "endif ()"
  )
endforeach ()
if (PKGS)
  list (APPEND DEPENDS_CONFIG "unset (_depwarn)")
endif ()
basis_join ("${DEPENDS_CONFIG}" "\n" DEPENDS_CONFIG)

# ============================================================================
# build tree configuration settings
# ============================================================================

if (BUILD_CONFIG_SETTINGS)
  set (INSTALL_PREFIX_CONFIG "${PROJECT_BINARY_DIR}")
  if (BUILD_EXAMPLE)
    set (EXAMPLE_DIR_CONFIG "${PROJECT_EXAMPLE_DIR}")
  else ()
    set (EXAMPLE_DIR_CONFIG)
  endif ()
  set (INCLUDE_DIR_CONFIG "${BINARY_INCLUDE_DIR};${PROJECT_INCLUDE_DIR}")
  set (LIBRARY_DIR_CONFIG "${BINARY_LIBRARY_DIR}")
  set (PYTHONPATH_CONFIG  "${BINARY_PYTHON_LIBRARY_DIR}")
  set (JYTHONPATH_CONFIG  "${BINARY_JYTHON_LIBRARY_DIR}")
  set (PERL5LIB_CONFIG    "${BINARY_PERL_LIBRARY_DIR}")
  set (MATLABPATH_CONFIG  "${BINARY_MATLAB_LIBRARY_DIR}")
  set (BASHPATH_CONFIG    "${BINARY_BASH_LIBRARY_DIR}")
  set (MODULES_DIR_CONFIG "${BINARY_LIBCONF_DIR}")
  return ()
endif ()

# ============================================================================
# installation configuration settings
# ============================================================================

basis_get_relative_path (INSTALL_PREFIX_CONFIG "${CMAKE_INSTALL_PREFIX}/${INSTALL_CONFIG_DIR}" "${CMAKE_INSTALL_PREFIX}")

## @brief Installation prefix.
set (INSTALL_PREFIX_CONFIG "\${CMAKE_CURRENT_LIST_DIR}/${INSTALL_PREFIX_CONFIG}")
## @brief Directory of example files.
if (BUILD_EXAMPLE)
  set (EXAMPLE_DIR_CONFIG "\${\${NS}INSTALL_PREFIX}/${INSTALL_EXAMPLE_DIR}")
else ()
  set (EXAMPLE_DIR_CONFIG)
endif ()
## @brief Include directories.
set (INCLUDE_DIR_CONFIG "\${\${NS}INSTALL_PREFIX}/${INSTALL_INCLUDE_DIR}")
## @brief Directory where libraries are located.
set (LIBRARY_DIR_CONFIG "\${\${NS}INSTALL_PREFIX}/${INSTALL_LIBRARY_DIR}")
## @brief Directory of Python modules.
set (PYTHONPATH_CONFIG "\${\${NS}INSTALL_PREFIX}/${INSTALL_PYTHON_LIBRARY_DIR}")
## @brief Directory of Jython modules.
set (JYTHONPATH_CONFIG "\${\${NS}INSTALL_PREFIX}/${INSTALL_JYTHON_LIBRARY_DIR}")
## @brief Directory of Perl modules.
set (PERL5LIB_CONFIG "\${\${NS}INSTALL_PREFIX}/${INSTALL_PERL_LIBRARY_DIR}")
## @brief Directory of MATLAB modules.
set (MATLABPATH_CONFIG "\${\${NS}INSTALL_PREFIX}/${INSTALL_MATLAB_LIBRARY_DIR}")
## @brief Directory of Bash modules.
set (BASHPATH_CONFIG "\${\${NS}INSTALL_PREFIX}/${INSTALL_BASH_LIBRARY_DIR}")
## @brief Directory of CMake package configuration files of project modules.
set (MODULES_DIR_CONFIG "\${\${NS}INSTALL_PREFIX}/${INSTALL_CONFIG_DIR}")
