# ============================================================================
# Copyright (c) 2011-2012 University of Pennsylvania
# Copyright (c) 2013-2016 Andreas Schuh
# All rights reserved.
#
# See COPYING file for license information or visit
# https://cmake-basis.github.io/download.html#license
# ============================================================================

##############################################################################
# @file  GenerateConfig.cmake
# @brief Generates package configuration files.
#
# This CMake script configures the \<package\>Config.cmake et al. files,
# once for the build tree and once for the install tree. Variables with a
# _CONFIG suffix are replaced in the default template file by either the
# value for the build or the install tree, respectively.
#
# If present, this script includes the @c PROJECT_CONFIG_DIR/ConfigBuild.cmake
# and/or @c PROJECT_CONFIG_DIR/ConfigInstall.cmake file before configuring the
# Config.cmake.in template. If a file @c PROJECT_CONFIG_DIR/Config.cmake.in
# exists, it is used as template. Otherwise, the default template file is used.
#
# Similarly, if the file @c PROJECT_CONFIG_DIR/ConfigVersion.cmake.in exists,
# it is used as template for the \<package\>ConfigVersion.cmake file. The same
# applies to ConfigUse.cmake.in.
#
# @ingroup CMakeTools
##############################################################################

basis_sanitize_for_regex (CMAKE_CURRENT_LIST_DIR_RE "${CMAKE_CURRENT_LIST_DIR}")

# ============================================================================
# names of output files
# ============================================================================

# Attention: This has to be done before configuring any files such that these
#            variables can be used by the template files.

## @addtogroup CMakeUtilities
#  @{


## @brief Package name.
set (CONFIG_PREFIX "${PROJECT_PACKAGE_CONFIG_PREFIX}")
## @brief Name of the CMake package configuration file.
set (CONFIG_FILE "${CONFIG_PREFIX}Config.cmake")
## @brief Name of the CMake package version file.
set (VERSION_FILE "${CONFIG_PREFIX}ConfigVersion.cmake")
## @brief Name of the CMake package use file.
set (USE_FILE "${CONFIG_PREFIX}Use.cmake")
## @brief Name of the CMake target exports file.
set (EXPORTS_FILE "${CONFIG_PREFIX}${BASIS_EXPORT_SUFFIX}.cmake")
## @brief Name of the CMake target exports file for custom targets.
set (CUSTOM_EXPORTS_FILE "${CONFIG_PREFIX}Custom${BASIS_EXPORT_SUFFIX}.cmake")


## @}
# end of Doxygen group

# ============================================================================
# export build targets
# ============================================================================

if (BASIS_EXPORT_ENABLED)
  basis_export_targets (
    FILE        "${EXPORTS_FILE}"
    CUSTOM_FILE "${CUSTOM_EXPORTS_FILE}"
  )
endif ()

# ============================================================================
# namespace
# ============================================================================

# code used at top of package configuration and use files to set package
# namespace prefix used for configuration variables
set (BASIS_NS
"# prefix used for variable names
set (NS \"${PROJECT_CONFIG_PREFIX}_\")

# allow caller to change namespace - used by projects with modules
if (\${NS}CONFIG_PREFIX)
  set (NS \"\${\${NS}CONFIG_PREFIX}_\")
endif ()"
)

# ============================================================================
# project configuration file
# ============================================================================

# ----------------------------------------------------------------------------
# choose template

if (EXISTS "${PROJECT_CONFIG_DIR}/Config.cmake.in")
  set (TEMPLATE "${PROJECT_CONFIG_DIR}/Config.cmake.in")
elseif (PROJECT_IS_MODULE OR PROJECT_IS_SUBMODULE OR PROJECT_IS_SUBPROJECT)
  set (TEMPLATE "${CMAKE_CURRENT_LIST_DIR}/ModuleConfig.cmake.in")
else ()
  set (TEMPLATE "${CMAKE_CURRENT_LIST_DIR}/Config.cmake.in")
endif ()

# ----------------------------------------------------------------------------
# provide code of BASIS config file as variable

if (NOT TEMPLATE MATCHES "^${CMAKE_CURRENT_LIST_DIR_RE}/")
  if (PROJECT_IS_MODULE OR PROJECT_IS_SUBMODULE OR PROJECT_IS_SUBPROJECT)
    file (READ "${CMAKE_CURRENT_LIST_DIR}/ModuleConfig.cmake.in" BASIS_TEMPLATE)
  else ()
    file (READ "${CMAKE_CURRENT_LIST_DIR}/Config.cmake.in"       BASIS_TEMPLATE)
  endif ()
  # remove file header
  string (REGEX REPLACE "^########.*########" "" BASIS_TEMPLATE "${BASIS_TEMPLATE}")
  string (STRIP "${BASIS_TEMPLATE}" BASIS_TEMPLATE)
else ()
  set (BASIS_TEMPLATE "")
endif ()

# ----------------------------------------------------------------------------
# build tree related configuration

set (BUILD_CONFIG_SETTINGS 1)
include ("${CMAKE_CURRENT_LIST_DIR}/BasisConfigSettings.cmake")
include ("${PROJECT_CONFIG_DIR}/ConfigSettings.cmake" OPTIONAL)

# ----------------------------------------------------------------------------
# configure project configuration file for build tree

string (CONFIGURE "${BASIS_TEMPLATE}" BASIS_CONFIG @ONLY)
configure_file ("${TEMPLATE}" "${BINARY_LIBCONF_DIR}/${CONFIG_FILE}" @ONLY)

if (NOT BASIS_BUILD_ONLY)

  # --------------------------------------------------------------------------
  # install tree related configuration

  set (BUILD_CONFIG_SETTINGS 0)
  include ("${CMAKE_CURRENT_LIST_DIR}/BasisConfigSettings.cmake")
  include ("${PROJECT_CONFIG_DIR}/ConfigSettings.cmake" OPTIONAL)

  # --------------------------------------------------------------------------
  # configure project configuration file for install tree
  
  string (CONFIGURE "${BASIS_TEMPLATE}" BASIS_CONFIG @ONLY)
  configure_file ("${TEMPLATE}" "${BINARY_CONFIG_DIR}/${CONFIG_FILE}" @ONLY)
  
  # --------------------------------------------------------------------------
  # install project configuration file
  
  install (
    FILES       "${BINARY_CONFIG_DIR}/${CONFIG_FILE}"
    DESTINATION "${INSTALL_CONFIG_DIR}"
  )

endif ()

# ============================================================================
# project version file
# ============================================================================
if (NOT PROJECT_IS_SUBMODULE)

  # --------------------------------------------------------------------------
  # choose template
  
  if (EXISTS "${PROJECT_CONFIG_DIR}/ConfigVersion.cmake.in")
    set (TEMPLATE "${PROJECT_CONFIG_DIR}/ConfigVersion.cmake.in")
  else ()
    set (TEMPLATE "${CMAKE_CURRENT_LIST_DIR}/ConfigVersion.cmake.in")
  endif ()
  
  # --------------------------------------------------------------------------
  # configure project configuration version file
  
  configure_file ("${TEMPLATE}" "${BINARY_LIBCONF_DIR}/${VERSION_FILE}" @ONLY)
  
  # --------------------------------------------------------------------------
  # install project configuration version file
  
  if (NOT BASIS_BUILD_ONLY)
    install (
      FILES       "${BINARY_LIBCONF_DIR}/${VERSION_FILE}"
      DESTINATION "${INSTALL_CONFIG_DIR}"
    )
  endif ()

endif ()
# ============================================================================
# project use file
# ============================================================================

# ----------------------------------------------------------------------------
# choose template

if (EXISTS "${PROJECT_CONFIG_DIR}/ConfigUse.cmake.in")
  set (TEMPLATE "${PROJECT_CONFIG_DIR}/ConfigUse.cmake.in")
elseif (EXISTS "${PROJECT_CONFIG_DIR}/Use.cmake.in")
  # backwards compatibility to version <= 0.1.8 of BASIS
  set (TEMPLATE "${PROJECT_CONFIG_DIR}/Use.cmake.in")
elseif (PROJECT_IS_MODULE OR PROJECT_IS_SUBMODULE OR PROJECT_IS_SUBPROJECT)
  set (TEMPLATE "${CMAKE_CURRENT_LIST_DIR}/ModuleConfigUse.cmake.in")
else ()
  set (TEMPLATE "${CMAKE_CURRENT_LIST_DIR}/ConfigUse.cmake.in")
endif ()

# ----------------------------------------------------------------------------
# provide code of BASIS use file as variable
if (NOT TEMPLATE MATCHES "^${CMAKE_CURRENT_LIST_DIR_RE}/")
  if (PROJECT_IS_MODULE OR PROJECT_IS_SUBMODULE OR PROJECT_IS_SUBPROJECT)
    file (READ "${CMAKE_CURRENT_LIST_DIR}/ModuleConfigUse.cmake.in" BASIS_USE)
  else ()
    file (READ "${CMAKE_CURRENT_LIST_DIR}/ConfigUse.cmake.in"       BASIS_USE)
  endif ()
  # remove file header
  string (REGEX REPLACE "^########.*########" "" BASIS_USE "${BASIS_USE}")
  string (STRIP "${BASIS_USE}" BASIS_USE)
else ()
  set (BASIS_USE "")
endif ()

# ----------------------------------------------------------------------------
# configure project use file

string (CONFIGURE "${BASIS_USE}" BASIS_USE @ONLY)
configure_file ("${TEMPLATE}" "${BINARY_LIBCONF_DIR}/${USE_FILE}" @ONLY)

# ----------------------------------------------------------------------------
# install project use file

if (NOT BASIS_BUILD_ONLY)
  install (
    FILES       "${BINARY_LIBCONF_DIR}/${USE_FILE}"
    DESTINATION "${INSTALL_CONFIG_DIR}"
  )
endif ()

unset (BASIS_NS)
unset (BASIS_TEMPLATE)
unset (BASIS_CONFIG)
unset (BASIS_USE)
unset (CMAKE_CURRENT_LIST_DIR_RE)
