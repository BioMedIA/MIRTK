# ============================================================================
# Copyright (c) 2011-2012 University of Pennsylvania
# Copyright (c) 2013-2016 Andreas Schuh
# All rights reserved.
#
# See COPYING file for license information or visit
# https://cmake-basis.github.io/download.html#license
# ============================================================================

##############################################################################
# @file  DirectoriesSettings.cmake
# @brief CMake variables of project directories.
#
# This file configures the project directory structure as defined by the
# Filesystem Hierarchy Standard for BASIS packages.
#
# @sa https://cmake-basis.github.io/standard/fhs/
#
# The project must follow the directory structure as defined by the
# <tt>PROJECT_&lt;*&gt;_DIR</tt> variables.
#
# Ideally, when changing the name of one of these directories, only the
# directory structure of the template needs to be updated. The BASIS CMake
# functions should not be required to change as they are supposed to use these
# variables instead of the actual names. Any change of the project directory
# structure has to be made with care, however, and backwards compatibility to
# previous releases of BASIS shall be maintained. Consider the use of the
# @c TEMPLATE_VERSION if required.
#
# @note The documentation of the variables can be found in Directories.cmake.in.
##############################################################################


# ============================================================================
# root directories of top-level project
# ============================================================================

if (NOT PROJECT_IS_MODULE)
  set (TOPLEVEL_PROJECT_SOURCE_DIR "${PROJECT_SOURCE_DIR}")
  set (TOPLEVEL_PROJECT_BINARY_DIR "${PROJECT_BINARY_DIR}")
endif ()

# ============================================================================
# local variables
# ============================================================================

if (BUNDLE_NAME AND NOT BUNDLE_NAME MATCHES "${PROJECT_PACKAGE_NAME_RE}")
  set (_BUNDLE "/${BUNDLE_NAME}")
else ()
  set (_BUNDLE)
endif ()
if (PROJECT_PACKAGE_VENDOR)
  set (_VENDOR "/${PROJECT_PACKAGE_VENDOR}")
else ()
  set (_VENDOR)
endif ()
set (_PACKAGE "/${PROJECT_PACKAGE_NAME}")
if (PROJECT_IS_SUBPROJECT)
  set (_MODULE "/${PROJECT_NAME}")
else ()
  set (_MODULE)
endif ()
if (UNIX)
  string (TOLOWER "${_VENDOR}"  _VENDOR)
  string (TOLOWER "${_PACKAGE}" _PACKAGE)
  string (TOLOWER "${_MODULE}"  _MODULE)
  string (TOLOWER "${_BUNDLE}"  _BUNDLE)
endif ()

# ============================================================================
# directories of site packages of script interpreters
# ============================================================================

# ----------------------------------------------------------------------------
set (PYTHON_SITELIB)
if (PYTHON_EXECUTABLE)
  execute_process (
    COMMAND "${PYTHON_EXECUTABLE}" "${BASIS_MODULE_PATH}/get_python_lib.py"
    RESULT_VARIABLE _RV
    OUTPUT_VARIABLE PYTHON_SITELIB
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
  )
  if (NOT _RV EQUAL 0)
    set (PYTHON_SITELIB)
  endif ()
endif ()
# ----------------------------------------------------------------------------
set (JYTHON_SITELIB)
if (JYTHON_EXECUTABLE)
  execute_process (
    COMMAND "${JYTHON_EXECUTABLE}" "${BASIS_MODULE_PATH}/get_python_lib.py"
    RESULT_VARIABLE _RV
    OUTPUT_VARIABLE JYTHON_SITELIB
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
  )
  if (NOT _RV EQUAL 0)
    set (JYTHON_SITELIB)
  endif ()
endif ()
# ----------------------------------------------------------------------------
# Perl
find_package (PerlLibs QUIET)
if (NOT PerlLibs_FOUND)
  set (PERL_SITELIB)
endif ()

# ============================================================================
# testing tree
# ============================================================================

set (TESTING_OUTPUT_DIR  "${TOPLEVEL_PROJECT_BINARY_DIR}/Testing/Temporary${_MODULE}")
set (TESTING_RUNTIME_DIR "${TOPLEVEL_PROJECT_BINARY_DIR}/Testing/bin${_MODULE}")
set (TESTING_LIBEXEC_DIR "${TOPLEVEL_PROJECT_BINARY_DIR}/Testing/lib${_MODULE}")
set (TESTING_LIBRARY_DIR "${TOPLEVEL_PROJECT_BINARY_DIR}/Testing/lib${_MODULE}")
set (TESTING_ARCHIVE_DIR "${TOPLEVEL_PROJECT_BINARY_DIR}/Testing/lib${_MODULE}")

foreach (_L IN ITEMS python jython perl matlab bash)
  string (TOUPPER "${_L}" _U)
  set (TESTING_${_U}_LIBRARY_DIR "${TOPLEVEL_PROJECT_BINARY_DIR}/Testing/lib/${_L}")
endforeach ()

# ============================================================================
# build tree
# ============================================================================

# set directories corresponding to the source tree directories
foreach (_P CODE CONFIG DATA DOC EXAMPLE MODULES TESTING)
  basis_get_relative_path (_D "${PROJECT_SOURCE_DIR}" "${PROJECT_${_P}_DIR}")
  set (BINARY_${_P}_DIR "${PROJECT_BINARY_DIR}/${_D}")
endforeach ()

set (BINARY_INCLUDE_DIR "${TOPLEVEL_PROJECT_BINARY_DIR}/include")
set (BINARY_RUNTIME_DIR "${TOPLEVEL_PROJECT_BINARY_DIR}/bin")
set (BINARY_LIBEXEC_DIR "${TOPLEVEL_PROJECT_BINARY_DIR}/lib${_MODULE}")
set (BINARY_LIBRARY_DIR "${TOPLEVEL_PROJECT_BINARY_DIR}/lib${_MODULE}")
set (BINARY_ARCHIVE_DIR "${TOPLEVEL_PROJECT_BINARY_DIR}/lib${_MODULE}")

if (WIN32)
  set (BINARY_LIBCONF_DIR "${TOPLEVEL_PROJECT_BINARY_DIR}/CMake")
else ()
  set (BINARY_LIBCONF_DIR "${TOPLEVEL_PROJECT_BINARY_DIR}/lib/cmake${_PACKAGE}")
endif ()

foreach (_L IN ITEMS python jython perl matlab bash)
  string (TOUPPER "${_L}" _U)
  set (BINARY_${_U}_LIBRARY_DIR "${TOPLEVEL_PROJECT_BINARY_DIR}/lib/${_L}")
endforeach ()

# set default CMake variables which are, however, not used by BASIS
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${BINARY_RUNTIME_DIR}")
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY "${BINARY_LIBRARY_DIR}")
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${BINARY_ARCHIVE_DIR}")

# ============================================================================
# install tree
# ============================================================================

# Attention: In order for CPack to work correctly, the destination paths have
#            to be given relative to CMAKE_INSTALL_PREFIX. Therefore, this
#            prefix must be excluded from the following paths!

# ----------------------------------------------------------------------------
# default installation prefix
string (REGEX REPLACE "[\\/]+$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")
# change default installation prefix used by CMake
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT OR NOT CMAKE_INSTALL_PREFIX)
  # <ProgramFilesDir>/<Vendor>/<Package>[-<version>]
  if (WIN32)
    get_filename_component (CMAKE_INSTALL_PREFIX "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\Windows\\CurrentVersion;ProgramFilesDir]" ABSOLUTE)
    if (NOT CMAKE_INSTALL_PREFIX OR CMAKE_INSTALL_PREFIX MATCHES "/registry")
      set (CMAKE_INSTALL_PREFIX "C:/Program Files")
    endif ()
    set (CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}${_VENDOR}${_PACKAGE}")
    if (NOT PROJECT_VERSION MATCHES "^0\\.0\\.0$")
      set (CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}-${PROJECT_VERSION}")
    endif ()
  # /opt/<vendor>/<package>[-<version>]
  else ()
    set (CMAKE_INSTALL_PREFIX "/opt${_VENDOR}${_PACKAGE}")
    if (NOT PROJECT_VERSION MATCHES "^0\\.0\\.0$")
      set (CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}-${PROJECT_VERSION}")
    endif ()
  endif ()
endif ()
set (CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}" CACHE PATH "Installation prefix." FORCE)

# ----------------------------------------------------------------------------
# installation scheme - non-cached, can be preset using -D option of CMake
set (BASIS_INSTALL_SCHEME "default" CACHE STRING "default, opt, usr, win")
set_property(CACHE BASIS_INSTALL_SCHEME PROPERTY STRINGS default opt usr win)
mark_as_advanced (BASIS_INSTALL_SCHEME)

if (BASIS_INSTALL_SCHEME MATCHES "default")
  string (TOLOWER "${CMAKE_INSTALL_PREFIX}" _CMAKE_INSTALL_PREFIX_L)
  basis_sanitize_for_regex (_BUNDLE_NAME_RE "${BUNDLE_NAME}")
  string (TOLOWER "{_BUNDLE_NAME_RE}" _BUNDLE_NAME_RE_L)
  string (TOUPPER "{_BUNDLE_NAME_RE}" _BUNDLE_NAME_RE_U)
  if (WIN32)
    set (BASIS_INSTALL_SCHEME win)
  elseif (NOT _BUNDLE AND _CMAKE_INSTALL_PREFIX_L MATCHES "/(.*[_-])?(${PROJECT_PACKAGE_NAME_RE}|${PROJECT_PACKAGE_NAME_RE_L}|${PROJECT_PACKAGE_NAME_RE_U})[_-]?") # e.g. /opt/<package>[-<version>]
    set (BASIS_INSTALL_SCHEME opt)
  elseif (_BUNDLE AND _CMAKE_INSTALL_PREFIX_L MATCHES "/(.*[_-])?(${_BUNDLE_NAME_RE}|${_BUNDLE_NAME_RE_L}|${_BUNDLE_NAME_RE_U})[_-]?") # e.g. /opt/<bundle>[-<version>]
    set (BASIS_INSTALL_SCHEME opt)
  else ()
    set (BASIS_INSTALL_SCHEME usr)
  endif ()
  unset (_CMAKE_INSTALL_PREFIX_L)
  unset (_BUNDLE_NAME_RE)
  unset (_BUNDLE_NAME_RE_L)
  unset (_BUNDLE_NAME_RE_U)
endif ()

if (NOT BASIS_INSTALL_SCHEME MATCHES "^(opt|usr|win|bundle)$")
  message (FATAL_ERROR "Invalid BASIS_INSTALL_SCHEME! Valid values are 'default', 'opt', 'usr', 'win'.")
endif ()

# ----------------------------------------------------------------------------
# installation directories
if (BASIS_INSTALL_SCHEME MATCHES "win") # e.g., CMAKE_INSTALL_PREFIX := <ProgramFilesDir>/<Vendor>/<Package>

  # --------------------------------------------------------------------------
  # bundled dependency
  if (_BUNDLE)
    # package configuration
    set (INSTALL_CONFIG_DIR  "CMake")
    # executables
    set (INSTALL_RUNTIME_DIR "Lib${_PACKAGE}${_MODULE}")
    set (INSTALL_LIBEXEC_DIR "Lib${_PACKAGE}${_MODULE}")
    # libraries
    set (INSTALL_INCLUDE_DIR "Include")
    set (INSTALL_LIBRARY_DIR "Lib${_PACKAGE}${_MODULE}")
    set (INSTALL_ARCHIVE_DIR "Lib${_PACKAGE}${_MODULE}")
    # shared data
    set (INSTALL_SHARE_DIR   "Share${_PACKAGE}${_MODULE}")
    set (INSTALL_DATA_DIR    "Data${_PACKAGE}${_MODULE}")
    set (INSTALL_EXAMPLE_DIR "Example${_PACKAGE}${_MODULE}")
    # documentation
    set (INSTALL_DOC_DIR     "Doc${_PACKAGE}${_MODULE}")
    set (INSTALL_MAN_DIR)
    set (INSTALL_TEXINFO_DIR)
  # --------------------------------------------------------------------------
  # main package
  else ()
    # package configuration
    set (INSTALL_CONFIG_DIR  "CMake")
    # executables
    set (INSTALL_RUNTIME_DIR "Bin")
    set (INSTALL_LIBEXEC_DIR "Lib${_MODULE}")
    # libraries
    set (INSTALL_INCLUDE_DIR "Include")
    set (INSTALL_LIBRARY_DIR "Lib${_MODULE}")
    set (INSTALL_ARCHIVE_DIR "Lib${_MODULE}")
    # shared data
    set (INSTALL_SHARE_DIR   "Share${_MODULE}")
    set (INSTALL_DATA_DIR    "Data${_MODULE}")
    set (INSTALL_EXAMPLE_DIR "Example${_MODULE}")
    # documentation
    set (INSTALL_DOC_DIR     "Doc${_MODULE}")
    set (INSTALL_MAN_DIR)
    set (INSTALL_TEXINFO_DIR)
  endif ()

elseif (BASIS_INSTALL_SCHEME MATCHES "usr") # e.g., CMAKE_INSTALL_PREFIX := /usr/local

  # --------------------------------------------------------------------------
  # bundled dependency
  if (_BUNDLE)
    # package configuration
    set (INSTALL_CONFIG_DIR  "lib/cmake${_BUNDLE}")
    # executables
    set (INSTALL_RUNTIME_DIR "lib${_BUNDLE}${_PACKAGE}${_MODULE}")
    set (INSTALL_LIBEXEC_DIR "lib${_BUNDLE}${_PACKAGE}${_MODULE}")
    # libraries
    set (INSTALL_INCLUDE_DIR "include")
    set (INSTALL_LIBRARY_DIR "lib${_BUNDLE}${_PACKAGE}${_MODULE}")
    set (INSTALL_ARCHIVE_DIR "lib${_BUNDLE}${_PACKAGE}${_MODULE}")
    # shared data
    set (INSTALL_SHARE_DIR   "share${_BUNDLE}${_PACKAGE}${_MODULE}")
    set (INSTALL_DATA_DIR    "share${_BUNDLE}${_PACKAGE}${_MODULE}/data")
    set (INSTALL_EXAMPLE_DIR "share${_BUNDLE}${_PACKAGE}${_MODULE}/example")
    # documentation
    set (INSTALL_DOC_DIR     "doc${_BUNDLE}${_PACKAGE}${_MODULE}")
    set (INSTALL_MAN_DIR     "share${_BUNDLE}${_PACKAGE}${_MODULE}/man")
    set (INSTALL_TEXINFO_DIR "share${_BUNDLE}${_PACKAGE}${_MODULE}/info")
  # --------------------------------------------------------------------------
  # main package
  else ()
    # package configuration
    set (INSTALL_CONFIG_DIR  "lib/cmake${_PACKAGE}")
    # executables
    set (INSTALL_RUNTIME_DIR "bin")
    set (INSTALL_LIBEXEC_DIR "lib${_PACKAGE}${_MODULE}")
    # libraries
    set (INSTALL_INCLUDE_DIR "include")
    set (INSTALL_LIBRARY_DIR "lib${_PACKAGE}${_MODULE}")
    set (INSTALL_ARCHIVE_DIR "lib${_PACKAGE}${_MODULE}")
    # shared data
    set (INSTALL_SHARE_DIR   "share${_PACKAGE}${_MODULE}")
    set (INSTALL_DATA_DIR    "share${_PACKAGE}${_MODULE}/data")
    set (INSTALL_EXAMPLE_DIR "share${_PACKAGE}${_MODULE}/example")
    # documentation
    set (INSTALL_DOC_DIR     "doc${_PACKAGE}${_MODULE}")
    set (INSTALL_MAN_DIR     "share/man")
    set (INSTALL_TEXINFO_DIR "share/info")
  endif ()

else () # e.g., CMAKE_INSTALL_PREFIX := /opt/<vendor>/<package>

  # --------------------------------------------------------------------------
  # bundled dependency
  if (_BUNDLE)
    # package configuration
    set (INSTALL_CONFIG_DIR "lib/cmake${_BUNDLE}")
    # executables
    set (INSTALL_RUNTIME_DIR "lib${_PACKAGE}${_MODULE}")
    set (INSTALL_LIBEXEC_DIR "lib${_PACKAGE}${_MODULE}")
    # libraries
    set (INSTALL_INCLUDE_DIR "include")
    set (INSTALL_LIBRARY_DIR "lib${_PACKAGE}${_MODULE}")
    set (INSTALL_ARCHIVE_DIR "lib${_PACKAGE}${_MODULE}")
    # shared data
    set (INSTALL_SHARE_DIR   "share${_PACKAGE}${_MODULE}")
    set (INSTALL_DATA_DIR    "share${_PACKAGE}${_MODULE}/data")
    set (INSTALL_EXAMPLE_DIR "share${_PACKAGE}${_MODULE}/example")
    # documentation
    set (INSTALL_DOC_DIR     "doc${_PACKAGE}${_MODULE}")
    set (INSTALL_MAN_DIR     "share${_PACKAGE}${_MODULE}/man")
    set (INSTALL_TEXINFO_DIR "share${_PACKAGE}${_MODULE}/info")
  else ()
    # package configuration
    set (INSTALL_CONFIG_DIR  "lib/cmake${_PACKAGE}")
    # executables
    set (INSTALL_RUNTIME_DIR "bin")
    set (INSTALL_LIBEXEC_DIR "lib${_MODULE}")
    # libraries
    set (INSTALL_INCLUDE_DIR "include")
    set (INSTALL_LIBRARY_DIR "lib${_MODULE}")
    set (INSTALL_ARCHIVE_DIR "lib${_MODULE}")
    # shared data
    set (INSTALL_SHARE_DIR   "share${_MODULE}")
    set (INSTALL_DATA_DIR    "share${_MODULE}/data")
    set (INSTALL_EXAMPLE_DIR "share${_MODULE}/example")
    # documentation
    set (INSTALL_DOC_DIR     "doc${_MODULE}")
    set (INSTALL_MAN_DIR     "man")
    set (INSTALL_TEXINFO_DIR "info")
  endif ()

endif ()

# ----------------------------------------------------------------------------
# private script libraries
#
# The modules of script libraries which are only intended for use by this
# package itself are installed within the package own installation
# prefix/subdirectories.
if (BASIS_INSTALL_SCHEME MATCHES "win")

  foreach (_L IN ITEMS Python Jython Perl Matlab Bash)
    string (TOUPPER "${_L}" _U)
    if (BASIS_COMPILE_SCRIPTS)
      if (_U MATCHES "PERL")
        set (INSTALL_${_U}_LIBRARY_DIR "Lib/Perl5")
      elseif (NOT _U MATCHES "MATLAB|BASH" AND ${_U}_VERSION_MAJOR AND DEFINED ${_U}_VERSION_MINOR)
        set (INSTALL_${_U}_LIBRARY_DIR "Lib/${_L}${${_U}_VERSION_MAJOR}.${${_U}_VERSION_MINOR}")
      else ()
        set (INSTALL_${_U}_LIBRARY_DIR "Lib/${_L}")
      endif ()
    else ()
      set (INSTALL_${_U}_LIBRARY_DIR "Lib/${_L}")
    endif ()
  endforeach ()

elseif (BASIS_INSTALL_SCHEME MATCHES "usr")

  if (_BUNDLE)
    set (_P "${_BUNDLE}")
  else ()
    set (_P "${_PACKAGE}")
  endif ()

  foreach (_L IN ITEMS python jython perl matlab bash)
    string (TOUPPER "${_L}" _U)
    if (BASIS_COMPILE_SCRIPTS)
      if (_U MATCHES "PERL")
        set (INSTALL_${_U}_LIBRARY_DIR "lib${_P}/perl5")
      elseif (NOT _U MATCHES "MATLAB|BASH" AND ${_U}_VERSION_MAJOR AND DEFINED ${_U}_VERSION_MINOR)
        set (INSTALL_${_U}_LIBRARY_DIR "lib${_P}/${_L}${${_U}_VERSION_MAJOR}.${${_U}_VERSION_MINOR}")
      else ()
        set (INSTALL_${_U}_LIBRARY_DIR "lib${_P}/${_L}")
      endif ()
    else ()
      set (INSTALL_${_U}_LIBRARY_DIR "lib${_P}/${_L}")
    endif ()
  endforeach ()
 
else () # opt

  foreach (_L IN ITEMS python jython perl matlab bash)
    string (TOUPPER "${_L}" _U)
    if (BASIS_COMPILE_SCRIPTS)
      if (_U MATCHES "PERL")
        set (INSTALL_${_U}_LIBRARY_DIR "lib/perl5")
      elseif (NOT _U MATCHES "MATLAB|BASH" AND ${_U}_VERSION_MAJOR AND DEFINED ${_U}_VERSION_MINOR)
        set (INSTALL_${_U}_LIBRARY_DIR "lib/${_L}${${_U}_VERSION_MAJOR}.${${_U}_VERSION_MINOR}")
      else ()
        set (INSTALL_${_U}_LIBRARY_DIR "lib/${_L}")
      endif ()
    else ()
      set (INSTALL_${_U}_LIBRARY_DIR "lib/${_L}")
    endif ()
  endforeach ()

endif ()

# ----------------------------------------------------------------------------
# public script libraries
#
# The modules of script libraries which are intended for use by external packages
# are installed in the respective installation directories of the particular
# interpreter. For example, in case of Python, the public Python modules are
# installed in the site-packages directory of the found PYTHON_EXECUTABLE.
# In particular the modules in the PROJECT_LIBRARY_DIR are intended for use
# by external packages. Other modules added using the basis_add_script_library()
# and basis_add_script() CMake functions are by default considered to be intended
# for internal use by the other modules and executable scripts.
#
# Note: For those interpreters of scripting languages which by themselves do
#       not define a common installation directory for site packages, the
#       installation directory for public modules may be identical to the
#       one for private modules. Moreover, the user has the option to disable
#       the installation of public modules in the system default site directories
#       in order to prevent the installation of files outside the CMAKE_INSTALL_PREFIX.

# reset directories if BASIS_INSTALL_SITE_PACKAGES option has been changed
if (DEFINED _BASIS_INSTALL_SITE_PACKAGES)
  set (_RESET FALSE)
  if (BASIS_INSTALL_SITE_PACKAGES AND NOT _BASIS_INSTALL_SITE_PACKAGES)
    set (_RESET TRUE)
  elseif (NOT BASIS_INSTALL_SITE_PACKAGES AND _BASIS_INSTALL_SITE_PACKAGES)
    set (_RESET TRUE)
  endif ()
  if (_RESET)
    foreach (_L IN ITEMS PYTHON JYTHON PERL)
      # do not reset if BASIS_INSTALL_SITE_PACKAGES is OFF and path is already relative
      if (IS_ABSOLUTE "${INSTALL_${_L}_SITE_DIR}" OR BASIS_INSTALL_SITE_PACKAGES)
        basis_update_value (INSTALL_${_L}_SITE_DIR)
      endif ()
    endforeach ()
  endif ()
  unset (_RESET)
endif ()
set (_BASIS_INSTALL_SITE_PACKAGES "${BASIS_INSTALL_SITE_PACKAGES}" CACHE INTERNAL "Previous value of BASIS_INSTALL_SITE_PACKAGES option." FORCE)

# try to determine default installation directories
if (BASIS_INSTALL_SITE_PACKAGES)
  # Python
  if (NOT INSTALL_PYTHON_SITE_DIR)
    set (INSTALL_PYTHON_SITE_DIR "${PYTHON_SITELIB}")
  endif ()
  # Jython
  if (NOT INSTALL_JYTHON_SITE_DIR)
    set (INSTALL_JYTHON_SITE_DIR "${JYTHON_SITELIB}")
  endif ()
  # Perl
  if (NOT INSTALL_PERL_SITE_DIR)
    set (INSTALL_PERL_SITE_DIR "${PERL_SITELIB}")
  endif ()
endif ()

# if it failed to determine the default installation directories by executing some
# code or command, use the directories used for private libraries instead
foreach (_U IN ITEMS PYTHON JYTHON PERL MATLAB BASH)
  if (NOT INSTALL_${_U}_SITE_DIR)
    set (INSTALL_${_U}_SITE_DIR "${INSTALL_${_U}_LIBRARY_DIR}")
  endif ()
endforeach ()

# cache directories - also so users can edit them
foreach (_L IN ITEMS Python Jython Perl MATLAB Bash)
  string (TOUPPER "${_L}" _U)
  set (INSTALL_${_U}_SITE_DIR "${INSTALL_${_U}_SITE_DIR}" CACHE PATH "Installation directory for public ${_L} modules." FORCE)
  mark_as_advanced (INSTALL_${_U}_SITE_DIR)
endforeach ()

# ============================================================================
# top-level references
# ============================================================================

if (NOT PROJECT_IS_MODULE)
  # source tree
  foreach (_D CODE CONFIG DATA DOC EXAMPLE INCLUDE MODULES TESTING)
    set (TOPLEVEL_PROJECT_${_D}_DIR "${PROJECT_${_D}_DIR}")
  endforeach ()
  # build tree
  foreach (_D CODE CONFIG DATA DOC EXAMPLE INCLUDE MODULES TESTING RUNTIME LIBEXEC LIBRARY ARCHIVE LIBCONF)
    set (TOPLEVEL_BINARY_${_D}_DIR "${BINARY_${_D}_DIR}")
  endforeach ()
  foreach (_L IN ITEMS PYTHON JYTHON PERL MATLAB BASH)
    set (TOPLEVEL_BINARY_${_L}_LIBRARY_DIR "${BINARY_${_L}_LIBRARY_DIR}")
  endforeach ()
  # installation
  foreach (_D IN ITEMS CONFIG INCLUDE RUNTIME LIBEXEC LIBRARY ARCHIVE DATA DOC EXAMPLE SHARE)
    set (TOPLEVEL_INSTALL_${_D}_DIR "${INSTALL_${_D}_DIR}")
  endforeach ()
  foreach (_L IN ITEMS PYTHON JYTHON PERL MATLAB BASH)
    set (TOPLEVEL_INSTALL_${_L}_LIBRARY_DIR "${INSTALL_${_L}_LIBRARY_DIR}")
  endforeach ()
endif ()

# ============================================================================
# clean up
# ============================================================================

unset (_D)
unset (_L)
unset (_U)
unset (_P)
unset (_RV)
unset (_VENDOR)
unset (_PACKAGE)
unset (_MODULE)
unset (_DEFAULT_SCHEME)


## @}
# end of Doxygen group
