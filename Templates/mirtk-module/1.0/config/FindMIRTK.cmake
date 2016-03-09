#.rst:
# FindMIRTK
# ---------
#
# Find the Medical Image Registration ToolKit (MIRTK).
#
# This module looks for a MIRTK installation using find_package in CONFIG mode.
# 
# This module does not itself define any :prop_tgt:`IMPORTED` target as these
# are defined by the ``MIRTKTargets.cmake`` file included by ``MIRTKConfig.cmake``
# if MIRTK has been found. For example, the Common library of the MIRTK is
# imported as target ``mirtk::LibCommon``. Similar for other modules.
# When a required MIRTK module was not found, the ``MIRTKConfig.cmake`` will
# report an error with the names of the required, optional, and missing modules.
#
# Result variables::
#
#   MIRTK_FOUND            - True if headers and required components were found.
#   MIRTK_<Module>_FOUND   - True if requested module was found.
#   MIRTK_VERSION          - Version of found MIRTK libraries.
#   MIRTK_VERSION_MAJOR    - Major version number of found MIRTK libraries.
#   MIRTK_VERSION_MINOR    - Minor version number of found MIRTK libraries.
#   MIRTK_VERSION_PATCH    - Patch number of found MIRTK libraries.
#   MIRTK_INCLUDE_DIRS     - MIRTK include directories.
#   MIRTK_LIBRARY_DIRS     - Link directories for MIRTK libraries.
#   MIRTK_MODULES          - List of official MIRTK modules.
#   MIRTK_MODULES_ENABLED  - List of installed MIRTK modules.
#   MIRTK_MODULES_FOUND    - List of requested and found MIRTK modules.
#   MIRTK_MODULES_NOTFOUND - List of not found optional MIRTK modules.
#
# By default, this module reads hints about search paths from variables::
#
#   DEPENDS_MIRTK_DIR - Either installation root or MIRTKConfig.cmake directory.
#   MIRTK_DIR         - Directory containing the MIRTKConfig.cmake file.
#   MIRTK_ROOT        - Root directory of MIRTK installation.
#
# This module considers the common ``MIRTK_DIR`` and ``MIRTK_ROOT`` CMake or environment
# variables to initialize the ``DEPENDS_MIRTK_DIR`` cache entry. The ``DEPENDS_MIRTK_DIR``
# is the non-internal cache entry visible in the CMake GUI. It is marked advanced
# when the MIRTK was found, and non-advanced otherwise. This variable can be set
# to either the installation prefix of MIRTK, i.e., the top-level directory, or the
# directory containing the ``MIRTKConfig.cmake`` file. It therefore is a hybrid of
# ``MIRTK_ROOT`` and ``MIRTK_DIR`` and replaces these. The common DEPENDS prefix
# for cache entries used to set the location of dependencies allows the grouping
# of these variables in the CMake GUI. This is a feature of the CMake BASIS
# basis_find_package command. As this command is not available without having found
# a MIRTK installation before, this module can be used to replicate a subset of
# the basis_find_package functionality for finding MIRTK.
#
# A custom name or prefix for the cache path entry can be set using the variables::
#
#   MIRTK_CACHE_PATH_NAME   - Name of user search path cache entry.
#   MIRTK_CACHE_PATH_PREFIX - Prefix of user search path cache entry prepended to "MIRTK_DIR".
#   MIRTK_CACHE_PATH_CMAKE  - If true, use CMake's default MIRTK_DIR cache entry.
#
# For example, to use ``MIRTK_DIR`` as cache variable::
#
#   set(MIRTK_CACHE_PATH_CMAKE TRUE)
#   find_package(MIRTK REQUIRED)
#
# or::
#
#   set(MIRTK_CACHE_PATH_NAME MIRTK_DIR)
#   find_package(MIRTK REQUIRED)
#
# To use ``MIRTK_ROOT`` as cache variable::
#
#   set(MIRTK_CACHE_PATH_NAME MIRTK_ROOT)
#   find_package(MIRTK REQUIRED)
#
# To use ``SYSTEM_MIRTK_DIR`` as cache variable::
#
#   set(MIRTK_CACHE_PATH_PREFIX SYSTEM_)
#   find_package(MIRTK REQUIRED)

#=============================================================================
# Copyright 2016 Imperial College London
# Copyright 2016 Andreas Schuh
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
# 
# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
# 
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.
# 
# * Neither the names of Kitware, Inc., the Insight Software Consortium,
#   nor the names of their contributors may be used to endorse or promote
#   products derived from this software without specific prior written
#   permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#=============================================================================

# Set MIRTKConfig.cmake directory from installation prefix path
function (_mirtk_root_to_config_dir OUT IN)
  if (IN)
    if (WIN32)
      set(${OUT} "${IN}/CMake" PARENT_SCOPE)
    else ()
      set(${OUT} "${IN}/lib/cmake/mirtk" PARENT_SCOPE)
    endif ()
  else ()
    set(${OUT} "NOTFOUND" PARENT_SCOPE)
  endif ()
endfunction ()

# Set installation prefix path from MIRTKConfig.cmake directory
function (_mirtk_config_to_root_dir OUT IN)
  if (IN)
    if (WIN32)
      string(REGEX REPLACE "/+CMake/*$" "" _prefix "${IN}")
    else ()
      string(REGEX REPLACE "/+lib/+cmake/+mirtk/*$" "" _prefix "${IN}")
    endif ()
  else ()
    set(_prefix "NOTFOUND")
  endif ()
  set(${OUT} "${_prefix}" PARENT_SCOPE)
endfunction ()

# Name of cache variable to use for user hint of MIRTK location
if (MIRTK_CACHE_PATH_NAME)
  set(_cache "${MIRTK_CACHE_PATH_NAME}")
elseif (MIRTK_CACHE_PATH_PREFIX)
  set(_cache "${MIRTK_CACHE_PATH_PREFIX}MIRTK_DIR")
elseif (MIRTK_CACHE_PATH_CMAKE)
  set(_cache MIRTK_DIR)
else ()
  set(_cache DEPENDS_MIRTK_DIR)
endif ()

# Add cache entry for user hint of MIRTK search location
if ("^${_cache}$" STREQUAL "^MIRTK_DIR$")

  # Use default MIRTK_DIR cache entry used by find_package in CONFIG mode
  set (MIRTK_DIR "${MIRTK_DIR}" CACHE PATH "Directory containing MIRTKConfig.cmake file.")

elseif ("^${_cache}$" STREQUAL "^MIRTK_ROOT$")

  # Use MIRTK_ROOT as cache entry
  set (MIRTK_ROOT "${MIRTK_ROOT}" CACHE PATH "Installation prefix/root directory of MIRTK.")

  # Set MIRTK_ROOT from -DMIRTK_DIR=<path> and mark the latter as INTERNAL
  if (MIRTK_DIR AND (NOT DEFINED _MIRTK_DIR OR (DEFINED _MIRTK_DIR AND NOT "^${MIRTK_DIR}$" STREQUAL "^${_MIRTK_DIR}$")))
    _mirtk_config_to_root_dir(_prefix "${MIRTK_DIR}")
    set_property(CACHE MIRTK_ROOT PROPERTY VALUE "${_prefix}")
  endif ()
  get_property(_cached CACHE MIRTK_DIR PROPERTY TYPE SET)
  if (_cached)
    set_property(CACHE MIRTK_DIR PROPERTY TYPE INTERNAL)
  endif ()

  # Set MIRTK_ROOT from environment variable
  #
  # Note that the MIRTK_DIR environment variable is used by CMake's find_package
  # as installation prefix search path rather than the location of the MIRTKConfig.cmake!
  # The MIRTK_DIR environment variable is thus equivalent to the MIRTK_ROOT CMake variable.
  if (NOT MIRTK_ROOT)
    foreach (_dir IN ITEMS "$ENV{MIRTK_ROOT}" "$ENV{MIRTK_DIR}")
      if (_dir)
        set_property(CACHE ${_cache} PROPERTY VALUE "${_dir}")
        break()
      endif ()
    endforeach ()
  endif ()

  # Set MIRTK_DIR based on MIRTK_ROOT
  _mirtk_root_to_config_dir(_config_dir "${MIRTK_ROOT}")
  set(MIRTK_DIR "${_config_dir}" CACHE INTERNAL "Directory containing MIRTKConfig.cmake file." FORCE)

else ()

  # Use more user friendly hybrid DEPENDS_MIRTK_DIR cache variable which allows grouping
  # of DEPENDS paths or custom named cache entry, but still consider MIRTK_ROOT and
  # MIRTK_DIR as more common alternatives set in the user shell environment or on the
  # CMake command line. The DEPENDS_<Package>_DIR is what CMake BASIS uses by default.
  set (${_cache} "${${_cache}}" CACHE PATH "Installation prefix of MIRTK or directory containing MIRTKConfig.cmake file.")

  # Override DEPENDS_MIRTK_DIR by alternative search path variable value if these
  # were specified on the command line using the -D option. Note that these variables
  # cannot be set in the CMake GUI because their type is changed here to INTERNAL.
  # This has two reasons, firstly to not have duplicate variables with different
  # names for the same purpose, and secondly to be able to recognize when their
  # value is changed using the -D command line option of the cmake command.
  if (MIRTK_DIR AND (NOT DEFINED _MIRTK_DIR OR (DEFINED _MIRTK_DIR AND NOT "^${MIRTK_DIR}$" STREQUAL "^${_MIRTK_DIR}$")))
    _mirtk_config_to_root_dir(_prefix "${MIRTK_DIR}")
    set_property(CACHE ${_cache} PROPERTY VALUE "${_prefix}")
  endif ()
  if (MIRTK_ROOT AND (NOT DEFINED _MIRTK_ROOT OR (DEFINED _MIRTK_ROOT AND NOT "^${MIRTK_ROOT}$" STREQUAL "^${_MIRTK_ROOT}$")))
    set_property(CACHE ${_cache} PROPERTY VALUE "${MIRTK_ROOT}")
  endif ()

  # Mark alternatives as internal cache entries
  foreach (_var IN ITEMS MIRTK_DIR MIRTK_ROOT)
    get_property(_cached CACHE ${_var} PROPERTY TYPE SET)
    if (_cached)
      set_property(CACHE ${_var} PROPERTY TYPE INTERNAL)
    endif ()
  endforeach ()

  # If still not set, use common environment variables to set DEPENDS_MIRTK_DIR
  if (NOT ${_cache})
    foreach (_dir IN ITEMS "$ENV{MIRTK_DIR}" "$ENV{MIRTK_ROOT}")
      if (_dir)
        set_property(CACHE ${_cache} PROPERTY VALUE "${_dir}")
        break()
      endif ()
    endforeach ()
  endif ()

  # Allow DEPENDS_MIRTK_DIR to be set to either the root directory...
  if (${_cache})
    list (INSERT CMAKE_PREFIX_PATH 0 "${${_cache}}")
  endif ()
  # ...or the directory containing the MIRTKConfig.cmake file
  set(MIRTK_DIR "${${_cache}}" CACHE INTERNAL "Directory containing MIRTKConfig.cmake file." FORCE)
endif ()

# Look for MIRTK installation
if (NOT MIRTK_FIND_QUIETLY)
  set(_msg "Looking for MIRTK")
  if (MIRTK_FIND_VERSION)
    set(_msg "${_msg} ${MIRTK_FIND_VERSION}")
  endif ()
  if (MIRTK_FIND_COMPONENTS)
    set(_msg "${_msg} [${MIRTK_FIND_COMPONENTS}]")
  endif ()
  message(STATUS "${_msg}...")
endif ()

set(_argv)
if (MIRTK_FIND_VERSION)
  list(APPEND _argv ${MIRTK_FIND_VERSION})
endif ()
if (MIRTK_FIND_VERSION_EXACT)
  list(APPEND _argv EXACT)
endif ()
set(_comps)
set(_comps_opt)
foreach (_comp IN LISTS MIRTK_FIND_COMPONENTS)
  if (MIRTK_FIND_REQUIRED_${_component})
    list(APPEND _comps ${_component})
  else ()
    list(APPEND _comps_opt ${_component})
  endif ()
endforeach ()
if (_comps)
  list(APPEND _argv COMPONENTS ${_comps})
endif ()
if (_comps_opt)
  list(APPEND _argv OPTIONAL_COMPONENTS ${_comps_opt})
endif ()

find_package(MIRTK ${_argv} CONFIG QUIET)

if (NOT MIRTK_FIND_QUIETLY)
  if (MIRTK_FOUND)
    message(STATUS "Looking for MIRTK... - found v${MIRTK_VERSION}")
  else ()
    message(STATUS "Looking for MIRTK... - not found")
  endif ()
endif ()

# Make internal search path cache entries consistent with non-internal cache entry
if (MIRTK_FOUND)
  mark_as_advanced(FORCE ${_cache})
  _mirtk_config_to_root_dir(_prefix "${MIRTK_DIR}")
else ()
  mark_as_advanced(CLEAR ${_cache})
  set(_prefix NOTFOUND)
endif ()
if (NOT "^${_cache}$" STREQUAL "^MIRTK_DIR$")
  set_property(CACHE MIRTK_DIR PROPERTY TYPE INTERNAL)
  set(_MIRTK_DIR "${MIRTK_DIR}" CACHE INTERNAL "Previous MIRTK_DIR value" FORCE)
endif ()
get_property(_cached CACHE MIRTK_ROOT PROPERTY TYPE SET)
if (_cached)
  set_property(CACHE MIRTK_ROOT PROPERTY VALUE "${_prefix}")
else ()
  set(MIRTK_ROOT "${_prefix}")
endif ()
if (NOT _cache MATCHES "^MIRTK_(DIR|ROOT)$")
  set_property(CACHE ${_cache} PROPERTY VALUE "${_prefix}")
endif ()

# Make internal cache copies of alternative search path variables
# so we can detect when a new value was specified using -D option
foreach (_var IN ITEMS MIRTK_DIR MIRTK_ROOT)
  if (NOT "^${_cache}$" STREQUAL "^${_var}$")
    get_property(_cached CACHE ${_var} PROPERTY TYPE SET)
    if (_cached)
      set(_${_var} "${${_var}}" CACHE INTERNAL "Previous value of ${_var} after last find_package(MIRTK)" FORCE)
    endif ()
  endif ()
endforeach ()

# Raise fatal error when MIRTK required but not found
if (MIRTK_FIND_REQUIRED AND NOT MIRTK_FOUND)
  set(_msg "Could not find MIRTK! Please ensure that it is installed"
           " in a standard system location or set ${_cache}")
  if ("^${_cache}$" STREQUAL "^MIRTK_DIR$")
    set(_msg "${_msg} to the directory containing the MIRTKConfig.cmake file.")
  elseif ("^${_cache}$" STREQUAL "^MIRTK_ROOT$")
    set(_msg "${_msg} to the installation prefix of MIRTK, i.e., the root directory.")
  else ()
    set(_msg "${_msg} either to the installation prefix of MIRTK, i.e., the root directory,"
                    " or the directory containing the MIRTKConfig.cmake file.")
  endif ()
  string(REPLACE ";" "" _msg "${_msg}")
  message(FATAL_ERROR "${_msg}")
endif ()

# Unset local variables
unset(_prefix)
unset(_config_dir)
unset(_cache)
unset(_argv)
unset(_comp)
unset(_comps)
unset(_comps_opt)
unset(_cached)
unset(_dir)
unset(_var)
unset(_msg)
