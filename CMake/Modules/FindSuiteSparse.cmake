#.rst:
# FindSuiteSparse
# ---------------
#
# Find the SuiteSparse libraries.
# http://faculty.cse.tamu.edu/davis/suitesparse.html
#
# This module reads hints about search locations from variables:
#
# ::
#
#    SuiteSparse_ROOT  - Root directory of SuiteSparse installation.
#                        Can be an environment variable instead. It is
#                        derived from the found SuiteSparse_INCLUDE_DIR
#                        if unset before find_package(SuiteSparse).
#                        Can also be an environment variable.
#
# This module considers the following CMake variables set by find_package:
#
# ::
#
#    SuiteSparse_FIND_COMPONENTS    - Names of requested libraries:
#                                       AMD, CAMD, COLAMD, CCOLAMD, CHOLMOD - ordering methods
#                                       SPQR     - multifrontal QR factorization
#                                       UMFPACK  - multifrontal LU factorization
#                                       CXSparse - sparse Cholesky factorization
#                                       KLU, BTF - sparse LU factorization
#                                       LDL      - sparse LDL factorization
#                                       RBio     - MATLAB toolbox for reading/writing
#                                                  sparse matrices in Rutherford/Boeing format
#    SuiteSparse_FIND_REQUIRED_<C>  - Whether SuiteSparse library <C> is required.
#                                     SuiteSparse is considered to be not found when at least
#                                     one required library or its include path is missing.
#    SuiteSparse_FIND_REQUIRED      - Raise FATAL_ERROR when required components not found.
#    SuiteSparse_FIND_QUIETLY       - Suppress all other (status) messages.
#
# This modules sets the following variables with the aggregated results of all
# found components, including copies of these variables with the all uppercase
# SUITESPARSE_ prefix.
#
# ::
#
#    SuiteSparse_FOUND          - Whether all required components were found.
#    SuiteSparse_VERSION        - Version for use in VERSION_LESS et al. comparisons.
#    SuiteSparse_VERSION_MAJOR  - Major library version number.
#    SuiteSparse_VERSION_MINOR  - Minor library version number.
#    SuiteSparse_VERSION_PATCH  - Minor library version number.
#    SuiteSparse_VERSION_STRING - Version string for output messages.
#    SuiteSparse_INCLUDE_DIR    - Include path of SuiteSparse_config.h (version >= 4) or UFconfig.h (version < 4)
#    SuiteSparse_LIBRARY_DIR    - Directory containing "config" library file.
#    SuiteSparse_INCLUDE_DIRS   - Include path of found libraries and their dependencies.
#    SuiteSparse_LIBRARIES      - Path of found libraries and their link dependencies.
#
#    SUITESPARSE_MAIN_VERSION   - Same as SuiteSparse_VERSION_MAJOR.
#    SUITESPARSE_SUB_VERSION    - Same as SuiteSparse_VERSION_MINOR.
#    SUITESPARSE_SUBSUB_VERSION - Same as SuiteSaprse_VERSION_PATCH.
#
# Additionally, the following variables are set for each found component, including
# the always looked for "config" component:
#
# ::
#
#    [SuiteSpares_]<C>_INCLUDE_DIR  - Include path of library header.
#    [SuiteSpares_]<C>_INCLUDE_DIRS - Include path of library header including dependencies.
#    [SuiteSpares_]<C>_LIBRARY      - Path of library file of component.
#    [SuiteSpares_]<C>_LIBRARIES    - Path of library file including link dependencies.
#
# where <C> is the name of the component with case indentical to the
# SuiteSparse_FIND_COMPONENTS argument in the prefixed version, and in all uppercase
# letters without the prefix, except for the "config" component.
#
# For example, SuiteSparse_CXSparse_LIBRARIES and CXSPARSE_LIBRARIES are the variables
# both containing the list of libraries for the CXSparse component. The uncached non-prefixed
# all uppercase variables are provided to be compatible with other FindSuiteSparse.cmake
# modules and to accomodate projects which use SuiteSparse libraries that were previously
# distributed separately such as, for example, UMFPACK. In this case, a custom
# FindUMFPACK.cmake module can be replaced by:
#
# ::
#
#    find_package(SuiteSparse COMPONENTS UMFPACK)
#
# and variables UMFPACK_FOUND, UMFPACK_INCLUDE_DIRS, and UMFPACK_LIBRARIES can be used
# to refer to the found UMFPACK installation. These variables would otherwise be set by
# a custom FindUMFPACK.cmake module.
#
# To debug the discovery of SuiteSparse libraries, set `SuiteSparse_DEBUG` to
# a true value before the find_package call. This Find module will then output
# additional diagnostic messages which help to identify problems finding a
# SuiteSparse installation.

#=============================================================================
# Copyright 2016 Andreas Schuh <andreas.schuh.84@gmail.com>
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distribute this file outside of CMake, substitute the full
#  License text for the above reference.)

# ------------------------------------------------------------------------------
# Print debug message when debugging of this Find module is enabled
macro (_suitesparse_debug_message msg)
  if (SuiteSparse_DEBUG)
    message("** FindSuiteSparse: ${msg}")
  endif ()
endmacro ()

# ------------------------------------------------------------------------------
# Get SuiteSparse component include path suffixes
function (_suitesparse_get_include_suffixes varname component)
  string(TOLOWER "${component}" component_lower)
  set(
    ${varname}
      "${component}"
      "include/${component}"
      "${component}/Include"
      "${component_lower}"
      "include/${component_lower}"
      "suitesparse"
      "include"
      "src"
    PARENT_SCOPE
  )
endfunction ()

# ------------------------------------------------------------------------------
# Get SuiteSparse component library path suffixes
macro (_suitesparse_get_library_suffixes varname component)
  set(
    ${varname}
      "${component}"
      "${component}/Lib"
      "lib/x86_64-linux-gnu"
      "lib"
      "lib64"
      "lib32"
  )
endmacro ()

# ------------------------------------------------------------------------------
# Remove duplicate include paths from list
macro (_suitesparse_remove_duplicate_paths varname)
  if (${varname})
    list(REMOVE_DUPLICATES ${varname})
  endif ()
endmacro ()

# ------------------------------------------------------------------------------
# Remove duplicate libraries from list
macro (_suitesparse_remove_duplicate_libraries varname)
  if (${varname})
    list(REVERSE ${varname})
    list(REMOVE_DUPLICATES ${varname})
    list(REVERSE ${varname})
  endif ()
endmacro ()

# ------------------------------------------------------------------------------
# List of components to look for, including inter-dependencies
# Note: Order of components must be suitable for static linking!
set(SuiteSparse_COMPONENTS)
foreach (_suitesparse_component IN LISTS SuiteSparse_FIND_COMPONENTS)
  if (NOT "^${_suitesparse_component}$" STREQUAL "^config$")
    list(APPEND SuiteSparse_COMPONENTS ${_suitesparse_component})
    if ("^${_suitesparse_component}$" STREQUAL "^UMFPACK$")
      list(APPEND SuiteSparse_COMPONENTS CHOLMOD COLAMD AMD)
    elseif ("^${_suitesparse_component}$" STREQUAL "^CHOLMOD$")
      list(APPEND SuiteSparse_COMPONENTS CCOLAMD COLAMD CAMD AMD)
    endif ()
  endif ()
endforeach ()
_suitesparse_remove_duplicate_libraries(SuiteSparse_COMPONENTS)

# ------------------------------------------------------------------------------
# Status message
if (NOT SuiteSparse_FIND_QUIETLY)
  set(_suitesparse_find_status "Looking for SuiteSparse")
  if (SuiteSparse_FIND_COMPONENTS)
    set(_suitesparse_find_status "${_suitesparse_find_status} [${SuiteSparse_COMPONENTS}]")
  endif ()
  if (NOT SuiteSparse_FIND_REQUIRED)
    set(_suitesparse_find_status "${_suitesparse_find_status} (optional)")
  endif ()
  message(STATUS "${_suitesparse_find_status}...")
endif ()

# ------------------------------------------------------------------------------
# Find "config" and set SuiteSparse_ROOT (if unset) used to find other components
if (NOT SuiteSparse_ROOT)
  file(TO_CMAKE_PATH "$ENV{SuiteSparse_ROOT}" SuiteSparse_ROOT)
endif ()

_suitesparse_debug_message("Hints:")
_suitesparse_debug_message("- SuiteSparse_ROOT = ${SuiteSparse_ROOT}")

_suitesparse_get_include_suffixes(_suitesparse_include_suffixes SuiteSparse_config)
if (SuiteSparse_ROOT)
  find_path(SuiteSparse_config_INCLUDE_DIR
    NAMES SuiteSparse_config.h UFconfig.h
    PATHS "${SuiteSparse_ROOT}"
    PATH_SUFFIXES ${_suitesparse_include_suffixes}
    NO_DEFAULT_PATH
  )
else ()
  find_path(SuiteSparse_config_INCLUDE_DIR
    NAMES SuiteSparse_config.h UFconfig.h
    PATH_SUFFIXES ${_suitesparse_include_suffixes}
  )
  if (SuiteSparse_config_INCLUDE_DIR)
    foreach (_suitesparse_include_suffix IN LISTS _suitesparse_include_suffixes)
      if (SuiteSparse_config_INCLUDE_DIR MATCHES "^(.*)/+${_suitesparse_include_suffix}/*$")
        string(REGEX REPLACE "/[iI]nclude" "" SuiteSparse_ROOT "${CMAKE_MATCH_1}")
        break()
      endif ()
    endforeach ()
  endif ()
endif ()
mark_as_advanced(SuiteSparse_config_INCLUDE_DIR)

_suitesparse_get_library_suffixes(_suitesparse_library_suffixes SuiteSparse_config)
find_library(SuiteSparse_config_LIBRARY
  NAMES suitesparseconfig
  PATHS "${SuiteSparse_ROOT}"
  PATH_SUFFIXES ${_suitesparse_library_suffixes}
  NO_DEFAULT_PATH
)
mark_as_advanced(SuiteSparse_config_LIBRARY)

if (SuiteSparse_config_INCLUDE_DIR AND SuiteSparse_config_LIBRARY)
  set(SuiteSparse_config_FOUND 1)
  set(SuiteSparse_config_INCLUDE_DIRS ${SuiteSparse_config_INCLUDE_DIR})
  set(SuiteSparse_config_LIBRARIES    ${SuiteSparse_config_LIBRARY})
  set(SuiteSparse_INCLUDE_DIR         ${SuiteSparse_config_INCLUDE_DIR})
  get_filename_component(SuiteSparse_LIBRARY_DIR "${SuiteSparse_config_LIBRARY}" PATH)
else ()
  set(SuiteSparse_config_FOUND 0)
  set(SuiteSparse_config_INCLUDE_DIRS)
  set(SuiteSparse_config_LIBRARIES)
  set(SuiteSparse_INCLUDE_DIR)
  set(SuiteSparse_LIBRARY_DIR)
endif ()

_suitesparse_debug_message("Config (FOUND=${SuiteSparse_config_FOUND}):")
_suitesparse_debug_message("- SuiteSparse_ROOT               = ${SuiteSparse_ROOT}")
_suitesparse_debug_message("- SuiteSparse_INCLUDE_DIR        = ${SuiteSparse_INCLUDE_DIR}")
_suitesparse_debug_message("- SuiteSparse_LIBRARY_DIR        = ${SuiteSparse_LIBRARY_DIR}")
_suitesparse_debug_message("- SuiteSparse_config_INCLUDE_DIR = ${SuiteSparse_config_INCLUDE_DIR}")
_suitesparse_debug_message("- SuiteSparse_config_LIBRARY     = ${SuiteSparse_config_LIBRARY}")

# ------------------------------------------------------------------------------
# Extract version information
if (SuiteSparse_config_FOUND)
  if (EXISTS "${SuiteSparse_config_INCLUDE_DIR}/SuiteSparse_config.h")
    set(_suitesparse_config_file "${SuiteSparse_config_INCLUDE_DIR}/SuiteSparse_config.h")
  else ()
    set(_suitesparse_config_file "${SuiteSparse_config_INCLUDE_DIR}/UFconfig.h")
  endif ()
  file(READ "${_suitesparse_config_file}" _suitesparse_config)
  if (_suitesparse_config MATCHES "#define SUITESPARSE_MAIN_VERSION ([0-9]+)")
    set(SuiteSparse_VERSION_MAJOR "${CMAKE_MATCH_1}")
    if (_suitesparse_config MATCHES "#define SUITESPARSE_SUB_VERSION ([0-9]+)")
      set(SuiteSparse_VERSION_MINOR "${CMAKE_MATCH_1}")
      if (_suitesparse_config MATCHES "#define SUITESPARSE_SUBSUB_VERSION ([0-9]+)")
        set(SuiteSparse_VERSION_PATCH "${CMAKE_MATCH_1}")
      else ()
        set(SuiteSparse_VERSION_PATCH 0)
      endif ()
    else ()
      set(SuiteSparse_VERSION_MINOR 0)
    endif ()
    set(SuiteSparse_VERSION
      "${SuiteSparse_VERSION_MAJOR}.${SuiteSparse_VERSION_MINOR}.${SuiteSparse_VERSION_PATCH}"
    )
    if (SuiteSparse_VERSION_PATCH)
      set(SuiteSparse_VERSION_STRING ${SuiteSparse_VERSION})
    else ()
      set(SuiteSparse_VERSION_STRING ${SuiteSparse_VERSION_MAJOR}.${SuiteSparse_VERSION_MINOR})
    endif ()
  else ()
    set(SuiteSparse_VERSION)
    set(SuiteSparse_VERSION_MAJOR)
    set(SuiteSparse_VERSION_MINOR)
    set(SuiteSparse_VERSION_PATCH)
    set(SuiteSparse_VERSION_STRING)
  endif ()

  _suitesparse_debug_message("Version information from ${_suitesparse_config_file}")
  _suitesparse_debug_message("- SuiteSparse_VERSION_STRING = ${SuiteSparse_VERSION_STRING}")
  _suitesparse_debug_message("- SuiteSparse_VERSION_MAJOR  = ${SuiteSparse_VERSION_MAJOR}")
  _suitesparse_debug_message("- SuiteSparse_VERSION_MINOR  = ${SuiteSparse_VERSION_MINOR}")
  _suitesparse_debug_message("- SuiteSparse_VERSION_PATCH  = ${SuiteSparse_VERSION_PATCH}")

  unset(_suitesparse_config)
  unset(_suitesparse_config_file)
endif ()

# ------------------------------------------------------------------------------
# Look for other requested components
foreach (_suitesparse_component IN LISTS SuiteSparse_COMPONENTS)
  if (NOT "^${_suitesparse_component}$" STREQUAL "^config$")
    string(TOLOWER "${_suitesparse_component}" _suitesparse_component_lower)

    # Find component include path
    _suitesparse_get_include_suffixes(_suitesparse_include_suffixes ${_suitesparse_component})
    if ("^${_suitesparse_component}$" STREQUAL "^SPQR$")
      set(_suitesparse_component_header SuiteSparseQR.hpp)
    else ()
      set(_suitesparse_component_header ${_suitesparse_component_lower}.h)
    endif ()

    find_path(SuiteSparse_${_suitesparse_component}_INCLUDE_DIR
      NAMES ${_suitesparse_component_header}
      HINTS "${SuiteSparse_INCLUDE_DIR}"
      PATHS "${SuiteSparse_ROOT}"
      PATH_SUFFIXES ${_suitesparse_include_suffixes}
      NO_DEFAULT_PATH
    )
    mark_as_advanced(SuiteSparse_${_suitesparse_component}_INCLUDE_DIR)

    # Find component library file
    set(SuiteSparse_${_suitesparse_component}_LIBRARIES)
    _suitesparse_get_library_suffixes(_suitesparse_library_suffixes ${_suitesparse_component})
    find_library(SuiteSparse_${_suitesparse_component}_LIBRARY
      NAMES ${_suitesparse_component_lower}
      HINTS "${SuiteSparse_LIBRARY_DIR}"
      PATHS "${SuiteSparse_ROOT}"
      PATH_SUFFIXES ${_suitesparse_library_suffixes}
      NO_DEFAULT_PATH
    )
    mark_as_advanced(SuiteSparse_${_suitesparse_component}_LIBRARY)

    # Mark component as either found or not and init aggregates
    if (SuiteSparse_${_suitesparse_component}_INCLUDE_DIR AND
        SuiteSparse_${_suitesparse_component}_LIBRARY)
      set(SuiteSparse_${_suitesparse_component}_FOUND 1)
      set(SuiteSparse_${_suitesparse_component}_INCLUDE_DIRS
        ${SuiteSparse_${_suitesparse_component}_INCLUDE_DIR}
      )
      set(SuiteSparse_${_suitesparse_component}_LIBRARIES
        ${SuiteSparse_${_suitesparse_component}_LIBRARY}
      )
    else ()
      set(SuiteSparse_${_suitesparse_component}_FOUND 0)
      set(SuiteSparse_${_suitesparse_component}_INCLUDE_DIRS)
      set(SuiteSparse_${_suitesparse_component}_LIBRARIES)
    endif ()

  endif ()
endforeach ()

# ------------------------------------------------------------------------------
# Add inter-library link dependencies

# Library suitesparseconfig version >= 4 requires rt library for timing by default
# when compiled on Unix, but not on OS X (which does not have librt)
if (UNIX AND NOT APPLE)
  find_library(SuiteSparse_RT_LIBRARY NAMES rt)
  mark_as_advanced(SuiteSparse_RT_LIBRARY)
  if (SuiteSparse_RT_LIBRARY)
    list(APPEND SuiteSparse_config_LIBRARIES ${SuiteSparse_RT_LIBRARY})
  endif ()
endif ()

# CHOLMOD: required [AMD, COLAMD], optional [CAMD, CCOLAMD]
if (SuiteSparse_CHOLMOD_FOUND)
  if (SuiteSparse_AMD_FOUND)
    list(APPEND SuiteSparse_CHOLMOD_INCLUDE_DIRS ${SuiteSparse_AMD_INCLUDE_DIRS})
    list(APPEND SuiteSparse_CHOLMOD_LIBRARIES    ${SuiteSparse_AMD_LIBRARY})
  else ()
    message(WARNING "SuiteSparse CHOLMOD library requires AMD library which was not found")
    set(SuiteSparse_CHOLMOD_FOUND "CHOLMOD requires AMD and COLAMD-NOTFOUND")
  endif ()
  if (SuiteSparse_COLAMD_FOUND)
    list(APPEND SuiteSparse_CHOLMOD_INCLUDE_DIRS ${SuiteSparse_COLAMD_INCLUDE_DIRS})
    list(APPEND SuiteSparse_CHOLMOD_LIBRARIES    ${SuiteSparse_COLAMD_LIBRARY})
  else ()
    message(WARNING "SuiteSparse CHOLMOD library requires COLAMD library which was not found")
    set(SuiteSparse_CHOLMOD_FOUND "CHOLMOD requires AMD and COLAMD-NOTFOUND")
  endif ()
  if (SuiteSparse_CAMD_FOUND)
    list(APPEND SuiteSparse_CHOLMOD_INCLUDE_DIRS ${SuiteSparse_CAMD_INCLUDE_DIRS})
    list(APPEND SuiteSparse_CHOLMOD_LIBRARIES    ${SuiteSparse_CAMD_LIBRARY})
  endif ()
  if (SuiteSparse_CCOLAMD_FOUND)
    list(APPEND SuiteSparse_CHOLMOD_INCLUDE_DIRS ${SuiteSparse_CCOLAMD_INCLUDE_DIRS})
    list(APPEND SuiteSparse_CHOLMOD_LIBRARIES    ${SuiteSparse_CCOLAMD_LIBRARY})
  endif ()
endif()

# UMFPACK: required [], optional [CHOLMOD, AMD]
if (SuiteSparse_UMFPACK_FOUND)
  if (SuiteSparse_CHOLMOD_FOUND)
    list(APPEND SuiteSparse_UMFPACK_INCLUDE_DIRS ${SuiteSparse_CHOLMOD_INCLUDE_DIRS})
    list(APPEND SuiteSparse_UMFPACK_LIBRARIES    ${SuiteSparse_CHOLMOD_LIBRARIES})
  elseif (SuiteSparse_AMD_FOUND)
    list(APPEND SuiteSparse_UMFPACK_INCLUDE_DIRS ${SuiteSparse_AMD_INCLUDE_DIRS})
    list(APPEND SuiteSparse_UMFPACK_LIBRARIES    ${SuiteSparse_AMD_LIBRARIES})
  endif()
endif()

# CXSparse: optional [TBB]
if (SuiteSparse_CXSparse_FOUND)
  find_package(TBB COMPONENTS tbb malloc QUIET)
  if (TBB_FOUND)
    list(APPEND SuiteSparse_${_suitesparse_component}_LIBRARIES ${TBB_LIBRARIES})
  endif ()
endif ()

# Find external link dependencies, but only when SuiteSparse "config" found
if (SuiteSparse_config_FOUND)
  # Add "config" to all components except CXSparse
  foreach (_suitesparse_component IN LISTS SuiteSparse_COMPONENTS)
    if (NOT "^{_suitesparse_component}$" STREQUAL "^CXSparse$")
      list(APPEND SuiteSparse_${_suitesparse_component}_INCLUDE_DIRS
        ${SuiteSparse_config_INCLUDE_DIRS}
      )
      list(APPEND SuiteSparse_${_suitesparse_component}_LIBRARIES
        ${SuiteSparse_config_LIBRARIES}
      )
    endif ()
  endforeach ()
  # BLAS (required)
  find_package(BLAS QUIET)
  if (BLAS_FOUND)
    foreach (_suitesparse_component IN LISTS SuiteSparse_COMPONENTS)
      list(APPEND SuiteSparse_${_suitesparse_component}_LIBRARIES ${BLAS_LIBRARIES})
    endforeach ()
  endif ()
  # METIS (optional)
  find_library(SuiteSparse_METIS_LIBRARY NAMES metis)
  mark_as_advanced(SuiteSparse_METIS_LIBRARY)
  if (SuiteSparse_METIS_LIBRARY)
    foreach (_suitesparse_component IN LISTS SuiteSparse_COMPONENTS)
      list(APPEND SuiteSparse_${_suitesparse_component}_LIBRARIES ${SuiteSparse_METIS_LIBRARY})
    endforeach ()
  endif ()
endif ()

# Remove duplicates
foreach (_suitesparse_component IN LISTS SuiteSparse_COMPONENTS)
  _suitesparse_remove_duplicate_paths(SuiteSparse_${_suitesparse_component}_INCLUDE_DIRS)
  _suitesparse_remove_duplicate_libraries(SuiteSparse_${_suitesparse_component}_LIBRARIES)

  _suitesparse_debug_message("${_suitesparse_component} (FOUND=${SuiteSparse_${_suitesparse_component}_FOUND}):")
  _suitesparse_debug_message("- SuiteSparse_${_suitesparse_component}_INCLUDE_DIR  = ${SuiteSparse_${_suitesparse_component}_INCLUDE_DIR}")
  _suitesparse_debug_message("- SuiteSparse_${_suitesparse_component}_INCLUDE_DIRS = [${SuiteSparse_${_suitesparse_component}_INCLUDE_DIRS}]")
  _suitesparse_debug_message("- SuiteSparse_${_suitesparse_component}_LIBRARY      = ${SuiteSparse_${_suitesparse_component}_LIBRARY}")
  _suitesparse_debug_message("- SuiteSparse_${_suitesparse_component}_LIBRARIES    = [${SuiteSparse_${_suitesparse_component}_LIBRARIES}]")
endforeach ()

# ------------------------------------------------------------------------------
# Aggregate found libraries
if (SuiteSparse_config_FOUND)
  set(SuiteSparse_INCLUDE_DIRS ${SuiteSparse_config_INCLUDE_DIRS})
  set(SuiteSparse_LIBRARIES    ${SuiteSparse_config_LIBRARIES})
  foreach (_suitesparse_component IN LISTS SuiteSparse_COMPONENTS)
    if (SuiteSparse_${_suitesparse_component}_FOUND)
      list(APPEND SuiteSparse_INCLUDE_DIRS SuiteSparse_${_suitesparse_component}_INCLUDE_DIRS)
      list(APPEND SuiteSparse_LIBRARIES    SuiteSparse_${_suitesparse_component}_LIBRARIES)
    endif ()
  endforeach ()
  _suitesparse_remove_duplicate_paths(SuiteSparse_INCLUDE_DIRS)
  _suitesparse_remove_duplicate_libraries(SuiteSparse_LIBRARIES)
else ()
  set(SuiteSparse_INCLUDE_DIRS)
  set(SuiteSparse_LIBRARIES)
endif ()

# ------------------------------------------------------------------------------
# Handle QUIET, REQUIRED, and [EXACT] VERSION arguments and set SuiteSparse_FOUND
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  SuiteSparse
  FOUND_VAR
    SuiteSparse_FOUND
  VERSION_VAR
    SuiteSparse_VERSION
  REQUIRED_VARS
    BLAS_FOUND
    SuiteSparse_config_INCLUDE_DIR
    SuiteSparse_config_LIBRARY
  HANDLE_COMPONENTS
)

if (NOT SuiteSparse_FIND_QUIETLY)
  if (SuiteSparse_FOUND)
    message(STATUS "${_suitesparse_find_status}... - found v${SuiteSparse_VERSION_STRING}")
  else ()
    message(STATUS "${_suitesparse_find_status}... - not found")
  endif ()
endif ()

# ------------------------------------------------------------------------------
# Set uncached variables without "SuiteSparse_" prefix for components which
# previously may have been installed separate from SuiteSparse and to be
# compatible with other FindSuiteSparse.cmake modules (e.g., Ceres solver project)
set(_suitesparse_components       AMD CAMD COLAMD CCOLAMD CHOLMOD UMFPACK SPQR CXSparse)
set(_suitesparse_varname_suffixes FOUND INCLUDE_DIR INCLUDE_DIRS LIBRARY LIBRARIES)
foreach (_suitesparse_component IN LISTS _suitesparse_components)
  if (SuiteSparse_${_suitesparse_component}_FOUND)
    foreach (_suitesparse_varname_suffix IN LISTS _suitesparse_varname_suffixes)
      string(TOUPPER "${_suitesparse_component}" _suitesparse_component_upper)
      set(${_suitesparse_component_upper}_${_suitesparse_varname_suffix}
        ${SuiteSparse_${_suitesparse_component}_${_suitesparse_varname_suffix}}
      )
    endforeach ()
  endif ()
endforeach ()
if (SuiteSparse_SPQR_FOUND)
  foreach (_suitesparse_varname_suffix IN LISTS _suitesparse_varname_suffixes)
    set(SUITESPARSEQR_${_suitesparse_varname_suffix}
      ${SuiteSparse_SPQR_${_suitesparse_varname_suffix}}
    )
  endforeach ()
endif ()

# ------------------------------------------------------------------------------
# Compatibility with other FindSuiteSparse.cmake modules
set(SUITESPARSE_FOUND          ${SuiteSparse_FOUND})
set(SUITESPARSE_INCLUDE_DIR    ${SuiteSparse_INCLUDE_DIR})
set(SUITESPARSE_INCLUDE_DIRS   ${SuiteSparse_INCLUDE_DIRS})
set(SUITESPARSE_LIBRARY_DIR    ${SuiteSparse_LIBRARY_DIR})
set(SUITESPARSE_LIBRARIES      ${SuiteSparse_LIBRARIES})
set(SUITESPARSE_VERSION        ${SuiteSparse_VERSION})
set(SUITESPARSE_VERSION_MAJOR  ${SuiteSparse_VERSION_MAJOR})
set(SUITESPARSE_VERSION_MINOR  ${SuiteSparse_VERSION_MINOR})
set(SUITESPARSE_VERSION_PATCH  ${SuiteSparse_VERSION_PATCH})
set(SUITESPARSE_MAIN_VERSION   ${SuiteSparse_VERSION_MAJOR})
set(SUITESPARSE_SUB_VERSION    ${SuiteSparse_VERSION_MINOR})
set(SUITESPARSE_SUBSUB_VERSION ${SuiteSparse_VERSION_PATCH})

# ------------------------------------------------------------------------------
# Unset local auxiliary variables
unset(_suitesparse_components)
unset(_suitesparse_component)
unset(_suitesparse_component_lower)
unset(_suitesparse_component_upper)
unset(_suitesparse_component_index)
unset(_suitesparse_component_header)
unset(_suitesparse_include_suffixes)
unset(_suitesparse_library_suffixes)
unset(_suitesparse_varname_suffixes)
unset(_suitesparse_varname_suffix)
unset(_suitesparse_find_status)
