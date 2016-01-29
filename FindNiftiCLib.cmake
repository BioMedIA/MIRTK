# ============================================================================
# Copyright (c) 2011-2012 University of Pennsylvania
# Copyright (c) 2013-2014 Andreas Schuh
# All rights reserved.
#
# See COPYING file for license information or visit
# http://opensource.andreasschuh.com/cmake-basis/download.html#license
# ============================================================================

##############################################################################
# @file  FindNiftiCLib.cmake
# @brief Find nifticlib package.
#
# @par Input variables:
# <table border="0">
#   <tr>
#     @tp @b NiftiCLib_DIR @endtp
#     <td>The nifticlib package files are searched under the specified root
#         directory. If they are not found there, the default search paths
#         are considered. This variable can also be set as environment variable.</td>
#   </tr>
#   <tr>
#     @tp @b NIFTICLIB_DIR @endtp
#     <td>Alternative environment variable for @p NiftiCLib_DIR.</td>
#   </tr>
#   <tr>
#     @tp @b NiftiCLib_USE_STATIC_LIB @endtp
#     <td>Forces this module to search for the static library. Otherwise,
#         the shared library is preferred.</td>
#   </tr>
# </table>
#
# @par Output variables:
# <table border="0">
#   <tr>
#     @tp  @b NiftiCLib_FOUND @endtp
#     <td>Whether the nifticlib package was found and the following CMake
#         variables are valid.</td>
#   </tr>
#   <tr>
#     @tp @b NiftiCLib_INCLUDE_DIR @endtp
#     <td>Cached include directory/ies.</td>
#   </tr>
#   <tr>
#     @tp @b NiftiCLib_INCLUDE_DIRS @endtp
#     <td>Alias for @p NiftiCLib_INCLUDE_DIR (not cached).</td>
#   </tr>
#   <tr>
#     @tp @b NiftiCLib_INCLUDES @endtp
#     <td>Alias for @p NiftiCLib_INCLUDE_DIR (not cached).</td>
#   </tr>
#   <tr>
#     @tp @b NiftiCLib_LIBRARY @endtp
#     <td>Path of @c niftiio library.</td>
#   </tr>
#   <tr>
#     @tp @b NiftiCLib_LIB @endtp
#     <td>Alias for @p NiftiCLib_LIBRARY (not cached).</td>
#   </tr>
#   <tr>
#     @tp @b NiftiCLib_LIBRARIES @endtp
#     <td>Path of @c niftiio library and prerequisite libraries.</td>
#   </tr>
# </table>
#
# @par Imported targets:
# <table border="0">
#   <tr>
#     @tp @b niftiio @endtp
#     <td>The library target of the @c nifticlib library.</td>
#   </tr>
# </table>
#
# @ingroup CMakeFindModules
##############################################################################

# ----------------------------------------------------------------------------
# initialize search
if (NOT NiftiCLib_DIR)
  if (NOT "$ENV{NIFTICLIB_DIR}" STREQUAL "")
    set (NiftiCLib_DIR "$ENV{NIFTICLIB_DIR}" CACHE PATH "Installation prefix for NiftiCLib." FORCE)
  else ()
    set (NiftiCLib_DIR "$ENV{NiftiCLib_DIR}" CACHE PATH "Installation prefix for NiftiCLib." FORCE)
  endif ()
endif ()

set (NiftiCLib_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})

if (NiftiCLib_USE_STATIC_LIB)
  if (WIN32)
    set (CMAKE_FIND_LIBRARY_SUFFIXES .lib)
  else ()
    set (CMAKE_FIND_LIBRARY_SUFFIXES .a)
  endif()
else ()
  if (WIN32)
    set (CMAKE_FIND_LIBRARY_SUFFIXES .dll .lib)
  elseif(APPLE)
    set (CMAKE_FIND_LIBRARY_SUFFIXES .dylib .a)
  else ()
    set (CMAKE_FIND_LIBRARY_SUFFIXES .so .a)
  endif()
endif ()

# ----------------------------------------------------------------------------
# find paths/files
if (NiftiCLib_DIR)

  find_path (
    NiftiCLib_INCLUDE_DIR
      NAMES         nifti/nifti1_io.h
      HINTS         ${NiftiCLib_DIR}
      PATH_SUFFIXES "include"
      DOC           "path to directory containing nifti1_io.h file."
      NO_DEFAULT_PATH
  )

  find_library (
    NiftiCLib_LIBRARY
      NAMES         niftiio
      HINTS         ${NiftiCLib_DIR}
      PATH_SUFFIXES lib
      DOC           "Path of niftiio library"
      NO_DEFAULT_PATH
  )

  find_library (
    NiftiCLib_znz_LIBRARY
      NAMES znz
      HINTS ENV LD_LIBRARY_PATH
      DOC   "Path of znz library"
  )

else ()

  find_path (
    NiftiCLib_INCLUDE_DIR
      NAMES         nifti/nifti1_io.h
      HINTS         ENV C_INCLUDE_PATH ENV CXX_INCLUDE_PATH
      DOC           "path to directory containing nifti1_io.h file."
  )

  find_library (
    NiftiCLib_LIBRARY
      NAMES niftiio
      HINTS ENV LD_LIBRARY_PATH
      DOC   "Path of niftiio library"
  )

  find_library (
    NiftiCLib_znz_LIBRARY
      NAMES znz
      HINTS ENV LD_LIBRARY_PATH
      DOC   "Path of znz library"
  )

endif ()

mark_as_advanced (NiftiCLib_INCLUDE_DIR)
mark_as_advanced (NiftiCLib_LIBRARY)
mark_as_advanced (NiftiCLib_znz_LIBRARY)

# ----------------------------------------------------------------------------
# prerequisites
if (NiftiCLib_USE_STATIC_LIB OR NiftiCLib_znz_LIBRARY MATCHES "\\.a$")
  find_package (ZLIB REQUIRED)
endif ()

set (NiftiCLib_LIBRARIES "${ZLIB_LIBRARIES}")
if (NiftiCLib_znz_LIBRARY)
  list (APPEND NiftiCLib_LIBRARIES "${NiftiCLib_znz_LIBRARY}")
endif ()
if (NiftiCLib_LIBRARY)
  list (APPEND NiftiCLib_LIBRARIES "${NiftiCLib_LIBRARY}")
endif ()

# ----------------------------------------------------------------------------
# import targets
if (NiftiCLib_znz_LIBRARY)
  if (NiftiCLib_USE_STATIC_LIB OR NiftiCLib_znz_LIBRARY MATCHES "\\.a$")
    add_library (niftiznz STATIC IMPORTED)
  else ()
    add_library (niftiznz SHARED IMPORTED)
  endif ()
  set_target_properties (
    niftiznz
    PROPERTIES
      IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
      IMPORTED_LOCATION                 "${NiftiCLib_znz_LIBRARY}"
      IMPORTED_LINK_INTERFACE_LIBRARIES "${ZLIB_LIBRARIES}"
  )
endif ()

if (NiftiCLib_LIBRARY)
  if (NiftiCLib_USE_STATIC_LIB OR NiftiCLib_LIBRARY MATCHES "\\.a$")
    add_library (niftiio STATIC IMPORTED)
  else ()
    add_library (niftiio SHARED IMPORTED)
  endif ()
  set_target_properties (
    niftiio
    PROPERTIES
      IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
      IMPORTED_LOCATION                 "${NiftiCLib_LIBRARY}"
  )
  if (TARGET niftiznz)
    set_target_properties (niftiio PROPERTIES IMPORTED_LINK_INTERFACE_LIBRARIES niftiznz)
  endif ()
endif ()

# ----------------------------------------------------------------------------
# aliases / backwards compatibility
if (NiftiCLib_INCLUDE_DIR)
  set (NiftiCLib_INCLUDE_DIRS "${NiftiCLib_INCLUDE_DIR}")
  if (NOT NiftiCLib_INCLUDE_DIR MATCHES "/nifti/?$")
    list (APPEND NiftiCLib_INCLUDE_DIRS "${NiftiCLib_INCLUDE_DIR}/nifti")
  endif ()
  set (NiftiCLib_INCLUDES "${NiftiCLib_INCLUDE_DIRS}")
endif ()

if (NiftiCLib_LIBRARY)
  set (NiftiCLib_LIB "${NiftiCLib_LIBRARY}")
endif ()

# ----------------------------------------------------------------------------
# reset CMake variables
set (CMAKE_FIND_LIBRARY_SUFFIXES ${NiftiCLib_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES})

# ----------------------------------------------------------------------------
# handle the QUIETLY and REQUIRED arguments and set *_FOUND to TRUE
# if all listed variables are found or TRUE
include (FindPackageHandleStandardArgs)

find_package_handle_standard_args (
  NiftiCLib
  REQUIRED_VARS
    NiftiCLib_INCLUDE_DIR
    NiftiCLib_LIBRARY
    NiftiCLib_znz_LIBRARY
)

set (NiftiCLib_FOUND ${NIFTICLIB_FOUND})

# ----------------------------------------------------------------------------
# set NiftiCLib_DIR
if (NOT NiftiCLib_DIR AND NiftiCLib_FOUND)
  string (REGEX REPLACE "include(/nifti)?/?" "" NiftiCLib_PREFIX "${NiftiCLib_INCLUDE_DIR}")
  set (NiftiCLib_DIR "${NiftiCLib_PREFIX}" CACHE PATH "Installation prefix for NiftiCLib." FORCE)
  unset (NiftiCLib_PREFIX)
endif ()
