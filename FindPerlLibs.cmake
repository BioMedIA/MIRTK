# ============================================================================
# Copyright (c) 2011-2012 University of Pennsylvania
# Copyright (c) 2013-2016 Andreas Schuh
# All rights reserved.
#
# See COPYING file for license information or visit
# http://opensource.andreasschuh.com/cmake-basis/download.html#license
# ============================================================================

##############################################################################
# @file  FindPerlLibs.cmake
# @brief Find Perl libraries. Fixes issue with CMake's default FindPerlLibs.
#
# @sa http://www.cmake.org/pipermail/cmake/2008-July/022638.html
#
# @ingroup CMakeFindModules
##############################################################################

# ----------------------------------------------------------------------------
# use CMake's FindPerlLibs.cmake module
set (PerlLibs_FIND_REQUIRED_BACKUP "${PerlLibs_FIND_REQUIRED}")
set (PerlLibs_FIND_QUIETLY_BACKUP  "${PerlLibs_FIND_QUIETLY}")
set (PerlLibs_FIND_REQUIRED FALSE)
set (PerlLibs_FIND_QUIETLY  TRUE)

include ("${CMAKE_ROOT}/Modules/FindPerlLibs.cmake")

set (PerlLibs_FIND_REQUIRED "${PerlLibs_FIND_REQUIRED_BACKUP}")
set (PerlLibs_FIND_QUIETLY  "${PerlLibs_FIND_QUIETLY_BACKUP}")
unset (PerlLibs_FIND_REQUIRED_BACKUP)
unset (PerlLibs_FIND_QUIETLY_BACKUP)

# ----------------------------------------------------------------------------
# try to fix any issues
if (PERL_EXECUTABLE)
  execute_process (
    COMMAND "${PERL_EXECUTABLE}" -MConfig -e "print \$Config{version}"
    OUTPUT_VARIABLE PERL_OUTPUT
    RESULT_VARIABLE PERL_RETURN_VALUE
  )
  if (PERL_RETURN_VALUE EQUAL 0)
    set (PERL_VERSION_STRING ${PERL_OUTPUT})
  endif ()
  if (PERL_VERSION_STRING MATCHES "([0-9]+)\\.([0-9]+)\\.([0-9]+)")
    set (PERL_VERSION_MAJOR "${CMAKE_MATCH_1}")
    set (PERL_VERSION_MINOR "${CMAKE_MATCH_2}")
    set (PERL_VERSION_PATCH "${CMAKE_MATCH_3}")
  else ()
    message (WARNING "Perl interpreter version string has unexpected format: ${PERL_VERSION_STRING}")
  endif ()

  # try to fix failure in PERL_INCLUDE_PATH
  if (PERL_INCLUDE_PATH MATCHES ".*-NOTFOUND")
    execute_process (
      COMMAND "${PERL_EXECUTABLE}" -MConfig -e "print \$Config{archlibexp}"
      OUTPUT_VARIABLE PERL_OUTPUT
      RESULT_VARIABLE PERL_RETURN_VALUE
    )
    if (NOT PERL_RETURN_VALUE)
      find_path (PERL_INCLUDE_PATH perl.h "${PERL_OUTPUT}/CORE")
    endif ()
  endif ()

  # try to fix failure in PERL_LIBRARY
  IF (PERL_LIBRARY MATCHES ".*-NOTFOUND")
    execute_process (
      COMMAND "${PERL_EXECUTABLE}" -MConfig -e "print \$Config{libperl}"
      OUTPUT_VARIABLE PERL_OUTPUT
      RESULT_VARIABLE PERL_RETURN_VALUE
    )
    if (NOT PERL_RETURN_VALUE)
      find_library (PERL_LIBRARY NAMES "${PERL_OUTPUT}" PATHS "${PERL_INCLUDE_PATH}")
    endif ()
  endif ()

  # unset local variables
  unset (PERL_OUTPUT)
  unset (PERL_RETURN_VALUE)
endif ()

# ----------------------------------------------------------------------------
# handle the QUIETLY and REQUIRED arguments and set *_FOUND to TRUE
# if all listed variables are found or TRUE
include (FindPackageHandleStandardArgs)

find_package_handle_standard_args (
  PerlLibs
  REQUIRED_VARS
    PERL_LIBRARY
    PERL_INCLUDE_PATH
  VERSION_VAR
    PERL_VERSION_STRING
)

if (NOT DEFINED PerlLibs_FOUND AND DEFINED PERLLIBS_FOUND)
  set (PerlLibs_FOUND "${PERLLIBS_FOUND}")
endif ()

# ----------------------------------------------------------------------------
# map names of variables from upper-case prefix to case-sensitive prefix
#
# This is normally automatically done by CMake BASIS find_package,
# but because "PerlLibs" sets variables also for "Perl" by *including*
# the FindPerl module, we have to do it here explicitly at least for
# those variables which are set by the FindPerl module.
if (NOT DEFINED Perl_FOUND AND DEFINED PERL_FOUND)
  set (Perl_FOUND "${PERL_FOUND}")
endif ()

if (Perl_FOUND)
  set (Perl_EXECUTABLE     "${PERL_EXECUTABLE}")
  set (Perl_VERSION_MAJOR  "${PERL_VERSION_MAJOR}")
  set (Perl_VERSION_MINOR  "${PERL_VERSION_MINOR}")
  set (Perl_VERSION_PATCH  "${PERL_VERSION_PATCH}")
  set (Perl_VERSION_STRING "${PERL_VERSION_STRING}")
endif ()

if (PerlLibs_FOUND)
  set (PerlLibs_LIBRARY        "${PERL_LIBRARY}")
  set (PerlLibs_INCLUDE_DIR    "${PERL_INCLUDE_PATH}")
  set (PerlLibs_VERSION_MAJOR  "${PERL_VERSION_MAJOR}")
  set (PerlLibs_VERSION_MINOR  "${PERL_VERSION_MINOR}")
  set (PerlLibs_VERSION_PATCH  "${PERL_VERSION_PATCH}")
  set (PerlLibs_VERSION_STRING "${PERL_VERSION_STRING}")
endif ()
