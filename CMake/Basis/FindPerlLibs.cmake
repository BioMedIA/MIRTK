# ============================================================================
# Copyright (c) 2011-2012 University of Pennsylvania
# Copyright (c) 2013-2014 Andreas Schuh
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
set (CMAKE_MODULE_PATH_BACKUP "${CMAKE_MODULE_PATH}")
set (PerlLIbs_FIND_REQUIRED_BACKUP "${PerlLibs_FIND_REQUIRED}")
set (PerlLIbs_FIND_QUIETLY_BACKUP "${PerlLibs_FIND_QUIETLY}")
set (CMAKE_MODULE_PATH)
set (PerlLibs_FIND_REQUIRED FALSE)
set (PerlLibs_FIND_QUIETLY  TRUE)

find_package (PerlLibs)

set (CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH_BACKUP}")
set (PerlLIbs_FIND_REQUIRED "${PerlLibs_FIND_REQUIRED_BACKUP}")
set (PerlLIbs_FIND_QUIETLY "${PerlLibs_FIND_QUIETLY_BACKUP}")
set (PerlLibs_FIND_REQUIRED_BACKUP)
set (PerlLibs_FIND_QUIETLY_BACKUP)
set (CMAKE_MODULE_PATH_BACKUP)

# ----------------------------------------------------------------------------
# try to fix any issues
if (PERL_EXECUTABLE)
  execute_process (
    COMMAND "${PERL_EXECUTABLE}" -MConfig -e "print \$Config{version}"
    OUTPUT_VARIABLE PERL_OUTPUT
    RESULT_VARIABLE PERL_RETURN_VALUE
  )
  if (PERL_RETURN_VALUE EQUAL 0)
    set (PERL_VERSION ${PERL_OUTPUT})
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
    PERL_VERSION
)

set (PerlLibs_FOUND "${PERLLIBS_FOUND}")
