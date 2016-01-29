# ============================================================================
# Copyright (c) 2011-2012 University of Pennsylvania
# Copyright (c) 2013-2014 Andreas Schuh
# All rights reserved.
#
# See COPYING file for license information or visit
# http://opensource.andreasschuh.com/cmake-basis/download.html#license
# ============================================================================

##############################################################################
# @file  FindBASH.cmake
# @brief Find BASH interpreter.
#
# Sets the CMake variables @c BASH_FOUND, @c BASH_EXECUTABLE,
# @c BASH_VERSION_STRING, @c BASH_VERSION_MAJOR, @c BASH_VERSION_MINOR, and
# @c BASH_VERSION_PATCH.
#
# @ingroup CMakeFindModules
##############################################################################

# ----------------------------------------------------------------------------
# find BASH executable
find_program (BASH_EXECUTABLE bash)
mark_as_advanced (BASH_EXECUTABLE)

# ----------------------------------------------------------------------------
# get version of found BASH executable
if (BASH_EXECUTABLE)
  execute_process (COMMAND "${BASH_EXECUTABLE}" --version OUTPUT_VARIABLE _BASH_STDOUT ERROR_VARIABLE _BASH_STDERR)
  if (_BASH_STDOUT MATCHES "version ([0-9]+)\\.([0-9]+)\\.([0-9]+)")
    set (BASH_VERSION_MAJOR "${CMAKE_MATCH_1}")
    set (BASH_VERSION_MINOR "${CMAKE_MATCH_2}")
    set (BASH_VERSION_PATCH "${CMAKE_MATCH_3}")
    set (BASH_VERSION_STRING "${BASH_VERSION_MAJOR}.${BASH_VERSION_MINOR}.${BASH_VERSION_PATCH}")
  else ()
    message (WARNING "Failed to determine version of Bash interpreter (${BASH_EXECUTABLE})! Error:\n${_BASH_STDERR}")
  endif ()
  unset (_BASH_STDOUT)
  unset (_BASH_STDERR)
endif ()

# ----------------------------------------------------------------------------
# handle the QUIETLY and REQUIRED arguments and set *_FOUND to TRUE
# if all listed variables are found or TRUE
include (FindPackageHandleStandardArgs)

find_package_handle_standard_args (
  BASH
  REQUIRED_VARS
    BASH_EXECUTABLE
)
