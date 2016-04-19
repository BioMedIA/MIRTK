# ============================================================================
# Copyright (c) 2011-2012 University of Pennsylvania
# Copyright (c) 2013-2016 Andreas Schuh
# All rights reserved.
#
# See COPYING file for license information or visit
# https://cmake-basis.github.io/download.html#license
# ============================================================================

##############################################################################
# @file  CheckPublicHeaders.cmake
# @brief CMake script used to check whether public headers were added/removed.
#
# This script removes the deprecated public headers from the build tree and
# then throws a fatal error if public header files were added or removed to
# the project.
#
# @ingroup CMakeUtilities
##############################################################################

# ----------------------------------------------------------------------------
# check arguments
if (NOT CMAKE_FILE)
  message (FATAL_ERROR "Missing argument CMAKE_FILE!")
endif ()

if (NOT REFERENCE_FILE)
  message (FATAL_ERROR "Missing argument REFERENCE_FILE!")
endif ()

if (NOT BINARY_INCLUDE_DIR)
  message (FATAL_ERROR "Missing argument BINARY_INCLUDE_DIR!")
endif ()

if (NOT PROJECT_INCLUDE_DIRS)
  message (FATAL_ERROR "Missing argument PROJECT_INCLUDE_DIRS!")
endif ()

if (NOT VARIABLE_NAME)
  set (VARIABLE_NAME "PUBLIC_HEADERS")
endif ()

# ----------------------------------------------------------------------------
# compare files
execute_process (
  COMMAND ${CMAKE_COMMAND} -E compare_files "${CMAKE_FILE}" "${REFERENCE_FILE}"
  RESULT_VARIABLE EXIT_CODE
  OUTPUT_QUIET
  ERROR_QUIET
)

if (EXIT_CODE EQUAL 0)
  set (CMAKE_FILE_DIFFERS FALSE)
else ()
  set (CMAKE_FILE_DIFFERS TRUE)
endif ()

# ----------------------------------------------------------------------------
# remove obsolete public headers
if (CMAKE_FILE_DIFFERS)
  unset (${VARIABLE_NAME})
  include ("${CMAKE_FILE}")
  set (_HEADERS "${${VARIABLE_NAME}}")
  unset (${VARIABLE_NAME})
  include ("${REFERENCE_FILE}")
  foreach (H IN LISTS _HEADERS)
    list (FIND ${VARIABLE_NAME} "${H}" IDX)
    if (IDX EQUAL -1)
      # TODO: this hard coded /include path may break for custom include directories
      string (REGEX REPLACE "^.*/include/" "" H "${H}")
      string (REGEX REPLACE "\\.in$"       "" H "${H}")
      file (REMOVE "${BINARY_INCLUDE_DIR}/${H}")
    endif ()
  endforeach ()
endif ()

# ----------------------------------------------------------------------------
# remove files if different
if (CMAKE_FILE_DIFFERS)
  file (REMOVE "${CMAKE_FILE}")
  file (REMOVE "${REFERENCE_FILE}")
  file (REMOVE "${REFERENCE_FILE}.update")
endif ()

# ----------------------------------------------------------------------------
# fatal error if files added/removed
if (CMAKE_FILE_DIFFERS AND ERRORMSG)
  message (FATAL_ERROR "${ERRORMSG}")
endif ()
