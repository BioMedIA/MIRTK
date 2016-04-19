# ============================================================================
# Copyright (c) 2011-2012 University of Pennsylvania
# Copyright (c) 2013-2016 Andreas Schuh
# All rights reserved.
#
# See COPYING file for license information or visit
# https://cmake-basis.github.io/download.html#license
# ============================================================================

##############################################################################
# @file  PostprocessChangeLog.cmake
# @brief Script used to postprocess ChangeLog generated from commit history.
#
# Usage: cmake -DCHANGELOG_FILE=file -DINPUTFORMAT=[SVN|SVN2CL|GIT]
#                   -P PostprocessChangeLog.cmake
#
# @ingroup CMakeTools
##############################################################################

# ----------------------------------------------------------------------------
# check required arguments
if (NOT CHANGELOG_FILE)
  message (FATAL_ERROR "Missing CHANGELOG_FILE argument!")
endif ()

if (NOT INPUTFORMAT)
  set (INPUTFORMAT SVN)
endif ()

# ----------------------------------------------------------------------------
# read change log
file (READ "${CHANGELOG_FILE}" CHANGELOG)

# ----------------------------------------------------------------------------
# git log
if (INPUTFORMAT MATCHES "GIT")

  # remove git-svn-id entries from commit message body
  string (REGEX REPLACE "[ \n\r\t]*git-svn-id:[ \n\r]*[^@]*@[0-9]+[ \n\t]+[-0-9a-z]*" "" CHANGELOG "${CHANGELOG}")
  # group entries of same date and same author
  string (REGEX MATCHALL "[^\n]+(\n|$)" LINES "${CHANGELOG}")
  # clear changelog
  set (CHANGELOG)
  # process changelog line-by-line and leave out duplicate date and author lines
  set (PREV)
  foreach (LINE IN LISTS LINES)
    string (REGEX REPLACE "[\n]$" "" LINE "${LINE}")
    # Hack: For some reason the regular expression used above to split the
    #       change log file into lines produces lines with a newline followed
    #       by a semicolon. This could be a bug in CMake as well.
    string (REGEX REPLACE "[\n];" "\n" LINE "${LINE}")
    # a line with a date and author marks the beginning of a log entry
    if (LINE MATCHES "^[0-9][0-9][0-9][0-9]-[0-9][0-9]?-[0-9][0-9]? [a-zA-Z ]+$")
      if (NOT PREV OR NOT LINE STREQUAL PREV)
        if (PREV)
          set (CHANGELOG "${CHANGELOG}\n")
        endif ()
        set (CHANGELOG "${CHANGELOG}${LINE}\n\n")
        set (PREV "${LINE}")
      endif ()
    else ()
      set (CHANGELOG "${CHANGELOG}${LINE}\n")
    endif ()
  endforeach ()

endif ()

# ----------------------------------------------------------------------------
# write change log
file (WRITE "${CHANGELOG_FILE}" "${CHANGELOG}")
