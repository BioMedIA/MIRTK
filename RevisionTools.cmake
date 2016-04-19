# ============================================================================
# Copyright (c) 2011-2012 University of Pennsylvania
# Copyright (c) 2013-2016 Andreas Schuh
# All rights reserved.
#
# See COPYING file for license information or visit
# https://cmake-basis.github.io/download.html#license
# ============================================================================

##############################################################################
# @file  RevisionTools.cmake
# @brief CMake functions and macros related to revision control systems.
#
# @ingroup CMakeTools
##############################################################################

if (__BASIS_REVISIONTOOLS_INCLUDED)
  return ()
else ()
  set (__BASIS_REVISIONTOOLS_INCLUDED TRUE)
endif ()


# ============================================================================
# required commands
# ============================================================================

find_package (Subversion QUIET)
find_package (Git        QUIET)


## @addtogroup CMakeUtilities
#  @{


# ============================================================================
# Subversion
# ============================================================================

# ----------------------------------------------------------------------------
## @brief Get current revision of file or directory.
#
# @param [in]  URL  Absolute path to directory or file. May also be a URL to the
#                   directory or file in the repository. A leading "file://" is
#                   automatically removed such that the svn command treats it as a
#                   local path.
# @param [out] REV  The revision number of URL. If URL is not under revision
#                   control or Subversion_SVN_EXECUTABLE is invalid, "0" is returned.
#
# @returns Sets @p REV to the revision of the working copy/repository
#          at URL @p URL.
function (basis_svn_get_revision URL REV)
  set (OUT 0)
  if (Subversion_SVN_EXECUTABLE)
    # remove "file://" from URL
    string (REGEX REPLACE "file://" "" TMP "${URL}")
    # retrieve SVN info
    execute_process (
      COMMAND         "${Subversion_SVN_EXECUTABLE}" info "${TMP}"
      OUTPUT_VARIABLE OUT
      ERROR_QUIET
      OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if (BASIS_DEBUG)
      message ("** basis_svn_get_revision()")
      message ("**   svn info: ${OUT}")
    endif ()
    # extract revision
    if (OUT MATCHES "^(.*\n)?Revision: ([^\n]+).*" AND NOT CMAKE_MATCH_2 STREQUAL "")
      set (OUT "${CMAKE_MATCH_2}")
    else ()
      set (OUT 0)
    endif ()
  endif ()
  # return
  set ("${REV}" "${OUT}" PARENT_SCOPE)
endfunction ()

# ----------------------------------------------------------------------------
## @brief Get revision number when directory or file was last changed.
#
# @param [in]  URL  Absolute path to directory or file. May also be a URL to the
#                   directory or file in the repository. A leading "file://" is
#                   automatically removed such that the svn command treats it as a
#                   local path.
# @param [out] REV  Revision number when URL was last modified. If URL is not
#                   under Subversion control or Subversion_SVN_EXECUTABLE is invalid,
#                   "0" is returned.
#
# @returns Sets @p REV to revision number at which the working copy/repository
#          specified by the URL @p URL was last modified.
function (basis_svn_get_last_changed_revision URL REV)
  set (OUT 0)
  if (Subversion_SVN_EXECUTABLE)
    # remove "file://" from URL
    string (REGEX REPLACE "file://" "" TMP "${URL}")
    # retrieve SVN info
    execute_process (
      COMMAND         "${Subversion_SVN_EXECUTABLE}" info "${TMP}"
      OUTPUT_VARIABLE OUT
      ERROR_QUIET
      OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if (BASIS_DEBUG)
      message ("** basis_svn_get_revision()")
      message ("**   svn info: ${OUT}")
    endif ()
    # extract last changed revision
    if (OUT MATCHES "^(.*\n)?Last Changed Rev: ([^\n]+).*" AND NOT CMAKE_MATCH_2 STREQUAL "")
      set (OUT "${CMAKE_MATCH_2}")
    else ()
      set (OUT 0)
    endif ()
  endif ()
  # return
  set ("${REV}" "${OUT}" PARENT_SCOPE)
endfunction ()

# ----------------------------------------------------------------------------
## @brief Get status of revision controlled file.
#
# @param [in]  URL    Absolute path to directory or file. May also be a URL to
#                     the directory or file in the repository.
#                     A leading "file://" will be removed such that the svn
#                     command treats it as a local path.
# @param [out] STATUS The status of URL as returned by 'svn status'.
#                     If the local directory or file is unmodified, an
#                     empty string is returned. An empty string is also
#                     returned when Subversion_SVN_EXECUTABLE is invalid.
#
# @returns Sets @p STATUS to the output of the <tt>svn info</tt> command.
function (basis_svn_status URL STATUS)
  if (Subversion_SVN_EXECUTABLE)
    # remove "file://" from URL
    string (REGEX REPLACE "file://" "" TMP "${URL}")

    # retrieve SVN status of URL
    execute_process (
      COMMAND         "${Subversion_SVN_EXECUTABLE}" status "${TMP}"
      OUTPUT_VARIABLE OUT
      ERROR_QUIET
    )

    # return
    set ("${STATUS}" "${OUT}" PARENT_SCOPE)
  else ()
    set ("${STATUS}" "" PARENT_SCOPE)
  endif ()
endfunction ()

# ============================================================================
# Git
# ============================================================================

# ----------------------------------------------------------------------------
# @brief Determine whether or not a given directory is a Git repository
function (basis_is_git_repository FLAG DIR)
  if (GITCOMMAND AND NOT GIT_EXECUTABLE)
    set (GIT_EXECUTABLE GITCOMMAND)
  endif ()
  if (GIT_EXECUTABLE)
    execute_process (
      COMMAND "${GIT_EXECUTABLE}" rev-parse
      WORKING_DIRECTORY "${DIR}"
      RESULT_VARIABLE RETVAL
      OUTPUT_QUIET
      ERROR_QUIET
    )
  else ()
    set (RETVAL 1)
  endif ()
  if (RETVAL EQUAL 0)
    set (${FLAG} "TRUE" PARENT_SCOPE)
  else ()
    set (${FLAG} "FALSE" PARENT_SCOPE)
  endif ()
endfunction ()

# ----------------------------------------------------------------------------
## @brief Get HEAD commit SHA of file or directory.
#
# @param [in]  URL  Absolute path to repository directory or single file.
# @param [out] REV  The short commit SHA when URL was last modified. If URL
#                   is not under Git control or GIT_EXECUTABLE is invalid,
#                   "0" is returned.
# @param [in]  ARGN Length of commit SHA to return.
#
# @returns Sets @p REV either to the HEAD commit SHA of the repository at
#          directory @p URL or the last commit which modified the file.
function (basis_git_get_revision URL REV)
  if (ARGC GREATER 3)
    message (FATAL_ERROR "basis_git_get_revision: Too many arguments")
  endif ()
  set (OUT 0)
  if (GITCOMMAND AND NOT GIT_EXECUTABLE)
    set (GIT_EXECUTABLE GITCOMMAND)
  endif ()
  if (GIT_EXECUTABLE)
    # remove "file://" from URL
    string (REGEX REPLACE "file://" "" DIR "${URL}")
    # retrieve Git commit SHA of HEAD
    if (IS_DIRECTORY "${DIR}")
      execute_process (
        COMMAND "${GIT_EXECUTABLE}" rev-parse HEAD
        WORKING_DIRECTORY "${DIR}"
        OUTPUT_VARIABLE OUT
        ERROR_QUIET
        OUTPUT_STRIP_TRAILING_WHITESPACE
      )
      if (BASIS_DEBUG)
        message ("** basis_git_get_revision()")
        message ("**   DIR: ${DIR}")
        message ("**   OUT: ${OUT}")
      endif ()
    # retrieve Git commit SHA when file was last modified
    else ()
      basis_get_filename_component(OBJ "${DIR}" NAME)
      basis_get_filename_component(DIR "${DIR}" PATH)
      execute_process (
        COMMAND "${GIT_EXECUTABLE}" log -n 1 -- "${OBJ}"
        WORKING_DIRECTORY "${DIR}"
        OUTPUT_VARIABLE OUT
        ERROR_QUIET
        OUTPUT_STRIP_TRAILING_WHITESPACE
      )
      if (BASIS_DEBUG)
        message ("** basis_git_get_revision()")
        message ("**   DIR: ${DIR}")
        message ("**   OBJ: ${OBJ}")
        message ("**   OUT: ${OUT}")
      endif ()
      # extract commit SHA
      if (OUT MATCHES "commit ([0-9a-f]+).*" AND NOT CMAKE_MATCH_1 STREQUAL "")
        set (OUT "${CMAKE_MATCH_1}")
      else ()
        set (OUT 0)
      endif ()
    endif ()
  endif ()
  if (ARGC EQUAL 3)
    string (SUBSTRING "${OUT}" 0 "${ARGV2}" OUT)
  endif ()
  # return
  set ("${REV}" "${OUT}" PARENT_SCOPE)
endfunction ()


# ============================================================================
# Meta
# ============================================================================

# ----------------------------------------------------------------------------
## @brief Get revision of file or directory.
#
# @param [in]  URL  Absolute path to directory or single file.
# @param [out] REV  Revision number when directory / repository or file was
#                   last modified, "0" is returned when no known revision
#                   control system is used or revision command not found.
#
# @returns Revision when directory of file was last modified.
function (basis_get_revision URL REV)
  # remove "file://" from URL
  string (REGEX REPLACE "file://" "" DIR "${URL}")
  # get directory path
  if (NOT IS_DIRECTORY "${URL}")
    basis_get_filename_component (DIR "${URL}" PATH)
  else ()
    set (DIR "${URL}")
  endif ()
  # check if directory is part of a Git repository
  basis_is_git_repository (IS_GIT_REPOSITORY "${DIR}")
  if (IS_GIT_REPOSITORY)
    basis_git_get_revision ("${URL}" OUT 7)
  else ()
    basis_svn_get_last_changed_revision ("${URL}" OUT)
  endif ()
  # return revision
  set ("${REV}" "${OUT}" PARENT_SCOPE)
endfunction ()


## @}
# end of Doxygen group
