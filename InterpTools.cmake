# ============================================================================
# Copyright (c) 2011-2012 University of Pennsylvania
# Copyright (c) 2013-2016 Andreas Schuh
# All rights reserved.
#
# See COPYING file for license information or visit
# https://cmake-basis.github.io/download.html#license
# ============================================================================

##############################################################################
# @file  InterpTools.cmake
# @brief Script interpreter tools.
#
# @ingroup CMakeTools
##############################################################################

if (__BASIS_INTERPTOOLS_INCLUDED)
  return ()
else ()
  set (__BASIS_INTERPTOOLS_INCLUDED TRUE)
endif ()


## @addtogroup CMakeUtilities
# @{


# ----------------------------------------------------------------------------
## @brief Get version of Python interpreter.
#
# @param [out] ARGV1 If given, the named variable is set to the version string
#                    of the Python interpreter. Otherwise, the variables
#                    @c PYTHON_VERSION_STRING, @c PYTHON_VERSION_MAJOR,
#                    @c PYTHON_VERSION_MINOR, and @c PYTHON_VERSION_PATCH are
#                    set in the scope of the caller.
function (basis_get_python_version)
  if (PYTHON_EXECUTABLE)
    execute_process(
      COMMAND "${PYTHON_EXECUTABLE}" -E -c "import sys; sys.stdout.write(';'.join([str(x) for x in sys.version_info[:3]]))"
      OUTPUT_VARIABLE VERSION
      RESULT_VARIABLE RETVAL
      ERROR_QUIET
    )
    if (RETVAL EQUAL 0)
      string (REPLACE ";" "." VERSION_STRING "${VERSION}")
      list (GET VERSION 0 VERSION_MAJOR)
      list (GET VERSION 1 VERSION_MINOR)
      list (GET VERSION 2 VERSION_PATCH)
      if (VERSION_PATCH EQUAL 0)
        string (REGEX REPLACE "\\.0$" "" VERSION_STRING "${VERSION_STRING}")
      endif()
    else ()
      # sys.version predates sys.version_info
      execute_process (
        COMMAND "${PYTHON_EXECUTABLE}" -E -c "import sys; sys.stdout.write(sys.version)"
        OUTPUT_VARIABLE VERSION
        RESULT_VARIABLE RETVAL
        ERROR_QUIET
      )
      if (RETVAL EQUAL 0)
        string (REGEX REPLACE " .*" "" VERSION_STRING "${VERSION}")
        string (REGEX REPLACE "^([0-9]+)\\.[0-9]+.*" "\\1" VERSION_MAJOR "${VERSION_STRING}")
        string (REGEX REPLACE "^[0-9]+\\.([0-9])+.*" "\\1" VERSION_MINOR "${VERSION_STRING}")
        if (VERSION_STRING MATCHES "^[0-9]+\\.[0-9]+\\.[0-9]+.*")
          string (REGEX REPLACE "^[0-9]+\\.[0-9]+\\.([0-9]+).*" "\\1" VERSION_PATCH "${VERSION_STRING}")
        else ()
          set (VERSION_PATCH "0")
        endif()
      else ()
        # sys.version was first documented for Python 1.5
        set (VERSION_STRING "1.4")
        set (VERSION_MAJOR  "1")
        set (VERSION_MINOR  "4")
        set (VERSION_PATCH  "0")
      endif ()
    endif ()
  else ()
    set (VERSION_STRING "0.0")
    set (VERSION_MAJOR  "0")
    set (VERSION_MINOR  "0")
    set (VERSION_PATCH  "0")
  endif ()
  if (ARGC EQUAL 1)
    set (${ARGV0} "${VERSION_STRING}" PARENT_SCOPE)
  else ()
    set (PYTHON_VERSION_STRING "${VERSION_STRING}" PARENT_SCOPE)
    set (PYTHON_VERSION_MAJOR  "${VERSION_MAJOR}"  PARENT_SCOPE)
    set (PYTHON_VERSION_MINOR  "${VERSION_MINOR}"  PARENT_SCOPE)
    set (PYTHON_VERSION_PATCH  "${VERSION_PATCH}"  PARENT_SCOPE)
  endif ()
endfunction ()

# ----------------------------------------------------------------------------
## @brief Get version of Jython interpreter.
#
# @param [out] ARGV1 If given, the named variable is set to the version string
#                    of the Jython interpreter. Otherwise, the variables
#                    @c JYTHON_VERSION_STRING, @c JYTHON_VERSION_MAJOR,
#                    @c JYTHON_VERSION_MINOR, and @c JYTHON_VERSION_PATCH are
#                    set in the scope of the caller.
function (basis_get_jython_version)
  if (JYTHON_EXECUTABLE)
    execute_process(
      COMMAND "${JYTHON_EXECUTABLE}" -c "import sys; sys.stdout.write(';'.join([str(x) for x in sys.version_info[:3]]))"
      OUTPUT_VARIABLE VERSION
      RESULT_VARIABLE RETVAL
      ERROR_QUIET
    )
    if (RETVAL EQUAL 0)
      string (REPLACE ";" "." VERSION_STRING "${VERSION}")
      list (GET VERSION 0 VERSION_MAJOR)
      list (GET VERSION 1 VERSION_MINOR)
      list (GET VERSION 2 VERSION_PATCH)
      if (VERSION_PATCH EQUAL 0)
        string (REGEX REPLACE "\\.0$" "" VERSION_STRING "${VERSION_STRING}")
      endif()
    else ()
      # sys.version predates sys.version_info
      execute_process (
        COMMAND "${JYTHON_EXECUTABLE}" -c "import sys; sys.stdout.write(sys.version)"
        OUTPUT_VARIABLE VERSION
        RESULT_VARIABLE RETVAL
        ERROR_QUIET
      )
      if (RETVAL EQUAL 0)
        string (REGEX REPLACE " .*" "" VERSION_STRING "${VERSION}")
        string (REGEX REPLACE "^([0-9]+)\\.[0-9]+.*" "\\1" VERSION_MAJOR "${VERSION_STRING}")
        string (REGEX REPLACE "^[0-9]+\\.([0-9])+.*" "\\1" VERSION_MINOR "${VERSION_STRING}")
        if (VERSION_STRING MATCHES "^[0-9]+\\.[0-9]+\\.[0-9]+.*")
          string (REGEX REPLACE "^[0-9]+\\.[0-9]+\\.([0-9]+).*" "\\1" VERSION_PATCH "${VERSION_STRING}")
        else ()
          set (VERSION_PATCH "0")
        endif()
      else ()
        set (VERSION_STRING "0.0")
        set (VERSION_MAJOR  "0")
        set (VERSION_MINOR  "0")
        set (VERSION_PATCH  "0")
      endif ()
    endif ()
  else ()
    set (VERSION_STRING "0.0")
    set (VERSION_MAJOR  "0")
    set (VERSION_MINOR  "0")
    set (VERSION_PATCH  "0")
  endif ()
  if (ARGC EQUAL 1)
    set (${ARGV0} "${VERSION_STRING}" PARENT_SCOPE)
  else ()
    set (JYTHON_VERSION_STRING "${VERSION_STRING}" PARENT_SCOPE)
    set (JYTHON_VERSION_MAJOR  "${VERSION_MAJOR}"  PARENT_SCOPE)
    set (JYTHON_VERSION_MINOR  "${VERSION_MINOR}"  PARENT_SCOPE)
    set (JYTHON_VERSION_PATCH  "${VERSION_PATCH}"  PARENT_SCOPE)
  endif ()
endfunction ()

# ----------------------------------------------------------------------------
## @brief Get version of Perl interpreter.
#
# @param [out] ARGV1 If given, the named variable is set to the version string
#                    of the Perl interpreter. Otherwise, the variables
#                    @c PERL_VERSION_STRING, @c PERL_VERSION_MAJOR,
#                    @c PERL_VERSION_MINOR, and @c PERL_VERSION_PATCH are
#                    set in the scope of the caller.
function (basis_get_perl_version)
  if (PERL_EXECUTABLE)
    execute_process (COMMAND "${PERL_EXECUTABLE}" --version OUTPUT_VARIABLE VERSION)
  else ()
    set (VERSION)
  endif ()
  if (VERSION MATCHES "[( ]v([0-9]+)\\.([0-9]+)\\.([0-9]+)[ )]")
    set (VERSION_MAJOR "${CMAKE_MATCH_1}")
    set (VERSION_MINOR "${CMAKE_MATCH_2}")
    set (VERSION_PATCH "${CMAKE_MATCH_3}")
    set (VERSION_STRING "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")
  else ()
    set (VERSION_STRING "0.0")
    set (VERSION_MAJOR  "0")
    set (VERSION_MINOR  "0")
    set (VERSION_PATCH  "0")
  endif ()
  if (ARGC EQUAL 1)
    set (${ARGV0} "${VERSION_STRING}" PARENT_SCOPE)
  else ()
    set (PERL_VERSION_STRING "${VERSION_STRING}" PARENT_SCOPE)
    set (PERL_VERSION_MAJOR  "${VERSION_MAJOR}"  PARENT_SCOPE)
    set (PERL_VERSION_MINOR  "${VERSION_MINOR}"  PARENT_SCOPE)
    set (PERL_VERSION_PATCH  "${VERSION_PATCH}"  PARENT_SCOPE)
  endif ()
endfunction ()

# ----------------------------------------------------------------------------
## @brief Get version of Bash interpreter.
#
# @param [out] ARGV1 If given, the named variable is set to the version string
#                    of the Bash interpreter. Otherwise, the variables
#                    @c BASH_VERSION_STRING, @c BASH_VERSION_MAJOR,
#                    @c BASH_VERSION_MINOR, and @c BASH_VERSION_PATCH are
#                    set in the scope of the caller.
function (basis_get_bash_version)
  if (BASH_EXECUTABLE)
    execute_process (COMMAND "${BASH_EXECUTABLE}" --version OUTPUT_VARIABLE VERSION)
  else ()
    set (VERSION)
  endif ()
  if (VERSION MATCHES "version ([0-9]+)\\.([0-9]+)\\.([0-9]+)")
    set (VERSION_MAJOR "${CMAKE_MATCH_1}")
    set (VERSION_MINOR "${CMAKE_MATCH_2}")
    set (VERSION_PATCH "${CMAKE_MATCH_3}")
    set (VERSION_STRING "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")
  else ()
    set (VERSION_STRING "0.0")
    set (VERSION_MAJOR  "0")
    set (VERSION_MINOR  "0")
    set (VERSION_PATCH  "0")
  endif ()
  if (ARGC EQUAL 1)
    set (${ARGV0} "${VERSION_STRING}" PARENT_SCOPE)
  else ()
    set (BASH_VERSION_STRING "${VERSION_STRING}" PARENT_SCOPE)
    set (BASH_VERSION_MAJOR  "${VERSION_MAJOR}"  PARENT_SCOPE)
    set (BASH_VERSION_MINOR  "${VERSION_MINOR}"  PARENT_SCOPE)
    set (BASH_VERSION_PATCH  "${VERSION_PATCH}"  PARENT_SCOPE)
  endif ()
endfunction ()


## @}
# end of Doxygen group
