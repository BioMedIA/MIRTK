# ============================================================================
# Copyright (c) 2011-2012 University of Pennsylvania
# Copyright (c) 2013-2016 Andreas Schuh
# All rights reserved.
#
# See COPYING file for license information or visit
# https://cmake-basis.github.io/download.html#license
# ============================================================================

##############################################################################
# @file  ExecuteProcess.cmake
# @brief Execute process using CMake script mode.
#
# This CMake script can be used as argument for the -P option of cmake, when
# another command shall be executed by CMake, for example, as custom build
# command. The advantage of using this script is that all options of the
# CMake command execute_process() can be used, i.e., a timeout can be
# specified.
#
# The arguments of the execute_process() command have to be specified via
# the -D option on the command line of cmake before the -P \<this script\>
# option is given. The name of the CMake variables must be equal the
# name of the arguments to the execute_process() command.
#
# Arguments of execute_process() which are considered:
#
# - @b COMMAND
# - @b WORKING_DIRECTORY
# - @b TIMEOUT
# - @b OUTPUT_FILE
# - @b ERROR_FILE
# - @b OUTPUT_QUIET
# - @b ERROR_QUIET
# - @b OUTPUT_STRIP_TRAILING_WHITESPACE
# - @b ERROR_STRIP_TRAILING_WHITESPACE
#
# Additionally, matching expressions (separated by ';') to identify error messages
# in the output streams stdout and stderr can be specified by the input argument
# ERROR_EXPRESSION. When the output of the executed command matches one of
# the error expressions, a fatal error message is displayed causing CMake to
# return the exit code 1.
#
# Setting VERBOSE to true enables verbose output messages.
#
# When the input argument LOG_ARGS evaluates to true, the values of COMMAND,
# WORKING_DIRECTORY, and TIMEOUT are added to the top of the output files
# specified by OUTPUT_FILE and ERROR_FILE.
#
# The arguments ARGS and ARGS_FILE can be used to specify (additional) command
# arguments. The content of the text file ARGS_FILE is read when it this file
# exists. Separate lines of this file are considered single arguments.
# The arguments specified by ARGS and ARGS_FILE are concatenated where the
# arguments given by ARGS follow after the ones read from the ARGS_FILE.
# All occurences of the string 'ARGS' in the COMMAND are replaced by these
# arguments. If no such string is present, the arguments are simply passed
# to the execute_process() command as its ARGS argument.
# The argument ARGS_SEPARATOR specifies the separator used to separate the
# arguments given by ARGS and ARGS_FILE when the 'ARGS' string in COMMAND
# is replaced. By default, it is set to ';'.
#
# Example:
# @code
# cmake -DCOMMAND='ls;-l' -DWORKING_DIRECTORY='/' -DTIMEOUT=60
#       -P ExecuteProcess.cmake
# @endcode
#
# The output of the executed process can further be searched for error expressions
# specified by the ERROR_EXPRESSION variable. If the process output matches
# this expression, a fatal CMake error is raised.
#
# Certain errors may be temporary such as in particular a license checkout
# error of the MATLAB Compiler. Such errors can be filtered out using the
# RETRY_EXPRESSION variable. If such error is detected, this script sleeps for
# RETRY_DELAY seconds and then executes the process again. This is done
# at maximum RETRY_ATTEMPTS times. If the retry attempts are exceeded, a
# fatal CMake error is raised instead.
#
# @sa http://www.cmake.org/cmake/help/cmake2.6docs.html#command:execute_process
#
# @ingroup CMakeUtilities
##############################################################################

if (POLICY CMP0053)
  cmake_policy (SET CMP0053 NEW)
endif ()

# ----------------------------------------------------------------------------
# unset environment variables that may cause problems otherwise
unset (ENV{PYTHONHOME})
unset (ENV{PYTHONPATH})

# ----------------------------------------------------------------------------
# initialize arguments
set (CONFIGURED_COMMAND "@BUILD_CMD@")

string (REPLACE "#" "@" AT_BUILD_CMD_AT "#BUILD_CMD#")
if (CONFIGURED_COMMAND AND NOT CONFIGURED_COMMAND STREQUAL "${AT_BUILD_CMD_AT}")
  set (COMMAND "${CONFIGURED_COMMAND}")
elseif (NOT COMMAND)
  message (FATAL_ERROR "No command specified for execute_process(): use \"-DCOMMAND=cmd\"")
endif ()

if (NOT ARGS_SEPARATOR)
  set (ARGS_SEPARATOR ";")
endif ()

if (ARGS_FILE)
  include ("${ARGS_FILE}" OPTIONAL)

  if (ARGS AND DEFINED ${ARGS})
    set (ARGS "${${ARGS}}")
  else ()
    set (ARGS "")
  endif ()
endif ()

if ("${COMMAND}" MATCHES "ARGS")
  string (REPLACE ";" "${ARGS_SEPARATOR}" TMP "${ARGS}")
  string (REPLACE "ARGS" "${ARGS_SEPARATOR}${TMP}" COMMAND "${COMMAND}")
  set (ARGS)
endif ()

set (EXECUTE_PROCESS_ARGS "COMMAND" "${COMMAND}")

if (ARGS)
  list (APPEND EXECUTE_PROCESS_ARGS "ARGS" "${ARGS}")
endif ()

list (APPEND EXECUTE_PROCESS_ARGS "RESULT_VARIABLE" "RETVAL")

if (TIMEOUT)
  list (APPEND EXECUTE_PROCESS_ARGS "TIMEOUT" "${TIMEOUT}")
endif ()

if (WORKING_DIRECTORY)
  list (APPEND EXECUTE_PROCESS_ARGS "WORKING_DIRECTORY" "${WORKING_DIRECTORY}")
endif ()

if (OUTPUT_FILE)
  list (APPEND EXECUTE_PROCESS_ARGS "OUTPUT_FILE" "${OUTPUT_FILE}")
endif ()

if (ERROR_FILE)
  list (APPEND EXECUTE_PROCESS_ARGS "ERROR_FILE" "${ERROR_FILE}")
endif ()

if (OUTPUT_QUIET)
  list (APPEND EXECUTE_PROCESS_ARGS "OUTPUT_QUIET")
endif ()

if (ERROR_QUIET)
  list (APPEND EXECUTE_PROCESS_ARGS "ERROR_QUIET")
endif ()

if (OUTPUT_STRIP_TRAILING_WHITESPACE)
  list (APPEND EXECUTE_PROCESS_ARGS "OUTPUT_STRIP_TRAILING_WHITESPACE")
endif ()

if (ERROR_STRIP_TRAILING_WHITESPACE)
  list (APPEND EXECUTE_PROCESS_ARGS "ERROR_STRIP_TRAILING_WHITESPACE")
endif ()

if (NOT OUTPUT_FILE)
  list (APPEND EXECUTE_PROCESS_ARGS "OUTPUT_VARIABLE" "STDOUT")
endif ()

if (NOT ERROR_FILE)
  list (APPEND EXECUTE_PROCESS_ARGS "ERROR_VARIABLE"  "STDERR")
endif ()

if (NOT RETRY_ATTEMPTS)
  set (RETRY_ATTEMPTS 0)
endif ()

if (NOT RETRY_DELAY)
  set (RETRY_DELAY 60)
endif ()

# --------------------------------------------------------------------------
# verbose message of command to execute
set (CMD)
foreach (ARG ${COMMAND})
  if (CMD)
    set (CMD "${CMD} ")
  endif ()
  if (ARG MATCHES " ")
    set (CMD "${CMD}\"${ARG}\"")
  else ()
    set (CMD "${CMD}${ARG}")
  endif ()
endforeach ()

if (VERBOSE)
  message ("${CMD}")
endif ()

# ----------------------------------------------------------------------------
# execute command
set (RETRY TRUE) # execute process at least once
while (RETRY)
  # --------------------------------------------------------------------------
  # whether retry is required is decided upon after process execution
  set (RETRY FALSE)

  # --------------------------------------------------------------------------
  # execute process
  execute_process (${EXECUTE_PROCESS_ARGS})

  # --------------------------------------------------------------------------
  # read in output from log files
  if (OUTPUT_FILE)
    file (READ "${OUTPUT_FILE}" STDOUT)
  endif ()

  if (ERROR_FILE)
    file (READ "${ERROR_FILE}" STDERR)
  endif ()

  # --------------------------------------------------------------------------
  # parse output for recoverable errors
  foreach (EXPR IN LISTS RETRY_EXPRESSION)
    if (STDOUT MATCHES "${EXPR}" OR STDERR MATCHES "${EXPR}")
      if (RETRY_ATTEMPTS GREATER 0)
        message ("Process output matches \"${EXPR}\", retry in ${RETRY_DELAY} seconds")
        # retry
        math (EXPR RETRY_ATTEMPTS "${RETRY_ATTEMPTS} - 1")
        set (RETRY TRUE)
      else ()
        # no retries left
        set (RETVAL 1)
      endif ()
      break ()
    endif ()
  endforeach ()

  # --------------------------------------------------------------------------
  # sleep for given amount of seconds before retrying to execute process
  if (RETRY)
    # use sleep command if available, i.e., on Unix and also on Windows
    # if the Windows 2003 Resource Kit is installed
    find_program (SLEEP sleep)
    if (SLEEP)
      execute_process (
        COMMAND "${SLEEP}" ${RETRY_DELAY}
        TIMEOUT ${RETRY_DELAY}
        ERROR_QUIET OUTPUT_QUIET
      )
    else ()
      # work-around using ping command if sleep command not available, i.e.,
      # on Windows where the Windows 2003 Resource Kit is not installed
      # See http://malektips.com/dos0017.html
      find_program (PING ping)
      if (WIN32)
        execute_process (
          COMMAND ${PING} 127.0.0.1 -n ${RETRY_DELAY} -w 1000
          TIMEOUT ${RETRY_DELAY}
          ERROR_QUIET OUTPUT_QUIET
        )
      else ()
        # usually not required as sleep command should be available
        execute_process (
          COMMAND ${PING} 127.0.0.1 -c ${RETRY_DELAY} -W 1000
          TIMEOUT ${RETRY_DELAY}
          ERROR_QUIET OUTPUT_QUIET
        )
      endif ()
    else ()
      message (WARNING "Cannot delay retries as neither sleep nor ping command is available!")
    endif ()
  endif ()
endwhile ()

# ----------------------------------------------------------------------------
# parse output for errors
foreach (EXPR IN LISTS ERROR_EXPRESSION)
  if (STDOUT MATCHES "${EXPR}" OR STDERR MATCHES "${EXPR}")
    set (RETVAL 1)
    break ()
  endif ()
endforeach ()

# ----------------------------------------------------------------------------
# give some hints for known error messages
if (NOT RETVAL EQUAL 0)
  set (HINTS)
  # MATLAB Compiler
  if (CMD MATCHES "mcc")
    if (STDERR MATCHES "The file.*appears to be a MEX-file.[ ]+It shadows the M-file.*but will not execute properly at runtime, as it does not export a function.*named 'mexFunction.'")
      set (HINTS "${HINTS}

  Note: The error that a MEX-file would shadow a M-file can be caused by a shared library that
        is required by the MEX-file but not found in the search path of the dynamic loader.
")
    endif ()
  endif ()
  # append hints to output
  if ("${STDOUT}" STREQUAL "${STDERR}")
    set (STDERR "${STDERR}${HINTS}")
    set (STDOUT "${STDOUT}${HINTS}")
  else ()
    set (STDERR "${STDERR}${HINTS}")
  endif ()
endif ()

# ----------------------------------------------------------------------------
# prepand command to log file
if (LOG_ARGS)
  if (OUTPUT_FILE)
    set (TMP "Command: ${CMD}\n\nWorking directory: ${WORKING_DIRECTORY}\n\n")
    set (TMP "${TMP}Timeout: ${TIMEOUT}\n\nOutput:\n\n${STDOUT}")
    file (WRITE "${OUTPUT_FILE}" "${TMP}")
    unset (TMP)
  endif ()

  if (ERROR_FILE AND NOT "${ERROR_FILE}" STREQUAL "${OUTPUT_FILE}")
    set (TMP "Command: ${CMD}\n\nWorking directory: ${WORKING_DIRECTORY}\n\n")
    set (TMP "${TMP}Timeout: ${TIMEOUT}\n\nOutput:\n\n${STDERR}")
    file (WRITE "${ERROR_FILE}" "${TMP}")
    unset (TMP)
  endif ()
endif ()

# ----------------------------------------------------------------------------
# print error message (and exit with exit code 1) on error
if (NOT RETVAL EQUAL 0)
  
  # print error
  if ("${STDOUT}" STREQUAL "${STDERR}")
    message (
      FATAL_ERROR "
Command: ${CMD}
Working directory: ${WORKING_DIRECTORY}
Timeout: ${TIMEOUT}
Output:
${STDOUT}")
  else ()
    message (
      FATAL_ERROR "
Command: ${CMD}
Working directory: ${WORKING_DIRECTORY}
Timeout: ${TIMEOUT}
Output (stdout):
${STDOUT}
Output (stderr):
${STDERR}")
  endif ()
endif ()
