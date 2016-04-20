# ============================================================================
# Copyright (c) 2011-2012 University of Pennsylvania
# Copyright (c) 2013-2016 Andreas Schuh
# All rights reserved.
#
# See COPYING file for license information or visit
# https://cmake-basis.github.io/download.html#license
# ============================================================================

##############################################################################
# @file  MatlabTools.cmake
# @brief Enables use of MATLAB Compiler and build of MEX-files.
#
# @ingroup CMakeTools
##############################################################################

if (__BASIS_MATLABTOOLS_INCLUDED)
  return ()
else ()
  set (__BASIS_MATLABTOOLS_INCLUDED TRUE)
endif ()


# ============================================================================
# modules
# ============================================================================

# Note: Required because generate_matlab_executable.cmake uses this module.

include (CMakeParseArguments)
include ("${CMAKE_CURRENT_LIST_DIR}/CommonTools.cmake")
include ("${CMAKE_CURRENT_LIST_DIR}/UtilitiesTools.cmake")

# ============================================================================
# options
# ============================================================================

## @addtogroup BasisSettings
#  @{

## @brief Enable/Disable compilation of MATLAB sources if the MATLAB Compiler is available.
option (BASIS_COMPILE_MATLAB "Enable compilation of MATLAB sources if MATLAB Compiler (mcc) is available." ON)
mark_as_advanced (BASIS_COMPILE_MATLAB)

# ----------------------------------------------------------------------------
## @brief Add global MATLAB MEX-script options to CMake cache
#
# @par MATLAB MEX-script options:
# <table border="0">
#   <tr>
#     @tp @b BASIS_MEX_TIMEOUT @endtp
#     <td>Timeout for MEX script execution.</td>
#   </tr>
#   <tr>
#     @tp @b BASIS_MEX_FLAGS @endtp
#     <td>Compile flags used to build MEX-files using the MEX script.</td>
#   </tr>
# </table>
macro(basis_add_mex_options)
  if (MATLAB_MEX_EXECUTABLE)
    set (BASIS_MEX_TIMEOUT "600" CACHE STRING "Timeout for MEX script execution")
    set (BASIS_MEX_FLAGS   ""    CACHE STRING "Common MEX switches (separated by ' '; use '\\' to mask ' ').")
    mark_as_advanced (BASIS_MEX_FLAGS)
    mark_as_advanced (BASIS_MEX_TIMEOUT)
  endif ()
endmacro ()

# ----------------------------------------------------------------------------
## @brief Add global MATLAB Compiler (mcc) options to CMake cache
#
# @par MATLAB Compiler options:
# <table border="0">
#   <tr>
#     @tp @b BASIS_MCC_MATLAB_MODE @endtp
#     <td>Enable/Disable invocation of MATLAB Compiler in MATLAB mode.</td>
#   </tr>
#   <tr>
#     @tp @b BASIS_MCC_FLAGS @endtp
#     <td>Compile flags used to build MATLAB Compiler targets.</td>
#   </tr>
#   <tr>
#     @tp @b BASIS_MCC_TIMEOUT @endtp
#     <td>Timeout for building MATLAB Compiler targets.</td>
#   </tr>
#   <tr>
#     @tp @b BASIS_MCC_RETRY_ATTEMPTS @endtp
#     <td>Maximum number of retries on MATLAB Compiler license checkout.</td>
#   </tr>
#   <tr>
#     @tp @b BASIS_MCC_RETRY_DELAY @endtp
#     <td>Delay between retries to build MATLAB Compiler compiled targets on license checkout errors.</td>
#   </tr>
# </table>
macro(basis_add_mcc_options)
  option (
    BASIS_MCC_MATLAB_MODE
    "Prefer MATLAB mode over standalone mode to invoke MATLAB Compiler to release MCC licence ASAP."
    OFF # using MATLAB mode is preferred when the license is shared
        # among users as it releases the license immediately once done
  )
  set (
    BASIS_MCC_FLAGS
      "-R -singleCompThread"
    CACHE STRING
      "Common MATLAB Compiler flags (separated by ' '; use '\\' to mask ' ')."
  )
  set (BASIS_MCC_TIMEOUT        "1800" CACHE STRING "Timeout for MATLAB Compiler execution")
  set (BASIS_MCC_RETRY_ATTEMPTS "4"    CACHE STRING "Maximum number of retries on MATLAB Compiler license checkout error.")
  set (BASIS_MCC_RETRY_DELAY    "30"   CACHE STRING "Delay between retries to build MATLAB Compiler compiled targets on license checkout error.")
  mark_as_advanced (BASIS_MCC_MATLAB_MODE)
  mark_as_advanced (BASIS_MCC_FLAGS)
  mark_as_advanced (BASIS_MCC_TIMEOUT)
  mark_as_advanced (BASIS_MCC_RETRY_ATTEMPTS)
  mark_as_advanced (BASIS_MCC_RETRY_DELAY)
endmacro ()

## @}
# end of Doxygen group

# ============================================================================
# utilities
# ============================================================================

# ----------------------------------------------------------------------------
## @brief Determine version of MATLAB installation.
#
# @param [out] VERSION Value returned by the "version" command of MATLAB or
#                      an empty string if execution of MATLAB failed.
#
# @returns Sets the variable named by @p VERSION to the full MATLAB version.
#
# @ingroup CMakeUtilities
function (basis_get_full_matlab_version VERSION)
  if (NOT MATLAB_EXECUTABLE)
    set (${VERSION} "" PARENT_SCOPE)
    return ()
  endif ()
  set (WORKING_DIR "${CMAKE_BINARY_DIR}/CMakeFiles")
  set (OUTPUT_FILE "${WORKING_DIR}/MatlabVersion.txt")
  # read MATLAB version from existing output file
  set (_MATLAB_VERSION)
  if (EXISTS "${OUTPUT_FILE}")
    file (READ "${OUTPUT_FILE}" LINES)
    string (REGEX REPLACE "\n"    ";" LINES "${LINES}")
    string (REGEX REPLACE "^;|;$" ""  LINES "${LINES}")
    list (LENGTH LINES NLINES)
    if (NLINES EQUAL 2)
      list (GET LINES 0 _MATLAB_EXECUTABLE)
      list (GET LINES 1 _MATLAB_VERSION)
      basis_sanitize_for_regex (RE "${MATLAB_EXECUTABLE}")
      if (NOT _MATLAB_EXECUTABLE MATCHES "^${RE}$")
        set (_MATLAB_VERSION)
      endif ()
      unset (RE)
    endif ()
  endif ()
  # run matlab command to write return value of "version" command to text file
  if (NOT _MATLAB_VERSION)
    message (STATUS "Determining MATLAB version...")
    set (CMD "${MATLAB_EXECUTABLE}" -nodesktop -nosplash -singleCompThread)
    if (WIN32)
      list (APPEND CMD -automation)
    else ()
      list (APPEND CMD -nojvm)
    endif ()
    file (WRITE "${WORKING_DIR}/basis_get_full_matlab_version.m" 
"% DO NOT EDIT. Automatically created by BASIS (basis_get_full_matlab_version).
fid = fopen ('${OUTPUT_FILE}', 'w')
if fid == -1, fprintf(2, '??? Error: Failed to open file ${OUTPUT_FILE} for writing!'), quit force, end
fprintf (fid, '${MATLAB_EXECUTABLE}\\n%s\\n', version)
fclose (fid)
quit force
"
    )
    execute_process (
      COMMAND           ${CMD} -r "cd('${WORKING_DIR}');basis_get_full_matlab_version;"
      WORKING_DIRECTORY "${WORKING_DIR}"
      RESULT_VARIABLE   RETVAL
      TIMEOUT           600 # MATLAB startup can be *very* slow the first time
      ERROR_VARIABLE    STDERR
      OUTPUT_QUIET
    )
    if (NOT RETVAL EQUAL 0 OR STDERR MATCHES "\\?\\?\\? Error")
      set (${VERSION} "" PARENT_SCOPE)
      if (RETVAL MATCHES "timeout")
        set (REASON ": ${RETVAL}")
      elseif (NOT RETVAL EQUAL 0)
        set (REASON " with exit code: ${RETVAL}")
      else ()
        set (REASON ": Failed to open file ${OUTPUT_FILE} for writing")
      endif ()
      message (STATUS "Determining MATLAB version... - failed${REASON}")
      return ()
    endif ()
    # wait until MATLAB process terminated entirely and wrote the (buffered?) file
    set (nsleep_count 0)
    set (nsleep_max  30)
    while (NOT EXISTS "${OUTPUT_FILE}")
      math (EXPR nsleep_count "${nsleep_count}+1")
      if (nsleep_count GREATER nsleep_max)
        message (STATUS "Determining MATLAB version... - failed: File ${OUTPUT_FILE} still non-existent after ${nsleep_max}s of successful MATLAB termination")
        return ()
      endif ()
      if (WIN32)
        execute_process (COMMAND ping 1.1.1.1 -n 1 -w 1000 OUTPUT_QUIET)
      else ()
        execute_process (COMMAND sleep 1 OUTPUT_QUIET)
      endif ()
    endwhile ()
    # read MATLAB version from text file
    file (READ "${OUTPUT_FILE}" LINES)
    string (REGEX REPLACE "\n"    ";" LINES "${LINES}")
    string (REGEX REPLACE "^;|;$" ""  LINES "${LINES}")
    list (LENGTH LINES NLINES)
    if (NLINES EQUAL 2)
      list (GET LINES 1 _MATLAB_VERSION)
    else ()
      set (${VERSION} "" PARENT_SCOPE)
      message (STATUS "Determining MATLAB version... - failed")
      return ()
    endif ()
    if (BASIS_VERBOSE)
      message (STATUS "Determining MATLAB version... - done: ${_MATLAB_VERSION}")
    else ()
      message (STATUS "Determining MATLAB version... - done")
    endif ()
  endif ()
  # return
  set (${VERSION} "${_MATLAB_VERSION}" PARENT_SCOPE)
endfunction ()

# ----------------------------------------------------------------------------
## @brief Get version of MATLAB installation.
#
# @param [out] ARGV1 If given, the named variable is set to the version string
#                    ("<major>.<minor>.<patch>") of the MATLAB installation.
#                    Otherwise, the variables @c MATLAB_VERSION_STRING,
#                    @c MATLAB_VERSION_MAJOR, @c MATLAB_VERSION_MINOR,
#                    @c MATLAB_VERSION_PATCH, and @c MATLAB_RELEASE are set
#                    in the scope of the caller.
#
# @ingroup CMakeUtilities
function (basis_get_matlab_version)
  if (ARGC GREATER 1)
    message (FATAL_ERROR "basis_get_matlab_version(): Too many arguments!")
  endif ()
  basis_get_full_matlab_version (VERSION)
  if (VERSION MATCHES "^([0-9]+)\\.([0-9]+)\\.([0-9]+)")
    set (VERSION_STRING "${CMAKE_MATCH_0}")
    set (VERSION_MAJOR  "${CMAKE_MATCH_1}")
    set (VERSION_MINOR  "${CMAKE_MATCH_2}")
    set (VERSION_PATCH  "${CMAKE_MATCH_3}")
  else ()
    set (VERSION_STRING "0.0")
    set (VERSION_MAJOR  "0")
    set (VERSION_MINOR  "0")
    set (VERSION_PATCH  "0")
  endif ()
  if (ARGC EQUAL 1)
    set (${ARGV0} "${VERSION_STRING}" PARENT_SCOPE)
  else ()
    set (MATLAB_VERSION_STRING "${VERSION_STRING}" PARENT_SCOPE)
    set (MATLAB_VERSION_MAJOR  "${VERSION_MAJOR}"  PARENT_SCOPE)
    set (MATLAB_VERSION_MINOR  "${VERSION_MINOR}"  PARENT_SCOPE)
    set (MATLAB_VERSION_PATCH  "${VERSION_PATCH}"  PARENT_SCOPE)
    if (VERSION MATCHES ".*\\\((.+)\\\)")
      set (MATLAB_RELEASE "${CMAKE_MATCH_1}" PARENT_SCOPE)
    else ()
      set (MATLAB_RELEASE "" PARENT_SCOPE)
    endif ()
  endif ()
endfunction ()

# ----------------------------------------------------------------------------
## @brief Get release version of MATLAB installation.
#
# @param [out] ARGV1 If given, the named variable is set to the release string
#                    of the MATLAB installation, e.g., "R2009b". Otherwise,
#                    the variable @c MATLAB_RELEASE is set in the scope of the
#                    caller.
#
# @ingroup CMakeUtilities
function (basis_get_matlab_release)
  if (ARGC GREATER 1)
    message (FATAL_ERROR "basis_get_matlab_release(): Too many arguments!")
  endif ()
  basis_get_full_matlab_version (VERSION)
  if (VERSION MATCHES ".*\\\((.+)\\\)")
    set (RELEASE "${CMAKE_MATCH_1}")
  else ()
    set (RELEASE "")
  endif ()
  if (ARGC EQUAL 1)
    set (${ARGV0} "${RELEASE}" PARENT_SCOPE)
  else ()
    set (MATLAB_RELEASE "${RELEASE}")
  endif ()
endfunction ()

# ----------------------------------------------------------------------------
## @brief Determine extension of MEX-files for this architecture.
#
# @param [out] ARGN The first argument ARGV0 is set to the extension of
#                   MEX-files (excluding '.'). If the CMake variable MEX_EXT
#                   is set, its value is returned. Otherwise, this function
#                   tries to determine it from the system information.
#                   If the extension could not be determined, an empty string
#                   is returned. If no argument is given, the extension is
#                   cached as the variable MEX_EXT.
#
# @returns Sets the variable named by the first argument to the
#          platform-specific extension of MEX-files.
#
# @ingroup CMakeUtilities
function (basis_mexext)
  # default return value
  set (MEXEXT "${MEX_EXT}")
  # use MEXEXT if possible
  if (NOT MEXEXT AND MATLAB_MEXEXT_EXECUTABLE)
    execute_process (
      COMMAND         "${MATLAB_MEXEXT_EXECUTABLE}"
      RESULT_VARIABLE RETVAL
      OUTPUT_VARIABLE MEXEXT
      ERROR_QUIET
      OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if (RETVAL)
      set (MEXEXT "")
    endif ()
  endif ()
  # otherwise, determine extension given CMake variables describing the system
  if (NOT MEXEXT)
    if (CMAKE_SYSTEM_NAME MATCHES "Linux")
      if (CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64" OR
          CMAKE_SIZEOF_VOID_P MATCHES 8)
        set (MEXEXT "mexa64")
      elseif (CMAKE_SYSTEM_PROCESSOR MATCHES "x86" OR
              CMAKE_SYSTEM_PROCESSOR MATCHES "i686")
        set (MEXEXT "mexglx")
      endif ()
    elseif (CMAKE_SYSTEM_NAME MATCHES "Windows")
      if (CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64" OR
          CMAKE_SIZEOF_VOID_P    MATCHES 8)
        set (MEXEXT "mexw64")
      elseif (CMAKE_SYSTEM_PROCESSOR MATCHES "x86" OR
              CMAKE_SYSTEM_PROCESSOR MATCHES "i686")
        set (MEXEXT "mexw32")
      endif ()
    elseif (CMAKE_SYSTEM_NAME MATCHES "Darwin")
      if (CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64" OR
          CMAKE_SIZEOF_VOID_P    MATCHES 8)
        set (MEXEXT "mexmaci64")
      else ()
        set (MEXEXT "mexmaci")
      endif ()
    elseif (CMAKE_SYSTEM_NAME MATCHES "SunOS")
      set (MEXEXT "mexs64")
    endif ()
  endif ()
  # return value
  if (ARGC GREATER 0)
    set ("${ARGV0}" "${MEXEXT}" PARENT_SCOPE)
  else ()
    if (NOT DEFINED MEX_EXT)
      set (MARKIT 1)
    else ()
      set (MARKIT 0)
    endif ()
    set (MEX_EXT "${MEXEXT}" CACHE STRING "The extension of MEX-files for this architecture." FORCE)
    if (MARKIT)
      mark_as_advanced (MEX_EXT)
    endif ()
  endif ()
endfunction ()

# ----------------------------------------------------------------------------
## @brief This function writes a MATLAB M-file with addpath() statements.
#
# This function writes an MATLAB M-file into the top directory of the build
# tree which contains an addpath() statement for each directory that was added
# via basis_include_directories().
#
# @returns Creates file add_\<project\>_paths.m in the current binary directory.
#
# @ingroup CMakeUtilities
function (basis_create_addpaths_mfile)
  basis_get_project_property (INCLUDE_DIRS PROPERTY PROJECT_INCLUDE_DIRS)
  basis_write_addpaths_mfile ("${CMAKE_CURRENT_BINARY_DIR}/add_${PROJECT_NAME_L}_paths.m" ${INCLUDE_DIRS})
endfunction ()

# ----------------------------------------------------------------------------
## @brief This function writes a MATLAB M-file with addpath() statements.
#
# @param [in] MFILE Name of M-file.
# @param [in] ARGN  The remaining arguments are the paths which should be added
#                   to the search path of MATLAB when this M-file is being
#                   executed. If the option APPEND is given, the paths are
#                   appended to the specified M-file. Otherwise, any existing
#                   file will be overwritten. The given directory paths can
#                   be relative, in which case they are interpreted relative
#                   to the location of the written M-file using the mfilename()
#                   function of MATLAB.
#
# @ingroup CMakeUtilities
macro (basis_write_addpaths_mfile MFILE)
  CMAKE_PARSE_ARGUMENTS (ARGN "APPEND" "" "" ${ARGN})
  if (NOT ARGN_APPEND)
    file (WRITE "${MFILE}" "% DO NOT edit. This file is automatically generated by BASIS.
[mfiledir, ~, ~, ~] = fileparts(mfilename('fullpath'));\n")
  endif ()
  foreach (P IN LISTS ARGN_UNPARSED_ARGUMENTS)
    if (P MATCHES "^\\.?$")
      file (APPEND "${MFILE}" "addpath(mfiledir);\n")
    elseif (IS_ABSOLUTE "${P}")
      file (APPEND "${MFILE}" "addpath('${P}');\n")
    else ()
      file (APPEND "${MFILE}" "addpath([mfiledir '/${P}']);\n")
    endif ()
  endforeach ()
endmacro ()

# ----------------------------------------------------------------------------
## @brief Generate MATLAB wrapper executable.
#
# This function writes a Bash script on Unix or a Windows Command script on
# Windows platforms which execute the specified MATLAB command using the -r
# option of the matlab executable and the -nodesktop and -nosplash options.
# It is used by the build scripts generated by the basis_build_mcc_target()
# in order to build an executable from MATLAB source files without the use
# of the MATLAB Compiler. In this case, the MATLAB source files are simply
# copied to the installation directory and the wrapper script written by
# this function used to execute the main function with the command-line
# arguments passed on to this executable.
#
# @param [in] OUTPUT_FILE    Name of the output executable file.
# @param [in] ARGN           The remaining options
# @par
# <table border=0>
#   <tr>
#     @tp @b DESTINATION dir @endtp
#     <td>Installation destination. (default: directory of @c OUTPUT_FILE)</td>
#   </tr>
#   <tr>
#     @tp @b COMMAND name @endtp
#     <td>Name of the MATLAB command to execute, i.e.,
#         the name of the main function.</td>
#   </tr>
#   <tr>
#     @tp @b STARTUP mfile @endtp
#     <td>Absolute path of a startup M-file.</td>
#   </tr>
#   <tr>
#     @tp @b MATLABPATH dir1[ dir2...]
#     <td>List of directories to be added to the MATLAB search path.</td>
#   </tr>
#   <tr>
#     @tp @b OPTIONS opt1[ opt2...]
#     <td>Additional options to pass on to the <tt>matlab</tt> executable.</td>
#   </tr>
# </table>
function (basis_generate_matlab_executable OUTPUT_FILE)
  CMAKE_PARSE_ARGUMENTS (ARGN "" "COMMAND;STARTUP;DESTINATION" "OPTIONS;MATLABPATH" ${ARGN})
  if (NOT MATLAB_EXECUTABLE)
    set (MATLAB_EXECUTABLE matlab)
  endif ()
  if (NOT OUTPUT_FILE)
    message ("basis_generate_matlab_executable(): Missing OUTPUT_FILE argument!")
  endif ()
  if (NOT ARGN_DESTINATION)
    get_filename_component (ARGN_DESTINATION "${OUTPUT_FILE}" PATH)
  endif ()
  # source path
  set (MATLABPATH)
  foreach (P IN LISTS ARGN_MATLABPATH)
    if (P MATCHES "^\\.?$")
      set (P "$__DIR__")
    elseif (NOT IS_ABSOLUTE "${P}")
      set (P "$__DIR__/${P}")
    endif ()
    list (APPEND MATLABPATH "${P}")
  endforeach ()
  if (MATLABPATH)
    list (REMOVE_DUPLICATES MATLABPATH)
  endif ()
  # startup script
  if (ARGN_STARTUP)
    get_filename_component (STARTUP_COMMAND "${ARGN_STARTUP}" NAME_WE)
    get_filename_component (STARTUP_DIR     "${ARGN_STARTUP}" PATH)
    get_filename_component (STARTUP_PKG     "${STARTUP_DIR}"  NAME)
    if (STARTUP_PKG MATCHES "^\\+")
      get_filename_component (STARTUP_DIR "${STARTUP_DIR}" PATH)
      string (REGEX REPLACE "^\\+" "" STARTUP_PKG "${STARTUP_PKG}")
      set (STARTUP_COMMAND "${STARTUP_PKG}.${STARTUP_COMMAND}")
    endif ()
    if (IS_ABSOLUTE "${STARTUP_DIR}")
      file (RELATIVE_PATH STARTUP_DIR "${ARGN_DESTINATION}" "${STARTUP_DIR}")
    endif ()
    if (STARTUP_DIR)
      set (STARTUP_DIR "$__DIR__/${STARTUP_DIR}")
    else ()
      set (STARTUP_DIR "$__DIR__")
    endif ()
    list (FIND MATLABPATH "${STARTUP_DIR}" IDX)
    if (IDX EQUAL -1)
      set (STARTUP_CODE ", addpath('${STARTUP_DIR}')")
    else ()
      set (STARTUP_CODE)
    endif ()
    set (STARTUP_CODE "${STARTUP_CODE}, ${STARTUP_COMMAND}")
  else ()
    set (STARTUP_CODE)
  endif ()
  # write wrapper executable
  if (MATLABPATH)
    basis_list_to_delimited_string (MATLABPATH "', '" NOAUTOQUOTE ${MATLABPATH})
    set (MATLABPATH ", addpath('${MATLABPATH}', '-begin')")
  else ()
    set (MATLABPATH)
  endif ()
  file (WRITE "${OUTPUT_FILE}"
    # note that Bash variables within the script are denoted by $var
    # instead of ${var} to prevent CMake from substituting these patterns
    "#! /bin/bash

readonly __DIR__=\"${BASIS_BASH___DIR__}\"

errlog=
finish()
{
    local status=0
    if [[ -n \"$errlog\" ]]; then
        grep '??? Error' \"$errlog\" &> /dev/null
        [[ $? -ne 0 ]] || status=1
        /bin/rm \"$errlog\"
    fi
    exit $status
}

if [[ -d \"$TMPDIR\" ]]; then
    tmpdir=$TMPDIR
else
    tmpdir=/tmp
fi

errlog=`mktemp \"$tmpdir/${ARGN_COMMAND}-log.XXXXXX\"`
[[ $? -eq 0 ]] || {
    echo \"Failed to create temporary log file in '$tmpdir'!\" 1>&2
    exit 1
}

args=
while [[ $# -gt 0 ]]; do
  [[ -z \"$args\" ]] || args=\"$args, \"
  args=\"$args'$1'\"
  shift
done

echo 'Launching MATLAB to execute ${ARGN_COMMAND} function...'
trap finish EXIT # DO NOT install trap earlier !
'${MATLAB_EXECUTABLE}' -nodesktop -nosplash ${ARGN_OPTIONS} \\
    -r \"try${MATLABPATH}${STARTUP_CODE}, ${ARGN_COMMAND}($args), catch err, fprintf(2, ['??? Error executing ${ARGN_COMMAND}\\n' err.message '\\n']), end, quit force\" \\
    2> >(tee \"$errlog\" >&2)"
  ) # end of file(WRITE) command
  if (UNIX)
    execute_process (COMMAND /bin/chmod +x "${OUTPUT_FILE}")
  endif ()
endfunction ()

# ============================================================================
# MEX-file target
# ============================================================================

# ----------------------------------------------------------------------------
## @brief Add MEX-file target.
#
# @note This function should not be used directly. Instead, it is called
#       by basis_add_library() if the (detected) programming language
#       of the given source code files is @c CXX (i.e., C/C++) and the @c MEX
#       type option is given.
#
# This function is used to add a shared library target which is built
# using the MATLAB MEX script (mex).
#
# By default, the BASIS C++ utilities library is added as link dependency.
# If none of the BASIS C++ utilities are used by this target, the option
# NO_BASIS_UTILITIES can be given. To enable this option by default, set the
# variable @c BASIS_UTILITIES to @c FALSE, best in the <tt>Settings.cmake</tt>
# file located in the @c PROJECT_CONFIG_DIR (add such file if missing).
# If the use of the BASIS C++ utilities is disabled by default, the
# @c USE_BASIS_UTILITIES option can be used to enable them for this target
# only. Note that the utilities library is a static library and thus the linker
# would simply not include any of the BASIS utility functions in the final
# binary file if not used. The only advantage of setting @c BASIS_UTILITIES to
# @c FALSE or to always specify @c NO_BASIS_UTILITIES if no target uses the
# utilities is that the BASIS utilities library will not be build in this case.
#
# A custom CMake build target with the following properties is added by this
# function to the build system. These properties are used by
# basis_build_mex_target() to generate a build script written in CMake
# code which is executed by a custom CMake command. Before the invokation of
# basis_build_mex_target(), the target properties can be modified using
# basis_set_target_properties().
#
# @note Custom BASIS build targets are finalized by BASIS using basis_project_end(),
#       i.e., the end of the root CMake configuration file of the (sub-)project.
#
# @par Properties on script library targets
# <table border=0>
#   <tr>
#     @tp @b MFILE file @endtp
#     <td>MATLAB source file with function prototype and documentation of MEX-file.
#         (default: none)</td>
#   </tr>
#   <tr>
#     @tp @b PREFIX prefix @endtp
#     <td>Output prefix of build MEX-file such as package name
#         (the prefix must include the leading + and trailing /).</td>
#   </tr>
# </table>
#
# @attention Properties documented as read-only must not be modified.
#
# An install command for the added library target is added by this function
# as well. The MEX-file will be installed as part of the specified @p COMPONENT
# in the @c INSTALL_LIBRARY_DIR on Unix and @c INSTALL_RUNTIME_DIR on Windows.
#
# @param [in] TARGET_NAME Name of build target.
# @param [in] ARGN        The remaining arguments are parsed and the following
#                         arguments extracted. All unparsed arguments are treated
#                         as the source files of the MEX-file.
# @par
# <table border="0">
#   <tr>
#     @tp @b COMPONENT name @endtp
#     <td>Name of installation component as part of which this MEX-file is being
#         installed if the @c LIBRARY_INSTALL_DIRECTORY property is not "none".
#         (default: @c BASIS_LIBRARY_COMPONENT)</td>
#   </tr>
#   <tr>
#     @tp @b [NO]EXPORT @endtp
#     <td>Whether to export this target. (default: @c TRUE)</td>
#   </tr>
#   <tr>
#     @tp @b NO_BASIS_UTILITIES @endtp
#     <td>Specify that the BASIS utilities are not used by this MEX-file and
#         hence no link dependency on the BASIS utilities shall be added.
#         (default: @c NOT BASIS_UTILITIES)</td>
#   </tr>
#   <tr>
#     @tp @b USE_BASIS_UTILITIES @endtp
#     <td>Specify that the BASIS utilities are used and required by this MEX-file
#         and hence a link dependency on the BASIS utilities must be added.
#         (default: @c BASIS_UTILITIES)</td>
#   </tr>
#   <tr>
#     @tp @b FINAL @endtp
#     <td>Finalize custom targets immediately. Any following target property changes
#         will have no effect. When this option is used, the custom target which
#         executes the custom build command is added in the current working directory.
#         Otherwise it will be added in the top-level source directory of the project.
#         Which with the Visual Studio generators adds the corresponding Visual Studio
#         Project files directly to the top-level build directory. This can be avoided
#         using this option or calling basis_finalize_targets() at the end of each
#         CMakeLists.txt file.</td>
#   </tr>
# </table>
#
# @returns Adds custom target to build MEX-file using the MEX script.
#
# @sa basis_add_library()
#
# @ingroup CMakeUtilities
function (basis_add_mex_file TARGET_NAME)
  # check target name
  basis_check_target_name ("${TARGET_NAME}")
  basis_make_target_uid (TARGET_UID "${TARGET_NAME}")
  message (STATUS "Adding MEX-file ${TARGET_UID}...")
  # required commands available ?
  if (NOT MATLAB_MEX_EXECUTABLE)
    message (FATAL_ERROR "MATLAB MEX script (mex) not found! It is required to build target ${TARGET_UID}."
                         " Forgot to add MATLAB{mex} as dependency? Otherwise, set MATLAB_MEX_EXECUTABLE manually and try again.")
  endif ()
  basis_add_mex_options()
  # parse arguments
  CMAKE_PARSE_ARGUMENTS (
    ARGN
      "USE_BASIS_UTILITIES;NO_BASIS_UTILITIES;EXPORT;NOEXPORT;FINAL"
      "COMPONENT;DESTINATION"
      ""
    ${ARGN}
  )
  set (SOURCES ${ARGN_UNPARSED_ARGUMENTS})
  basis_set_flag (ARGN EXPORT ${BASIS_EXPORT_DEFAULT})
  if (ARGN_USE_BASIS_UTILITIES AND ARGN_NO_BASIS_UTILITIES)
    message (FATAL_ERROR "Target ${TARGET_UID}: Options USE_BASIS_UTILITIES and NO_BASIS_UTILITIES are mutually exclusive!")
  endif ()
  if (ARGN_USE_BASIS_UTILITIES)
    set (USES_BASIS_UTILITIES TRUE)
  elseif (ARGN_NO_BASIS_UTILITIES)
    set (USES_BASIS_UTILITIES FALSE)
  else ()
    set (USES_BASIS_UTILITIES ${BASIS_UTILITIES})
  endif ()
  basis_mexext (MEXEXT)
  # IS_TEST flag
  basis_sanitize_for_regex (RE "${PROJECT_TESTING_DIR}")
  if (CMAKE_CURRENT_SOURCE_DIR MATCHES "^${RE}")
    set (IS_TEST TRUE)
  else ()
    set (IS_TEST FALSE)
  endif ()
  # installation component
  if (NOT ARGN_COMPONENT)
    set (ARGN_COMPONENT "${BASIS_LIBRARY_COMPONENT}")
  endif ()
  if (NOT ARGN_COMPONENT)
    set (ARGN_COMPONENT "Unspecified")
  endif ()
  # installation directory
  if (ARGN_DESTINATION)
    if (ARGN_DESTINATION MATCHES "^[nN][oO][nN][eE]$")
      set (ARGN_DESTINATION)
    elseif (IS_ABSOLUTE "${ARGN_DESTINATION}")
      file (RELATIVE_PATH ARGN_DESTINATION "${CMAKE_INSTALL_PREFIX}" "${ARGN_DESTINATION}")
    endif ()
  else ()
    set (ARGN_DESTINATION "${INSTALL_MATLAB_LIBRARY_DIR}")
  endif ()
  # configure (.in) source files
  basis_configure_sources (SOURCES ${SOURCES})
  # add custom target
  add_custom_target (${TARGET_UID} ALL SOURCES ${SOURCES})
  get_directory_property (INCLUDE_DIRS INCLUDE_DIRECTORIES)
  get_directory_property (LINK_DIRS    LINK_DIRECTORIES)
  if (MATLAB_LIBRARY_DIR)
    list (INSERT LINK_DIRS 0 "${MATLAB_LIBRARY_DIR}")
  endif ()
  set_target_properties (
    ${TARGET_UID}
    PROPERTIES
      LANGUAGE                  "CXX"
      BASIS_TYPE                MEX
      BASIS_UTILITIES           ${USES_BASIS_UTILITIES}
      BASIS_INCLUDE_DIRECTORIES "${INCLUDE_DIRS}"
      BASIS_LINK_DIRECTORIES    "${LINK_DIRS}"
      BUILD_DIRECTORY           "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${TARGET_UID}"
      SOURCE_DIRECTORY          "${CMAKE_CURRENT_SOURCE_DIR}"
      BINARY_DIRECTORY          "${CMAKE_CURRENT_BINARY_DIR}"
      LIBRARY_OUTPUT_DIRECTORY  "${BINARY_MATLAB_LIBRARY_DIR}"
      LIBRARY_INSTALL_DIRECTORY "${ARGN_DESTINATION}"
      LIBRARY_COMPONENT         "${ARGN_COMPONENT}"
      COMPILE_FLAGS             "${BASIS_MEX_FLAGS}"
      LINK_FLAGS                ""
      LINK_DEPENDS              "${LINK_DEPENDS}"
      PREFIX                    ""
      OUTPUT_NAME               ""
      SUFFIX                    ".${MEXEXT}"
      MFILE                     ""
      TEST                      ${IS_TEST}
      EXPORT                    ${EXPORT}
  )
  # link to BASIS utilities
  if (USES_BASIS_UTILITIES)
    basis_target_link_libraries (.${TARGET_UID} basis)
  endif ()
  # finalize target
  if (ARGN_FINAL)
    basis_finalize_targets (${TARGET_UID})
  endif ()
  # add target to list of targets
  basis_set_project_property (APPEND PROPERTY TARGETS "${TARGET_UID}")
  message (STATUS "Adding MEX-file ${TARGET_UID}... - done")
endfunction ()

# ============================================================================
# MATLAB Compiler target
# ============================================================================

# ----------------------------------------------------------------------------
## @brief Add MATLAB Compiler target.
#
# @note This function should not be used directly. Instead, it is called
#       by either basis_add_executable() or basis_add_library() if the
#       (detected) programming language of the given source code files is
#       @c MATLAB.
#
# This function is used to add an executable or shared library target which is
# built using the MATLAB Compiler (MCC).
#
# A custom CMake build target with the following properties is added by this
# function to the build system. These properties are used by
# basis_build_mcc_target() to generate a build script written in CMake
# code which is executed by a custom CMake command. Before the invokation of
# basis_build_mcc_target(), the target properties can be modified using
# basis_set_target_properties().
#
# @note Custom BASIS build targets are finalized by BASIS using basis_project_end(),
#       i.e., the end of the root CMake configuration file of the (sub-)project.
#
# @par Properties on MATLAB Compiler targets
# <table border=0>
#   <tr><td>TODO</td></tr>
# </table>
#
# An install command for the added executable or library target is added by
# this function as well. The executable will be installed as part of the
# @p RUNTIME_COMPONENT in the directory @c INSTALL_RUNTIME_DIR. The runtime
# library will be installed as part of the @p RUNTIME_COMPONENT in the directory
# @c INSTALL_LIBRARY_DIR on Unix and @c INSTALL_RUNTIME_DIR on Windows.
# Static/import libraries will be installed as part of the @p LIBRARY_COMPONENT
# in the directory @c INSTALL_ARCHIVE_DIR.
#
# @note If this function is used within the @c PROJECT_TESTING_DIR, the built
#       executable is output to the @c BINARY_TESTING_DIR directory tree instead.
#       Moreover, no installation rules are added. Test executables are further
#       not exported, regardless of the value of the @c EXPORT property.
#
# @param [in] TARGET_NAME Name of build target.
# @param [in] ARGN        The remaining arguments are parsed and the following
#                         arguments extracted. All unparsed arguments are treated
#                         as the MATLAB or C/C++ source files, respectively.
# @par
# <table border="0">
#   <tr>
#     @tp <b>EXECUTABLE</b>|<b>LIBEXEC</b>|<b>SHARED</b> @endtp
#     <td>Type of the MATLAB Compiler target which can be either a stand-alone
#         executable, an auxiliary executable, or a shared library.
#         (default: @c EXECUTABLE)</td>
#   </tr>
#   <tr>
#     @tp @b COMPONENT name @endtp
#     <td>Name of component as part of which this executable or library will be
#         installed if the @c RUNTIME_INSTALL_DIRECTORY or @c LIBRARY_INSTALL_DIRECTORY
#         property is not "none". Used only if @p RUNTIME_COMPONENT or
#         @p LIBRARY_COMPONENT not specified.
#         (default: see @p RUNTIME_COMPONENT and @p LIBRARY_COMPONENT arguments)</td>
#   </tr>
#   <tr>
#     @tp @b DESTINATION dir @endtp
#     <td>Installation directory for executable or runtime and library component
#         of shared library relative to @c CMAKE_INSTALL_PREFIX. Used only if
#         @p RUNTIME_DESTINATION or @p LIBRARY_DESTINATION not specified.
#         If "none" (case-insensitive) is given as argument, no default installation
#         rules are added. (default: see @p RUNTIME_DESTINATION and
#         @p LIBRARY_DESTINATION arguments)</td>
#   </tr>
#   <tr>
#     @tp @b LIBRARY_COMPONENT name @endtp
#     <td>Name of component as part of which import/static library will be intalled
#         if a shared library is build and the @c LIBRARY_INSTALL_DIRECTORY property is
#         not "none". (default: @c COMPONENT if specified or @c BASIS_LIBRARY_COMPONENT
#         otherwise)</td>
#   </tr>
#   <tr>
#     @tp @b LIBRARY_DESTINATION dir @endtp
#     <td>Installation directory of the library component relative to
#         @c CMAKE_INSTALL_PREFIX. If "none" (case-insensitive) is given as argument or
#         an executable is build, no installation rule for the library component is added.
#         (default: @c INSTALL_ARCHIVE_DIR)</td>
#   </tr>
#   <tr>
#     @tp @b HEADER_DESTINATION dir @endtp
#     <td>Installation directory of the library header file relative to
#         @c INSTALL_INCLUDE_DIR. If "none" (case-insensitive) is given as argument or
#         an executable is build, no installation rule for the library header file is added.
#         (default: @c INSTALL_INCLUDE_DIR)</td>
#   </tr>
#   <tr>
#     @tp @b RUNTIME_COMPONENT name @endtp
#     <td>Name of component as part of which executable or runtime library, respectively,
#         will be installed if the @c RUNTIME_INSTALL_DIRECTORY property is not "none".
#         (default: @c COMPONENT if specified or @c BASIS_RUNTIME_COMPONENT otherwise)</td>
#   </tr>
#   <tr>
#     @tp @b RUNTIME_DESTINATION dir @endtp
#     <td>Installation directory of the executable or runtime component of the shared library
#         relative to @c CMAKE_INSTALL_PREFIX. If "none" (case-insensitive) is given as argument,
#         no installation rule for the runtime library is added.
#         (default: @c INSTALL_LIBRARY_DIR for shared libraries on Unix or
#         @c INSTALL_RUNTIME_DIR otherwise)</td>
#   </tr>
#   <tr>
#     @tp @b [NO]EXPORT @endtp
#     <td>Whether to export this target. (default: @c TRUE)</td>
#   </tr>
#   <tr>
#     @tp @b NO_BASIS_UTILITIES @endtp
#     <td>Specify that the BASIS utilities are not used by this executable or shared library
#         and hence no link dependency on the BASIS utilities shall be added.
#         (default: @c NOT BASIS_UTILITIES)</td>
#   </tr>
#   <tr>
#     @tp @b USE_BASIS_UTILITIES @endtp
#     <td>Specify that the BASIS utilities are used and required by this executable
#         or shared library, respectively, and hence a link dependency on the BASIS utilities
#         must be added.
#         (default: @c BASIS_UTILITIES)</td>
#   </tr>
#   <tr>
#     @tp @b FINAL @endtp
#     <td>Finalize custom targets immediately. Any following target property changes
#         will have no effect. When this option is used, the custom target which
#         executes the custom build command is added in the current working directory.
#         Otherwise it will be added in the top-level source directory of the project.
#         Which with the Visual Studio generators adds the corresponding Visual Studio
#         Project files directly to the top-level build directory. This can be avoided
#         using this option or calling basis_finalize_targets() at the end of each
#         CMakeLists.txt file.</td>
#   </tr>
# </table>
#
# @todo Consider NO_BASIS_UTILITIES and USE_BASIS_UTILITIES options after the BASIS
#       utilities for MATLAB have been implemented.
#
# @returns Adds custom target which builds depending on the @p BASIS_TYPE property
#          either an executable or a shared library using the MATLAB Compiler.
#
# @sa basis_add_executable()
# @sa basis_add_library()
#
# @ingroup CMakeUtilities
function (basis_add_mcc_target TARGET_NAME)
  # check target name
  basis_check_target_name ("${TARGET_NAME}")
  basis_make_target_uid (TARGET_UID "${TARGET_NAME}")
  # parse arguments
  CMAKE_PARSE_ARGUMENTS (
    ARGN
      "SHARED;EXECUTABLE;LIBEXEC;USE_BASIS_UTILITIES;NO_BASIS_UTILITIES;EXPORT;NOEXPORT;FINAL"
      "COMPONENT;RUNTIME_COMPONENT;LIBRARY_COMPONENT;DESTINATION;RUNTIME_DESTINATION;LIBRARY_DESTINATION;HEADER_DESTINATION"
      ""
    ${ARGN}
  )
  set (SOURCES "${ARGN_UNPARSED_ARGUMENTS}")
  basis_set_flag (ARGN EXPORT ${BASIS_EXPORT_DEFAULT})
  if (ARGN_USE_BASIS_UTILITIES AND ARGN_NO_BASIS_UTILITIES)
    message (FATAL_ERROR "Target ${TARGET_UID}: Options USE_BASIS_UTILITIES and NO_BASIS_UTILITIES are mutually exclusive!")
  endif ()
  if (ARGN_USE_BASIS_UTILITIES)
    set (USES_BASIS_UTILITIES TRUE)
  elseif (ARGN_NO_BASIS_UTILITIES)
    set (USES_BASIS_UTILITIES FALSE)
  else ()
    set (USES_BASIS_UTILITIES ${BASIS_UTILITIES})
  endif ()
  if (ARGN_SHARED AND (ARGN_EXECUTABLE OR ARGN_LIBEXEC))
    message (FATAL_ERROR "Target ${TARGET_UID}: Options SHARED and EXECUTABLE or LIBEXEC are mutually exclusive!")
  endif ()
  if (ARGN_SHARED)
    set (TYPE LIBRARY)
  else ()
    set (TYPE EXECUTABLE)
  endif ()
  string (TOLOWER "${TYPE}" type)
  message (STATUS "Adding MATLAB ${type} ${TARGET_UID}...")
  # IS_TEST flag
  basis_sanitize_for_regex (RE "${PROJECT_TESTING_DIR}")
  if (CMAKE_CURRENT_SOURCE_DIR MATCHES "^${RE}")
    set (IS_TEST TRUE)
  else ()
    set (IS_TEST FALSE)
  endif ()
  # output directory
  if (IS_TEST)
    set (LIBRARY_OUTPUT_DIRECTORY "${TESTING_LIBRARY_DIR}")
    if (ARGN_LIBEXEC)
      set (RUNTIME_OUTPUT_DIRECTORY "${TESTING_LIBEXEC_DIR}")
    else ()
      set (RUNTIME_OUTPUT_DIRECTORY "${TESTING_RUNTIME_DIR}")
    endif ()
  else ()
    set (LIBRARY_OUTPUT_DIRECTORY "${BINARY_LIBRARY_DIR}")
    if (ARGN_LIBEXEC)
      set (RUNTIME_OUTPUT_DIRECTORY "${BINARY_LIBEXEC_DIR}")
    else ()
      set (RUNTIME_OUTPUT_DIRECTORY "${BINARY_RUNTIME_DIR}")
    endif ()
  endif ()
  # installation component
  if (ARGN_COMPONENT)
    if (NOT ARGN_LIBRARY_COMPONENT)
      set (ARGN_LIBRARY_COMPONENT "${ARGN_COMPONENT}")
    endif ()
    if (NOT ARGN_RUNTIME_COMPONENT)
      set (ARGN_RUNTIME_COMPONENT "${ARGN_COMPONENT}")
    endif ()
  endif ()
  if (NOT ARGN_RUNTIME_COMPONENT)
    set (ARGN_RUNTIME_COMPONENT "${BASIS_RUNTIME_COMPONENT}")
  endif ()
  if (NOT ARGN_RUNTIME_COMPONENT)
    set (ARGN_RUNTIME_COMPONENT "Unspecified")
  endif ()
  if (NOT ARGN_LIBRARY_COMPONENT)
    set (ARGN_LIBRARY_COMPONENT "${BASIS_LIBRARY_COMPONENT}")
  endif ()
  if (NOT ARGN_LIBRARY_COMPONENT)
    set (ARGN_LIBRARY_COMPONENT "Unspecified")
  endif ()
  # installation directories
  if (ARGN_DESTINATION)
    if (NOT ARGN_RUNTIME_DESTINATION)
      set (ARGN_RUNTIME_DESTINATION "${ARGN_DESTINATION}")
    endif ()
    if (NOT ARGN_LIBRARY_DESTINATION)
      set (ARGN_LIBRARY_DESTINATION "${ARGN_DESTINATION}")
    endif ()
    if (NOT ARGN_HEADER_DESTINATION)
      set (ARGN_HEADER_DESTINATION "${ARGN_DESTINATION}")
    endif ()
  endif ()
  if (NOT ARGN_RUNTIME_DESTINATION AND NOT IS_TEST)
    if (ARGN_LIBEXEC)
      set (ARGN_RUNTIME_DESTINATION "${INSTALL_LIBEXEC_DIR}")
    else ()
      set (ARGN_RUNTIME_DESTINATION "${INSTALL_RUNTIME_DIR}")
    endif ()
  endif ()
  if (NOT ARGN_LIBRARY_DESTINATION AND NOT IS_TEST)
    set (ARGN_LIBRARY_DESTINATION "${INSTALL_LIBRARY_DIR}")
  endif ()
  if (NOT ARGN_HEADER_DESTINATION AND NOT IS_TEST)
    set (ARGN_HEADER_DESTINATION ".")
  endif ()
  if (ARGN_RUNTIME_DESTINATION MATCHES "^[nN][oO][nN][eE]$")
    set (ARGN_RUNTIME_DESTINATION)
  endif ()
  if (ARGN_LIBRARY_DESTINATION MATCHES "^[nN][oO][nN][eE]$")
    set (ARGN_LIBRARY_DESTINATION)
  endif ()
  if (ARGN_HEADER_DESTINATION MATCHES "^[nN][oO][nN][eE]$")
    set (ARGN_HEADER_DESTINATION)
  elseif (NOT IS_ABSOLUTE ARGN_HEADER_DESTINATION)
    set (ARGN_HEADER_DESTINATION "${BINARY_INCLUDE_DIR}/${ARGN_HEADER_DESTINATION}")
  endif ()
  # whether to compile and compilation flags (for mcc)
  if ("^${TYPE}$" STREQUAL "^LIBRARY$")
    set (COMPILE TRUE)
  elseif (BASIS_COMPILE_MATLAB)
    if (MATLAB_MCC_EXECUTABLE)
      set (COMPILE TRUE)
    else ()
      set (COMPILE FALSE)
      message (WARNING "MATLAB Compiler not found. Will generate a wrapper script for target"
                       " ${TARGET_UID} which executes the MATLAB code using the -r option of"
                       " the MATLAB interpreter. It is recommended to compile the MATLAB code"
                       " using the MATLAB Compiler if possible, however. Therefore, make sure"
                       " that the MATLAB Compiler is available and check the value of the"
                       " advanced MATLAB_MCC_EXECUTABLE variable in CMake."
                       "\nMake sure to include MATLAB{mcc} as project dependency.")
    endif ()
  else ()
    set (COMPILE FALSE)
  endif ()
  if (WIN32 AND NOT COMPILE)
    # TODO implement generation of Windows Command on Windows
    set (CONTACT)
    if (PROJECT_CONTACT)
      set (CONTACT "\n\nYou may further want to contact ${PROJECT_CONTACT} in order to ask"
                   " for a binary distribution package which contains pre-build binaries"
                   " created using the MATLAB Compiler and download the MATLAB Compiler"
                   " Runtime only if no MATLAB Compiler license is available to you.")
      basis_list_to_string (CONTACT ${CONTACT})
    endif ()
    if (NOT BASIS_COMPILE_MATLAB)
      basis_update_value (BASIS_COMPILE_MATLAB ON)
    endif ()
    message (FATAL_ERROR "The optional generation of a Windows Command which executes"
                         " the MATLAB code using the -r option of the MATLAB interpreter"
                         " as an alternative to the build of the MATLAB sources using"
                         " the MATLAB Compiler is not yet implemented. You will have"
                         " to obtain a MATLAB Compiler license and set the advanced"
                         " MATLAB_MCC_EXECUTABLE variable in CMake or use this package"
                         " on a Unix system instead.${CONTACT}")
  endif ()
  if (COMPILE)
    if (NOT MATLAB_MCC_EXECUTABLE)
      message (FATAL_ERROR "MATLAB Compiler not found! It is required to build target ${TARGET_UID}."
                           " Ensure that MATLAB{mcc} is declared as project dependency"
                           " and check the setting of MATLAB_DIR and/or MATLAB_MCC_EXECUTABLE.")
    endif ()
    basis_add_mcc_options()
  else ()
    if (NOT MATLAB_EXECUTABLE)
      message (FATAL_ERROR "MATLAB not found! It is required to build target ${TARGET_UID}."
                           " Ensure that MATLAB{matlab} is declared as project dependency"
                           " and check the setting of MATLAB_DIR and/or MATLAB_EXECUTABLE.")
    endif ()
  endif ()
  if ("^${TYPE}$" STREQUAL "^EXECUTABLE$")
    set (COMPILE_FLAGS "${BASIS_MCC_FLAGS}")
  else ()
    set (COMPILE_FLAGS "")
  endif ()
  # output file name prefix/suffix
  if ("^${TYPE}$" STREQUAL "^EXECUTABLE$")
    set (PREFIX)
    if (WIN32 AND "^${COMPILE_FLAGS}$" STREQUAL "^NOMCC$")
      set (SUFFIX ".cmd")
    else ()
      set (SUFFIX)
    endif ()
  else ()
    if (WIN32)
      set (PREFIX)
      set (SUFFIX .lib) # link library file extension
    elseif (APPLE)
      set (PREFIX lib)
      set (SUFFIX .dylib)
    else ()
      set (PREFIX lib)
      set (SUFFIX .so)
    endif ()
  endif ()
  # configure (.in) source files
  basis_configure_sources (SOURCES ${SOURCES})
  # add custom target
  add_custom_target (${TARGET_UID} ALL SOURCES ${SOURCES})
  basis_get_target_name (OUTPUT_NAME "${TARGET_UID}")
  get_directory_property (INCLUDE_DIRS INCLUDE_DIRECTORIES)
  get_directory_property (LINK_DIRS    LINK_DIRECTORIES)
  set_target_properties (
    ${TARGET_UID}
    PROPERTIES
      LANGUAGE                  "MATLAB"
      BASIS_TYPE                "MCC_${TYPE}"
      BASIS_UTILITIES           FALSE # TODO Implement utilities for MATLAB
      BASIS_INCLUDE_DIRECTORIES "${INCLUDE_DIRS}"
      BASIS_LINK_DIRECTORIES    "${LINK_DIRS}"
      BUILD_DIRECTORY           "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${TARGET_UID}"
      SOURCE_DIRECTORY          "${CMAKE_CURRENT_SOURCE_DIR}"
      BINARY_DIRECTORY          "${CMAKE_CURRENT_BINARY_DIR}"
      LIBRARY_OUTPUT_DIRECTORY  "${LIBRARY_OUTPUT_DIRECTORY}"
      LIBRARY_INSTALL_DIRECTORY "${ARGN_LIBRARY_DESTINATION}"
      LIBRARY_HEADER_DIRECTORY  "${ARGN_HEADER_DESTINATION}"
      LIBRARY_COMPONENT         "${ARGN_LIBRARY_COMPONENT}"
      RUNTIME_OUTPUT_DIRECTORY  "${RUNTIME_OUTPUT_DIRECTORY}"
      RUNTIME_INSTALL_DIRECTORY "${ARGN_RUNTIME_DESTINATION}"
      RUNTIME_COMPONENT         "${ARGN_RUNTIME_COMPONENT}"
      PREFIX                    "${PREFIX}"
      OUTPUT_NAME               "${OUTPUT_NAME}"
      SUFFIX                    "${SUFFIX}"
      COMPILE_FLAGS             "${COMPILE_FLAGS}"
      COMPILE                   "${COMPILE}"
      LINK_DEPENDS              ""
      EXPORT                    ${EXPORT}
      LIBEXEC                   ${ARGN_LIBEXEC}
      TEST                      ${IS_TEST}
  )
  # finalize target
  if (ARGN_FINAL)
    basis_finalize_targets (${TARGET_UID})
  endif ()
  # add target to list of targets
  basis_set_project_property (APPEND PROPERTY TARGETS "${TARGET_UID}")
  message (STATUS "Adding MATLAB ${type} ${TARGET_UID}... - done")
endfunction ()

# ============================================================================
# custom build commands
# ============================================================================

# ----------------------------------------------------------------------------
## @brief Add custom command for build of MEX-file.
#
# This function is called by basis_finalize_targets() which in turn is called
# by basis_project_end(), i.e., the end of the root CMake configuration file
# of the (sub-)project.
#
# @param [in] TARGET_UID Name/UID of custom target added by basis_add_mex_file().
#
# @sa basis_add_mex_file()
#
# @ingroup CMakeUtilities
function (basis_build_mex_file TARGET_UID)
  # does this target exist ?
  basis_get_target_uid (TARGET_UID "${TARGET_UID}")
  if (NOT TARGET "${TARGET_UID}")
    message (FATAL_ERROR "Unknown build target: ${TARGET_UID}")
  endif ()
  if (BASIS_VERBOSE)
    message (STATUS "Adding build command for target ${TARGET_UID}...")
  endif ()
  # get target properties
  basis_get_target_name (TARGET_NAME ${TARGET_UID})
  set (
    PROPERTIES
      BASIS_TYPE
      BASIS_UTILITIES
      BASIS_INCLUDE_DIRECTORIES
      BASIS_LINK_DIRECTORIES
      BUILD_DIRECTORY
      SOURCE_DIRECTORY
      BINARY_DIRECTORY
      LIBRARY_OUTPUT_DIRECTORY
      LIBRARY_INSTALL_DIRECTORY
      LIBRARY_COMPONENT
      PREFIX
      OUTPUT_NAME
      SUFFIX
      COMPILE_FLAGS
      LINK_DEPENDS
      LINK_FLAGS
      MFILE
      EXPORT
      SOURCES
  )
  get_target_property (IS_TEST ${TARGET_UID} TEST)
  foreach (PROPERTY ${PROPERTIES})
    get_target_property (${PROPERTY} ${TARGET_UID} ${PROPERTY})
    if (NOT ${PROPERTY})
      set (${PROPERTY})
    endif ()
  endforeach ()
  # sanity check of property values
  if (NOT BASIS_TYPE MATCHES "^MEX$")
    message (FATAL_ERROR "Target ${TARGET_UID}: Invalid BASIS_TYPE: ${BASIS_TYPE}")
  endif ()
  list (GET SOURCES 0 BUILD_DIR) # CMake <3.1 stores path to internal build directory here
  if (BUILD_DIR MATCHES "CMakeFiles")
    list (REMOVE_AT SOURCES 0)
  endif ()
  set (BUILD_DIR "${BUILD_DIRECTORY}.dir")
  if (NOT SOURCES)
    message (FATAL_ERROR "Target ${TARGET_UID}: Empty SOURCES list!"
                         " Have you accidentally modified this read-only property or"
                         " is your (newer) CMake version not compatible with BASIS?")
  endif ()
  if (NOT LIBRARY_COMPONENT)
    set (LIBRARY_COMPONENT "Unspecified")
  endif ()
  if (MFILE)
    if (NOT IS_ABSOLUTE "${MFILE}")
      set (MFILE "${SOURCE_DIRECTORY}/${MFILE}")
    endif ()
    if (NOT EXISTS "${MFILE}")
      message (FATAL_ERROR "M-file ${MFILE} of MEX-file target ${TARGET_UID} does not exist!")
    endif ()
  endif ()
  # output name
  if (NOT OUTPUT_NAME)
    set (OUTPUT_NAME "${TARGET_NAME}")
  endif ()
  if (SUFFIX)
    set (OUTPUT_NAME "${OUTPUT_NAME}${SUFFIX}")
  endif ()
  string (REGEX REPLACE "/+" "/" PREFIX "${PREFIX}")
  string (REGEX REPLACE "/$" ""  PREFIX "${PREFIX}")
  if (PREFIX AND NOT PREFIX MATCHES "^/")
    set (PREFIX "/${PREFIX}")
  endif ()
  # initialize dependencies of custom build command
  set (DEPENDS ${SOURCES})
  # get list of libraries to link to
  set (LINK_LIBS)
  foreach (LIB ${LINK_DEPENDS})
    basis_get_target_uid (UID "${LIB}")
    if (TARGET ${UID})
      basis_get_target_location (LIB_FILE ${UID} ABSOLUTE)
      string (REPLACE "<${BASIS_GE_CONFIG}>" "{CONFIG}" LIB_FILE "${LIB_FILE}")
      list (APPEND DEPENDS ${UID})
    else ()
      set (LIB_FILE "${LIB}")
    endif ()
    list (APPEND LINK_LIBS "${LIB_FILE}")
  endforeach ()
  get_filename_component (OUTPUT_NAME_WE "${OUTPUT_NAME}" NAME_WE)
  # decompose user supplied MEX switches
  macro (extract VAR)
    string (REGEX REPLACE "${VAR}=\"([^\"]+)\"|${VAR}=([^\" ])*" "" COMPILE_FLAGS "${COMPILE_FLAGS}")
    if (CMAKE_MATCH_1)
      set (${VAR} "${CMAKE_MATCH_1}")
    elseif (CMAKE_MATCH_2)
      set (${VAR} "${CMAKE_MATCH_2}")
    else ()
      set (${VAR})
    endif ()
  endmacro ()
  if (UNIX)
    extract (CC)
    extract (CFLAGS)
    extract (CXX)
    extract (CXXFLAGS)
    extract (CLIBS)
    extract (CXXLIBS)
    extract (LD)
    extract (LDXX)
    extract (LDFLAGS)
    extract (LDCXXFLAGS)
    if (LINK_FLAGS)
      set (LDFLAGS "${LDFLAGS} ${LINK_FLAGS}")
    endif ()
    # set defaults for not provided options
    if (NOT CC)
      set (CC "${CMAKE_C_COMPILER}")
    endif ()
    if (NOT CFLAGS)
      set (CFLAGS "${CMAKE_C_FLAGS}")
    endif ()
    if (NOT CFLAGS MATCHES "( |^)-fPIC( |$)")
      set (CFLAGS "-fPIC ${CFLAGS}")
    endif ()
    if (NOT CXX)
      set (CXX "${CMAKE_CXX_COMPILER}")
    endif ()
    if (NOT CXXFLAGS)
      set (CXXFLAGS "${CMAKE_CXX_FLAGS}")
    endif ()
    if (NOT CXXFLAGS MATCHES "( |^)-fPIC( |$)")
      set (CXXFLAGS "-fPIC ${CXXFLAGS}")
    endif ()
    if (NOT LD)
      set (LD "${CMAKE_CXX_COMPILER}") # do not use CMAKE_LINKER here
    endif ()
    if (NOT LDFLAGS)
      set (LDFLAGS "\$LDFLAGS ${CMAKE_SHARED_LINKER_FLAGS}")
    endif ()
    # We chose to use CLIBS and CXXLIBS instead of the -L and -l switches
    # to add also link libraries added via basis_target_link_libraries()
    # because the MEX script will not use these arguments if CLIBS or CXXLIBS
    # is set. Moreover, the -l switch can only be used to link to a shared
    # library and not a static one (on UNIX).
    #foreach (LIB ${LINK_LIBS})
    #  if (LIB MATCHES "[/\\\.]")
    #    set (CXXLIBS "${CXXLIBS} ${LIB}")
    #  endif ()
    #endforeach ()
  endif ()
  # get remaining switches
  basis_string_to_list (MEX_USER_ARGS "${COMPILE_FLAGS}")
  # assemble MEX switches
  set (MEX_ARGS)
  if (UNIX)
    list (APPEND MEX_ARGS "CC=${CC}" "CFLAGS=${CFLAGS}")           # C compiler and flags
    if (CLIBS)
      list (APPEND MEX_ARGS "CLIBS=${CLIBS}")                      # C link libraries
    endif ()
    list (APPEND MEX_ARGS "CXX=${CXX}" "CXXFLAGS=${CXXFLAGS}")     # C++ compiler and flags
    if (CXXLIBS)
      list (APPEND MEX_ARGS "CXXLIBS=${CXXLIBS}")                  # C++ link libraries
    endif ()
    if (LD)
      list (APPEND MEX_ARGS "LD=${LD}")                            # C linker
    endif ()
    if (LDFLAGS)
      list (APPEND MEX_ARGS "LDFLAGS=${LDFLAGS}")                  # C link flags
    endif ()
    if (LDCXX)
      list (APPEND MEX_ARGS "LDCXX=${LDCXX}")                      # C++ linker
    endif ()
    if (LDCXXFLAGS)
      list (APPEND MEX_ARGS "LDCXXFLAGS=${LDCXXFLAGS}")            # C++ link flags
    endif ()
  endif ()
  list (APPEND MEX_ARGS "-outdir" "${BUILD_DIR}")                # output directory
  list (APPEND MEX_ARGS "-output" "${OUTPUT_NAME_WE}")           # output name (w/o extension)
  foreach (INCLUDE_PATH ${BASIS_INCLUDE_DIRECTORIES})            # include directories
    list (FIND MEX_ARGS "-I${INCLUDE_PATH}" IDX)                 # as specified via
    if (INCLUDE_PATH AND IDX EQUAL -1)                           # basis_include_directories()
      list (APPEND MEX_ARGS "-I${INCLUDE_PATH}")
    endif ()
  endforeach ()
  set (MEX_LIBPATH)
  set (MEX_LIBS)
  foreach (LINK_DIR ${BASIS_LINK_DIRECTORIES})                    # link directories
    if (WIN32)
      string (REPLACE "/" "\\" LINK_DIR "${LINK_DIR}")
      set (LINK_DIR "/LIBPATH:\\\"${LINK_DIR}\\\"")
    else ()
      set (LINK_DIR "-L${LINK_DIR}")
    endif ()
    list (APPEND MEX_LIBPATH "${LINK_DIR}")                      # as specified via basis_link_directories()
  endforeach ()
  foreach (LIBRARY ${LINK_LIBS})                                # link libraries
    get_filename_component (LINK_DIR "${LIBRARY}" PATH)         # as specified via basis_target_link_libraries()
    get_filename_component (LINK_LIB "${LIBRARY}" NAME)
    string (REGEX REPLACE "\\.(so|dylib)(\\.[0-9]+)*$" "" LINK_LIB "${LINK_LIB}")
    if (CMAKE_SYSTEM_NAME MATCHES "Darwin")
      # cf. https://github.com/cmake-basis/BASIS/issues/443
      string (REGEX REPLACE "\\.a(\\.[0-9]+)*$" "" LINK_LIB "${LINK_LIB}")
    endif ()
    string (REGEX REPLACE "^-l" "" LINK_LIB "${LINK_LIB}")
    if (LINK_DIR)
      if (WIN32)
        string (REPLACE "/" "\\" LINK_DIR "${LINK_DIR}")
        set (LINK_DIR "/LIBPATH:\\\"${LINK_DIR}\\\"")
      else ()
        set (LINK_DIR "-L${LINK_DIR}")
      endif ()
      list (APPEND MEX_LIBPATH "${LINK_DIR}")
    endif ()
    if (WIN32)
      if (NOT LINK_LIB MATCHES "\\.lib$")
        set (LINK_LIB "${LINK_LIB}.lib")
      endif ()
    else ()
      string (REGEX REPLACE "^lib" "" LINK_LIB "${LINK_LIB}")
      set (LINK_LIB "-l${LINK_LIB}")
    endif ()
    list (APPEND MEX_LIBS "${LINK_LIB}")
  endforeach ()
  if (MEX_LIBPATH)
    list (REMOVE_DUPLICATES MEX_LIBPATH)
  endif ()
  # do not remove duplicate entries in MEX_LIBS which may be needed
  # to resolve (cyclic) dependencies between statically linked libraries
  # (cf. https://github.com/cmake-basis/BASIS/issues/444)
  if (MEX_LIBPATH OR MEX_LIBS)
    if (WIN32)
      basis_list_to_delimited_string (MEX_LIBPATH " " NOAUTOQUOTE ${MEX_LIBPATH})
      basis_list_to_delimited_string (MEX_LIBS    " " ${MEX_LIBS})
      list (APPEND MEX_ARGS "LINKFLAGS#$LINKFLAGS ${MEX_LIBPATH} ${MEX_LIBS}")
    else ()
      list (APPEND MEX_ARGS ${MEX_LIBPATH} ${MEX_LIBS})
    endif ()
  endif ()
  # other user switches 
  list (APPEND MEX_ARGS ${MEX_USER_ARGS})
  # source files
  list (APPEND MEX_ARGS ${SOURCES})
  # build command for invocation of MEX script
  set (BUILD_CMD     "${MATLAB_MEX_EXECUTABLE}" -v ${MEX_ARGS})
  set (BUILD_SCRIPT  "${BUILD_DIR}/build.cmake")
  set (BUILD_LOG     "${BUILD_DIR}/build.log")
  set (BUILD_OUTPUT  "${LIBRARY_OUTPUT_DIRECTORY}${PREFIX}/${OUTPUT_NAME}")
  set (BUILD_OUTPUTS "${BUILD_OUTPUT}")
  if (MFILE)
    set (BUILD_MFILE "${LIBRARY_OUTPUT_DIRECTORY}${PREFIX}/${OUTPUT_NAME_WE}.m")
    list (APPEND BUILD_OUTPUTS "${BUILD_MFILE}")
  else ()
    set (BUILD_MFILE)
  endif ()
  # configure build script
  configure_file ("${BASIS_SCRIPT_EXECUTE_PROCESS}" "${BUILD_SCRIPT}" @ONLY)
  # relative paths used for comments of commands
  file (RELATIVE_PATH REL "${CMAKE_BINARY_DIR}" "${BUILD_OUTPUT}")
  # add custom command to build executable using MEX script
  add_custom_command (
    OUTPUT "${BUILD_OUTPUT}"
    # rebuild when input sources were modified
    DEPENDS "${BUILD_SCRIPT}" "${CMAKE_CURRENT_LIST_FILE}" ${DEPENDS}
    # invoke MEX script, wrapping the command in CMake execute_process()
    # command allows for inspection of command output for error messages
    # and specification of timeout
    COMMAND "${CMAKE_COMMAND}"
            "-DCONFIG=$<${BASIS_GE_CONFIG}>"
            "-DWORKING_DIRECTORY=${BUILD_DIR}"
            "-DTIMEOUT=${BASIS_MEX_TIMEOUT}"
            "-DERROR_EXPRESSION=[E|e]rror:"
            "-DOUTPUT_FILE=${BUILD_LOG}"
            "-DERROR_FILE=${BUILD_LOG}"
            "-DVERBOSE=OFF"
            "-DLOG_ARGS=ON"
            "-P" "${BUILD_SCRIPT}"
    # post-build command
    COMMAND "${CMAKE_COMMAND}" -E copy   "${BUILD_DIR}/${OUTPUT_NAME}" "${BUILD_OUTPUT}"
    COMMAND "${CMAKE_COMMAND}" -E remove "${BUILD_DIR}/${OUTPUT_NAME}"
    # inform user where build log can be found
    COMMAND "${CMAKE_COMMAND}" -E echo "Build log written to ${BUILD_LOG}"
    # comment
    COMMENT "Building MEX-file ${REL}..."
    VERBATIM
  )
  if (BUILD_MFILE)
    add_custom_command (
      OUTPUT  "${BUILD_MFILE}"
      DEPENDS "${MFILE}"
      COMMAND "${CMAKE_COMMAND}" -E copy "${MFILE}" "${BUILD_MFILE}"
      COMMENT "Copying M-file of ${REL}..."
    )
  endif ()
  # add custom target
  add_custom_target (_${TARGET_UID} DEPENDS ${BUILD_OUTPUTS} SOURCES ${SOURCES})
  add_dependencies (${TARGET_UID} _${TARGET_UID})
  # cleanup on "make clean"
  set_property (
    DIRECTORY
    APPEND PROPERTY
      ADDITIONAL_MAKE_CLEAN_FILES
        "${BUILD_DIR}/${OUTPUT_NAME}"
        "${BUILD_OUTPUTS}"
        "${BUILD_LOG}"
  )
  # export target
  if (EXPORT)
    basis_add_custom_export_target (${TARGET_UID} "${IS_TEST}")
  endif ()
  # install MEX-file
  if (LIBRARY_INSTALL_DIRECTORY)
    install (
      FILES       ${BUILD_OUTPUTS}
      DESTINATION "${LIBRARY_INSTALL_DIRECTORY}${PREFIX}"
      COMPONENT   "${LIBRARY_COMPONENT}"
    )
  endif ()
  if (BASIS_VERBOSE)
    message (STATUS "Adding build command for target ${TARGET_UID}... - done")
  endif ()
endfunction ()

# ----------------------------------------------------------------------------
## @brief Add custom command for build of MATLAB Compiler target.
#
# This function is called by basis_finalize_targets() which in turn is called
# by basis_project_end(), i.e., the end of the root CMake configuration file
# of the (sub-)project.
#
# @param [in] TARGET_UID Name/UID of custom target added by basis_add_mcc_target().
#
# @sa basis_add_mcc_target()
#
# @ingroup CMakeUtilities
function (basis_build_mcc_target TARGET_UID)
  # does this target exist ?
  basis_get_target_uid (TARGET_UID "${TARGET_UID}")
  if (NOT TARGET "${TARGET_UID}")
    message (FATAL_ERROR "Unknown target ${TARGET_UID}!")
  endif ()
  if (BASIS_VERBOSE)
    message (STATUS "Adding build command for target ${TARGET_UID}...")
  endif ()
  # get target properties
  basis_get_target_name (TARGET_NAME ${TARGET_UID})
  set (
    PROPERTIES
      BASIS_TYPE
      BASIS_UTILITIES
      BASIS_INCLUDE_DIRECTORIES
      BASIS_LINK_DIRECTORIES
      BUILD_DIRECTORY
      SOURCE_DIRECTORY
      BINARY_DIRECTORY
      LIBRARY_OUTPUT_DIRECTORY
      LIBRARY_INSTALL_DIRECTORY
      LIBRARY_HEADER_DIRECTORY
      LIBRARY_COMPONENT
      RUNTIME_OUTPUT_DIRECTORY
      RUNTIME_INSTALL_DIRECTORY
      RUNTIME_COMPONENT
      PREFIX
      OUTPUT_NAME
      SUFFIX
      SOURCES
      COMPILE_FLAGS
      COMPILE
      LINK_DEPENDS
      EXPORT
  )
  get_target_property (IS_TEST ${TARGET_UID} TEST)
  foreach (PROPERTY ${PROPERTIES})
    get_target_property (${PROPERTY} ${TARGET_UID} ${PROPERTY})
    if (NOT ${PROPERTY})
      set (${PROPERTY})
    endif ()
  endforeach ()
  # sanity checks of property values
  set (EXECUTABLE FALSE)
  set (LIBEXEC    FALSE)
  set (LIBRARY    FALSE)
  if (BASIS_TYPE MATCHES "^MCC_(EXECUTABLE|LIBEXEC|LIBRARY)$")
    set (${CMAKE_MATCH_1} TRUE)
    if (LIBEXEC)
      set (EXECUTABLE TRUE)
    endif ()
  else ()
    message (FATAL_ERROR "Target ${TARGET_UID}: Invalid BASIS_TYPE: ${BASIS_TYPE}")
  endif ()
  list (GET SOURCES 0 BUILD_DIR) # CMake <3.1 stores path to internal build directory here
  if (BUILD_DIR MATCHES "CMakeFiles")
    list (REMOVE_AT SOURCES 0)
  endif ()
  set (BUILD_DIR "${BUILD_DIRECTORY}.dir")
  if (NOT SOURCES)
    message (FATAL_ERROR "Target ${TARGET_UID}: Empty SOURCES list!"
                         " Have you accidentally modified this read-only property or"
                         " is your (newer) CMake version not compatible with BASIS?")
  endif ()
  list (GET SOURCES 0 MAIN_SOURCE) # entry point
  if (NOT RUNTIME_COMPONENT)
    set (RUNTIME_COMPONENT "Unspecified")
  endif ()
  if (NOT LIBRARY_COMPONENT)
    set (LIBRARY_COMPONENT "Unspecified")
  endif ()
  # output name
  if (NOT OUTPUT_NAME)
    set (OUTPUT_NAME "${TARGET_NAME}")
  endif ()
  if (PREFIX)
    set (OUTPUT_NAME "${PREFIX}${OUTPUT_NAME}")
  endif ()
  if (SUFFIX)
    set (OUTPUT_NAME "${OUTPUT_NAME}${SUFFIX}")
  endif ()
  get_filename_component (OUTPUT_NAME_WE "${OUTPUT_NAME}" NAME_WE)
  # MCC only allows alpha-numeric characters and underscores
  # TODO: Figure out how to build a shared library without this restriction
  #       (cf. https://github.com/cmake-basis/BASIS/issues/410).
  if (NOT OUTPUT_NAME MATCHES "[a-zA-Z][_a-zA-Z0-9]*")
    message (FATAL_ERROR "Target ${TARGET_UID} has invalid output name ${OUTPUT_NAME}."
                         "MCC only allows alpha-numeric characters and underscores. "
                         "See GitHub issue #410 for updates on this restriction at\n"
                         "https://github.com/cmake-basis/BASIS/issues/410")
  endif ()
  set (MCC_OUTPUT_NAME    "${OUTPUT_NAME}")
  set (MCC_OUTPUT_NAME_WE "${OUTPUT_NAME_WE}")
  #string (REGEX REPLACE "\\+|-" "_" MCC_OUTPUT_NAME "${OUTPUT_NAME}")
  #get_filename_component (MCC_OUTPUT_NAME_WE "${MCC_OUTPUT_NAME}" NAME_WE)
  # initialize dependencies of custom build command
  set (DEPENDS ${SOURCES})
  # build output file and comment
  file (RELATIVE_PATH REL "${CMAKE_BINARY_DIR}" "${BUILD_DIR}/${OUTPUT_NAME}")
  if (LIBRARY)
    set (BUILD_OUTPUT "${LIBRARY_OUTPUT_DIRECTORY}/${OUTPUT_NAME}")
    set (BUILD_COMMENT "Building MATLAB library ${REL}...")
  else ()
    set (BUILD_OUTPUT "${RUNTIME_OUTPUT_DIRECTORY}/${OUTPUT_NAME}")
    set (BUILD_COMMENT "Building MATLAB executable ${REL}...")
  endif ()
  # --------------------------------------------------------------------------
  # assemble build command for build of executable wrapper script
  if (EXECUTABLE AND NOT COMPILE)
    # used to recognize source files which are located in the build tree
    basis_sanitize_for_regex (BINARY_CODE_DIR_RE "${BINARY_CODE_DIR}")
    # main MATLAB function and search path
    get_filename_component (MATLAB_COMMAND "${MAIN_SOURCE}" NAME_WE)
    get_filename_component (SOURCE_DIR     "${MAIN_SOURCE}" PATH)
    get_filename_component (SOURCE_PACKAGE "${SOURCE_DIR}"       NAME)
    if (SOURCE_PACKAGE MATCHES "^\\+")
      get_filename_component (SOURCE_DIR "${SOURCE_DIR}" PATH)
      set (MATLAB_COMMAND "${SOURCE_PACKAGE}.${MATLAB_COMMAND}")
      string (REGEX REPLACE "^\\+" "" MATLAB_COMMAND "${MATLAB_COMMAND}")
    else ()
      set (SOURCE_PACKAGE "${MATLAB_COMMAND}")
    endif ()
    basis_get_relative_path (DIR "${PROJECT_SOURCE_DIR}" "${SOURCE_DIR}")
    set (BINARY_DIR "${PROJECT_BINARY_DIR}/${DIR}") # location of configured sources
    # output file
    set (OUTPUT_FILE "${BUILD_OUTPUT}")
    get_filename_component (OUTPUT_DIR "${OUTPUT_FILE}" PATH)
    # installation
    set (INSTALL_FILE       "${BUILD_DIR}/${OUTPUT_NAME}")                          # file to be installed
    set (INSTALL_DIR        "${CMAKE_INSTALL_PREFIX}/${RUNTIME_INSTALL_DIRECTORY}") # location of installed wrapper
    set (INSTALL_SOURCE_DIR "${INSTALL_MATLAB_LIBRARY_DIR}")                        # location of installed MATLAB sources
    if (NOT SOURCE_PACKAGE MATCHES "^\\+")
      set (INSTALL_SOURCE_DIR "${INSTALL_SOURCE_DIR}/${SOURCE_PACKAGE}")
    endif ()
    # startup file
    if (SOURCE_PACKAGE MATCHES "^\\+")
      set (BUILD_STARTUP_FILE   "${BINARY_DIR}/${SOURCE_PACKAGE}/startup.m")
      set (INSTALL_STARTUP_FILE "${BUILD_DIR}/startup.m")
      set (INSTALL_STARTUP_DIR  "${CMAKE_INSTALL_PREFIX}/${INSTALL_SOURCE_DIR}/${SOURCE_PACKAGE}")
    else ()
      set (BUILD_STARTUP_FILE   "${BINARY_DIR}/startup.m")
      set (INSTALL_STARTUP_FILE "${BUILD_DIR}/startup.m")
      set (INSTALL_STARTUP_DIR  "${CMAKE_INSTALL_PREFIX}/${INSTALL_SOURCE_DIR}")
    endif ()
    get_filename_component (BUILD_STARTUP_DIR "${BUILD_STARTUP_FILE}" PATH)
    list (APPEND BUILD_OUTPUT "${BUILD_STARTUP_FILE}")
    list (APPEND BUILD_OUTPUT "${INSTALL_STARTUP_FILE}")
    # MATLAB search path
    #
    # The following paths are written to the startup M-file such
    # that they can conveniently be added to the search path within an
    # interactive MATLAB sesssion. The path to this startup file is
    # added to the search path by the wrapper executable first, and
    # then this startup file is evaluated which adds all additional
    # paths. If no startup file is specified, all paths are added
    # to the search path on the command-line in the wrapper script.
    set (BUILD_MATLABPATH)
    set (INSTALL_MATLABPATH)
    # if any source file was configured and hence is located in the
    # build tree instead of the source tree, add corresponding build
    # tree path to BUILD_MATLABPATH as well
    if (SOURCES MATCHES "^${BINARY_CODE_DIR_RE}")
      file (RELATIVE_PATH REL "${BUILD_STARTUP_DIR}" "${BINARY_DIR}")
      if (REL)
        list (APPEND BUILD_MATLABPATH "${REL}")
      else ()
        list (APPEND BUILD_MATLABPATH ".")
      endif ()
    endif ()
    list (APPEND BUILD_MATLABPATH "${SOURCE_DIR}")
    # The following is not required because startup script is located
    # in the same directory as the other sources and basis_generate_matlab_executable()
    # adds this path automatically in order to run the startup script.
    # Moreover, if anyone wants to run the startup script, they have to
    # add this path manually or make it the current directory first.
    #file (RELATIVE_PATH REL "${INSTALL_STARTUP_DIR}" "${CMAKE_INSTALL_PREFIX}/${INSTALL_SOURCE_DIR}")
    #if (REL)
    #  list (APPEND INSTALL_MATLABPATH "${REL}")
    #else ()
    #  list (APPEND INSTALL_MATLABPATH ".")
    #endif ()
    # link dependencies, i.e., MEX-files
    foreach (LINK_DEPEND ${LINK_DEPENDS})
      basis_get_target_uid (UID "${LINK_DEPEND}")
      if (TARGET ${UID})
        basis_get_target_location (LINK_DEPEND ${UID} ABSOLUTE)
        if (LINK_DEPEND MATCHES "\\.mex")
          get_filename_component (LINK_PATH "${LINK_DEPEND}" PATH)
          list (APPEND BUILD_MATLABPATH "${LINK_PATH}")
          list (APPEND DEPENDS ${UID})
        endif ()
        basis_get_target_location (LINK_DEPEND ${UID} POST_INSTALL)
        basis_get_target_property (BUNDLED     ${UID} BUNDLED)
        basis_get_target_property (IMPORTED    ${UID} IMPORTED)
        if (NOT IMPORTED OR BUNDLED)
          file (RELATIVE_PATH REL "${INSTALL_STARTUP_DIR}" "${LINK_DEPEND}")
          if (REL)
            set (LINK_DEPEND "${REL}")
          else ()
            set (LINK_DEPEND ".")
          endif ()
        endif ()
        if (LINK_DEPEND MATCHES "\\.mex")
          get_filename_component (LINK_PATH "${LINK_DEPEND}" PATH)
          list (APPEND INSTALL_MATLABPATH "${LINK_PATH}")
        endif ()
      elseif (IS_ABSOLUTE "${LINK_DEPEND}")
        if (IS_DIRECTORY "${LINK_DEPEND}")
          list (APPEND BUILD_MATLABPATH   "${LINK_DEPEND}")
          list (APPEND INSTALL_MATLABPATH "${LINK_DEPEND}")
        elseif (EXISTS "${LINK_DEPEND}" AND LINK_DEPEND MATCHES "\\.mex")
          get_filename_component (LINK_PATH "${LINK_DEPEND}" PATH)
          list (APPEND BUILD_MATLABPATH   "${LINK_PATH}")
          list (APPEND INSTALL_MATLABPATH "${LINK_PATH}")
        endif ()
      endif ()
    endforeach ()
    # cleanup directory paths
    string (REGEX REPLACE "/+"     "/"    BUILD_MATLABPATH   "${BUILD_MATLABPATH}")
    string (REGEX REPLACE "/+"     "/"    INSTALL_MATLABPATH "${INSTALL_MATLABPATH}")
    string (REGEX REPLACE "/(;|$)" "\\1"  BUILD_MATLABPATH   "${BUILD_MATLABPATH}")
    string (REGEX REPLACE "/(;|$)" "\\1"  INSTALL_MATLABPATH "${INSTALL_MATLABPATH}")
    # remove duplicates
    if (BUILD_MATLABPATH)
      list (REMOVE_DUPLICATES BUILD_MATLABPATH)
    endif ()
    if (INSTALL_MATLABPATH)
      list (REMOVE_DUPLICATES INSTALL_MATLABPATH)
    endif ()
    # configure build script
    set (BUILD_SCRIPT "${BUILD_DIR}/build.cmake")
    configure_file ("${BASIS_MODULE_PATH}/generate_matlab_executable.cmake.in" "${BUILD_SCRIPT}" @ONLY)
    # add custom command to build wrapper executable
    add_custom_command (
      OUTPUT          ${BUILD_OUTPUT}
      # rebuild when input sources were modified
      MAIN_DEPENDENCY "${MAIN_SOURCE}"
      DEPENDS         "${BUILD_SCRIPT}" "${CMAKE_CURRENT_LIST_FILE}" ${DEPENDS}
      # invoke MATLAB Compiler in either MATLAB or standalone mode
      # wrapping command in CMake execute_process () command allows for inspection
      # parsing of command output for error messages and specification of timeout
      COMMAND         "${CMAKE_COMMAND}" "-P" "${BUILD_SCRIPT}"
      # comment
      COMMENT         "${BUILD_COMMENT}"
    )
    # install source files - preserving relative paths in SOURCE_DIR
    foreach (SOURCE IN LISTS SOURCES)
      get_filename_component  (REL "${SOURCE}" PATH)
      if (SOURCE MATCHES "^${BINARY_CODE_DIR_RE}")
        basis_get_relative_path (REL "${BINARY_DIR}" "${REL}")
      else ()
        basis_get_relative_path (REL "${SOURCE_DIR}" "${REL}")
      endif ()
      if (REL MATCHES "^\\.?$|^\\.\\./")
        install (
          FILES       "${SOURCE}"
          DESTINATION "${INSTALL_SOURCE_DIR}"
          COMPONENT   "${RUNTIME_COMPONENT}"
        )
      else ()
        install (
          FILES       "${SOURCE}"
          DESTINATION "${INSTALL_SOURCE_DIR}/${REL}"
          COMPONENT   "${RUNTIME_COMPONENT}"
        )
      endif ()
    endforeach ()
  # --------------------------------------------------------------------------
  # assemble build command for build using MATLAB Compiler
  else ()
    set (INSTALL_FILE "${BUILD_OUTPUT}") # file to be installed
    # get list of libraries to link to (e.g., MEX-file)
    set (LINK_LIBS)
    foreach (LIB ${LINK_DEPENDS})
      basis_get_target_uid (UID "${LIB}")
      if (TARGET ${UID})
        basis_get_target_location (LIB_FILE ${UID} ABSOLUTE)
        string (REPLACE "<${BASIS_GE_CONFIG}>" "{CONFIG}" LIB_FILE "${LIB_FILE}")
        list (APPEND DEPENDS ${UID})
      else ()
        set (LIB_FILE "${LIB}")
      endif ()
      list (APPEND LINK_LIBS "${LIB_FILE}")
    endforeach ()
    # MATLAB search path
    foreach (P ${SOURCES})
      get_filename_component (P "${P}" PATH)
      list (APPEND MATLABPATH "${P}")
    endforeach ()
    list (APPEND MATLABPATH ${BASIS_INCLUDE_DIRECTORIES})
    # MATLAB Compiler arguments
    string (REGEX REPLACE " +" ";" MCC_ARGS "${COMPILE_FLAGS}") # user specified flags
    foreach (P IN LISTS MATLABPATH)                             # search path, -I options
      string (REGEX REPLACE "/(\\+|@).*$" "" P "${P}")          # remove package/class subdirectories
      list (FIND MCC_ARGS "${P}" IDX)                           # skip if already added
      if (IDX GREATER 0)
        math (EXPR IDX "${IDX} - 1")
        list (GET MCC_ARGS ${IDX} OPT)
        if (NOT OPT MATCHES "^-I$")
          set (IDX -1)
        endif ()
      endif ()
      if (EXISTS "${P}" AND IDX EQUAL -1)
        list (APPEND MCC_ARGS -I "${P}")
      endif ()
    endforeach ()
    if (LIBRARY)
      list (APPEND MCC_ARGS -W "cpplib:${MCC_OUTPUT_NAME_WE}" -T link:lib) # build library
    else ()                                                     #      or
      list (APPEND MCC_ARGS -m -o "${MCC_OUTPUT_NAME_WE}")      # build standalone application
    endif ()
    list (APPEND MCC_ARGS -d "${BUILD_DIR}")                    # (temp) output directory
    list (APPEND MCC_ARGS ${SOURCES})                           # source M-files
    foreach (LIB ${LINK_LIBS})                                  # link libraries, e.g. MEX-files
      list (FIND MCC_ARGS "${LIB}" IDX)
      if (LIB AND IDX EQUAL -1)
        list (APPEND MCC_ARGS "-a" "${LIB}")
      endif ()
    endforeach ()
    # build command for invocation of MATLAB Compiler in standalone mode
    set (BUILD_CMD    "${MATLAB_MCC_EXECUTABLE}" ${MCC_USER_ARGS} ${MCC_ARGS})
    set (BUILD_LOG    "${BUILD_DIR}/build.log")
    set (BUILD_SCRIPT "${BUILD_DIR}/build.cmake")
    set (WORKING_DIR  "${SOURCE_DIRECTORY}")
    set (MATLAB_MODE  OFF)
    # build command for invocation of MATLAB Compiler in MATLAB mode
    if (BASIS_MCC_MATLAB_MODE)
      set (MATLAB_MODE ON)
      if (NOT MATLAB_EXECUTABLE)
        message (WARNING "MATLAB executable not found. It is required to build target ${TARGET_UID} in MATLAB mode."
                         " Either set the advanced option BASIS_MCC_MATLAB_MODE to OFF or add MATLAB{matlab} as dependency."
                         " Will build target ${TARGET_UID} in standalone mode instead.")
        set (MATLAB_MODE OFF)
      endif ()
      if (MATLAB_MODE)
        basis_list_to_delimited_string (ARGS "', '" NOAUTOQUOTE ${MCC_USER_ARGS} ${MCC_ARGS})
        set (
          BUILD_CMD
            "${MATLAB_EXECUTABLE}" # run MATLAB
            "-nosplash"            # do not display splash screen on start up
            "-nodesktop"           # run in command line mode
            "-nojvm"               # we do not need the Java Virtual Machine
            "-r" "try, mcc('-v', '${ARGS}'), catch err, fprintf(2, err.message), end, quit force"
        )
      endif ()
    endif ()
    # post-build command
    set (POST_BUILD_COMMAND)
    if (NOT "^${MCC_OUTPUT_NAME}$" STREQUAL "^${OUTPUT_NAME}$")
      list (APPEND POST_BUILD_COMMAND
        COMMAND "${CMAKE_COMMAND}" -E copy
                "${BUILD_DIR}/${MCC_OUTPUT_NAME}"
                "${BUILD_DIR}/${OUTPUT_NAME}"
        COMMAND "${CMAKE_COMMAND}" -E copy
                "${BUILD_DIR}/${MCC_OUTPUT_NAME_WE}.h"
                "${BUILD_DIR}/${OUTPUT_NAME_WE}.h"
      )
    endif ()
    if (LIBRARY)
      list (APPEND POST_BUILD_COMMAND
        COMMAND "${CMAKE_COMMAND}" -E copy
                "${BUILD_DIR}/${OUTPUT_NAME}"
                "${LIBRARY_OUTPUT_DIRECTORY}/${OUTPUT_NAME}"
        COMMAND "${CMAKE_COMMAND}" -E copy
                "${BUILD_DIR}/${OUTPUT_NAME_WE}.h"
                "${LIBRARY_HEADER_DIRECTORY}/${OUTPUT_NAME_WE}.h"
      )
      if (WIN32)
        list (APPEND POST_BUILD_COMMAND
          COMMAND "${CMAKE_COMMAND}" -E copy
                  "${BUILD_DIR}/${OUTPUT_NAME_WE}.dll"
                  "${LIBRARY_OUTPUT_DIRECTORY}/${OUTPUT_NAME_WE}.dll")
      endif ()
    else ()
      if (CMAKE_SYSTEM_NAME STREQUAL "Darwin")
        # TODO: This file should be regenerated if it is missing.
        file (WRITE "${BUILD_DIR}/${OUTPUT_NAME}" "#!/bin/bash\nexec $(dirname $BASH_SOURCE)/${OUTPUT_NAME_WE}.app/Contents/MacOS/${MCC_OUTPUT_NAME_WE}")
        execute_process (COMMAND chmod +x "${BUILD_DIR}/${OUTPUT_NAME}")
        list (APPEND POST_BUILD_COMMAND
          COMMAND "${CMAKE_COMMAND}" -E copy_directory
                  "${BUILD_DIR}/${OUTPUT_NAME_WE}.app"
                  "${RUNTIME_OUTPUT_DIRECTORY}/${OUTPUT_NAME_WE}.app"
          COMMAND "${CMAKE_COMMAND}" -E copy
                  "${BUILD_DIR}/${OUTPUT_NAME}"
                  "${RUNTIME_OUTPUT_DIRECTORY}/${OUTPUT_NAME}"
        )
      else ()
        list (APPEND POST_BUILD_COMMAND
          COMMAND "${CMAKE_COMMAND}" -E copy
                  "${BUILD_DIR}/${OUTPUT_NAME}"
                  "${RUNTIME_OUTPUT_DIRECTORY}/${OUTPUT_NAME}"
        )
      endif ()
    endif ()
    # configure build script
    configure_file ("${BASIS_SCRIPT_EXECUTE_PROCESS}" "${BUILD_SCRIPT}" @ONLY)
    # add custom command to build executable using MATLAB Compiler
    add_custom_command (
      OUTPUT ${BUILD_OUTPUT}
      # rebuild when input sources were modified
      MAIN_DEPENDENCY "${MAIN_SOURCE}"
      DEPENDS         ${DEPENDS}
      # invoke MATLAB Compiler in either MATLAB or standalone mode
      # wrapping command in CMake execute_process() command allows for inspection
      # of command output for error messages and specification of timeout
      COMMAND "${CMAKE_COMMAND}"
              "-DCONFIG=$<${BASIS_GE_CONFIG}>"
              "-DWORKING_DIRECTORY=${WORKING_DIR}"
              "-DTIMEOUT=${BASIS_MCC_TIMEOUT}"
              "-DRETRY_EXPRESSION=License checkout failed"
              "-DRETRY_ATTEMPTS=${BASIS_MCC_RETRY_ATTEMPTS}"
              "-DRETRY_DELAY=${BASIS_MCC_RETRY_DELAY}"
              "-DERROR_EXPRESSION=[E|e]rror:|Illegal output name"
              "-DOUTPUT_FILE=${BUILD_LOG}"
              "-DERROR_FILE=${BUILD_LOG}"
              "-DVERBOSE=OFF"
              "-DLOG_ARGS=ON"
              "-P" "${BUILD_SCRIPT}"
      # post build command(s)
      ${POST_BUILD_COMMAND}
      # inform user where build log can be found
      COMMAND "${CMAKE_COMMAND}" -E echo "Build log written to ${BUILD_LOG}"
      # comment
      COMMENT "${BUILD_COMMENT}"
      VERBATIM
    )
  endif ()
  # --------------------------------------------------------------------------
  # add custom target
  add_custom_target (_${TARGET_UID} DEPENDS ${BUILD_OUTPUT} SOURCES ${SOURCES})
  add_dependencies (${TARGET_UID} _${TARGET_UID})
  # cleanup on "make clean"
  set (ADDITIONAL_MAKE_CLEAN_FILES "${BUILD_OUTPUT}")
  if (COMPILE)
    list (APPEND ADDITIONAL_MAKE_CLEAN_FILES
      "${BUILD_DIR}/${MCC_OUTPUT_NAME_WE}.prj"
      "${BUILD_DIR}/mccExcludedFiles.log"
      "${BUILD_DIR}/mccBuild.log"
      "${BUILD_DIR}/readme.txt"
    )
    if (LIBRARY)
      list (APPEND ADDITIONAL_MAKE_CLEAN_FILES
        "${BUILD_DIR}/${OUTPUT_NAME_WE}.h"
        "${BUILD_DIR}/${MCC_OUTPUT_NAME_WE}.h"
        "${BUILD_DIR}/${MCC_OUTPUT_NAME_WE}.c"
        "${BUILD_DIR}/${MCC_OUTPUT_NAME_WE}.exports"
      )
    else ()
      list (APPEND ADDITIONAL_MAKE_CLEAN_FILES
        "${BUILD_DIR}/${OUTPUT_NAME}"
        "${BUILD_DIR}/${MCC_OUTPUT_NAME}"
        "${BUILD_DIR}/run_${MCC_OUTPUT_NAME_WE}.sh"
        "${BUILD_DIR}/${MCC_OUTPUT_NAME_WE}_main.c"
        "${BUILD_DIR}/${MCC_OUTPUT_NAME_WE}_mcc_component_data.c"
      )
    endif ()
  endif ()
  set_property (DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES "${ADDITIONAL_MAKE_CLEAN_FILES}")
  unset (ADDITIONAL_MAKE_CLEAN_FILES)
  # export target
  if (EXPORT)
    basis_add_custom_export_target (${TARGET_UID} "${IS_TEST}")
  endif ()
  # install executable or library
  if (LIBRARY)
    if (LIBRARY_HEADER_DIRECTORY)
      file (RELATIVE_PATH INCLUDE_DIR_SUFFIX "${BINARY_INCLUDE_DIR}" "${LIBRARY_HEADER_DIRECTORY}")
      if (NOT INCLUDE_DIR_SUFFIX MATCHES "^../")
        string (REGEX REPLACE "/$" "" INCLUDE_DIR_SUFFIX "${INCLUDE_DIR_SUFFIX}")
        if (INCLUDE_DIR_SUFFIX STREQUAL ".")
          set (INCLUDE_DIR_SUFFIX)
        else ()
          set (INCLUDE_DIR_SUFFIX "/${INCLUDE_DIR_SUFFIX}")
        endif ()
        install (
          FILES       "${BUILD_DIR}/${OUTPUT_NAME_WE}.h"
          DESTINATION "${INSTALL_INCLUDE_DIR}${INCLUDE_DIR_SUFFIX}"
          COMPONENT   "${LIBRARY_COMPONENT}"
        )
      endif ()
    endif ()
    if (LIBRARY_INSTALL_DIRECTORY)
      if (WIN32)
        install (
          FILES       "${BUILD_DIR}/${OUTPUT_NAME_WE}.dll"
          DESTINATION "${RUNTIME_INSTALL_DIRECTORY}"
          COMPONENT   "${RUNTIME_COMPONENT}"
        )
      endif ()
      install (
        FILES       "${INSTALL_FILE}"
        DESTINATION "${LIBRARY_INSTALL_DIRECTORY}"
        COMPONENT   "${LIBRARY_COMPONENT}"
      )
    endif ()
  else ()
    if (RUNTIME_INSTALL_DIRECTORY)
      if (CMAKE_SYSTEM_NAME STREQUAL "Darwin")
        install (
          DIRECTORY   "${BUILD_DIR}/${OUTPUT_NAME_WE}.app"
          DESTINATION "${RUNTIME_INSTALL_DIRECTORY}"
          COMPONENT   "${RUNTIME_COMPONENT}"
          USE_SOURCE_PERMISSIONS
        )
      endif ()
      install (
        PROGRAMS    "${INSTALL_FILE}"
        DESTINATION "${RUNTIME_INSTALL_DIRECTORY}"
        COMPONENT   "${RUNTIME_COMPONENT}"
      )
    endif ()
  endif ()
  # done
  if (BASIS_VERBOSE)
    message (STATUS "Adding build command for target ${TARGET_UID}... - done")
  endif ()
endfunction ()
