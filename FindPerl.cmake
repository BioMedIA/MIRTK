##############################################################################
# @file  FindPerl.cmake
# @brief Find Perl interpreter.
#
# @par Output variables:
# <table border="0">
#   <tr>
#     @tp @b Perl_FOUND @endtp
#     <td>Was the Python executable found.</td>
#   </tr>
#   <tr>
#     @tp @b PERL_FOUND @endtp
#     <td>Alias for @b Perl_FOUND for backwards compatibility.</td>
#   </tr>
#   <tr>
#     @tp @b PERL_EXECUTABLE @endtp
#     <td>Path to the Perl interpreter.</td>
#   </tr>
#   <tr>
#     @tp @b PERL_VERSION_STRING @endtp
#     <td>Perl version found e.g. 5.12.4.</td>
#   </tr>
#   <tr>
#     @tp @b PERL_VERSION_MAJOR @endtp
#     <td>Perl major version found e.g. 5.</td>
#   </tr>
#   <tr>
#     @tp @b PERL_VERSION_MINOR @endtp
#     <td>Perl minor version found e.g. 12.</td>
#   </tr>
#   <tr>
#     @tp @b PERL_VERSION_PATCH @endtp
#     <td>Perl patch version found e.g. 4.</td>
#   </tr>
# </table>
#
# @note This module has been copied from CMake 2.8.5 and modified to also
#       obtain the version information of the found Perl interpreter.
#
# @ingroup CMakeFindModules
##############################################################################

#=============================================================================
# Copyright 2001-2009 Kitware, Inc.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.
#
# * Neither the names of Kitware, Inc., the Insight Software Consortium,
#   nor the names of their contributors may be used to endorse or promote
#   products derived from this software without specific prior written
#   permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#=============================================================================

include (FindCygwin)

set (PERL_POSSIBLE_BIN_PATHS "${CYGWIN_INSTALL_PATH}/bin")

if (WIN32)
  get_filename_component (
    ActivePerl_CurrentVersion 
      "[HKEY_LOCAL_MACHINE\\SOFTWARE\\ActiveState\\ActivePerl;CurrentVersion]" 
    NAME
  )
  set (PERL_POSSIBLE_BIN_PATHS ${PERL_POSSIBLE_BIN_PATHS}
    "C:/Perl/bin" 
    "[HKEY_LOCAL_MACHINE\\SOFTWARE\\ActiveState\\ActivePerl\\${ActivePerl_CurrentVersion}]/bin"
  )
endif ()

find_program (PERL_EXECUTABLE NAMES perl PATHS ${PERL_POSSIBLE_BIN_PATHS})

if (PERL_EXECUTABLE)
  execute_process (COMMAND "${PERL_EXECUTABLE}" --version OUTPUT_VARIABLE _Perl_STDOUT ERROR_VARIABLE _Perl_STDERR)
  if (_Perl_STDOUT MATCHES "[( ]v([0-9]+)\\.([0-9]+)\\.([0-9]+)[ )]")
    set (PERL_VERSION_MAJOR "${CMAKE_MATCH_1}")
    set (PERL_VERSION_MINOR "${CMAKE_MATCH_2}")
    set (PERL_VERSION_PATCH "${CMAKE_MATCH_3}")
    set (PERL_VERSION_STRING "${PERL_VERSION_MAJOR}.${PERL_VERSION_MINOR}.${PERL_VERSION_PATCH}")
  else ()
    message (WARNING "Failed to determine version of Perl interpreter (${PERL_EXECUTABLE})! Error:\n${_Perl_STDERR}")
  endif ()
  unset (_Perl_STDOUT)
  unset (_Perl_STDERR)
endif ()

include (FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS (Perl DEFAULT_MSG PERL_EXECUTABLE)

mark_as_advanced (PERL_EXECUTABLE)
