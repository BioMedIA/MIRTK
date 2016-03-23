# ============================================================================
# CMake - Cross Platform Makefile Generator
# Copyright 2000-2009 Kitware, Inc., Insight Software Consortium
# All rights reserved.
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
#
# ----------------------------------------------------------------------------
#
# The above copyright and license notice applies to distributions of
# CMake in source and binary form.  Some source files contain additional
# notices of original copyright by their contributors; see each source
# for details.  Third-party software packages supplied with CMake under
# compatible licenses provide their own copyright notices documented in
# corresponding subdirectories.
#
# ----------------------------------------------------------------------------
#
# CMake was initially developed by Kitware with the following sponsorship:
#
# * National Library of Medicine at the National Institutes of Health
#   as part of the Insight Segmentation and Registration Toolkit (ITK).
#
# * US National Labs (Los Alamos, Livermore, Sandia) ASC Parallel
#   Visualization Initiative.
#
# * National Alliance for Medical Image Computing (NAMIC) is funded by the
#   National Institutes of Health through the NIH Roadmap for Medical Research,
#   Grant U54 EB005149.
#
# * Kitware, Inc.
# ============================================================================

##############################################################################
# @file  FindITK.cmake
# @brief Find an ITK installation or build tree.
#
# When ITK is found, the ITKConfig.cmake file is sourced to setup the
# location and configuration of ITK.  Please read this file, or
# ITKConfig.cmake.in from the ITK source tree for the full list of
# definitions.  Of particular interest is ITK_USE_FILE, a CMake source file
# that can be included to set the include directories, library directories,
# and preprocessor macros.  In addition to the variables read from
# ITKConfig.cmake, this find module also defines
#
# @par Output variables:
# <table border="0">
#   <tr>
#     @tp @b ITK_DIR @endtp
#     <td>The directory containing ITKConfig.cmake.  
#         This is either the root of the build tree, 
#         or the lib/InsightToolkit directory.  
#         This is the only cache entry.</td>
#   </tr>
#   <tr>
#     @tp @b ITK_FOUND @endtp
#     <td>Whether ITK was found.  If this is true, @c ITK_DIR is okay.</td>
#   </tr>
#   <tr>
#     @tp @b USE_ITK_FILE @endtp
#     <td>The full path to the <tt>UseITK.cmake</tt> file.  
#         This is provided for backward 
#         compatability. Use @c ITK_USE_FILE instead.</td>
#   </tr>
# </table>
#
# @ingroup CMakeFindModules
##############################################################################

# Use the Config mode of the find_package() command to find ITKConfig.
# If this succeeds (possibly because ITK_DIR is already set), the
# command will have already loaded ITKConfig.cmake and set ITK_FOUND.
IF(NOT ITK_FOUND)
  SET(_ITK_REQUIRED "")
  IF(ITK_FIND_REQUIRED)
    SET(_ITK_REQUIRED REQUIRED)
  ENDIF()
  SET(_ITK_QUIET "")
  IF(ITK_FIND_QUIETLY)
    SET(_ITK_QUIET QUIET)
  ENDIF()
  FIND_PACKAGE(ITK ${ITK_FIND_VERSION} ${_ITK_REQUIRED} ${_ITK_QUIET} NO_MODULE
    NAMES ITK InsightToolkit
    CONFIGS ITKConfig.cmake
    )
ENDIF()

IF(ITK_FOUND)
  # Set USE_ITK_FILE for backward-compatability.
  SET(USE_ITK_FILE ${ITK_USE_FILE})
ENDIF()
