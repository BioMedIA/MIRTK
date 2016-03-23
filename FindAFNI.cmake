# ============================================================================
# Copyright (c) 2011-2012 University of Pennsylvania
# Copyright (c) 2013-2014 Andreas Schuh
# All rights reserved.
#
# See COPYING file for license information or visit
# http://opensource.andreasschuh.com/cmake-basis/download.html#license
# ============================================================================

##############################################################################
# @file  FindAFNI.cmake
# @brief Find programs of the AFNI package.
#
# @sa http://afni.nimh.nih.gov/afni/
#
# @par Input variables:
# <table border="0">
#   <tr>
#     @tp @b AFNI_DIR @endtp
#     <td>The AFNI package files are searched under the specified root
#         directory. If they are not found there, the default search paths
#         are considered. This variable can also be set as environment variable.</td>
#   </tr>
#   <tr>
#     @tp @b AFNI_FIND_COMPONENTS @endtp
#     <td>List of components, i.e., AFNI programs, to look for.
#         Specify using COMPONENTS argument of find_package() command.</td>
#   </tr>
#   <tr>
#     @tp @b AFNI_FIND_OPTIONAL_COMPONENTS @endtp
#     <td>List of optional components, i.e., AFNI programs, to look for.
#         Specify using OPTIONAL_COMPONENTS argument of find_package() command.</td>
#   </tr>
# </table>
#
# @par Output variables:
# <table border="0">
#   <tr>
#     @tp @b AFNI_FOUND @endtp
#     <td>Whether all required programs of the AFNI package were found. If only
#         optional programs were searched, this variable is set to @c TRUE if
#         all named programs were found.</td>
#   </tr>
#   <tr>
#     @tp @b AFNI_&lt;tool&gt;_EXECUTABLE @endtp
#     <td>Absolute path of the corresponding found AFNI program, e.g., @c AFNI_3dcalc_EXECUTABLE.</td>
#   </tr>
# </table>
#
# @ingroup CMakeFindModules
##############################################################################

# ----------------------------------------------------------------------------
# initialize search
if (NOT AFNI_DIR)
  set (AFNI_DIR "$ENV{AFNI_DIR}" CACHE PATH "Installation prefix of AFNI." FORCE)
endif ()

if (NOT AFNI_FIND_COMPONENTS AND NOT AFNI_FIND_OPTIONAL_COMPONENTS)
  message (FATAL_ERROR "The FindAFNI.cmake module requires a list of AFNI programs to look for"
                       " specified using the COMPONENTS and/or OPTIONAL_COMPONENTS argument"
                       " of the find_package() command, e.g.:"
                       "\n"
                       "find_package(AFNI COMPONENTS 3dcalc)")
endif ()

# ----------------------------------------------------------------------------
# private helper macro to look for a particular required or optional AFNI program
macro (_AFNI_find_program NAME REQUIRED)
  if (AFNI_DIR)
    find_program (AFNI_${NAME}_EXECUTABLE NAMES ${NAME} HINTS ${AFNI_DIR} PATH_SUFFIXES bin NO_DEFAULT_PATH)
  else ()
    find_program (AFNI_${NAME}_EXECUTABLE NAMES ${NAME})
  endif ()
  if (NOT AFNI_${NAME}_EXECUTABLE)
    if (AFNI_FIND_COMPONENTS)
      # looking either only for required components or for both optional and
      # and required components; in this case, let AFNI_FOUND reflect only
      # whether all required components were found, but ignore the optional ones;
      # the caller can still check AFNI_<tool>_EXECUTABLE explicitly for these
      # optional components to see whether or not a particular AFNI programs was found
      if (REQUIRED)
        set (AFNI_FOUND FALSE)
      endif ()
    else ()
      # looking only for optional components anyway, so let AFNI_FOUND reflect
      # if all of these optional components were found instead
      set (AFNI_FOUND FALSE)
    endif ()
    if (REQUIRED)
      list (APPEND _AFNI_MISSING_COMPONENTS ${NAME})
    else ()
      list (APPEND _AFNI_MISSING_OPTIONAL_COMPONENTS ${NAME})
    endif ()
  endif ()
  mark_as_advanced (AFNI_${NAME}_EXECUTABLE)
endmacro ()

# ----------------------------------------------------------------------------
# find AFNI program(s)
set (AFNI_FOUND TRUE)

set (_AFNI_MISSING_COMPONENTS)
foreach (_AFNI_TOOL IN LISTS AFNI_FIND_COMPONENTS)
  _AFNI_find_program (${_AFNI_TOOL} TRUE)
endforeach ()

set (_AFNI_MISSING_OPTIONAL_COMPONENTS)
foreach (_AFNI_TOOL IN LISTS AFNI_FIND_OPTIONAL_COMPONENTS)
  _AFNI_find_program (${_AFNI_TOOL} FALSE)
endforeach ()

# ----------------------------------------------------------------------------
# handle the QUIETLY and REQUIRED arguments
set (_AFNI_HELP_MESSAGE "Please set AFNI_DIR to the directory containing these executables or specify the location of each not found executable using the advanced AFNI_<tool>_EXECUTABLE CMake variables.")

if (_AFNI_MISSING_COMPONENTS)
  message (FATAL_ERROR "Could NOT find the following required AFNI program(s):\n${_AFNI_MISSING_COMPONENTS}\n${_AFNI_HELP_MESSAGE}")
elseif (_AFNI_MISSING_OPTIONAL_COMPONENTS AND NOT AFNI_FIND_QUIETLY)
  message (WARNING "Could NOT find the following optional AFNI program(s):\n${_AFNI_MISSING_OPTIONAL_COMPONENTS}\n${_AFNI_HELP_MESSAGE}")
endif ()

# ----------------------------------------------------------------------------
# set AFNI_DIR
if (NOT AFNI_DIR)
  foreach (_AFNI_TOOL IN LISTS AFNI_FIND_COMPONENTS AFNI_FIND_OPTIONAL_COMPONENTS)
    if (AFNI_${_AFNI_TOOL}_EXECUTABLE)
      get_filename_component (AFNI_DIR "${AFNI_${_AFNI_TOOL}_EXECUTABLE}" PATH)
      string (REGEX REPLACE "/bin/?" "" AFNI_DIR "${AFNI_DIR}")
      set (AFNI_DIR "${AFNI_DIR}" CACHE PATH "Installation prefix of AFNI." FORCE)
      break ()
    endif ()
  endforeach ()
endif ()


unset (_AFNI_TOOL)
unset (_AFNI_HELP_MESSAGE)
unset (_AFNI_MISSING_COMPONENTS)
unset (_AFNI_MISSING_OPTIONAL_COMPONENTS)
