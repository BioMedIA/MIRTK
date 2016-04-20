#.rst
# FindPython
# ----------
#
# Find Python components such as the interpreter, C libraries, or Python modules.
#
# This module uses find_package to look for these components using the
# specialized Find modules named FindPythonInterp.cmake, FindPythonLibs.cmake,
# and FindPythonModules.cmake. A Python_FIND_COMPONENTS list entry that does
# not match either "Interp" or "Libs" is assumed to be the name of a Python
# module to be searched. See the aforementioned Find modules for more details.

#=============================================================================
# Copyright 2016 Andreas Schuh <andreas.schuh.84@gmail.com>
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distribute this file outside of CMake, substitute the full
#  License text for the above reference.)

# ------------------------------------------------------------------------------
# Components to look for
set (_Python_FIND_Interp  FALSE)
set (_Python_FIND_Libs    FALSE)
set (_Python_FIND_Modules FALSE)
set (_Python_FIND_Modules_OPTIONAL)
set (_Python_FIND_Modules_REQUIRED)

if (NOT Python_FIND_COMPONENTS)
  set (Python_FIND_COMPONENTS Interp)
endif ()

foreach (_Python_COMPONENT IN LISTS Python_FIND_COMPONENTS)
  if (_Python_COMPONENT MATCHES "^(Interp|Libs)$")
    set (_Python_FIND_${_Python_COMPONENT} TRUE)
  else ()
    set (_Python_FIND_Modules TRUE)
    if (Python_FIND_REQUIRED_${_Python_COMPONENT})
      list (APPEND _Python_FIND_Modules_REQUIRED ${_Python_COMPONENT})
    else ()
      list (APPEND _Python_FIND_Modules_OPTIONAL ${_Python_COMPONENT})
    endif ()
  endif ()
endforeach ()

# ------------------------------------------------------------------------------
# Verbose message
if (NOT Python_FIND_QUIETLY)
  set(_Python_FIND_STATUS "Looking for Python [${Python_FIND_COMPONENTS}]")
  if (NOT Python_FIND_REQUIRED)
    set(_Python_FIND_STATUS "${_Python_FIND_STATUS} (optional)")
  endif ()
  message(STATUS "${_Python_FIND_STATUS}...")
endif ()

# ------------------------------------------------------------------------------
# Look for Python components
set (_Python_REQUIRED_VARS)

if (_Python_FIND_Interp)
  find_package (PythonInterp QUIET MODULE)
  set (Python_Interp_FOUND ${PythonInterp_FOUND})
  list (APPEND _Python_REQUIRED_VARS PYTHON_EXECUTABLE)
endif ()

if (Python_FIND_Libs)
  find_package (PythonLibs QUIET MODULE)
  set (Python_Libs_FOUND ${PythonLibs_FOUND})
  list (APPEND _Python_REQUIRED_VARS PYTHON_LIBRARIES PYTHON_INCLUDE_DIRS)
endif ()

if (Python_FIND_Modules)
  find_package (PythonModules QUIET MODULE
    COMPONENTS ${Python_FIND_Modules_REQUIRED}
    OPTIONAL_COMPONENTS ${Python_FIND_Modules_OPTIONAL}
  )
  foreach (_Python_COMPONENT IN LISTS _Python_FIND_Modules_REQUIRED _Python_FIND_Modules_OPTIONAL)
    set (Python_${_Python_COMPONENT}_FOUND ${PythonModules_${_Python_COMPONENT}_FOUND})
  endforeach ()
endif ()

# ------------------------------------------------------------------------------
# Handle QUIET, REQUIRED, and [EXACT] VERSION arguments and set Python_FOUND
if (PYTHON_VERSION_STRING)
  set (_Python_VERSION ${PYTHON_VERSION_STRING})
else ()
  set (_Python_VERSION ${PYTHONLIBS_VERSION_STRING})
endif ()

if (_Python_REQUIRED_VARS)
  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(Python
    REQUIRED_VARS ${_Python_REQUIRED_VARS}
    VERSION_VAR   _Python_VERSION
    HANDLE_COMPONENTS
  )
else ()
  set (Python_FOUND 1)
  set (PYTHON_FOUND 1)
endif ()

# ------------------------------------------------------------------------------
# Verbose message
if (NOT Python_FIND_QUIETLY)
  if (Python_FOUND)
    if (_Python_VERSION)
      message(STATUS "${_Python_FIND_STATUS}... - found v${_Python_VERSION}")
    else ()
      message(STATUS "${_Python_FIND_STATUS}... - found")
    endif ()
  else ()
    message(STATUS "${_Python_FIND_STATUS}... - not found")
  endif ()
endif ()

# ------------------------------------------------------------------------------
# Unset local variables
unset (_Python_VERSION)
unset (_Python_FIND_STATUS)
unset (_Python_FIND_Interp)
unset (_Python_FIND_Libs)
unset (_Python_FIND_Modules)
unset (_Python_FIND_Modules_REQUIRED)
unset (_Python_FIND_Modules_OPTIONAL)
unset (_Python_COMPONENT)
unset (_Python_REQUIRED_VARS)
