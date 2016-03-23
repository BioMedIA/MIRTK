# ==============================================================================
# Medical Image Registration ToolKit (MIRTK)
#
# Copyright 2016 Imperial College London
# Copyright 2016 Andreas Schuh
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# ==============================================================================

if (COMMAND mirtk_project_begin)
  return()
endif ()

# ------------------------------------------------------------------------------
## Initialize MIRTK project (see basis_project_begin)
macro(mirtk_project_begin)
  # Override default BASIS settings and hide unused cache entries
  if (NOT PROJECT_IS_MODULE)
    set(BASIS_NAMESPACE_DELIMITER_CMAKE ::)
    set(BUILD_MODULES_BY_DEFAULT    OFF       CACHE INTERNAL "" FORCE)
    set(BASIS_BUILD_ONLY            OFF       CACHE INTERNAL "" FORCE)
    set(BASIS_ALL_DOC               OFF       CACHE INTERNAL "" FORCE)
    set(BASIS_COMPILE_MATLAB        ON        CACHE INTERNAL "" FORCE)
    set(BASIS_COMPILE_SCRIPTS       OFF       CACHE INTERNAL "" FORCE)
    set(BASIS_INSTALL_RPATH         ON        CACHE INTERNAL "" FORCE)
    set(BASIS_INSTALL_SCHEME        "default" CACHE INTERNAL "" FORCE)
    set(BASIS_INSTALL_SITE_PACKAGES OFF       CACHE INTERNAL "" FORCE)
    set(BASIS_REGISTER              ON        CACHE INTERNAL "" FORCE)
    set(BASIS_REVISION_INFO         ON        CACHE INTERNAL "" FORCE)
    set(BASIS_SUPERBUILD_MODULES    OFF       CACHE INTERNAL "" FORCE)
    foreach (lang IN ITEMS BASH CXX PYTHON PERL)
      set(BUILD_BASIS_UTILITIES_FOR_${lang} OFF CACHE INTERNAL "" FORCE)
    endforeach ()
    unset(lang)
  endif ()
  # mirtk_add_executable sets the INSTALL_RPATH itself
  set(BASIS_INSTALL_RPATH OFF)
  # Export of build targets
  if (PROJECT_IS_MODULE)
    set(BASIS_EXPORT_ENABLED OFF)
  else ()
    set(BASIS_EXPORT_ENABLED ON)
  endif ()
  set(BASIS_EXPORT_SUFFIX "Targets")
  # Use CMAKE_CONFIGURATION_TYPES as selection of CMAKE_BUILD_TYPE and hide it
  basis_is_cached(_CACHED_CONFIGURATION_TYPES CMAKE_CONFIGURATION_TYPES)
  if (_CACHED_CONFIGURATION_TYPES)
    basis_is_cached(_CACHED_BUILD_TYPE CMAKE_BUILD_TYPE)
    if (_CACHED_BUILD_TYPE AND CMAKE_CONFIGURATION_TYPES)
      set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS ${CMAKE_CONFIGURATION_TYPES})
    endif ()
    unset(_CACHED_BUILD_TYPE)
    set_property(CACHE CMAKE_CONFIGURATION_TYPES PROPERTY TYPE INTERNAL)
  endif ()
  unset(_CACHED_CONFIGURATION_TYPES)
  # CodeBlocks generators add this cache entry, hide it when set
  basis_is_cached(_CACHED_CODEBLOCKS_EXECUTABLE CMAKE_CODEBLOCKS_EXECUTABLE)
  if (_CACHED_CODEBLOCKS_EXECUTABLE)
    if (CMAKE_CODEBLOCKS_EXECUTABLE)
      mark_as_advanced(FORCE CMAKE_CODEBLOCKS_EXECUTABLE)
    else ()
      mark_as_advanced(CLEAR CMAKE_CODEBLOCKS_EXECUTABLE)
    endif ()
  endif ()
  # Choose between shared or static linkage (shared is recommended)
  if (NOT BUILD_SHARED_LIBS_DEFAULT)
    set (BUILD_SHARED_LIBS_DEFAULT ON)
  endif ()
  option(BUILD_SHARED_LIBS "Request build of shared libraries" ${BUILD_SHARED_LIBS_DEFAULT})
  mark_as_advanced(BUILD_SHARED_LIBS)
  # Enable profiling of program execution (cf. Common/include/mirtkProfiling.h)
  if (NOT WITH_PROFILING_DEFAULT)
    set (WITH_PROFILING_DEFAULT ON)
  endif ()
  option(WITH_PROFILING "Enable profiling of program execution" ${WITH_PROFILING_DEFAULT})
  if (WITH_PROFILING)
    add_definitions(-DMIRTK_WITH_PROFILING)
  endif ()
  # Start BASIS project
  basis_project_begin()
endmacro ()
