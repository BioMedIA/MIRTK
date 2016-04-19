# ============================================================================
# Copyright (c) 2011-2012 University of Pennsylvania
# Copyright (c) 2013-2016 Andreas Schuh
# All rights reserved.
#
# See COPYING file for license information or visit
# https://cmake-basis.github.io/download.html#license
# ============================================================================

##############################################################################
# @file  SlicerTools.cmake
# @brief Definition of tools for Slicer Extensions.
##############################################################################

# ----------------------------------------------------------------------------
# include guard
if (__BASIS_SLICERTOOLS_INCLUDED)
  return ()
else ()
  set (__BASIS_SLICERTOOLS_INCLUDED TRUE)
endif ()

# ============================================================================
# meta-data
# ============================================================================

# ----------------------------------------------------------------------------
## @brief Define project meta-data of Slicer module.
#
# This macro should be used instead of basis_project() for a Slicer module.
# It extends the considered meta-data by some additional variables that
# have to be set for a Slicer module and identifies this project as a Slicer
# module.
#
# @ingroup CMakeAPI
macro (basis_slicer_module)
  CMAKE_PARSE_ARGUMENTS (
    PROJECT
      "${BASIS_METADATA_LIST_SWITCH}"
      "${BASIS_METADATA_LIST_SINGLE};${BASIS_SLICER_METADATA_LIST_SINGLE}"
      "${BASIS_METADATA_LIST_MULTI};${BASIS_SLICER_METADATA_LIST_MULTI}"
    ${ARGN}
  )
  foreach (_L IN LISTS BASIS_SLICER_METADATA_LIST_MULTI)
    if (_L MATCHES "^(CATEGORY|CONTRIBUTORS)$")
      basis_list_to_delimited_string (PROJECT_${_L} ", " ${PROJECT_${_L}})
    else ()
      basis_list_to_string (PROJECT_${_L} ${PROJECT_${_L}})
    endif ()
  endforeach ()
  set (PROJECT_IS_SLICER_MODULE TRUE)
  basis_project_check_metadata ()
endmacro ()


## @addtogroup CMakeUtilities
# @{


# ============================================================================
# initialization
# ============================================================================

# ----------------------------------------------------------------------------
## @brief Initialize slicer module meta-data.
#
# At the moment, only one module, most often the top-level project can be a
# Slicer module, i.e., call basis_slicer_module() either in the
# BasisProject.cmake file of the top-level project or in the corresponding file
# of at most one of the project modules.
macro (basis_slicer_module_initialize)
  if (NOT PROJECT_IS_MODULE)
    set (_N 0)
    foreach (_M IN LISTS PROJECT_MODULES_ENABLED)
      if (${_M}_IS_SLICER_MODULE)
        math (EXPR _N "${_N} + 1")
      endif ()
    endforeach ()
    if (PROJECT_IS_SLICER_MODULE)
      if (_N GREATER 0)
        message (FATAL_ERROR "BASIS does not support projects with multiple Slicer modules!")
      endif ()
      set (_FIND_SLICER TRUE)
    elseif (_N GREATER 0)
      if (_N GREATER 1)
        message (FATAL_ERROR "BASIS does not support projects with multiple Slicer modules!")
      endif ()
      set (_FIND_SLICER TRUE)
    else ()
      set (_FIND_SLICER FALSE)
    endif ()
    # convert PROJECT_* or <Module>_* to MODULE_* variables
    if (PROJECT_IS_SLICER_MODULE)
      foreach (_D IN LISTS BASIS_SLICER_METADATA_LIST)
        if (DEFINED PROJECT_${_D})
          set (MODULE_${_D} "${PROJECT_${_D}}")
        endif ()
      endforeach ()
    else ()
      foreach (_M IN LISTS PROJECT_MODULES_ENABLED)
        if (${_M}_IS_SLICER_MODULE)
          foreach (_D IN LISTS BASIS_SLICER_METADATA_LIST)
            if (DEFINED ${_M}_${_D})
              set (MODULE_${_D} "${${_M}_${_D}}")
            endif ()
          endforeach ()
          break ()
        endif ()
      endforeach ()
    endif ()
    # find Slicer package
    if (_FIND_SLICER)
      # set metadata
      set (MODULE_DESCRIPTION   "${PROJECT_DESCRIPTION}")
      set (MODULE_NAME          "${PROJECT_NAME}")
      set (MODULE_README_FILE   "${PROJECT_SOURCE_DIR}/README.txt")
      set (MODULE_LICENSE_FILE  "${PROJECT_SOURCE_DIR}/COPYING.txt")
      set (MODULE_MAJOR_VERSION "${PROJECT_VERSION_MAJOR}")
      set (MODULE_MINOR_VERSION "${PROJECT_VERSION_MINOR}")
      set (MODULE_PATCH_VERSION "${PROJECT_VERSION_PATCH}")
      # skip project() command
      set (Slicer_SKIP_PROJECT_COMMAND        TRUE)
      set (Slicer_SKIP_USE_FILE_INCLUDE_CHECK TRUE)
      # work-around for bug 1901 in Slicer 4.1.0
      set (Slicer_SKIP_SlicerBlockModuleToExtensionMetadata TRUE)
      basis_slicer_module_to_extension_metadata ()
      # find Slicer
      basis_find_package (Slicer REQUIRED)
      basis_use_package  (Slicer REQUIRED)
    endif ()
    unset (_M)
    unset (_N)
    unset (_FIND_SLICER)
  endif ()
endmacro ()

# ============================================================================
# auxiliary functions
# ============================================================================

# ----------------------------------------------------------------------------
## @brief Copy MODULE_* metadata to EXTENSION_* metadata of Slicer Extension.
#
# This function is required to work-around bug 1901 in Slicer 4.1.0. It basically
# implements what the SlicerBlockModuleToExtensionMetadata.cmake module which
# can be found in the Extensions/CMake/ directory of the Slicer 4 source tree
# is implementing. The list of metadata has been copied from this particular
# CMake module of the 4.1.0 release of Slicer.
#
# @sa http://www.na-mic.org/Bug/view.php?id=1901
function (basis_slicer_module_to_extension_metadata)
  foreach (varname IN ITEMS ${BASIS_SLICER_METADATA_LIST} DESCRIPTION)
    if (DEFINED MODULE_${varname} AND NOT DEFINED EXTENSION_${varname})
      set (EXTENSION_${varname} ${MODULE_${varname}} PARENT_SCOPE)
    endif ()
  endforeach ()
endfunction ()


## @}
# Doxygen group CMakeUtilities
