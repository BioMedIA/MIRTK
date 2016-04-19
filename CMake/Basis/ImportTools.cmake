# ============================================================================
# Copyright (c) 2011-2012 University of Pennsylvania
# Copyright (c) 2013-2016 Andreas Schuh
# All rights reserved.
#
# See COPYING file for license information or visit
# https://cmake-basis.github.io/download.html#license
# ============================================================================

##############################################################################
# @file  ImportTools.cmake
# @brief Functions and macros for the import of targets.
#
# @ingroup CMakeTools
##############################################################################

if (__BASIS_IMPORTTOOLS_INCLUDED)
  return ()
else ()
  set (__BASIS_IMPORTTOOLS_INCLUDED TRUE)
endif ()


## @addtogroup CMakeUtilities
#  @{

# ----------------------------------------------------------------------------
## @brief Set target property.
#
# This function is overwritten by BASIS in order to update the information
# about imported build targets.
#
# @note Do not use this function in your CMakeLists.txt configuration files.
#       Use basis_set_target_properties() instead.
#
# @note Due to a bug in CMake (http://www.cmake.org/Bug/view.php?id=12303),
#       except of the first property given directly after the @c PROPERTIES keyword,
#       only properties listed in @c BASIS_PROPERTIES_ON_TARGETS can be set.
#
# @param [in] ARGN List of arguments for
#                  <a href="http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:set_target_properties">
#                  set_target_properties()</a>.
#
# @sa http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:set_target_properties
function (set_target_properties)
  # target names
  list (FIND ARGN "PROPERTIES" IDX)
  if (IDX EQUAL -1)
    message (FATAL_ERROR "Missing PROPERTIES argument!")
  elseif (IDX EQUAL 0)
    message (FATAL_ERROR "No targets specified!")
  endif ()
  set (INDICES)
  set (I 0)
  while (I LESS IDX)
    list (APPEND INDICES ${I})
    math (EXPR I "${I} + 1")
  endwhile ()
  list (GET ARGN ${INDICES} TARGETS)
  # remaining arguments are property value pairs
  list (REMOVE_AT ARGN ${INDICES} ${IDX})
  # set target properties
  #
  # Note: By iterating over the properties, the empty property values
  #       are correctly passed on to CMake's set_target_properties()
  #       command, while
  #       _set_target_properties(${TARGET_UIDS} PROPERTIES ${ARGN})
  #       (erroneously) discards the empty elements in ARGN.
  if (BASIS_DEBUG)
    message ("** set_target_properties:")
    message ("**   Target(s):  ${TARGETS}")
    message ("**   Properties: [${ARGN}]")
  endif ()
  list (LENGTH ARGN N)
  while (N GREATER 1)
    list (GET ARGN 0 PROPERTY)
    list (GET ARGN 1 VALUE)
    list (REMOVE_AT ARGN 0 1)
    list (LENGTH ARGN N)
    # The following loop is only required b/c CMake's ARGV and ARGN
    # lists do not support arguments which are themselves lists.
    # Therefore, we need a way to decide when the list of values for a
    # property is terminated. Hence, we only allow known properties
    # to be set, except for the first property where the name follows
    # directly after the PROPERTIES keyword.
    while (N GREATER 0)
      list (GET ARGN 0 ARG)
      if (ARG MATCHES "${BASIS_PROPERTIES_ON_TARGETS_RE}")
        break ()
      endif ()
      list (APPEND VALUE "${ARG}")
      list (REMOVE_AT ARGN 0)
      list (LENGTH ARGN N)
    endwhile ()
    if (BASIS_DEBUG)
      message ("**   -> ${PROPERTY} = [${VALUE}]")
    endif ()
    # check property name
    if (PROPERTY MATCHES "^$")
      message (FATAL_ERROR "Empty property name given!")
    # if property is related to the location of an imported target,
    # update corresponding project properties
    elseif (PROPERTY MATCHES "^IMPORTED_LOCATION")
      list (GET TARGETS 0 TARGET)
      basis_update_imported_location (${TARGET} ${PROPERTY} "${VALUE}")
    # if property is related to the type of an imported target,
    # update corresponding project properties
    elseif (PROPERTY MATCHES "^BASIS_TYPE$")
      list (GET TARGETS 0 TARGET)
      basis_update_imported_type (${TARGET} "${VALUE}")
    endif ()
    # set target property
    _set_target_properties (${TARGETS} PROPERTIES ${PROPERTY} "${VALUE}")
  endwhile ()
  # make sure that every property had a corresponding value
  if (NOT N EQUAL 0)
    message (FATAL_ERROR "No value given for target property ${ARGN}")
  endif ()
endfunction ()

# ----------------------------------------------------------------------------
## @brief Add imported target.
#
# Imported targets are only valid in the scope where they were imported.
# In order to be able to add the information of the imported executable targets
# to the ExecutableTargetInfo modules of the BASIS utilities which are configured
# during the finalization of the (top-level) project, the information of
# imported targets has to be stored in the global scope. Therefore, internal
# cache variables prefixed by the name of the project are used
# (see basis_set_project_property()):
#
# <table border="0">
#   <tr>
#     @tp @b IMPORTED_TARGETS @endtp
#     <td>List of imported targets.</td>
#   </tr>
#   <tr>
#     @tp @b IMPORTED_TYPES @endtp
#     <td>Types of imported targets.</td>
#   </tr>
#   <tr>
#     @tp @b IMPORTED_LOCATIONS @endtp
#     <td>Locations of imported target files.</td>
#   </tr>
#   <tr>
#     @tp @b IMPORTED_RANKS @endtp
#     <td>Rank of current imported locations. This rank value is used to decide
#         whether the current location takes precedence over another imported
#         location. For example, IMPORTED_LOCATION_&lt;a&gt;, may be preferred
#         over IMPORTED_LOCATION_&lt;b&gt;.
#   </tr>
# </table>
#
# @param [in] TARGET Name (UID) of the imported target.
# @param [in] TYPE   Type of the imported target.
#
# @sa basis_update_imported_location()
function (basis_add_imported_target TARGET TYPE)
  if (BASIS_DEBUG AND BASIS_VERBOSE)
    message ("** basis_add_imported_target:")
    message ("**   Target:   ${TARGET}")
    message ("**   Type:     ${TYPE}")
    message ("**   Bundled:  ${BUNDLE_PROJECT}")
  endif ()
  if (NOT TYPE STREQUAL "INTERFACE")
    if (BUNDLE_PROJECT)
      _set_target_properties (${TARGET} PROPERTIES BUNDLED TRUE)
    else ()
      _set_target_properties (${TARGET} PROPERTIES BUNDLED FALSE)
    endif ()
  endif ()
  # if target was added before
  basis_get_project_property (TARGETS PROPERTY IMPORTED_TARGETS)
  if (TARGETS)
    list (FIND TARGETS "${TARGET}" IDX)
    if (NOT IDX EQUAL -1)
      # do nothing
      return ()
    endif ()
  endif ()
  # otherwise, add it to the project properties
  basis_set_project_property (APPEND PROPERTY IMPORTED_TARGETS   "${TARGET}")
  basis_set_project_property (APPEND PROPERTY IMPORTED_TYPES     "${TYPE}")
  basis_set_project_property (APPEND PROPERTY IMPORTED_LOCATIONS "NOTFOUND")
  basis_set_project_property (APPEND PROPERTY IMPORTED_RANKS     10)
endfunction ()

# ----------------------------------------------------------------------------
## @brief Update location of imported target.
#
# @param [in] TARGET     Name (UID) of the imported target.
# @param [in] PROPERTY   Target location property. Either IMPORTED_LOCATION
#                        or IMPORTED_LOCATION_&lt;config&gt;, where &lt;config&gt;
#                        is one of the imported build configurations.
#                        This argument is used to decide whether to keep
#                        the current target information or to replace it
#                        by the new one.
# @param [in] LOCATION   Location of imported target.
function (basis_update_imported_location TARGET PROPERTY LOCATION)
  if (BASIS_DEBUG)
    message ("** basis_update_imported_location:")
    message ("**   Target:   ${TARGET}")
    message ("**   Location: ${LOCATION}")
  endif ()
  # get index of imported target
  basis_get_project_property (TARGETS PROPERTY IMPORTED_TARGETS)
  list (FIND TARGETS "${TARGET}" IDX)
  if (IDX EQUAL -1)
    # imported targets have to be added via basis_add_imported_target() first
    # otherwise, ignore target here and do not update the non-existent information
    return ()
  endif ()
  # get current information of target
  basis_get_project_property (TYPES     PROPERTY IMPORTED_TYPES)
  basis_get_project_property (LOCATIONS PROPERTY IMPORTED_LOCATIONS)
  basis_get_project_property (RANKS     PROPERTY IMPORTED_RANKS)
  list (GET TYPES ${IDX} TYPE)
  list (GET RANKS ${IDX} CURRENT_RANK)
  # decide whether current information shall be overwritten
  if (CMAKE_BUILD_TYPE)
    string (TOUPPER "${CMAKE_BUILD_TYPE}" C)
  else ()
    set (C "NOCONFIG")
  endif ()
  set (
    RANKING
      # first pick
      "IMPORTED_LOCATION_${C}"    # 0) prefer location corresponding to current configuration
      "IMPORTED_LOCATION"         # 1) then use non-configuration specific location
      "IMPORTED_LOCATION_RELEASE" # 2) otherwise use RELEASE version if available
      # 3) last pick, use first imported executable
  )
  list (FIND RANKING "${PROPERTY}" RANK)
  if (RANK EQUAL -1)
    set (RANK 3)
  endif ()
  # bail out if current information shall be kept
  if (NOT "${RANK}" LESS "${CURRENT_RANK}")
    return ()
  endif ()
  # remove current information
  list (REMOVE_AT TYPES     ${IDX})
  list (REMOVE_AT LOCATIONS ${IDX})
  list (REMOVE_AT RANKS     ${IDX})
  # add imported information
  list (LENGTH TYPES N)
  if (IDX LESS N)
    list (INSERT TYPES     ${IDX} "${TYPE}")
    list (INSERT LOCATIONS ${IDX} "${LOCATION}")
    list (INSERT RANKS     ${IDX} "${RANK}")
  else ()
    list (APPEND TYPES     "${TYPE}")
    list (APPEND LOCATIONS "${LOCATION}")
    list (APPEND RANKS     "${RANK}")
  endif ()
  # update project properties
  basis_set_project_property (PROPERTY IMPORTED_TYPES     "${TYPES}")
  basis_set_project_property (PROPERTY IMPORTED_LOCATIONS "${LOCATIONS}")
  basis_set_project_property (PROPERTY IMPORTED_RANKS     "${RANKS}")
endfunction ()

# ----------------------------------------------------------------------------
## @brief Update type of imported target.
#
# This function is in particular called in basis_set_target_properties()
# if the BASIS_TYPE property of custom BASIS targets is set after the
# imported target was added with the initial type UNKNOWN.
#
# @param [in] TARGET Name (UID) of the imported target.
# @param [in] TYPE   Type of imported target.
function (basis_update_imported_type TARGET TYPE)
  # get index of imported target
  basis_get_project_property (TARGETS PROPERTY IMPORTED_TARGETS)
  list (FIND TARGETS "${TARGET}" IDX)
  if (IDX EQUAL -1)
    # imported targets have to be added via basis_add_imported_target() first
    # otherwise, ignore target here and do not update the non-existent information
    return ()
  endif ()
  # get current type of imported target
  basis_get_project_property (TYPES PROPERTY IMPORTED_TYPES)
  list (GET TYPES ${IDX} CURRENT_TYPE)
  # bail out if current type shall be kept
  if (NOT CURRENT_TYPE MATCHES "^UNKNOWN$")
    return ()
  endif ()
  # replace current type
  list (REMOVE_AT TYPES ${IDX})
  list (LENGTH TYPES N)
  if (IDX LESS N)
    list (INSERT TYPES ${IDX} ${TYPE})
  else ()
    list (APPEND TYPES ${TYPE})
  endif ()
  # update project property
  basis_set_project_property (PROPERTY IMPORTED_TYPES "${TYPES}")
endfunction ()


## @}
# end of Doxygen group
