# ============================================================================
# Copyright (c) 2011-2012 University of Pennsylvania
# Copyright (c) 2013-2014 Carnegie Mellon University
# Copyright (c) 2013-2016 Andreas Schuh
# All rights reserved.
#
# See COPYING file for license information or visit
# https://cmake-basis.github.io/download.html#license
# ============================================================================

##############################################################################
# @file  CommonTools.cmake
# @brief Definition of common CMake functions.
#
# @ingroup CMakeTools
##############################################################################

if (__BASIS_COMMONTOOLS_INCLUDED)
  return ()
else ()
  set (__BASIS_COMMONTOOLS_INCLUDED TRUE)
endif ()


include (CMakeParseArguments)


## @addtogroup CMakeUtilities
#  @{


# ============================================================================
# find other packages
# ============================================================================

# ----------------------------------------------------------------------------
## @brief Overloaded find_package() command.
#
# This macro calls CMake's
# <a href="http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:find_package">
# find_package()</a> command and converts obsolete all uppercase "<PKG>_<VAR>"
# variables to case-sensitive "<Pkg>_<VAR>" variables.
# It further ensures that the global variables CMAKE_FIND_LIBRARY_SUFFIXES
# and CMAKE_FIND_EXECUTABLE_SUFFIX are reset to the values they had before
# the call to find_package(). This is required if the "Find<Pkg>.cmake" module
# has modified these variables, but not restored their initial value.
macro (find_package)
  if (BASIS_DEBUG)
    message ("** find_package(${ARGV})")
  endif ()
  # attention: find_package() can be recursive. Hence, use "stack" to keep
  #            track of library suffixes. Further note that we need to
  #            maintain a list of lists, which is not supported by CMake.
  list (APPEND _BASIS_FIND_LIBRARY_SUFFIXES "{${CMAKE_FIND_LIBRARY_SUFFIXES}}")
  list (APPEND _BASIS_FIND_EXECUTABLE_SUFFIX "${CMAKE_FIND_EXECUTABLE_SUFFIX}")
  _find_package(${ARGV})
  # map common uppercase <PKG>_* variables to case-preserving <Pkg>_*
  string (TOUPPER "${ARGV0}" _FP_ARGV0_U)
  foreach (_FP_VAR IN ITEMS FOUND DIR USE_FILE
                            VERSION VERSION_STRING
                            VERSION_MAJOR VERSION_MINOR VERSION_SUBMINOR VERSION_PATCH
                            MAJOR_VERSION MINOR_VERSION SUBMINOR_VERSION PATCH_VERSION
                            INCLUDE_DIR INCLUDE_DIRS INCLUDE_PATH
                            LIBRARY_DIR LIBRARY_DIRS LIBRARY_PATH
                            EXECUTABLE COMPILER CONVERTER)
    if (NOT ${ARGV0}_${_FP_VAR} AND DEFINED ${_FP_ARGV0_U}_${_FP_VAR})
      set (${ARGV0}_${_FP_VAR} "${${_FP_ARGV0_U}_${_FP_VAR}}")
    endif ()
  endforeach ()
  unset (_FP_VAR)
  unset (_FP_ARGV0_U)
  # restore CMAKE_FIND_LIBRARY_SUFFIXES
  string (REGEX REPLACE ";?{([^}]*)}$" "" _BASIS_FIND_LIBRARY_SUFFIXES "${_BASIS_FIND_LIBRARY_SUFFIXES}")
  set (CMAKE_FIND_LIBRARY_SUFFIXES "${CMAKE_MATCH_1}")
  # restore CMAKE_FIND_EXECUTABLE_SUFFIX
  list (LENGTH _BASIS_FIND_EXECUTABLE_SUFFIX _FP_LAST)
  if (_FP_LAST GREATER 0)
    math (EXPR _FP_LAST "${_FP_LAST} - 1")
    list (REMOVE_AT _BASIS_FIND_EXECUTABLE_SUFFIX ${_FP_LAST})
  endif ()
  unset (_FP_LAST)
endmacro ()

# ----------------------------------------------------------------------------
## @brief Tokenize dependency specification.
#
# This function parses a dependency specification such as
# "ITK-4.1{TestKernel,IO}" into the package name, i.e., ITK, the requested
# (minimum) package version(s), i.e., 4.1, and a list of package components, i.e.,
# TestKernel and IO. A valid dependency specification must specify the package
# name of the dependency (case-sensitive). The version and components
# specification are optional. Note that the components specification may
# be separated by an arbitrary number of whitespace characters including
# newlines. The same applies to the specification of the components themselves.
# This allows one to format the dependency specification as follows, for example:
# @code
# ITK {
#   TestKernel,
#   IO
# }
# @endcode
#
# VTK-7|6{}
#
# @param [in]  DEP Dependency specification, i.e., "<Pkg>[-<version>[|...]][{<Component1>[,...]}]".
# @param [out] PKG Package name.
# @param [out] VER Package version(s).
# @param [out] CMP List of package components.
function (basis_tokenize_dependency DEP PKG VER CMP)
  set (CMPS)
  if (DEP MATCHES "^([^ ]+)[ \\n\\t]*{([^}]*)}$")
    set (DEP "${CMAKE_MATCH_1}")
    string (REPLACE "," ";" COMPONENTS "${CMAKE_MATCH_2}")
    foreach (C IN LISTS COMPONENTS)
      string (STRIP "${C}" C)
      list (APPEND CMPS ${C})
    endforeach ()
  endif ()
  if (DEP MATCHES "^(.*)-([0-9]+(\\.[0-9]+)?(\\.[0-9]+)?(\\|[0-9]+(\\.[0-9]+)?(\\.[0-9]+)?)*)$")
    string (REPLACE "|" ";" CMAKE_MATCH_2 "${CMAKE_MATCH_2}")
    set (${PKG} "${CMAKE_MATCH_1}" PARENT_SCOPE)
    set (${VER} "${CMAKE_MATCH_2}" PARENT_SCOPE)
  else ()
    set (${PKG} "${DEP}" PARENT_SCOPE)
    set (${VER} ""       PARENT_SCOPE)
  endif ()
  set (${CMP} "${CMPS}" PARENT_SCOPE)
endfunction ()

## @brief Get installation prefix given path of <PKG>Config directory
function (_basis_config_to_prefix_dir PKG PKG_DIR PREFIX)
  if (PKG_DIR)
    if (APPLE)
      if (PKG_DIR MATCHES "^(.*)/([^/]+)\\.(framework|app)(/|$)")
        set (${PREFIX} "${CMAKE_MATCH_1}/${CMAKE_MATCH_2}.${CMAKE_MATCH_3}" PARENT_SCOPE)
        return ()
      endif ()
    endif ()
    string (REGEX REPLACE "/+(cmake|CMake)/*$" "" prefix "${PKG_DIR}")
    set (SUBDIR_RE "/+(lib[0-9]*(/[^/]+)?|share)(/+cmake)?")
    if (UNIX AND prefix MATCHES "${SUBDIR_RE}(/+[^/]+|/*$)")
      get_filename_component (subdir "${prefix}" NAME)
      string (TOLOWER "${subdir}" subdir)
      string (TOLOWER "${PKG}"    pkg)
      if ("^${subdir}$" STREQUAL "^${pkg}$")
        get_filename_component (prefix "${prefix}" PATH)
      endif ()
      string(REGEX REPLACE "${SUBDIR_RE}/*$" "" prefix "${prefix}")
    endif ()
  else ()
    set(prefix "NOTFOUND")
  endif ()
  set(${PREFIX} "${prefix}" PARENT_SCOPE)
endfunction ()

# ----------------------------------------------------------------------------
## @brief Find external software package or other project module.
#
# This macro replaces CMake's
# <a href="https://cmake.org/cmake/help/v2.8.12/cmake.html#command:find_package">
# find_package()</a> command and extends its functionality.
# In particular, if the given package name is the name of another module
# of this project (the top-level project), it ensures that this module is
# found instead of an external package.
#
# If the package is found, but only optionally used, i.e., the @c REQUIRED
# argument was not given to this macro, a <tt>WITH_&lt;PACKAGE&gt;</tt> option 
# is added by this macro which is by default @c OFF. This option can be set to
# @c ON by the user in order to require the discovery of the package.
# When this option is OFF, the (non-cached) <tt>USE_&lt;PACKAGE&gt;</tt> variable
# can be used by project developers to still request the discovery of the
# optional package, but no error is raised when the package is not found.
# This allows a project to use an optional dependency when an installation
# is found regardless of the <tt>WITH_&lt;PACKAGE&gt;</tt> option. Note
# that when <tt>USE_&lt;PACKAGE&gt;</tt> is defined, no <tt>WITH_&lt;PACKAGE&gt;</tt>
# entry is added by this macro to the cache.
#
# @param [in] PACKAGE Name of software package or other project module.
#                     Optionally, the package name can include a version
#                     specification as suffix which is separated from the
#                     package name using a dash (-), i.e., &lt;Package&gt;[-major[.minor[.patch[.tweak]]]].
#                     Multiple alternative versions have to be separated by
#                     a pipe character "|", the logical OR.
#                     If a version specification is given, it is passed on as
#                     @c version argument to CMake's
#                     <a href="https://cmake.org/cmake/help/v2.8.12/cmake.html#command:find_package">
#                     find_package()</a> command. The discovery of multiple
#                     alternative versions is only supported for the CONFIG
#                     mode of the find_package command.
# @param [in] ARGN    Additional arguments for
#                     <a href="https://cmake.org/cmake/help/v2.8.12/cmake.html#command:find_package">
#                     find_package()</a>.
#
# @retval &lt;PACKAGE&gt;_FOUND            Whether the given package was found.
# @retval &lt;PACKAGE&gt;_COMPONENTS_FOUND Names of found components.
#                                          Optional components are only included when
#                                          "<PACKAGE>_<COMPONENT>_FOUND" is set to @c TRUE
#                                          for each found component by the find_package call,
#                                          i.e., either the "Find<PACKAGE>" module or the
#                                          "<PACKAGE>Config" file.
#
# @sa https://cmake.org/cmake/help/v2.8.12/cmake.html#command:find_package
#
# @ingroup CMakeAPI
macro (basis_find_package PACKAGE)
  # Note that this MUST be a macro such that non-cached variables
  # set by find_package are set within the callers scope
  # ------------------------------------------------------------------------
  # parse arguments
  set (_BFP_OPTIONS
    QUIET
    REQUIRED
    MODULE
    NO_MODULE
    CONFIG
    NO_NOTFOUND_ERROR
  )
  set (_BFP_MULTI_ARGS
    COMPONENTS
    OPTIONAL_COMPONENTS
  )
  cmake_parse_arguments (
    _BFP_ARGN
      "${_BFP_OPTIONS}"
      ""
      "${_BFP_MULTI_ARGS}"
    ${ARGN}
  )
  # ------------------------------------------------------------------------
  # tokenize dependency specification
  basis_tokenize_dependency ("${PACKAGE}" PKG _BFP_VERSIONS _BFP_COMPONENTS)
  list (GET _BFP_ARGN_UNPARSED_ARGUMENTS 0 _BFP_VERSION)
  if (_BFP_VERSION MATCHES "^[0-9]+(\\.[0-9]+)*(\\|[0-9]+(\\.[0-9]+)*)*$")
    if (_BFP_VERSIONS)
      message (FATAL_ERROR "Cannot use both version specification as part of "
                           "package name and explicit version argument.")
    endif ()
    list (REMOVE_AT _BFP_ARGN_UNPARSED_ARGUMENTS 0)
    string(REPLACE "|" ";" _BFP_VERSIONS ${_BFP_VERSION})
  endif ()
  list (LENGTH _BFP_VERSIONS _BFP_VERSIONS_COUNT)
  if (_BFP_ARGN_MODULE AND _BFP_VERSIONS GREATER 1)
    message (FATAL_ERROR "Cannot look for multiple alternative package versions"
                         " in MODULE mode. The CONFIG|NO_MODULE mode of find_package"
                         " is used in this case. When MODULE mode is required"
                         " by package ${PKG}, only one version can be specified.")
  endif ()
  string (TOLOWER "${PKG}" PKG_L)
  string (TOUPPER "${PKG}" PKG_U)
  # ------------------------------------------------------------------------
  # set <PKG>_FIND_REQUIRED_<CMP>
  foreach (_BFP_CMP IN LISTS _BFP_COMPONENTS)
    set (${PKG}_FIND_REQUIRED_${_BFP_CMP} ${_BFP_ARGN_REQUIRED})
  endforeach ()
  foreach (_BFP_CMP IN LISTS _BFP_ARGN_COMPONENTS)
    set (${PKG}_FIND_REQUIRED_${_BFP_CMP} TRUE)
  endforeach ()
  foreach (_BFP_CMP IN LISTS _BFP_ARGN_OPTIONAL_COMPONENTS)
    set (${PKG}_FIND_REQUIRED_${_BFP_CMP} FALSE)
  endforeach ()
  list (APPEND _BFP_ARGN_COMPONENTS ${_BFP_COMPONENTS})
  if (_BFP_ARGN_COMPONENTS)
    list (REMOVE_DUPLICATES _BFP_ARGN_COMPONENTS)
  endif ()
  # ------------------------------------------------------------------------
  # prefix of package variable names set by Find<Pkg> or <Pkg>Config
  if (PKG MATCHES "^((P|J)ython)Interp$")
    string (TOUPPER "${CMAKE_MATCH_1}" _BFP_NS)
  else ()
    set (_BFP_NS "${PKG}")
  endif ()
  # ------------------------------------------------------------------------
  # some debugging output
  if (BASIS_DEBUG)
    set (_BFP_ARGS)
    if (_BFP_ARGN_REQUIRED)
      list (APPEND _BFP_ARGS REQUIRED)
    endif ()
    if (_BFP_ARGN_QUIET)
      list (APPEND _BFP_ARGS QUIET)
    endif ()
    if (_BFP_ARGN_MODULE)
      list (APPEND _BFP_ARGS MODULE)
    endif ()
    if (_BFP_ARGN_NO_MODULE OR _BFP_ARGN_CONFIG)
      list (APPEND _BFP_ARGS CONFIG)
    endif ()
    list (APPEND _BFP_ARGS ${_BFP_ARGN_UNPARSED_ARGUMENTS})
    set (_BFP_INFO "** basis_find_package()")
    set (_BFP_INFO "${_BFP_INFO}\n**     Package:    ${PKG}")
    if (_BFP_VERSIONS)
    set (_BFP_INFO "${_BFP_INFO}\n**     Versions:   ${_BFP_VERSIONS}")
    endif ()
    set (_BFP_INFO "${_BFP_INFO}\n**     Arguments:  [${_BFP_ARGS}]")
    if (_BFP_ARGN_COMPONENTS OR _BFP_ARGN_OPTIONAL_COMPONENTS)
    set (_BFP_INFO "${_BFP_INFO}\n**     Components: ")
      if (_BFP_ARGN_COMPONENTS AND _BFP_ARGN_OPTIONAL_COMPONENTS)
        set (_BFP_INFO "${_BFP_INFO}[${_BFP_ARGN_COMPONENTS}] and [${_BFP_ARGN_OPTIONAL_COMPONENTS}] (optional)")
      elseif (_BFP_ARGN_COMPONENTS)
        set (_BFP_INFO "${_BFP_INFO}[${_BFP_ARGN_COMPONENTS}]")
      else ()
        set (_BFP_INFO "${_BFP_INFO}[${_BFP_ARGN_OPTIONAL_COMPONENTS}] (optional)")
      endif ()
    endif ()
    message ("${_BFP_INFO}")
    unset (_BFP_INFO)
    unset (_BFP_ARGS)
  endif ()
  # ------------------------------------------------------------------------
  # find other modules of same project
  set (_BFP_IS_PROJECT FALSE)
  set (_BFP_IS_MODULE  FALSE)
  if (PROJECT_IS_MODULE)
    # allow modules to specify top-level project as dependency,
    # respectively, other modules as components of top-level project
    if ("^${PKG}$" STREQUAL "^${TOPLEVEL_PROJECT_NAME}$")
      if (_BFP_ARGN_COMPONENTS OR _BFP_ARGN_OPTIONAL_COMPONENTS)
        if (BASIS_DEBUG)
          message ("**     This is the top-level project. Components must be other modules of this project.")
        endif ()
        foreach (_BFP_CMP IN LISTS _BFP_ARGN_COMPONENTS _BFP_ARGN_OPTIONAL_COMPONENTS)
          list (FIND PROJECT_MODULES "${_BFP_CMP}" _BFP_CMPIDX)
          if (_BFP_CMPIDX EQUAL -1)
            message (FATAL_ERROR "Module ${PROJECT_NAME} has module ${_BFP_CMP} of top-level project ${PKG}"
                                 " as dependency, but no such module exists.")
          endif ()
          list (FIND PROJECT_MODULES_ENABLED "${_BFP_CMP}" _BFP_CMPIDX)
          if (_BFP_CMPIDX EQUAL -1)
            if (BASIS_DEBUG)
              message ("**     Identified it as disabled module of this project.")
            endif ()
            if (${PKG}_FIND_REQUIRED_${_BFP_CMP})
              if (_BFP_ARGN_REQUIRED)
                message (FATAL_ERROR
                  "Module ${PROJECT_NAME} requires module ${_BFP_CMP} of top-level project ${PKG},"
                  " but module ${_BFP_CMP} is not enabled."
                )
              elseif (NOT _BFP_ARGN_QUIET)
                message (STATUS
                  "Module ${PROJECT_NAME} optionally uses module ${_BFP_CMP} of top-level project ${PKG},"
                  " but module ${_BFP_CMP} is not enabled."
                )
              endif ()
            endif ()
            set (${PKG}_${_BFP_CMP}_FOUND FALSE)
          else ()
            if (BASIS_DEBUG)
              message ("**     Identified it as enabled module of this project.")
            endif ()
            include ("${BINARY_LIBCONF_DIR}/${TOPLEVEL_PROJECT_PACKAGE_CONFIG_PREFIX}${_BFP_CMP}Config.cmake")
            set (${PKG}_${_BFP_CMP}_FOUND TRUE)
          endif ()
        endforeach ()
      else ()
        if (BASIS_DEBUG)
          message ("**     This is the top-level project.")
        endif ()
      endif ()
      set (${PKG}_FOUND    TRUE)
      set (_BFP_IS_PROJECT TRUE)
    # look for other module of top-level project
    else ()
      list (FIND PROJECT_MODULES "${PKG}" _BFP_PKGIDX)
      if (NOT _BFP_PKGIDX EQUAL -1)
        set (_BFP_IS_MODULE TRUE)
        list (FIND PROJECT_MODULES_ENABLED "${PKG}" _BFP_PKGIDX)
        if (_BFP_PKGIDX EQUAL -1)
          if (BASIS_DEBUG)
            message ("**     Identified it as disable module of this project.")
          endif ()
          if (_BFP_ARGN_REQUIRED)
            message (FATAL_ERROR
              "Module ${PROJECT_NAME} requires module ${_BFP_CMP} of top-level project ${PKG},"
              " but module ${_BFP_CMP} is not enabled."
            )
          elseif (NOT _BFP_ARGN_QUIET)
            message (STATUS
              "Module ${PROJECT_NAME} optionally uses module ${_BFP_CMP} of top-level project ${PKG},"
              " but module ${_BFP_CMP} is not enabled."
            )
          endif ()
          set (${PKG}_FOUND FALSE)
        else ()
          if (BASIS_DEBUG)
            message ("**     Identified it as enabled module of this project.")
          endif ()
          include ("${BINARY_LIBCONF_DIR}/${TOPLEVEL_PROJECT_PACKAGE_CONFIG_PREFIX}${PKG}Config.cmake")
          set (${PKG}_FOUND TRUE)
        endif ()
      endif ()
    endif ()
  # --------------------------------------------------------------------------
  # find bundled packages
  else ()
    list (FIND BUNDLE_PROJECTS "${PKG}" _BFP_PKGIDX)
    if (NOT _BFP_PKGIDX EQUAL -1)
      if  (EXISTS "${CMAKE_INSTALL_PREFIX}/${INSTALL_CONFIG_DIR}/${PKG}Config.cmake")
        set (_BFP_CONFIG_FILE "${CMAKE_INSTALL_PREFIX}/${INSTALL_CONFIG_DIR}/${PKG}Config.cmake")
      else ()
        if (EXISTS "${CMAKE_INSTALL_PREFIX}/${INSTALL_CONFIG_DIR}/${PKG_L}-config.cmake")
          set (_BFP_CONFIG_FILE "${CMAKE_INSTALL_PREFIX}/${INSTALL_CONFIG_DIR}/${PKG_L}-config.cmake")
        else ()
          set (_BFP_CONFIG_FILE)
        endif ()
      endif ()
      if (_BFP_CONFIG_FILE)
        if (BASIS_DEBUG)
          message ("**     Identified it as other package of this bundle.")
        endif ()
        get_filename_component (_BFP_CONFIG_DIR "${_BFP_CONFIG_FILE}" PATH)
        _basis_config_to_prefix_dir (${PKG} "${_BFP_CONFIG_DIR}" _BFP_PREFIX)
        basis_set_or_update_value (DEPENDS_${PKG}_DIR "${_BFP_PREFIX}" PATH)
        include ("${_BFP_CONFIG_FILE}")
        set (${PKG}_FOUND TRUE)
        unset (_BFP_CONFIG_DIR)
      endif ()
      unset (_BFP_CONFIG_FILE)
    endif ()
  endif ()
  # --------------------------------------------------------------------------
  # otherwise, look for external package
  if (NOT _BFP_IS_PROJECT AND NOT _BFP_IS_MODULE)
    # ------------------------------------------------------------------------
    # provide option which allows users to request use of optional packages
    #
    # - WITH_<PKG>  == ON:  optional dependency is required
    # - WITH_<PKG>  == OFF:
    #   - USE_<PKG> == ON:  dependency is looked for and used if available
    #   - USE_<PKG> == OFF: dependency is ignored
    #
    # The default of the WITH_<PKG> option is OFF unless WITH_<PKG>_DEFAULT
    # is set in the Depends.cmake, the CMake command-line using the -D switch
    # or the root CMakeLists.txt file before the basis_project_begin call.
    # The (uncached) WITH_${PKG}_DEFAULT variable can be used by a project
    # to require optional dependencies by default, e.g., to enable optional
    # program features that depend on these external packages.
    if (NOT _BFP_ARGN_REQUIRED)
      if (NOT DEFINED WITH_${PKG}_DEFAULT)
        set (WITH_${PKG}_DEFAULT OFF)
      endif ()
      if (DEFINED USE_${PKG} AND NOT DEFINED WITH_${PKG})
        set (WITH_${PKG} ${WITH_${PKG}_DEFAULT})
      else ()
        option (WITH_${PKG} "Build with optional ${PKG} dependency" ${WITH_${PKG}_DEFAULT})
      endif ()
    endif ()
    # look for external package only if required, built with optional dependency
    # enabled by user (cf. WITH_<PKG> option above) or deprecated -DUSE_<PKG>=ON
    if (_BFP_ARGN_REQUIRED OR WITH_${PKG} OR USE_${PKG})
      # ----------------------------------------------------------------------
      # Use more user friendly hybrid DEPENDS_<PKG>_DIR cache variable which
      # allows grouping of DEPENDS paths cache entry, but still consider more
      # common variables named <PKG>_DIR, <PKG>_ROOT, <PKG>ROOT, or <PKG>_ROOT_DIR
      # set in the user shell environment or on the CMake command-line using -D.
      set (
        DEPENDS_${PKG}_DIR "${DEPENDS_${PKG}_DIR}" CACHE PATH
        "Top-level installation directory of ${PKG} or directory containing ${PKG}Config.cmake file."
      )
      if (DEPENDS_${PKG}_DIR)
        file (TO_CMAKE_PATH "${DEPENDS_${PKG}_DIR}" _BFP_PREFIX)
        if (NOT "^${DEPENDS_${PKG}_DIR}$" STREQUAL "^${_BFP_PREFIX}$")
          basis_update_value (DEPENDS_${PKG}_DIR "${_BFP_PREFIX}")
        endif ()
      endif ()
      # Names of considered alternative find search path variables excl. <PKG>_DIR
      set (_BFP_PKG_DIR_VARS
        ${PKG}_ROOT     ${PKG_U}_ROOT
        ${PKG}ROOT      ${PKG_U}ROOT
        ${PKG}_ROOT_DIR ${PKG_U}_ROOT_DIR
      )
      list (REMOVE_DUPLICATES _BFP_PKG_DIR_VARS)
      # Override DEPENDS_<PKG>_DIR by alternative search path variable value if these
      # were specified on the command line using the -D option. Note that these variables
      # cannot be set in the CMake GUI because their type is changed here to INTERNAL.
      # This has two reasons, firstly to not have duplicate variables with different
      # names for the same purpose, and secondly to be able to recognize when their
      # value is changed using the -D command line option of the cmake command.
      #
      # Order of precedence:
      # 1. <PKG>_DIR
      # 2. <PKG>_ROOT... CMake variable
      # 3. <PKG>_ROOT... environment variable
      foreach (_BFP_VAR IN LISTS _BFP_PKG_DIR_VARS)
        file (TO_CMAKE_PATH "${${_BFP_VAR}}" _BFP_PREFIX) # CMake (cache) variable
        if (_BFP_PREFIX)
          # first configure run or new value specified using -D option of cmake command
          if (NOT DEFINED _DEPENDS_${PKG}_DIR OR (DEFINED _DEPENDS_${PKG}_DIR AND NOT "^${_BFP_PREFIX}$" STREQUAL "^${_DEPENDS_${PKG}_DIR}$"))
            basis_update_value (DEPENDS_${PKG}_DIR "${_BFP_PREFIX}")
            break ()
          endif ()
        endif ()
      endforeach ()
      if (${PKG}_DIR) # find_package CONFIG mode variable
        # first configure run or new value specified using -D option of cmake command
        if (NOT DEFINED _${PKG}_DIR OR (DEFINED _${PKG}_DIR AND NOT "^${${PKG}_DIR}$" STREQUAL "^${_${PKG}_DIR}$"))
          file (TO_CMAKE_PATH "${${PKG}_DIR}" _BFP_PREFIX)
          _basis_config_to_prefix_dir(${PKG} "${_BFP_PREFIX}" _BFP_PREFIX)
          basis_update_value (DEPENDS_${PKG}_DIR "${_BFP_PREFIX}")
        endif ()
      endif ()
      # mark alternatives as internal cache entries
      foreach (_BFP_VAR IN LISTS _BFP_PKG_DIR_VARS)
        basis_is_cached (_BFP_CACHED ${_BFP_VAR})
        if (_BFP_CACHED)
          basis_update_type_of_variable (${_BFP_VAR} INTERNAL)
        endif ()
      endforeach ()
      # if still not set, use common environment variables to set DEPENDS_<PKG>_DIR
      if (NOT DEPENDS_${PKG}_DIR)
        foreach (_BFP_VAR IN LISTS _BFP_PKG_DIR_VARS)
          file (TO_CMAKE_PATH "$ENV{${_BFP_VAR}}" _BFP_PREFIX) # shell environment variable
          if (_BFP_PREFIX)
            basis_update_value (DEPENDS_${PKG}_DIR "${_BFP_PREFIX}")
            break()
          endif ()
        endforeach ()
      endif ()
      # circumvent issue with find_package interpreting <PKG>_DIR relative
      # to the current binary directory instead of the top-level directory
      if (DEPENDS_${PKG}_DIR AND NOT IS_ABSOLUTE "${DEPENDS_${PKG}_DIR}")
        get_filename_component (_BFP_PREFIX "${CMAKE_BINARY_DIR}/${DEPENDS_${PKG}_DIR}" ABSOLUTE)
        basis_update_value (DEPENDS_${PKG}_DIR "${_BFP_PREFIX}")
      endif ()
      # --------------------------------------------------------------------
      # reset _?<PKG>[_-]* variables when DEPENDS_<PKG>_DIR has changed
      if (_DEPENDS_${PKG}_DIR AND DEPENDS_${PKG}_DIR)
        if (NOT "^${_DEPENDS_${PKG}_DIR}$" STREQUAL "^${DEPENDS_${PKG}_DIR}$")
          get_cmake_property (_BFP_VARS VARIABLES)
          basis_sanitize_for_regex (PKG_RE     "${PKG}")
          basis_sanitize_for_regex (PKG_U_RE   "${PKG_U}")
          basis_sanitize_for_regex (_BFP_NS_RE "${_BFP_NS}")
          foreach (_BFP_VAR IN LISTS _BFP_VARS)
            if (_BFP_VAR MATCHES "^_?(${_BFP_NS_RE}|${PKG_RE}|${PKG_U_RE})[_-]")
              basis_is_cached (_BFP_CACHED ${_BFP_VAR})
              if (_BFP_CACHED)
                set_property (CACHE ${_BFP_VAR} PROPERTY VALUE "${_BFP_VAR}-NOTFOUND")
                set_property (CACHE ${_BFP_VAR} PROPERTY TYPE  INTERNAL)
              elseif (NOT _BFP_VAR MATCHES "^${PKG_RE}_FIND_")
                unset (${_BFP_VAR})
              endif ()
            endif ()
          endforeach ()
          unset (PKG_RE)
          unset (PKG_U_RE)
          unset (_BFP_NS_RE)
        endif ()
      endif ()
      # ----------------------------------------------------------------------
      # determine if additional components of found package should be discovered
      if (${PKG}_FOUND)
        set (_BFP_NO_FIND_PACKAGE TRUE)
        if (${PKG}_COMPONENTS_FOUND)
          # previously called with COMPONENTS
          set (_BFP_FIND_COMPONENTS)
          set (_BFP_FIND_OPTIONAL_COMPONENTS)
          if (_BFP_ARGN_COMPONENTS OR _BFP_ARGN_OPTIONAL_COMPONENTS)
            foreach (_BFP_CMP IN LISTS _BFP_ARGN_COMPONENTS _BFP_ARGN_OPTIONAL_COMPONENTS)
              list (FIND ${PKG}_COMPONENTS_FOUND "${_BFP_CMP}" _BFP_CMPIDX)
              if (_BFP_CMPIDX EQUAL -1)
                if (${PKG}_FIND_REQUIRED_${_BFP_CMP})
                  list (APPEND _BFP_FIND_COMPONENTS "${_BFP_CMP}")
                elseif (NOT DEFINED ${PKG}_${_BFP_CMP}_FOUND)
                  list (APPEND _BFP_FIND_OPTIONAL_COMPONENTS "${_BFP_CMP}")
                endif ()
              endif ()
            endforeach ()
          else ()
            # Not sure if "default" components were found when find_package
            # was previously invoked with the COMPONENTS argument, but
            # now without it. This depends on the Find<PKG> module.
            set (_BFP_NO_FIND_PACKAGE FALSE)
          endif ()
        else ()
          # previously called without COMPONENTS
          set (_BFP_FIND_COMPONENTS          ${_BFP_ARGN_COMPONENTS})
          set (_BFP_FIND_OPTIONAL_COMPONENTS ${_BFP_ARGN_OPTIONAL_COMPONENTS})
        endif ()
        if (_BFP_FIND_COMPONENTS OR _BFP_FIND_OPTIONAL_COMPONENTS)
          set (_BFP_NO_FIND_PACKAGE FALSE)
        endif ()
      else ()
        set (_BFP_NO_FIND_PACKAGE          FALSE)
        set (_BFP_FIND_COMPONENTS          ${_BFP_ARGN_COMPONENTS})
        set (_BFP_FIND_OPTIONAL_COMPONENTS ${_BFP_ARGN_OPTIONAL_COMPONENTS})
      endif ()
      # ----------------------------------------------------------------------
      # look for external package if not found or additional components needed
      if (NOT ${PKG}_FOUND AND NOT _BFP_ARGN_REQUIRED AND NOT WITH_${PKG} AND ((DEFINED USE_${PKG} AND NOT USE_${PKG}) OR NOT DEFINED USE_${PKG}))
        # skip if package is optional and user did not ask us to look for it
        # when package was found before, still perform steps below to ensure
        # that everything is set properly even when find_package was called
        mark_as_advanced (FORCE DEPENDS_${PKG}_DIR)
      else ()
        # status message with information what we are looking for
        if (_BFP_ARGN_QUIET)
          set (_BFP_STATUS)
        else ()
          set (_BFP_STATUS "Looking for ${PKG}")
          if (_BFP_VERSIONS)
            string (REPLACE ";" " or " _BFP_VERSION_STRING "${_BFP_VERSIONS}")
            set (_BFP_STATUS "${_BFP_STATUS} ${_BFP_VERSION_STRING}")
            unset (_BFP_VERSION_STRING)
          endif ()
          if (_BFP_FIND_COMPONENTS)
            set (_BFP_STATUS "${_BFP_STATUS} [${_BFP_FIND_COMPONENTS}]")
          endif ()
          if (_BFP_FIND_OPTIONAL_COMPONENTS)
            set (_BFP_STATUS "${_BFP_STATUS}, optional components [${_BFP_FIND_OPTIONAL_COMPONENTS}]")
          endif ()
          if (NOT _BFP_ARGN_REQUIRED)
            set (_BFP_STATUS "${_BFP_STATUS} (optional)")
          endif ()
          set (_BFP_STATUS "${_BFP_STATUS}...")
          message (STATUS "${_BFP_STATUS}")
        endif ()
        # make copy of previous <Pkg>_VERSION_STRING if already set which is used
        # as "last resort" when variable not set by the following find_package
        set (_BFP_VERSION_STRING "${${_BFP_NS}_VERSION_STRING}")
        unset (${_BFP_NS}_VERSION_STRING)
        # set internal <Pkg>_DIR used by find_package to locate <Pkg>Config
        set (${PKG}_DIR "${DEPENDS_${PKG}_DIR}" CACHE INTERNAL "Directory containing ${PKG}Config.cmake file." FORCE)
        # call find_package if not all components found yet
        set (_BFP_FOUND "${${PKG}_FOUND}") # used to decide what the intersection of
                                           # multiple find invocations for the same
                                           # package with different components will be
                                           # for the setting of <PKG>_FOUND
        if (NOT _BFP_NO_FIND_PACKAGE)
          # reset <PKG>_FOUND
          set (${PKG}_FOUND   FALSE)
          set (${PKG_U}_FOUND FALSE)
          # make copy of find_* search path variables
          foreach (_BFP_VAR IN ITEMS CMAKE_PREFIX_PATH CMAKE_PROGRAM_PATH)
            set (_BFP_${_BFP_VAR} "${${_BFP_VAR}}")
          endforeach ()
          # insert DEPENDS_<PKG>_DIR into find_* search paths
          if (IS_DIRECTORY "${DEPENDS_${PKG}_DIR}")
            # add directory to CMAKE_PROGRAM_PATH when the name
            # of the last subdirectory is "bin", "Bin", "sbin",
            # or "texbin" (i.e., MacTeX's "/Library/TeX/texbin" path)
            if (DEPENDS_${PKG}_DIR MATCHES "/([bB]in|sbin|texbin)/+$")
              list (INSERT CMAKE_PROGRAM_PATH 0 "${DEPENDS_${PKG}_DIR}")
              list (REMOVE_DUPLICATES CMAKE_PROGRAM_PATH)
            # add directory to CMAKE_PREFIX_PATH otherwise as users
            # tend to specify the installation prefix instead of the
            # actual directory containing the package configuration file
            else ()
              list (INSERT CMAKE_PREFIX_PATH 0 "${DEPENDS_${PKG}_DIR}")
              list (REMOVE_DUPLICATES CMAKE_PREFIX_PATH)
            endif ()
          endif ()
          # now look for the package using find_package
          set (_BFP_FIND_PACKAGE_ARGS ${_BFP_ARGN_UNPARSED_ARGUMENTS})
          if (_BFP_ARGN_QUIET OR "^${PKG}$" STREQUAL "^Boost$")
            list (APPEND _BFP_FIND_PACKAGE_ARGS "QUIET")
          endif ()
          if (_BFP_FIND_COMPONENTS OR _BFP_FIND_OPTIONAL_COMPONENTS)
            if (${PKG}_COMPONENTS_FOUND OR _BFP_FIND_COMPONENTS)
              list (APPEND _BFP_FIND_PACKAGE_ARGS "COMPONENTS" ${${PKG}_COMPONENTS_FOUND} ${_BFP_FIND_COMPONENTS})
            endif ()
            if (_BFP_FIND_OPTIONAL_COMPONENTS)
              list (APPEND _BFP_FIND_PACKAGE_ARGS "OPTIONAL_COMPONENTS" ${_BFP_FIND_OPTIONAL_COMPONENTS})
            endif ()
          endif ()
          if (${_BFP_VERSIONS_COUNT} GREATER 1)
            list (APPEND _BFP_FIND_PACKAGE_ARGS CONFIG)
            if (DEPENDS_${PKG}_DIR)
              foreach (_BFP_VERSION ${_BFP_VERSIONS})
                find_package (${PKG} ${_BFP_VERSION} ${_BFP_FIND_PACKAGE_ARGS} QUIET
                  PATHS ${DEPENDS_${PKG}_DIR} NO_DEFAULT_PATH
                )
                if (${PKG}_FOUND)
                  break ()
                endif ()
              endforeach ()
            endif ()
            if (NOT ${PKG}_FOUND)
              foreach (_BFP_VERSION ${_BFP_VERSIONS})
                find_package (${PKG} ${_BFP_VERSION} ${_BFP_FIND_PACKAGE_ARGS} QUIET)
                if (${PKG}_FOUND)
                  break ()
                endif ()
              endforeach ()
            endif ()
          else ()
            if (_BFP_ARGN_MODULE)
              list (APPEND _BFP_FIND_PACKAGE_ARGS MODULE)
            endif ()
            if (_BFP_ARGN_NO_MODULE OR _BFP_ARGN_CONFIG)
              list (APPEND _BFP_FIND_PACKAGE_ARGS CONFIG)
            endif ()
            find_package (${PKG} ${_BFP_VERSIONS} ${_BFP_FIND_PACKAGE_ARGS} QUIET)
          endif ()
          unset (_BFP_FIND_PACKAGE_ARGS)
          # restore find_* search path variables
          foreach (_BFP_VAR IN ITEMS CMAKE_PREFIX_PATH CMAKE_PROGRAM_PATH)
            set (${_BFP_VAR} "${_BFP_${_BFP_VAR}}")
            unset (_BFP_${_BFP_VAR})
          endforeach ()
          # ensure that <PKG>_DIR is still an internal cache entry
          basis_update_type_of_variable (${PKG}_DIR INTERNAL)
          # force reinclusion of package use file
          if (${PKG}_FOUND)
            if (${PKG}_USE_FILE_INCLUDED)
              set (${PKG}_USE_FILE_INCLUDED 0)
            endif ()
            if (BASIS_USE_${PKG}_INCLUDED)
              set (BASIS_USE_${PKG}_INCLUDED FALSE)
            endif ()
          endif ()
        endif () # NOT _BFP_NO_FIND_PACKAGE
        # set common <Pkg>_VERSION_STRING variable if possible and not set
        if (NOT DEFINED ${_BFP_NS}_VERSION_STRING)
          if (DEFINED ${_BFP_NS}_VERSION_MAJOR)
            set (${_BFP_NS}_VERSION_STRING ${${_BFP_NS}_VERSION_MAJOR})
            if (DEFINED ${_BFP_NS}_VERSION_MINOR)
              set (${_BFP_NS}_VERSION_STRING ${${_BFP_NS}_VERSION_STRING}.${${_BFP_NS}_VERSION_MINOR})
              if (DEFINED ${_BFP_NS}_VERSION_PATCH)
                set (${_BFP_NS}_VERSION_STRING ${${_BFP_NS}_VERSION_STRING}.${${_BFP_NS}_VERSION_PATCH})
              elseif (DEFINED ${_BFP_NS}_SUBMINOR_VERSION) # e.g., FindBoost
                set (${_BFP_NS}_VERSION_STRING ${${_BFP_NS}_VERSION_STRING}.${${_BFP_NS}_SUBMINOR_VERSION})
              endif ()
            endif ()
          elseif (DEFINED ${_BFP_NS}_MAJOR_VERSION)
            set (${_BFP_NS}_VERSION_STRING ${${_BFP_NS}_MAJOR_VERSION})
            if (DEFINED ${_BFP_NS}_MINOR_VERSION)
              set (${_BFP_NS}_VERSION_STRING ${${_BFP_NS}_VERSION_STRING}.${${_BFP_NS}_MINOR_VERSION})
              if (DEFINED ${PKG}_PATCH_VERSION)
                set (${_BFP_NS}_VERSION_STRING ${${_BFP_NS}_VERSION_STRING}.${${_BFP_NS}_PATCH_VERSION})
              elseif (DEFINED ${_BFP_NS}_SUBMINOR_VERSION) # e.g., FindBoost
                set (${_BFP_NS}_VERSION_STRING ${${_BFP_NS}_VERSION_STRING}.${${_BFP_NS}_SUBMINOR_VERSION})
              endif ()
            endif ()
          elseif (DEFINED ${_BFP_NS}_VERSION)
            set (${_BFP_NS}_VERSION_STRING ${${_BFP_NS}_VERSION})
          else ()
            set (${_BFP_NS}_VERSION_STRING "${_BFP_VERSION_STRING}")
          endif ()
        endif ()
        unset (_BFP_VERSION_STRING)
        # update DEPENDS_<PKG>_DIR from variables set by find_package
        if (${PKG}_FOUND)
          if (${PKG}_DIR AND IS_ABSOLUTE "${${PKG}_DIR}" AND
              (EXISTS "${${PKG}_DIR}/${PKG}Config.cmake" OR
               EXISTS "${${PKG}_DIR}/${PKG_L}-config.cmake"))
            _basis_config_to_prefix_dir(${PKG} "${${PKG}_DIR}" _BFP_PREFIX)
            basis_update_value (DEPENDS_${PKG}_DIR "${_BFP_PREFIX}")
          elseif (NOT DEPENDS_${PKG}_DIR)
            if (${_BFP_NS}_INCLUDE_DIR)
              list (GET ${_BFP_NS}_INCLUDE_DIR 0 _BFP_PREFIX)
              string (REGEX REPLACE "^(.*)/[iI]ncludes?(/.*)?$" "\\1" _BFP_PREFIX "${_BFP_PREFIX}")
            elseif (${_BFP_NS}_LIBRARY_DIR)
              list (GET ${_BFP_NS}_LIBRARY_DIR 0 _BFP_PREFIX)
              string (REGEX REPLACE "^(.*)/[lL]ib(s|exec|64)?(/.*)?$" "\\1" _BFP_PREFIX "${_BFP_PREFIX}")
            else ()
              set (_BFP_VARS
                ${_BFP_NS}_EXECUTABLE # e.g., FindBASH
                ${_BFP_NS}_COMPILER   # e.g., FindLATEX
                ${_BFP_NS}_CONVERTER  # e.g., FindLATEX
              )
              foreach (_BFP_VAR IN LISTS _BFP_VARS)
                if (${_BFP_VAR})
                  get_filename_component (_BFP_PREFIX "${${_BFP_VAR}}" PATH)
                  break ()
                endif ()
              endforeach ()
            endif ()
            basis_update_value (DEPENDS_${PKG}_DIR "${_BFP_PREFIX}")
          endif ()
        else ()
          basis_update_value (DEPENDS_${PKG}_DIR "NOTFOUND")
        endif ()
        # make internal copy of DEPENDS_<PKG>_DIR used to detect change
        set (_DEPENDS_${PKG}_DIR "${DEPENDS_${PKG}_DIR}" CACHE INTERNAL "(Previous) DEPENDS_${PKG}_DIR value." FORCE)
        # make internal search path cache entries consistent with DEPENDS_<PKG>_DIR
        foreach (_BFP_VAR IN LISTS _BFP_PKG_DIR_VARS)
          basis_is_cached (_BFP_CACHED ${_BFP_VAR})
          if (_BFP_CACHED)
            set (${_BFP_VAR} "${DEPENDS_${PKG}_DIR}" CACHE INTERNAL "Installation prefix of ${PKG}." FORCE)
          endif ()
        endforeach ()
        # make internal copy of <PKG>_DIR used to detect change via -D option
        #
        # Note: All other alternative variables such as <PKG>_ROOT are forced to
        #       be equal DEPENDS_<PKG>_DIR. Only <PKG>_DIR usually points to the
        #       <PKG>Config, while DEPENDS_<PKG>_DIR is the installation prefix.
        set (_${PKG}_DIR "${${PKG}_DIR}" CACHE INTERNAL "(Previous) ${PKG}_DIR value." FORCE)
        # status message with information about found package
        if (_BFP_STATUS)
          if (${PKG}_FOUND)
            if (_BFP_NO_FIND_PACKAGE)
              set (_BFP_STATUS "${_BFP_STATUS} - already found")
            else ()
              set (_BFP_STATUS "${_BFP_STATUS} - found")
            endif ()
            if ("^${PKG}$" STREQUAL "^Boost$")
              set (_BFP_STATUS "${_BFP_STATUS} v${${_BFP_NS}_MAJOR_VERSION}.${${_BFP_NS}_MINOR_VERSION}.${${_BFP_NS}_SUBMINOR_VERSION}")
            elseif (DEFINED ${_BFP_NS}_VERSION_STRING AND NOT ${_BFP_NS}_VERSION_STRING MATCHES "^(0(\\.0)?(\\.0)?)?$")
              set (_BFP_STATUS "${_BFP_STATUS} v${${_BFP_NS}_VERSION_STRING}")
            endif ()
            if (BASIS_VERBOSE AND DEPENDS_${PKG}_DIR)
              set (_BFP_STATUS "${_BFP_STATUS} at ${DEPENDS_${PKG}_DIR}")
            endif ()
          else ()
            set (_BFP_STATUS "${_BFP_STATUS} - not found")
          endif ()
          message (STATUS "${_BFP_STATUS}")
        endif ()
        unset (_BFP_STATUS)
        # raise error when a required package was not found
        if (NOT ${PKG}_FOUND AND NOT _BFP_ARGN_NO_NOTFOUND_ERROR AND (_BFP_ARGN_REQUIRED OR WITH_${PKG}))
          set (_BFP_ERROR)
          if (PROJECT_IS_MODULE)
            set (_BFP_ERROR "Module")
          else ()
            set (_BFP_ERROR "Project")
          endif ()
          set (_BFP_ERROR "${_BFP_ERROR} ${PROJECT_NAME}")
          if (_BFP_ARGN_REQUIRED)
            set (_BFP_ERROR "${_BFP_ERROR} requires ${PKG}")
          else ()
            set (_BFP_ERROR "${_BFP_ERROR} was requested to be build with ${PKG}")
          endif ()
          if (_BFP_VERSIONS)
            set (_BFP_ERROR "${_BFP_ERROR} version ${_BFP_VERSIONS}")
          endif ()
          set (_BFP_ERROR "${_BFP_ERROR}. Please ensure that the package is installed in a"
                          " standard system location or set DEPENDS_${PKG}_DIR to the"
                          " installation prefix (i.e., root directory of the installation).")
          if (NOT _BFP_ARGN_REQUIRED AND DEFINED WITH_${PKG})
            set (_BFP_ERROR "${_BFP_ERROR} To disable features which require this optional dependency,"
                            " set the WITH_${PKG} option to OFF and try again.")
          endif ()
          if (DEFINED ${PKG}_DIR)
            set (_BFP_ERROR "${_BFP_ERROR}\nThe DEPENDS_${PKG}_DIR variable can alternatively be set"
                            " to the directory containing a ${PKG}Config.cmake or ${PKG_L}-config.cmake"
                            " file. If no such file exists, contact either the developer of"
                            " this project or CMake BASIS to provide a Find${PKG}.cmake file.")
          endif ()
          basis_list_to_string(_BFP_ERROR ${_BFP_ERROR})
          message (FATAL_ERROR "\n${_BFP_ERROR}\n")
        endif ()
        # update list of found components
        if (${PKG}_FOUND)
          if (_BFP_FIND_COMPONENTS)
            list (APPEND ${PKG}_COMPONENTS_FOUND ${_BFP_FIND_COMPONENTS})
          endif ()
          foreach (_BFP_CMP IN LISTS _BFP_FIND_OPTIONAL_COMPONENTS)
            if (${PKG}_${_BFP_CMP}_FOUND)
              list (APPEND ${PKG}_COMPONENTS_FOUND ${_BFP_CMP})
            endif ()
          endforeach ()
          if (${PKG}_COMPONENTS_FOUND)
            list (REMOVE_DUPLICATES ${PKG}_COMPONENTS_FOUND)
          endif ()
        endif ()
        # if previously this package or components of it where found and the
        # re-discovery of the package or additional components is only optional,
        # set <PKG>_FOUND to TRUE again
        if (_BFP_FOUND AND NOT _BFP_ARGN_REQUIRED AND NOT WITH_${PKG})
          set (${PKG}_FOUND TRUE)
        endif ()
        unset (_BFP_FOUND)
      endif ()
    endif ()
  endif ()
  # --------------------------------------------------------------------------
  # unset locally used variables
  foreach (_BFP_CMP IN LISTS _BFP_ARGN_COMPONENTS _BFP_ARGN_OPTIONAL_COMPONENTS)
    unset (${PKG}_FIND_REQUIRED_${_BFP_CMP})
  endforeach ()
  foreach (_BFP_VAR IN LISTS _BFP_OPTIONS _BFP_MULTI_ARGS)
    unset (_BFP_ARGN_${_BFP_VAR})
  endforeach ()
  unset (_BFP_COMPONENTS)
  unset (_BFP_ARGN_UNPARSED_ARGUMENTS)
  unset (_BFP_OPTIONS)
  unset (_BFP_MULTI_ARGS)
  unset (_BFP_IS_MODULE)
  unset (_BFP_IS_PROJECT)
  unset (_BFP_VERSION)
  unset (_BFP_VERSIONS)
  unset (_BFP_VERSIONS_COUNT)
  unset (_BFP_PKGIDX)
  unset (_BFP_CMPIDX)
  unset (_BFP_CMP)
  unset (_BFP_VAR)
  unset (_BFP_VARS)
  unset (_BFP_CACHED)
  unset (_BFP_TYPE)
  unset (_BFP_PREFIX)
  unset (_BFP_PKG_DIR_VARS)
  unset (_BFP_FIND_COMPONENTS)
  unset (_BFP_FIND_OPTIONAL_COMPONENTS)
  unset (_BFP_NO_FIND_PACKAGE)
  unset (_BFP_NS)
  unset (PKG)
  unset (PKG_L)
  unset (PKG_U)
endmacro ()

# ----------------------------------------------------------------------------
## @brief Use found package.
#
# This macro includes the package's use file if the variable @c &lt;Pkg&gt;_USE_FILE
# is defined. Otherwise, it adds the include directories to the search path
# for include paths if possible. Therefore, the corresponding package
# configuration file has to set the proper CMake variables, i.e.,
# either @c &lt;Pkg&gt;_INCLUDES, @c &lt;Pkg&gt;_INCLUDE_DIRS, or @c &lt;Pkg&gt;_INCLUDE_DIR.
#
# If the given package name is the name of another module of this project
# (the top-level project), this function includes the use file of the specified
# module.
#
# @note As some packages still use all captial variables instead of ones
#       prefixed by a string that follows the same capitalization as the
#       package's name, this function also considers these if defined instead.
#       Hence, if @c &lt;PKG&gt;_INCLUDES is defined, but not @c &lt;Pkg&gt;_INCLUDES, it
#       is used in place of the latter.
#
# @note According to an email on the CMake mailing list, it is not a good idea
#       to use basis_link_directories() any more given that the arguments to
#       basis_target_link_libraries() are absolute paths to the library files.
#       Therefore, this code is commented and not used. It remains here as a
#       reminder only.
#
# @param [in] PACKAGE Name of other package. Optionally, the package name
#                     can include a version specification as suffix which
#                     is separated by the package name using a dash (-), i.e.,
#                     &lt;Package&gt;[-major[.minor[.patch[.tweak]]]].
#                     A version specification is simply ignored by this macro.
#
# @ingroup CMakeAPI
macro (basis_use_package PACKAGE)
  # tokenize package specification
  basis_tokenize_dependency ("${PACKAGE}" PKG VER CMPS)
  # use package
  foreach (A IN ITEMS "WORKAROUND FOR NOT BEING ABLE TO USE RETURN")
    if (BASIS_DEBUG)
      message ("** basis_use_package()")
      message ("**     Package: ${PKG}")
    endif ()
    if (PROJECT_IS_MODULE)
      # ignore BASIS as module dependency
      # important if BASIS itself is a project module
      if ("^${PKG}$" STREQUAL "^BASIS$")
        if (BASIS_DEBUG)
          message ("**     Ignoring BASIS dependency as it clearly is used already by the top-level project.")
        endif ()
        break ()
      # allow modules to specify top-level project as dependency
      elseif ("^${PKG}$" STREQUAL "^${TOPLEVEL_PROJECT_NAME}$")
        if (CMPS)
          if (BASIS_DEBUG)
            message ("**     These are other modules of the top-level project.")
          endif ()
          foreach (CMP IN LISTS CMPS)
            list (FIND PROJECT_MODULES "${CMP}" CMPIDX)
            if (CMPIDX EQUAL -1)
              message (FATAL_ERROR "Module ${PROJECT_NAME} has module ${CMP} of project ${TOPLEVEL_PROJECT_NAME}"
                                   " as dependency, but no such module exists.")
            endif ()
            list (FIND PROJECT_MODULES_ENABLED "${CMP}" CMPIDX)
            if (NOT CMPIDX EQUAL -1)
              if (BASIS_DEBUG)
                message ("**     Include package use file of module ${CMP}.")
              endif ()
              include ("${${${CMP}_CONFIG_PREFIX}_USE_FILE}")
              if (PROJECT_PACKAGE_NAME)
                add_definitions(-DHAVE_${PROJECT_PACKAGE_NAME}_${CMP})
              else ()
                add_definitions(-DHAVE_${TOPLEVEL_PROJECT_NAME}_${CMP})
              endif ()
            endif ()
          endforeach ()
          unset (CMPIDX)
          unset (CMP)
        else ()
          if (BASIS_DEBUG)
            message ("**     This is the top-level project.")
          endif ()
        endif ()
        break () # instead of return()
      # use other module of top-level project
      else ()
        list (FIND PROJECT_MODULES "${PKG}" PKGIDX)
        if (NOT PKGIDX EQUAL -1 AND ${PKG}_FOUND)
          if (BASIS_DEBUG)
            message ("**     Include package use file of other module.")
          endif ()
          include ("${${PKG}_USE_FILE}")
          add_definitions(-DHAVE_${PKG})
          unset (PKGIDX)
          break () # instead of return()
        endif ()
        unset (PKGIDX)
      endif ()
    endif ()
    # if this package is an external project, i.e., a project build as part
    # of the same superbuild as this project, set BUNDLE_PROJECT to TRUE.
    # it is used by (basis_)link_directories() and add_library() to mark
    # the imported link directories and target as belonging to the same
    # installation. this is in particular important for the RPATH settings.
    # whether this package is an external project or not, is decided by the
    # BUNDLE_PROJECTS variable which must be set using the -D option of
    # cmake to a list naming all the other packages which are part of the
    # superbuild.
    if (BUNDLE_PROJECTS)
      list (FIND BUNDLE_PROJECTS "${PKG}" PKGIDX)
      if (PKGIDX EQUAL -1)
        set (BUNDLE_PROJECT FALSE)
      else ()
        set (BUNDLE_PROJECT TRUE)
      endif ()
      unset (PKGIDX)
    endif ()
    # use external package
    if (${PKG}_FOUND)
      # use package only if basis_use_package() not invoked before
      if (BASIS_USE_${PKG}_INCLUDED)
        if (BASIS_DEBUG)
          message ("**     External package used before already.")
        endif ()
        break ()
      endif ()
      if (${PKG}_USE_FILE)
        if (BASIS_DEBUG)
          message ("**     Include package use file of external package.")
        endif ()
        if ("^${PKG}$" STREQUAL "^BASIS$")
          include ("${${PKG}_USE_FILE}" NO_POLICY_SCOPE)
        else ()
          include ("${${PKG}_USE_FILE}")
        endif ()
      else ()
        if (BASIS_DEBUG)
          message ("**     Use variables which were set by basis_find_package().")
        endif ()
        # OpenCV
        if ("^${PKG}$" STREQUAL "^OpenCV$")
          # the cv.h may be found as part of PerlLibs, the include path of
          # which is added at first by BASISConfig.cmake
          if (OpenCV_INCLUDE_DIRS)
            basis_include_directories (BEFORE ${OpenCV_INCLUDE_DIRS})
          elseif (OpenCV_INCLUDE_DIR)
            basis_include_directories (BEFORE ${OpenCV_INCLUDE_DIR})
          endif ()
        # generic
        else ()
          if (${PKG}_INCLUDE_DIRS)
            basis_include_directories (${${PKG}_INCLUDE_DIRS})
          elseif (${PKG}_INCLUDES)
            basis_include_directories (${${PKG}_INCLUDES})
          elseif (${PKG}_INCLUDE_PATH)
            basis_include_directories (${${PKG}_INCLUDE_PATH})
          elseif (${PKG}_INCLUDE_DIR)
            basis_include_directories (${${PKG}_INCLUDE_DIR})
          endif ()
        endif ()
      endif ()
      add_definitions(-DHAVE_${PKG})
      set (BASIS_USE_${PKG}_INCLUDED TRUE)
    elseif (ARGC GREATER 1 AND "^${ARGV1}$" STREQUAL "^REQUIRED$")
      if (BASIS_DEBUG)
        basis_dump_variables ("${PROJECT_BINARY_DIR}/VariablesAfterFind${PKG}.cmake")
      endif ()
      message (FATAL_ERROR "Package ${PACKAGE} not found!")
    endif ()
    # reset switch that identifies currently imported targets and link directories
    # as belonging to an external project which is part of the same superbuild
    set (BUNDLE_PROJECT FALSE)
  endforeach ()
endmacro ()

# ============================================================================
# basis_get_filename_component / basis_get_relative_path
# ============================================================================

# ----------------------------------------------------------------------------
## @brief Fixes CMake's
#         <a href="http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:get_filename_component">
#         get_filename_component()</a> command.
#
# The get_filename_component() command of CMake returns the entire portion
# after the first period (.) [including the period] as extension. However,
# only the component following the last period (.) [including the period]
# should be considered to be the extension.
#
# @note Consider the use of the basis_get_filename_component() macro as
#       an alias to emphasize that this function is different from CMake's
#       <a href="http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:get_filename_component">
#       get_filename_component()</a> command.
#
# @todo Fix issue http://public.kitware.com/Bug/view.php?id=15743 which
#       affects also basis_get_relative_path.
#
# @param [in,out] ARGN Arguments as accepted by get_filename_component().
#
# @returns Sets the variable named by the first argument to the requested
#          component of the given file path.
#
# @sa http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:get_filename_component
# @sa basis_get_filename_component()
function (get_filename_component)
  if (ARGC LESS 3)
    message (FATAL_ERROR "[basis_]get_filename_component(): At least three arguments required!")
  elseif (ARGC GREATER 4)
    message (FATAL_ERROR "[basis_]get_filename_component(): Too many arguments!")
  endif ()
  list (GET ARGN 0 VAR)
  list (GET ARGN 1 STR)
  list (GET ARGN 2 CMD)
  if (CMD MATCHES "^EXT")
    _get_filename_component (${VAR} "${STR}" ${CMD})
    string (REGEX MATCHALL "\\.[^.]*" PARTS "${${VAR}}")
    list (LENGTH PARTS LEN)
    if (LEN GREATER 1)
      math (EXPR LEN "${LEN} - 1")
      list (GET PARTS ${LEN} ${VAR})
    endif ()
  elseif (CMD MATCHES "NAME_WE")
    _get_filename_component (${VAR} "${STR}" NAME)
    string (REGEX REPLACE "\\.[^.]*$" "" ${VAR} ${${VAR}})
  else ()
    _get_filename_component (${VAR} "${STR}" ${CMD})
  endif ()
  if (ARGC EQUAL 4)
    if (NOT "^${ARGV3}$" STREQUAL "^CACHE$")
      message (FATAL_ERROR "[basis_]get_filename_component(): Invalid fourth argument: ${ARGV3}!")
    else ()
      set (${VAR} "${${VAR}}" CACHE STRING "")
    endif ()
  else ()
    set (${VAR} "${${VAR}}" PARENT_SCOPE)
  endif ()
endfunction ()

# ----------------------------------------------------------------------------
## @brief Alias for the overwritten get_filename_component() function.
#
# @sa get_filename_component()
#
# @ingroup CMakeAPI
macro (basis_get_filename_component)
  get_filename_component (${ARGN})
endmacro ()

# ----------------------------------------------------------------------------
## @brief Get path relative to a given base directory.
#
# Unlike the file(RELATIVE_PATH ...) command of CMake which if @p PATH and
# @p BASE are the same directory returns an empty string, this function
# returns a dot (.) in this case instead.
#
# @param [out] REL  @c PATH relative to @c BASE.
# @param [in]  BASE Path of base directory. If a relative path is given, it
#                   is made absolute using basis_get_filename_component()
#                   with ABSOLUTE as last argument.
# @param [in]  PATH Absolute or relative path. If a relative path is given
#                   it is made absolute using basis_get_filename_component()
#                   with ABSOLUTE as last argument.
#
# @returns Sets the variable named by the first argument to the relative path.
#
# @sa http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:file
#
# @ingroup CMakeAPI
function (basis_get_relative_path REL BASE PATH)
  if (BASE MATCHES "^$")
    message (FATAL_ERROR "Empty string given where (absolute) base directory path expected!")
  endif ()
  if (PATH MATCHES "^$")
    set (PATH ".")
  endif ()
  # Attention: http://public.kitware.com/Bug/view.php?id=15743
  basis_get_filename_component (PATH "${PATH}" ABSOLUTE)
  basis_get_filename_component (BASE "${BASE}" ABSOLUTE)
  if (NOT PATH)
    message (FATAL_ERROR "basis_get_relative_path(): No PATH given!")
  endif ()
  if (NOT BASE)
    message (FATAL_ERROR "basis_get_relative_path(): No BASE given!")
  endif ()
  file (RELATIVE_PATH P "${BASE}" "${PATH}")
  if ("${P}" STREQUAL "")
    set (P ".")
  endif ()
  set (${REL} "${P}" PARENT_SCOPE)
endfunction ()

## @brief Create a string from a list of variables indicating if they are defined and their values.
#
# Useful for debug and user errors, for example:
# @code
# set(VAR1 "I'm a string")
# set(VAR2 2)
# basis_variable_value_status(VAR_INFO_STRING VAR1 VAR2 VAR3)
# message(STATUS ${VAR_INFO_STRING})
# @endcode
#
# @param[out] VAR_INFO_STRING The output string variable that will set with the debug string.
# @param[in]  ARGN            List of variables to be put into a string along with their value.
function(basis_variable_value_status VAR_INFO_STRING)
  set (OUTPUT_STRING)
  foreach (VARIABLE_NAME IN ITEMS ${ARGN})
    if (DEFINED ${VARIABLE_NAME})
      set (OUTPUT_STRING "${OUTPUT_STRING}\n  variable name: ${VARIABLE_NAME}  value: ${${VARIABLE_NAME}}")
    else ()
      set (OUTPUT_STRING "${OUTPUT_STRING}\n  variable name: ${VARIABLE_NAME}  value is not defined")
    endif ()
  endforeach ()
  set (${VAR_INFO_STRING} ${OUTPUT_STRING} PARENT_SCOPE)
endfunction()

## @brief Checks for a list of variables required later in the script.
#
# Produces a clear error message explaining the problem and how to fix it if they are not present.
#
# @code
# basis_variable_check(
#    REQUIRED
#       LIBRARY1_INCLUDE_DIRS
#       LIBRARY2_INCLUDE_DIRS
#       LIBRARY2_LIBRARIES
#    OPTIONAL
#       LIBRARY3_INCLUDE_DIRS
#       LIBRARY3_LIBRARIES
#    OPTIONAL_PATH
#       
# )
# @endcode
#
# @param [in] ARGN This argument list is parsed and the following
#                  arguments are extracted.
# @par
# <table border="0">
#   <tr>
#     @tp @b REQUIRED var... @endtp
#     <td>List of variables that MUST be set to run this script correctly. 
#         Will produce a FATAL_ERROR message explaining which variables 
#         are misisng and exit the cmake script.</td>
#   </tr>
#   <tr>
#     @tp @b OPTIONAL var... @endtp
#     <td>List of variables need not be set to run this script correctly.</td>
#   </tr>
#   <tr>
#     @tp @b PATH_EXISTS var... @endtp
#     <td>List of path variables that MUST be set to a location that exists.</td>
#   </tr>
#   <tr>
#     @tp @b OPTIONAL_PATH_EXISTS var... @endtp
#     <td>List of path variables that are optional, but once set must be empty
#         or provide a path to location that exists.</td>
#   </tr>
# </table>
function(basis_variable_check)

  set(options ) # currently none
  set(oneValueArgs ) # currently none
  set(multiValueArgs REQUIRED OPTIONAL PATH_EXISTS OPTIONAL_PATH_EXISTS )
  cmake_parse_arguments(VARIABLE_CONFIGURATION "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )
  
  
  #-------------------------------------------
  # Create the error strings for missing REQUIRED variables
  set(MISSING_REQUIRED_VARIABLES)
  foreach(VARIABLE_NAME IN LISTS VARIABLE_CONFIGURATION_REQUIRED)
     if(NOT ${VARIABLE_NAME})
        list(APPEND MISSING_REQUIRED_VARIABLES ${VARIABLE_NAME})
     endif()
  endforeach()
  if(MISSING_REQUIRED_VARIABLES)
    basis_variable_value_status(MISSING_REQUIRED_VARIABLES_STATUS ${MISSING_REQUIRED_VARIABLES})
    set(MISSING_VARIABLES_ERROR "\nThe following variables are marked as REQUIRED but they are not set to a valid value. Please define the variables correctly in your cmake script or on the command line using -D. ${MISSING_REQUIRED_VARIABLES_STATUS}")
  endif(MISSING_REQUIRED_VARIABLES)
  
  
  #-------------------------------------------
  # Create and print the warning strings for missing OPTIONAL variables
  set(MISSING_OPTIONAL_VARIABLES)
  foreach(VARIABLE_NAME IN LISTS VARIABLE_CONFIGURATION_OPTIONAL)
     if(NOT ${VARIABLE_NAME})
        list(APPEND MISSING_OPTIONAL_VARIABLES ${VARIABLE_NAME})
     endif()
  endforeach()
  if(MISSING_OPTIONAL_VARIABLES)
    basis_variable_value_status(MISSING_OPTIONAL_VARIABLES_STATUS ${MISSING_OPTIONAL_VARIABLES})
    set(MISSING_VARIABLES_WARNING "\nThe following variables are marked as OPTIONAL but they are not set to a valid value. Please define the variables correctly in your cmake script or on the command line using -D. ${MISSING_OPTIONAL_VARIABLES_STATUS}")
    message(AUTHOR_WARNING "${MISSING_VARIABLES_WARNING}")
  endif(MISSING_OPTIONAL_VARIABLES)
  
  
  
  #-------------------------------------------
  # Create the error strings for missing or nonexistant REQUIRED PATH variables
  set(MISSING_PATH_EXISTS)
  foreach(VARIABLE_NAME IN LISTS VARIABLE_CONFIGURATION_PATH_EXISTS)
     if(NOT DEFINED ${VARIABLE_NAME} OR NOT EXISTS ${${VARIABLE_NAME}})
        list(APPEND MISSING_PATH_EXISTS ${VARIABLE_NAME})
     endif()
  endforeach()
  if(MISSING_PATH_EXISTS)
    basis_variable_value_status(MISSING_PATH_EXISTS_STATUS ${MISSING_PATH_EXISTS})
    set(PATH_EXISTS_ERROR "\nThe following PATH variables are marked as REQUIRED but they are not set to a valid location. Please define the variables correctly in your cmake script or on the command line using -D. ${MISSING_PATH_EXISTS_STATUS}")
  endif(MISSING_PATH_EXISTS)
  
  
  
  #-------------------------------------------
  # Create the error strings for missing or nonexistant OPTIONAL PATH variables
  set(MISSING_OPTIONAL_PATH_EXISTS)
  foreach(VARIABLE_NAME IN LISTS VARIABLE_CONFIGURATION_OPTIONAL_PATH_EXISTS)
     if(DEFINED ${VARIABLE_NAME} AND NOT ${${VARIABLE_NAME}} STREQUAL "" AND NOT EXISTS ${${VARIABLE_NAME}})
       # add VARIABLE_NAME to error list if a nonempty path is defined, but does not point to a real location
        list(APPEND MISSING_OPTIONAL_PATH_EXISTS ${VARIABLE_NAME})
     endif()
  endforeach()
  if(MISSING_OPTIONAL_PATH_EXISTS)
    basis_variable_value_status(MISSING_OPTIONAL_PATH_EXISTS_STATUS ${MISSING_OPTIONAL_PATH_EXISTS})
    #debug:
    #message(STATUS "MISSING_OPTIONAL_PATH_EXISTS: ${MISSING_OPTIONAL_PATH_EXISTS}")
    #message(STATUS "MISSING_OPTIONAL_PATH_EXISTS_STATUS: ${MISSING_OPTIONAL_PATH_EXISTS_STATUS}")
    set(OPTIONAL_PATH_EXISTS_ERROR "\nThe following PATH variables are marked as OPTIONAL but they are not set to a valid location. Please define the variables correctly in your cmake script or on the command line using -D. ${MISSING_OPTIONAL_PATH_EXISTS_STATUS}")
  endif(MISSING_OPTIONAL_PATH_EXISTS)
  
  #-------------------------------------------
  # Print all fatal errors at once, because the script will halt
  if(MISSING_VARIABLES_ERROR OR PATH_EXISTS_ERROR OR OPTIONAL_PATH_EXISTS_ERROR)
    message(FATAL_ERROR "${MISSING_VARIABLES_ERROR}${PATH_EXISTS_ERROR}${OPTIONAL_PATH_EXISTS_ERROR}\n")
  endif(MISSING_VARIABLES_ERROR OR PATH_EXISTS_ERROR OR OPTIONAL_PATH_EXISTS_ERROR)

endfunction(basis_variable_check)



# ============================================================================
# name / version
# ============================================================================

# ----------------------------------------------------------------------------
## @brief Convert string to lowercase only or mixed case.
#
# Strings in all uppercase or all lowercase are converted to all lowercase
# letters because these are usually used for acronymns. All other strings
# are returned unmodified with the one exception that the first letter has
# to be uppercase for mixed case strings.
#
# This function is in particular used to normalize the project name for use
# in installation directory paths and namespaces.
#
# @param [out] OUT String in CamelCase.
# @param [in]  STR String.
function (basis_normalize_name OUT STR)
  # strings in all uppercase or all lowercase such as acronymns are an
  # exception and shall be converted to all lowercase instead
  string (TOLOWER "${STR}" L)
  string (TOUPPER "${STR}" U)
  if ("${STR}" STREQUAL "${L}" OR "${STR}" STREQUAL "${U}")
    set (${OUT} "${L}" PARENT_SCOPE)
  # change first letter to uppercase
  else ()
    string (SUBSTRING "${U}"   0  1 A)
    string (SUBSTRING "${STR}" 1 -1 B)
    set (${OUT} "${A}${B}" PARENT_SCOPE)
  endif ()
endfunction ()

# ----------------------------------------------------------------------------
## @brief Extract version numbers from version string.
#
# @param [in]  VERSION Version string in the format "MAJOR[.MINOR[.PATCH]]".
# @param [out] MAJOR   Major version number if given or 0.
# @param [out] MINOR   Minor version number if given or 0.
# @param [out] PATCH   Patch number if given or 0.
#
# @returns See @c [out] parameters.
function (basis_version_numbers VERSION MAJOR MINOR PATCH)
  if (VERSION MATCHES "([0-9]+)(\\.[0-9]+)?(\\.[0-9]+)?(rc[1-9][0-9]*|[a-z]+)?")
    if (CMAKE_MATCH_1)
      set (VERSION_MAJOR ${CMAKE_MATCH_1})
    else ()
      set (VERSION_MAJOR 0)
    endif ()
    if (CMAKE_MATCH_2)
      set (VERSION_MINOR ${CMAKE_MATCH_2})
      string (REGEX REPLACE "^\\." "" VERSION_MINOR "${VERSION_MINOR}")
    else ()
      set (VERSION_MINOR 0)
    endif ()
    if (CMAKE_MATCH_3)
      set (VERSION_PATCH ${CMAKE_MATCH_3})
      string (REGEX REPLACE "^\\." "" VERSION_PATCH "${VERSION_PATCH}")
    else ()
      set (VERSION_PATCH 0)
    endif ()
  else ()
    set (VERSION_MAJOR 0)
    set (VERSION_MINOR 0)
    set (VERSION_PATCH 0)
  endif ()
  set ("${MAJOR}" "${VERSION_MAJOR}" PARENT_SCOPE)
  set ("${MINOR}" "${VERSION_MINOR}" PARENT_SCOPE)
  set ("${PATCH}" "${VERSION_PATCH}" PARENT_SCOPE)
endfunction ()

# ============================================================================
# set
# ============================================================================

# ----------------------------------------------------------------------------
## @brief Set flag given mutually exclusive
#         ARGN_&lt;FLAG&gt; and ARGN_NO&lt;FLAG&gt; function arguments.
#
# @param [in]  PREFIX  Prefix of function arguments. Set to the first argument
#                      of the CMAKE_PARSE_ARGUMENTS() command.
# @param [out] FLAG    Name of flag.
# @param [in]  DEFAULT Default flag value if neither <tt>ARGN_&lt;FLAG;gt;</tt>
#                      nor <tt>ARGN_NO&lt;FLAG;gt;</tt> evaluates to true.
macro (basis_set_flag PREFIX FLAG DEFAULT)
  if (${PREFIX}_${FLAG} AND ${PREFIX}_NO${FLAG})
    message (FATAL_ERROR "Options ${FLAG} and NO${FLAG} are mutually exclusive!")
  endif ()
  if (${PREFIX}_${FLAG})
    set (${FLAG} TRUE)
  elseif (${PREFIX}_NO${FLAG})
    set (${FLAG} FALSE)
  else ()
    set (${FLAG} ${DEFAULT})
  endif ()
endmacro ()

# ----------------------------------------------------------------------------
## @brief Determine if cache entry exists.
#
# @param [out] VAR   Name of boolean result variable.
# @param [in]  ENTRY Name of cache entry.
macro (basis_is_cached VAR ENTRY)
  if (DEFINED ${ENTRY})
    get_property (${VAR} CACHE ${ENTRY} PROPERTY TYPE SET)
  else ()
    set (${VAR} FALSE)
  endif ()
endmacro ()

# ----------------------------------------------------------------------------
## @brief Set type of variable.
#
# If the variable is cached, the type is updated, otherwise, a cache entry
# of the given type with the current value of the variable is added.
#
# @param [in] VAR  Name of variable.
# @param [in] TYPE Desired type of variable.
# @param [in] ARGN Optional DOC string used if variable was not cached before.
macro (basis_set_or_update_type VAR TYPE)
  basis_is_cached (_CACHED ${VAR})
  if (_CACHED)
    set_property (CACHE ${VAR} PROPERTY TYPE ${TYPE})
  else ()
    set (${VAR} "${${VAR}}" CACHE ${TYPE} "${ARGN}" FORCE)
  endif ()
  unset (_CACHED)
endmacro ()

# ----------------------------------------------------------------------------
## @brief Change type of cached variable.
#
# If the variable is not cached, nothing is done.
macro (basis_update_type_of_variable VAR TYPE)
  basis_is_cached (_CACHED ${VAR})
  if (_CACHED)
    set_property (CACHE ${VAR} PROPERTY TYPE ${TYPE})
  endif ()
  unset (_CACHED)
endmacro ()

# ----------------------------------------------------------------------------
## @brief Set variable value.
#
# If the variable is cached, this function will update the cache value,
# otherwise, it simply sets the CMake variable uncached to the given value(s).
macro (basis_set_or_update_value VAR)
  basis_is_cached (_CACHED ${VAR})
  if (_CACHED)
    if (ARGC GREATER 1)
      set_property (CACHE ${VAR} PROPERTY VALUE ${ARGN})
    else ()
      set (${VAR} "" CACHE INTERNAL "" FORCE)
    endif ()
  else ()
    set (${VAR} ${ARGN})
  endif ()
  unset (_CACHED)
endmacro ()

# ----------------------------------------------------------------------------
## @brief Update cache variable.
macro (basis_update_value VAR)
  basis_is_cached (_CACHED ${VAR})
  if (_CACHED)
    set_property (CACHE ${VAR} PROPERTY VALUE ${ARGN})
  endif ()
  unset (_CACHED)
endmacro ()

# ----------------------------------------------------------------------------
## @brief Set value of variable only if variable is not set already.
#
# @param [out] VAR  Name of variable.
# @param [in]  ARGN Arguments to set() command excluding variable name.
#
# @returns Sets @p VAR if its value was not valid before.
macro (basis_set_if_empty VAR)
  if (NOT ${VAR})
    set (${VAR} ${ARGN})
  endif ()
endmacro ()

# ----------------------------------------------------------------------------
## @brief Set value of variable only if variable is not defined yet.
#
# @param [out] VAR  Name of variable.
# @param [in]  ARGN Arguments to set() command excluding variable name.
#
# @returns Sets @p VAR if it was not defined before.
macro (basis_set_if_not_set VAR)
  if (NOT DEFINED "${VAR}")
    set ("${VAR}" ${ARGN})
  endif ()
endmacro ()

# ----------------------------------------------------------------------------
## @brief Set value of variable to either 0 or 1 based on option value
#
# This function can be used to convert option values from TRUE/ON to 1 and
# FALSE/OFF to 0 such that they can be used to configure a config.h.in header.
#
# @param [out] VAR Name of configuration header variable.
# @param [in]  OPT Value of CMake option.
macro (basis_set_config_option VAR OPT)
  if (${OPT})
    set ("${VAR}" 1)
  else ()
    set ("${VAR}" 0)
  endif ()
endmacro ()

# ----------------------------------------------------------------------------
## @brief Set path relative to script file.
#
# This function can be used in script configurations. It takes a variable
# name and a path as input arguments. If the given path is relative, it makes
# it first absolute using @c PROJECT_SOURCE_DIR. Then the path is made
# relative to the directory of the built script file. A CMake variable of the
# given name is set to the specified relative path. Optionally, a third
# argument, the path used for building the script for the install tree
# can be passed as well. If a relative path is given as this argument,
# it is made absolute by prefixing it with @c CMAKE_INSTALL_PREFIX instead.
#
# @note This function may only be used in script configurations such as
#       in particular the ScriptConfig.cmake.in file. It requires that the
#       variables @c __DIR__ and @c BUILD_INSTALL_SCRIPT are set properly.
#       These variables are set by the basis_configure_script() function.
#       Moreover, it makes use of the global @c CMAKE_INSTALL_PREFIX and
#       @c PROJECT_SOURCE_DIR variables.
#
# @param [out] VAR   Name of the variable.
# @param [in]  PATH  Path to directory or file.
# @param [in]  ARGV3 Path to directory or file inside install tree.
#                    If this argument is not given, PATH is used for both
#                    the build and install tree version of the script.
#
# @ingroup CMakeAPI
function (basis_set_script_path VAR PATH)
  if (NOT __DIR__)
    message (FATAL_ERROR "__DIR__ not set! Note that basis_set_script_path() may"
                         " only be used in script configurations (e.g., ScriptConfig.cmake.in).")
  endif ()
  if (ARGC GREATER 3)
    message (FATAL_ERROR "Too many arguments given for function basis_set_script_path()")
  endif ()
  if (ARGC EQUAL 3 AND BUILD_INSTALL_SCRIPT)
    set (PREFIX "${CMAKE_INSTALL_PREFIX}")
    set (PATH   "${ARGV2}")
  else ()
    set (PREFIX "${PROJECT_SOURCE_DIR}")
  endif ()
  if (NOT IS_ABSOLUTE "${PATH}")
    set (PATH "${PREFIX}/${PATH}")
  endif ()
  basis_get_relative_path (PATH "${__DIR__}" "${PATH}")
  if (NOT PATH)
    set (${VAR} "." PARENT_SCOPE)
  else ()
    string (REGEX REPLACE "/+$" "" PATH "${PATH}")
    set (${VAR} "${PATH}" PARENT_SCOPE)
  endif ()
endfunction ()

# ============================================================================
# set/get any property
# ============================================================================

# ----------------------------------------------------------------------------
## @brief Convert list into regular expression.
#
# This function is in particular used to convert a list of property names
# such as &lt;CONFIG&gt;_OUTPUT_NAME, e.g., the list @c BASIS_PROPERTIES_ON_TARGETS,
# into a regular expression which can be used in pattern matches.
#
# @param [out] REGEX Name of variable for resulting regular expression.
# @param [in]  ARGN  List of patterns which may contain placeholders in the
#                    form of "<this is a placeholder>". These are replaced
#                    by the regular expression "[^ ]+".
macro (basis_list_to_regex REGEX)
  string (REGEX REPLACE "<[^>]+>" "[^ ]+" ${REGEX} "${ARGN}")
  string (REGEX REPLACE ";" "|" ${REGEX} "${${REGEX}}")
  set (${REGEX} "^(${${REGEX}})$")
endmacro ()

# ----------------------------------------------------------------------------
## @brief Output current CMake variables to file.
function (basis_dump_variables RESULT_FILE)
  set (DUMP)
  get_cmake_property (VARIABLE_NAMES VARIABLES)
  foreach (V IN LISTS VARIABLE_NAMES)
    if (NOT V MATCHES "^_|^RESULT_FILE$|^ARGC$|^ARGV[0-9]?$|^ARGN_")
      set (VALUE "${${V}}")
      # sanitize value for use in set() command
      string (REPLACE "\\" "\\\\" VALUE "${VALUE}") # escape backspaces
      string (REPLACE "\"" "\\\"" VALUE "${VALUE}") # escape double quotes
      # Escape ${VAR} by \${VAR} such that CMake does not evaluate it.
      # Escape $STR{VAR} by \$STR{VAR} such that CMake does not report a
      # syntax error b/c it expects either ${VAR}, $ENV{VAR}, or $CACHE{VAR}.
      # Escape @VAR@ by \@VAR\@ such that CMake does not evaluate it.
      string (REGEX REPLACE "([^\\])\\\$([^ ]*){" "\\1\\\\\$\\2{" VALUE "${VALUE}")
      string (REGEX REPLACE "([^\\])\\\@([^ ]*)\@" "\\1\\\\\@\\2\\\\\@" VALUE "${VALUE}")
      # append variable to output file
      set (DUMP "${DUMP}set (${V} \"${VALUE}\")\n")
    endif ()
  endforeach ()
  file (WRITE "${RESULT_FILE}" "# CMake variables dump created by BASIS\n${DUMP}")
endfunction ()

# ----------------------------------------------------------------------------
## @brief Write CMake script file which sets the named variable to the
#         specified (list of) values.
function (basis_write_list FILENAME VARIABLE)
  file (WRITE "${FILENAME}" "# Automatically generated. Do not edit this file!\nset (${VARIABLE}\n")
  foreach (V IN LISTS ARGN)
    file (APPEND "${FILENAME}" "  \"${V}\"\n")
  endforeach ()
  file (APPEND "${FILENAME}" ")\n")
endfunction ()

# ----------------------------------------------------------------------------
## @brief Set a named property in a given scope.
#
# This function replaces CMake's
# <a href="http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:set_property">
# set_property()</a> command.
#
# @param [in] SCOPE The argument for the @p SCOPE parameter of
#                   <a href="http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:set_property">
#                   set_property()</a>.
# @param [in] ARGN  Arguments as accepted by.
#                   <a href="http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:set_property">
#                   set_property()</a>.
#
# @returns Sets the specified property.
#
# @sa http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:set_property
#
# @ingroup CMakeAPI
if (BASIS_USE_TARGET_UIDS)
  function (basis_set_property SCOPE)
    if (SCOPE MATCHES "^TARGET$|^TEST$")
      # map target/test names to UIDs
      list (LENGTH ARGN ARGN_LENGTH)
      if (ARGN_LENGTH EQUAL 0)
        message (FATAL_ERROR "basis_set_property(${SCOPE}): Expected arguments after SCOPE!")
      endif ()
      set (IDX 0)
      set (ARG)
      while (IDX LESS ARGN_LENGTH)
        list (GET ARGN ${IDX} ARG)
        if (ARG MATCHES "^APPEND$")
          math (EXPR IDX "${IDX} + 1")
          list (GET ARGN ${IDX} ARG)
          if (NOT ARG MATCHES "^PROPERTY$")
            message (FATAL_ERROR "basis_set_properties(${SCOPE}): Expected PROPERTY keyword after APPEND!")
          endif ()
          break ()
        elseif (ARG MATCHES "^PROPERTY$")
          break ()
        else ()
          if (SCOPE MATCHES "^TEST$")
            basis_get_test_uid (UID "${ARG}")
          else ()
            basis_get_target_uid (UID "${ARG}")
          endif ()
          list (INSERT ARGN ${IDX} "${UID}")
          math (EXPR IDX "${IDX} + 1")
          list (REMOVE_AT ARGN ${IDX}) # after insert to avoid index out of range
        endif ()
      endwhile ()
      if (IDX EQUAL ARGN_LENGTH)
        message (FATAL_ERROR "basis_set_properties(${SCOPE}): Missing PROPERTY keyword!")
      endif ()
      math (EXPR IDX "${IDX} + 1")
      list (GET ARGN ${IDX} ARG)
      # property name matches DEPENDS
      if (ARG MATCHES "DEPENDS")
        math (EXPR IDX "${IDX} + 1")
        while (IDX LESS ARGN_LENGTH)
          list (GET ARGN ${IDX} ARG)
          if (SCOPE MATCHES "^TEST$")
            basis_get_test_uid (UID "${ARG}")
          else ()
            basis_get_target_uid (UID "${ARG}")
          endif ()
          list (INSERT ARGN ${IDX} "${UID}")
          math (EXPR IDX "${IDX} + 1")
          list (REMOVE_AT ARGN ${IDX}) # after insert to avoid index out of range
        endwhile ()
      endif ()
    endif ()
    if (BASIS_DEBUG)
      message ("** basis_set_property():")
      message ("**   Scope:     ${SCOPE}")
      message ("**   Arguments: [${ARGN}]")
    endif ()
    set_property (${SCOPE} ${ARGN})
  endfunction ()
else ()
  macro (basis_set_property)
    set_property (${ARGV})
  endmacro ()
endif ()

# ----------------------------------------------------------------------------
## @brief Get a property.
#
# This function replaces CMake's
# <a href="http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:get_property">
# get_property()</a> command.
#
# @param [out] VAR     Property value.
# @param [in]  SCOPE   The argument for the @p SCOPE argument of
#                      <a href="http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:get_property">
#                      get_property()</a>.
# @param [in]  ELEMENT The argument for the @p ELEMENT argument of
#                      <a href="http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:get_property">
#                      get_property()</a>.
# @param [in]  ARGN    Arguments as accepted by
#                      <a href="http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:get_property">
#                      get_property()</a>.
#
# @returns Sets @p VAR to the value of the requested property.
#
# @sa http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:get_property
#
# @ingroup CMakeAPI
if (BASIS_USE_TARGET_UIDS)
  function (basis_get_property VAR SCOPE ELEMENT)
    if (SCOPE MATCHES "^TARGET$")
      basis_get_target_uid (ELEMENT "${ELEMENT}")
    elseif (SCOPE MATCHES "^TEST$")
      basis_get_test_uid (ELEMENT "${ELEMENT}")
    endif ()
    get_property (VALUE ${SCOPE} ${ELEMENT} ${ARGN})
    set ("${VAR}" "${VALUE}" PARENT_SCOPE)
  endfunction ()
else ()
  macro (basis_get_property)
    get_property (${ARGV})
  endmacro ()
endif ()

# ----------------------------------------------------------------------------
## @brief Set project-global property.
#
# Set property associated with current project/module. The property is in
# fact just a cached variable whose name is prefixed by the project's name.
function (basis_set_project_property)
  CMAKE_PARSE_ARGUMENTS (
    ARGN
      "APPEND"
      "PROJECT"
      "PROPERTY"
    ${ARGN}
  )

  if (NOT ARGN_PROJECT)
    set (ARGN_PROJECT "${PROJECT_NAME}")
  endif ()
  if (NOT ARGN_PROPERTY)
    message (FATAL_ERROR "Missing PROPERTY argument!")
  endif ()

  list (GET ARGN_PROPERTY 0 PROPERTY_NAME)
  list (REMOVE_AT ARGN_PROPERTY 0) # remove property name from values

  if (ARGN_APPEND)
    basis_get_project_property (CURRENT ${ARGN_PROJECT} ${PROPERTY_NAME})
    if (NOT "${CURRENT}" STREQUAL "")
      list (INSERT ARGN_PROPERTY 0 "${CURRENT}")
    endif ()
  endif ()

  set (
    ${ARGN_PROJECT}_${PROPERTY_NAME}
      "${ARGN_PROPERTY}"
    CACHE INTERNAL
      "Property ${PROPERTY_NAME} of project ${ARGN_PROJECT}."
    FORCE
  )
endfunction ()

# ----------------------------------------------------------------------------
## @brief Get project-global property value.
#
# Example:
# @code
# basis_get_project_property(TARGETS)
# basis_get_project_property(TARGETS ${PROJECT_NAME})
# basis_get_project_property(TARGETS ${PROJECT_NAME} TARGETS)
# basis_get_project_property(TARGETS PROPERTY TARGETS)
# @endcode
#
# @param [out] VARIABLE Name of result variable.
# @param [in]  ARGN     See the example uses. The optional second argument
#                       is either the name of the project similar to CMake's
#                       get_target_property() command or the keyword PROPERTY
#                       followed by the name of the property.
function (basis_get_project_property VARIABLE)
  if (ARGC GREATER 3)
    message (FATAL_ERROR "Too many arguments!")
  endif ()
  if (ARGC EQUAL 1)
    set (ARGN_PROJECT "${PROJECT_NAME}")
    set (ARGN_PROPERTY "${VARIABLE}")
  elseif (ARGC EQUAL 2)
    if (ARGV1 MATCHES "^PROPERTY$")
      message (FATAL_ERROR "Expected argument after PROPERTY keyword!")
    endif ()
    set (ARGN_PROJECT  "${ARGV1}")
    set (ARGN_PROPERTY "${VARIABLE}")
  else ()
    if (ARGV1 MATCHES "^PROPERTY$")
      set (ARGN_PROJECT "${PROJECT_NAME}")
    else ()
      set (ARGN_PROJECT  "${ARGV1}")
    endif ()
    set (ARGN_PROPERTY "${ARGV2}")
  endif ()
  set (${VARIABLE} "${${ARGN_PROJECT}_${ARGN_PROPERTY}}" PARENT_SCOPE)
endfunction ()

# ============================================================================
# list / string manipulations
# ============================================================================

# ----------------------------------------------------------------------------
## @brief Sanitize string variable for use in regular expression.
#
# @note This function may not work for all cases, but is used in particular
#       to sanitize project names, target names, namespace identifiers,...
#
#       This takes all of the dollar signs, and other special characters and
#       adds escape characters such as backslash as necessary.
#
# @param [out] OUT String that can be used in regular expression.
# @param [in]  STR String to sanitize.
macro (basis_sanitize_for_regex OUT STR)
  string (REGEX REPLACE "([.+*?^$])" "\\\\\\1" ${OUT} "${STR}")
endmacro ()

# ----------------------------------------------------------------------------
## @brief Concatenates all list elements into a single string.
#
# The list elements are concatenated without any delimiter in between.
# Use basis_list_to_delimited_string() to specify a delimiter such as a
# whitespace character or comma (,) as delimiter.
#
# @param [out] STR  Output string.
# @param [in]  ARGN Input list.
#
# @returns Sets @p STR to the resulting string.
#
# @sa basis_list_to_delimited_string()
function (basis_list_to_string STR)
  set (OUT)
  foreach (ELEM IN LISTS ARGN)
    set (OUT "${OUT}${ELEM}")
  endforeach ()
  set ("${STR}" "${OUT}" PARENT_SCOPE)
endfunction ()

# ----------------------------------------------------------------------------
## @brief Concatenates all list elements into a single delimited string.
#
# @param [out] STR   Output string.
# @param [in]  DELIM Delimiter used to separate list elements.
#                    Each element which contains the delimiter as substring
#                    is surrounded by double quotes (") in the output string.
# @param [in]  ARGN  Input list. If this list starts with the argument
#                    @c NOAUTOQUOTE, the automatic quoting of list elements
#                    which contain the delimiter is disabled.
#
# @returns Sets @p STR to the resulting string.
#
# @see basis_join
#
# @todo consider replacing basis_list_to_delimited_string with basis_join
function (basis_list_to_delimited_string STR DELIM)
  set (OUT)
  set (AUTOQUOTE TRUE)
  if (ARGN)
    list (GET ARGN 0 FIRST)
    if (FIRST MATCHES "^NOAUTOQUOTE$")
      list (REMOVE_AT ARGN 0)
      set (AUTOQUOTE FALSE)
    endif ()
  endif ()
  basis_sanitize_for_regex (DELIM_RE "${DELIM}")
  foreach (ELEM ${ARGN})
    if (OUT)
      set (OUT "${OUT}${DELIM}")
    endif ()
    if (AUTOQUOTE AND ELEM MATCHES "${DELIM_RE}")
      set (OUT "${OUT}\"${ELEM}\"")
    else ()
      set (OUT "${OUT}${ELEM}")
    endif ()
  endforeach ()
  set ("${STR}" "${OUT}" PARENT_SCOPE)
endfunction ()

# ----------------------------------------------------------------------------
## @brief Concatenates all list elements into a single delimited string.
#
# @param [in]  VALUES       Input list string.
# @param [in]  DELIMITER    Delimiter glue used to separate list elements.
#                           Each element which contains the delimiter as substring
#                           is surrounded by double quotes (") in the output string.
# @param [out] OUTPUT       Output string variable name.
#
# @code
#   set( letters "" "\;a" b c "d\;d" )
#   basis_join("${letters}" ":" output)
#   message("${output}") # :;a:b:c:d;d
# @endcode
#
# @returns Sets @p OUTPUT to the resulting string.
#
# @see basis_list_to_delimited_string
function(basis_join VALUES DELIMITER OUTPUT)
  string (REGEX REPLACE "([^\\]|^);" "\\1${DELIMITER}" _TMP_STR "${VALUES}")
  string (REGEX REPLACE "[\\](.)" "\\1" _TMP_STR "${_TMP_STR}") #fixes escaping
  set (${OUTPUT} "${_TMP_STR}" PARENT_SCOPE)
endfunction()

# ----------------------------------------------------------------------------
## @brief Splits a string at space characters into a list.
#
# @todo Probably this can be done in a better way...
#       Difficulty is, that string(REPLACE) does always replace all
#       occurrences. Therefore, we need a regular expression which matches
#       the entire string. More sophisticated regular expressions should do
#       a better job, though.
#
# @param [out] LST  Output list.
# @param [in]  STR  Input string.
#
# @returns Sets @p LST to the resulting CMake list.
function (basis_string_to_list LST STR)
  set (TMP "${STR}")
  set (OUT)
  # 1. extract elements such as "a string with spaces"
  while (TMP MATCHES "\"[^\"]*\"")
    string (REGEX REPLACE "^(.*)\"([^\"]*)\"(.*)$" "\\1\\3" TMP "${TMP}")
    if (OUT)
      set (OUT "${CMAKE_MATCH_2};${OUT}")
    else (OUT)
      set (OUT "${CMAKE_MATCH_2}")
    endif ()
  endwhile ()
  # 2. extract other elements separated by spaces (excluding first and last)
  while (TMP MATCHES " [^\" ]+ ")
    string (REGEX REPLACE "^(.*) ([^\" ]+) (.*)$" "\\1\\3" TMP "${TMP}")
    if (OUT)
      set (OUT "${CMAKE_MATCH_2};${OUT}")
    else (OUT)
      set (OUT "${CMAKE_MATCH_2}")
    endif ()
  endwhile ()
  # 3. extract first and last elements (if not done yet)
  if (TMP MATCHES "^[^\" ]+")
    set (OUT "${CMAKE_MATCH_0};${OUT}")
  endif ()
  if (NOT "${CMAKE_MATCH_0}" STREQUAL "${TMP}" AND TMP MATCHES "[^\" ]+$")
    set (OUT "${OUT};${CMAKE_MATCH_0}")
  endif ()
  # return resulting list
  set (${LST} "${OUT}" PARENT_SCOPE)
endfunction ()

# ----------------------------------------------------------------------------
## @brief Compare two lists.
#
# @param [out] RESULT Result of comparison.
# @param [in]  LIST1  Name of variable holding the first list.
# @param [in]  LIST2  Name of varaible holding the second list.
#
# @retval 0 The two lists are not identical.
# @retval 1 Both lists have identical elements (not necessarily in the same order).
macro (basis_compare_lists RESULT LIST1 LIST2)
  set (_L1 "${${LIST1}}")
  set (_L2 "${${LIST2}}")
  list (SORT _L1)
  list (SORT _L2)
  if ("${_L1}" STREQUAL "${_L2}")
    set (RESULT TRUE)
  else ()
    set (RESULT FALSE)
  endif ()
  unset (_L1)
  unset (_L2)
endmacro ()

# ============================================================================
# name <=> UID
# ============================================================================

# ----------------------------------------------------------------------------
## @brief Derive target name from source file name.
#
# @param [out] TARGET_NAME Target name.
# @param [in]  SOURCE_FILE Source file.
# @param [in]  ARGN        Third argument to get_filename_component().
#                          If not specified, the given path is only sanitized.
#
# @returns Target name derived from @p SOURCE_FILE.
function (basis_get_source_target_name TARGET_NAME SOURCE_FILE)
  # remove ".in" suffix from file name
  string (REGEX REPLACE "\\.in$" "" OUT "${SOURCE_FILE}")
  # get name component
  if (ARGC GREATER 2)
    get_filename_component (OUT "${OUT}" ${ARGV2})
  endif ()
  # replace special characters
  string (REGEX REPLACE "[./\\]" "_" OUT "${OUT}")
  # return
  set (${TARGET_NAME} "${OUT}" PARENT_SCOPE)
endfunction ()

# ----------------------------------------------------------------------------
## @brief Strip of top-level package name from target UID if present.
#
# If @c BASIS_USE_FULLY_QUALIFIED_TARGET_UID is @c ON, the top-level package
# name is always preserved and this operation does nothing.
#
# @param[in,out] TARGET_UID "Global" target name, i.e., actual CMake target name.
#
# @returns Sets @p TARGET_UID to the (stripped) UID.
function (basis_strip_target_uid TARGET_UID)
  if (NOT BASIS_USE_FULLY_QUALIFIED_UIDS)
    basis_sanitize_for_regex (RE "${TOPLEVEL_PROJECT_NAMESPACE_CMAKE}")
    string (REGEX REPLACE "^\\.*${RE}\\." "" UID "${${TARGET_UID}}")
    set ("${TARGET_UID}" "${UID}" PARENT_SCOPE)
  endif ()
endfunction ()

# ----------------------------------------------------------------------------
## @brief Make target UID from given target name.
#
# This function is intended for internal use only.
#
# @param [out] TARGET_UID  "Global" target name, i.e., actual CMake target name.
# @param [in]  TARGET_NAME Target name used as argument to BASIS CMake functions.
#
# @returns Sets @p TARGET_UID to the UID of the build target @p TARGET_NAME.
#
# @sa basis_make_target_uid()
function (_basis_make_target_uid TARGET_UID TARGET_NAME)
  if (TARGET_NAME MATCHES "^\\.+(.*)$")
    set (${TARGET_UID} "${CMAKE_MATCH_1}" PARENT_SCOPE)
  else ()
    set (UID "${PROJECT_NAMESPACE_CMAKE}.${TARGET_NAME}")
    basis_strip_target_uid (UID)
    set (${TARGET_UID} "${UID}" PARENT_SCOPE)
  endif ()
endfunction ()

# ----------------------------------------------------------------------------
## @brief Make target UID from given target name.
#
# This function is intended for use by the basis_add_*() functions only.
#
# If @c BASIS_USE_TARGET_UIDS is set to @c OFF, this operation
# always just sets the @p TARGET_UID to the given @p TARGET_NAME.
#
# @param [out] TARGET_UID  "Global" target name, i.e., actual CMake target name.
# @param [in]  TARGET_NAME Target name used as argument to BASIS CMake functions.
#
# @returns Sets @p TARGET_UID to the UID of the build target @p TARGET_NAME.
#
# @sa basis_get_target_uid()
if (BASIS_USE_TARGET_UIDS)
  function (basis_make_target_uid TARGET_UID TARGET_NAME)
    _basis_make_target_uid (UID "${TARGET_NAME}")
    set ("${TARGET_UID}" "${UID}" PARENT_SCOPE)
  endfunction ()
else ()
  function (basis_make_target_uid TARGET_UID TARGET_NAME)
    if (TARGET_NAME MATCHES "^\\.+(.*)$")
      set ("${TARGET_UID}" "${CMAKE_MATCH_1}" PARENT_SCOPE)
    else ()
      set ("${TARGET_UID}" "${TARGET_NAME}" PARENT_SCOPE)
    endif ()
  endfunction ()
endif ()

# ----------------------------------------------------------------------------
## @brief Get "global" target name, i.e., actual CMake target name.
#
# In order to ensure that CMake target names are unique across modules of
# a BASIS project, the target name given to the BASIS CMake functions is
# converted by basis_make_target_uid() into a so-called target UID which is
# used as actual CMake target name. This function can be used to get for a
# given target name or UID the closest match of a known target UID.
#
# The individual parts of the target UID, i.e, package name,
# module name, and target name are separated by a dot (.).
# If @c BASIS_USE_FULLY_QUALIFIED_UIDS is set to @c OFF, the common part of
# all target UIDs is removed by this function from the target UID.
# When the target is exported, however, this common part will be
# prefixed again. This is done by the basis_export_targets() function.
#
# Note that names of imported targets are not prefixed in any case.
#
# The counterpart basis_get_target_name() can be used to convert the target UID
# back to the target name without namespace prefix.
#
# If @c BASIS_USE_TARGET_UIDS is set to @c OFF, this operation
# always just sets the @p TARGET_UID to the given @p TARGET_NAME.
#
# @note At the moment, BASIS does not support modules which themselves have
#       modules again. This would require a more nested namespace hierarchy
#       and makes things unnecessarily complicated.
#
# @param [out] TARGET_UID  "Global" target name, i.e., actual CMake target name.
# @param [in]  TARGET_NAME Target name used as argument to BASIS CMake functions.
#
# @returns Sets @p TARGET_UID to the UID of the build target @p TARGET_NAME.
#
# @sa basis_get_target_name()
if (BASIS_USE_TARGET_UIDS)
  function (basis_get_target_uid TARGET_UID TARGET_NAME)
    # in case of a leading namespace separator, do not modify target name
    if (TARGET_NAME MATCHES "^\\.+(.*)$")
      set (UID "${CMAKE_MATCH_1}")
    # otherwise,
    else ()
      set (UID "${TARGET_NAME}")
      basis_sanitize_for_regex (BASE_RE "${TOPLEVEL_PROJECT_NAMESPACE_CMAKE}")
      # try prepending namespace or parts of it until target is known,
      # first assuming the simplified UIDs without the common prefix
      # of this package which applies to targets of this package
      if (NOT BASIS_USE_FULLY_QUALIFIED_UIDS AND NOT TARGET "${UID}")
        string (REGEX REPLACE "^${BASE_RE}\\." "" PREFIX "${PROJECT_NAMESPACE_CMAKE}")
        while (PREFIX)
          if (TARGET "${PREFIX}.${TARGET_NAME}")
            set (UID "${PREFIX}.${TARGET_NAME}")
            break ()
          else ()
            if (PREFIX MATCHES "(.*)\\.[^.]+")
              set (PREFIX "${CMAKE_MATCH_1}")
            else ()
              break ()
            endif ()
          endif ()
        endwhile ()
      endif ()
      # and then with the fully qualified UIDs for imported targets
      if (NOT TARGET "${UID}")
        set (PREFIX "${PROJECT_NAMESPACE_CMAKE}")
        while (PREFIX)
          if (TARGET "${PREFIX}.${TARGET_NAME}")
            set (UID "${PREFIX}.${TARGET_NAME}")
            break ()
          else ()
            if (PREFIX MATCHES "(.*)\\.[^.]+")
              set (PREFIX "${CMAKE_MATCH_1}")
            else ()
              break ()
            endif ()
          endif ()
        endwhile ()
      endif ()
    endif ()
    # strip off top-level namespace part (optional)
    if (NOT BASIS_USE_FULLY_QUALIFIED_UIDS)
      string (REGEX REPLACE "^${BASE_RE}\\." "" UID "${UID}")
    endif ()
    # return
    set ("${TARGET_UID}" "${UID}" PARENT_SCOPE)
  endfunction ()
else ()
  function (basis_get_target_uid TARGET_UID TARGET_NAME)
    if (TARGET_NAME MATCHES "^\\.+(.*)$")
      set (${TARGET_UID} "${CMAKE_MATCH_1}" PARENT_SCOPE)
    else ()
      set ("${TARGET_UID}" "${TARGET_NAME}" PARENT_SCOPE)
    endif ()
  endfunction ()
endif ()

# ----------------------------------------------------------------------------
## @brief Get fully-qualified target name.
#
# This function always returns a fully-qualified target UID, no matter if
# the option @c BASIS_USE_TARGET_UIDS or @c BASIS_USE_FULLY_QUALIFIED_UIDS
# is @c OFF. Note that if @c BASIS_USE_FULLY_QUALIFIED_UIDS is @c ON, the
# returned target UID may not be the actual name of a CMake target.
#
# @param [out] TARGET_UID  Fully-qualified target UID.
# @param [in]  TARGET_NAME Target name used as argument to BASIS CMake functions.
#
# @sa basis_get_target_uid()
function (basis_get_fully_qualified_target_uid TARGET_UID TARGET_NAME)
  if (BASIS_USE_TARGET_UIDS)
    basis_get_target_uid (UID "${TARGET_NAME}")
    if (NOT BASIS_USE_FULLY_QUALIFIED_UIDS)
      if (TARGET "${UID}")
        get_target_property (IMPORTED "${UID}" IMPORTED)
        if (NOT IMPORTED)
          set (UID "${TOPLEVEL_PROJECT_NAMESPACE_CMAKE}.${UID}")
        endif ()
      else ()
        set (UID "${TOPLEVEL_PROJECT_NAMESPACE_CMAKE}.${UID}")
      endif ()
    endif ()
  else ()
    if (TARGET "${TARGET_NAME}")
      get_target_property (IMPORTED "${TARGET_NAME}" IMPORTED)
    else ()
      set (IMPORTED FALSE)
    endif ()
    if (IMPORTED)
      set (UID "${TARGET_NAME}")
    else ()
      set (UID "${PROJECT_NAMESPACE_CMAKE}.${TARGET_NAME}")
    endif ()
  endif ()
  set (${TARGET_UID} "${UID}" PARENT_SCOPE)
endfunction ()

# ----------------------------------------------------------------------------
## @brief Get namespace of build target without check of UID.
#
# If @c BASIS_USE_TARGET_UIDS is set to @c OFF, this operation
# always just sets the @p TARGET_NS to an empty string.
#
# @param [out] TARGET_NS  Namespace part of target UID.
# @param [in]  TARGET_UID Target UID.
if (BASIS_USE_TARGET_UIDS)
  function (_basis_get_target_namespace TARGET_NS TARGET_UID)
    if (UID MATCHES "^(.*)\\.")
      set ("${TARGET_NS}" "${CMAKE_MATCH_1}" PARENT_SCOPE)
    else ()
      set ("${TARGET_NS}" "" PARENT_SCOPE)
    endif ()
  endfunction ()
else ()
  function (_basis_get_target_namespace TARGET_NS TARGET_UID)
    set ("${TARGET_NS}" "" PARENT_SCOPE)
  endfunction ()
endif ()

# ----------------------------------------------------------------------------
## @brief Get namespace of build target.
#
# If @c BASIS_USE_TARGET_UIDS is set to @c OFF, this operation
# always just sets the @p TARGET_NS to an empty string.
#
# @param [out] TARGET_NS  Namespace part of target UID.
# @param [in]  TARGET_UID Target UID/name.
if (BASIS_USE_TARGET_UIDS)
  function (basis_get_target_namespace TARGET_NS TARGET_UID)
    basis_get_fully_qualified_target_uid (UID "${TARGET_UID}")
    _basis_get_target_namespace (NS "${UID}")
    set ("${TARGET_NS}" "${NS}" PARENT_SCOPE)
  endfunction ()
else ()
  function (basis_get_target_namespace TARGET_NS TARGET_UID)
    set ("${TARGET_NS}" "" PARENT_SCOPE)
  endfunction ()
endif ()

# ----------------------------------------------------------------------------
## @brief Get "local" target name, i.e., BASIS target name without check of UID.
#
# If @c BASIS_USE_TARGET_UIDS is set to @c OFF, this operation
# always just sets the @p TARGET_NAME to the given @p TARGET_UID.
#
# @param [out] TARGET_NAME Target name used as argument to BASIS functions.
# @param [in]  TARGET_UID  "Global" target name, i.e., actual CMake target name.
#
# @returns Sets @p TARGET_NAME to the name of the build target with UID @p TARGET_UID.
#
# @sa basis_get_target_name(), basis_get_target_uid()
if (BASIS_USE_TARGET_UIDS)
  function (_basis_get_target_name TARGET_NAME TARGET_UID)
    # strip off namespace of current project
    basis_sanitize_for_regex (RE "${PROJECT_NAMESPACE_CMAKE}")
    string (REGEX REPLACE "^${RE}\\." "" NAME "${UID}")
    # return
    set (${TARGET_NAME} "${NAME}" PARENT_SCOPE)
  endfunction ()
else ()
  function (_basis_get_target_name TARGET_NAME TARGET_UID)
    set (${TARGET_NAME} "${TARGET_UID}" PARENT_SCOPE)
  endfunction ()
endif ()

# ----------------------------------------------------------------------------
## @brief Get "local" target name, i.e., BASIS target name.
#
# If @c BASIS_USE_TARGET_UIDS is set to @c OFF, this operation
# always just sets the @p TARGET_NAME to the given @p TARGET_UID.
#
# @param [out] TARGET_NAME Target name used as argument to BASIS functions.
# @param [in]  TARGET_UID  "Global" target name, i.e., actual CMake target name.
#
# @returns Sets @p TARGET_NAME to the name of the build target with UID @p TARGET_UID.
#
# @sa basis_get_target_uid()
if (BASIS_USE_TARGET_UIDS)
  function (basis_get_target_name TARGET_NAME TARGET_UID)
    basis_get_fully_qualified_target_uid (UID "${TARGET_UID}")
    _basis_get_target_name (NAME "${UID}")
    set (${TARGET_NAME} "${NAME}" PARENT_SCOPE)
  endfunction ()
else ()
  function (basis_get_target_name TARGET_NAME TARGET_UID)
    set (${TARGET_NAME} "${TARGET_UID}" PARENT_SCOPE)
  endfunction ()
endif ()

# ----------------------------------------------------------------------------
## @brief Checks whether a given name is a valid target name.
#
# Displays fatal error message when target name is invalid.
#
# @param [in] TARGET_NAME Desired target name.
#
# @returns Nothing.
function (basis_check_target_name TARGET_NAME)
  # reserved target name ?
  foreach (PATTERN IN LISTS BASIS_RESERVED_TARGET_NAMES)
    if ("^${TARGET_NAME}$" STREQUAL "^${PATTERN}$")
      message (FATAL_ERROR "Target name \"${TARGET_NAME}\" is reserved and cannot be used.")
    endif ()
  endforeach ()
  # invalid target name ?
  if (NOT TARGET_NAME MATCHES "^\\.?[a-zA-Z]([a-zA-Z0-9_+.]|-)*$|^__init__(_py)?$")
    message (FATAL_ERROR "Target name '${TARGET_NAME}' is invalid.\nChoose a target name"
                         " which only contains alphanumeric characters,"
                         " '_', '-', or '+', and starts with a letter."
                         " The only exception from this rule is __init__[_py] for"
                         " a __init__.py script.\n")
  endif ()

  # unique ?
  basis_get_target_uid (TARGET_UID "${TARGET_NAME}")
  if (TARGET "${TARGET_UID}")
    message (FATAL_ERROR "There exists already a target named ${TARGET_UID}."
                         " Target names must be unique.")
  endif ()
endfunction ()

# ----------------------------------------------------------------------------
## @brief Make test UID from given test name.
#
# This function is intended for use by the basis_add_test() only.
#
# If @c BASIS_USE_TARGET_UIDS is set to @c OFF, this operation
# always just sets the @p TEST_UID to the given @p TEST_NAME.
#
# @param [out] TEST_UID  "Global" test name, i.e., actual CTest test name.
# @param [in]  TEST_NAME Test name used as argument to BASIS CMake functions.
#
# @returns Sets @p TEST_UID to the UID of the test @p TEST_NAME.
#
# @sa basis_get_test_uid()
function (basis_make_test_uid TEST_UID TEST_NAME)
  basis_make_target_uid (UID "${TEST_NAME}")
  set ("${TEST_UID}" "${UID}" PARENT_SCOPE)
endfunction ()

# ----------------------------------------------------------------------------
## @brief Get "global" test name, i.e., actual CTest test name.
#
# In order to ensure that CTest test names are unique across BASIS projects,
# the test name used by a developer of a BASIS project is converted by this
# function into another test name which is used as actual CTest test name.
#
# The function basis_get_test_name() can be used to convert the unique test
# name, the test UID, back to the original test name passed to this function.
#
# If @c BASIS_USE_TARGET_UIDS is set to @c OFF, this operation
# always just sets the @p TEST_UID to the given @p TEST_NAME.
#
# @param [out] TEST_UID  "Global" test name, i.e., actual CTest test name.
# @param [in]  TEST_NAME Test name used as argument to BASIS CMake functions.
#
# @returns Sets @p TEST_UID to the UID of the test @p TEST_NAME.
#
# @sa basis_get_test_name()
if (BASIS_USE_TARGET_UIDS)
  function (basis_get_test_uid TEST_UID TEST_NAME)
    if (TEST_NAME MATCHES "^\\.+(.*)$")
      set ("${TEST_UID}" "${CMAKE_MATCH_1}" PARENT_SCOPE)
    else ()
      set (UID "${PROJECT_NAMESPACE_CMAKE}.${TEST_NAME}")
      basis_strip_target_uid (UID ${UID})
      set ("${TEST_UID}" "${UID}" PARENT_SCOPE)
    endif ()
  endfunction ()
else ()
  function (basis_get_test_uid TEST_UID TEST_NAME)
    if (TEST_NAME MATCHES "^\\.+(.*)$")
      set ("${TEST_UID}" "${CMAKE_MATCH_1}" PARENT_SCOPE)
    else ()
      set ("${TEST_UID}" "${TEST_NAME}" PARENT_SCOPE)
    endif ()
  endfunction ()
endif ()

# ----------------------------------------------------------------------------
## @brief Get "global" test name, i.e., actual CTest test name.
#
# This function always returns a fully-qualified test UID, no matter if
# the option @c BASIS_USE_FULLY_QUALIFIED_UIDS is @c OFF. Note that
# if this option is @c ON, the returned test UID may not be the
# actual name of a CMake test.
#
# If @c BASIS_USE_TARGET_UIDS is set to @c OFF, this operation
# always just sets the @p TEST_UID to the given @p TEST_NAME.
#
# @param [out] TEST_UID  Fully-qualified test UID.
# @param [in]  TEST_NAME Test name used as argument to BASIS CMake functions.
#
# @sa basis_get_test_uid()
if (BASIS_USE_TARGET_UIDS)
  function (basis_get_fully_qualified_test_uid TEST_UID TEST_NAME)
    if (TEST_NAME MATCHES "^\\.+(.*)$")
      set ("${TEST_UID}" "${CMAKE_MATCH_1}" PARENT_SCOPE)
    else ()
      set ("${TEST_UID}" "${TOPLEVEL_PROJECT_NAMESPACE_CMAKE}.${TEST_NAME}" PARENT_SCOPE)
    endif ()
  endfunction ()
else ()
  function (basis_get_fully_qualified_test_uid TEST_UID TEST_NAME)
    if (TEST_NAME MATCHES "^\\.+(.*)$")
      set ("${TEST_UID}" "${CMAKE_MATCH_1}" PARENT_SCOPE)
    else ()
      set ("${TEST_UID}" "${TEST_NAME}" PARENT_SCOPE)
    endif ()
  endfunction ()
endif ()

# ----------------------------------------------------------------------------
## @brief Get namespace of test.
#
# If @c BASIS_USE_TARGET_UIDS is set to @c OFF, this operation
# always just sets the @p TEST_NS to an empty string.
#
# @param [out] TEST_NS  Namespace part of test UID. If @p TEST_UID is
#                       no UID, i.e., does not contain a namespace part,
#                       the namespace of this project is returned.
# @param [in]  TEST_UID Test UID/name.
if (BASIS_USE_TARGET_UIDS)
  function (basis_get_test_namespace TEST_NS TEST_UID)
    if (TEST_UID MATCHES "^(.*)\\.")
      set (${TEST_NS} "${CMAKE_MATCH_1}" PARENT_SCOPE)
    else ()
      set (${TEST_NS} "" PARENT_SCOPE)
    endif ()
  endfunction ()
else ()
  function (basis_get_test_namespace TEST_NS TEST_UID)
    set ("${TEST_UID}" "" PARENT_SCOPE)
  endfunction ()
endif ()

# ----------------------------------------------------------------------------
## @brief Get "local" test name, i.e., BASIS test name.
#
# If @c BASIS_USE_TARGET_UIDS is set to @c OFF, this operation
# always just sets the @p TEST_NAME to the given @p TEST_UID.
#
# @param [out] TEST_NAME Test name used as argument to BASIS functions.
# @param [in]  TEST_UID  "Global" test name, i.e., actual CTest test name.
#
# @returns Sets @p TEST_NAME to the name of the test with UID @p TEST_UID.
#
# @sa basis_get_test_uid()
if (BASIS_USE_TARGET_UIDS)
  function (basis_get_test_name TEST_NAME TEST_UID)
    if (TEST_UID MATCHES "([^.]+)$")
      set (${TEST_NAME} "${CMAKE_MATCH_1}" PARENT_SCOPE)
    else ()
      set (${TEST_NAME} "" PARENT_SCOPE)
    endif ()
  endfunction ()
else ()
  function (basis_get_test_name TEST_NAME TEST_UID)
    set ("${TEST_NAME}" "${TEST_UID}" PARENT_SCOPE)
  endfunction ()
endif ()

# ----------------------------------------------------------------------------
## @brief Checks whether a given name is a valid test name.
#
# Displays fatal error message when test name is invalid.
#
# @param [in] TEST_NAME Desired test name.
#
# @returns Nothing.
function (basis_check_test_name TEST_NAME)
  # reserved test name ?
  foreach (PATTERN IN LISTS BASIS_RESERVED_TARGET_NAMES)
    if ("^${TARGET_NAME}$" STREQUAL "^${PATTERN}$")
      message (FATAL_ERROR "Test name \"${TARGET_NAME}\" is reserved and cannot be used.")
    endif ()
  endforeach ()
  # invalid test name ?
  if (NOT TEST_NAME MATCHES "^\\.?[a-zA-Z]([a-zA-Z0-9_+.]|-)*$")
    message (FATAL_ERROR "Test name ${TEST_NAME} is invalid.\nChoose a test name "
                         " which only contains alphanumeric characters,"
                         " '_', '-', or '+', and starts with a letter.\n")
  endif ()
endfunction ()

# ============================================================================
# common target tools
# ============================================================================

# ----------------------------------------------------------------------------
## @brief Whether a given target exists.
#
# This function should be used instead of the if(TARGET) command of CMake
# because target names are mapped by BASIS to target UIDs.
#
# @param [out] RESULT_VARIABLE Boolean result variable.
# @param [in]  TARGET_NAME     Name which to check whether it is a target.
#
# @sa basis_make_target_uid()
# @sa basis_get_target_uid()
macro (basis_exists_target RESULT_VARIABLE TARGET_NAME)
  basis_get_target_uid (_UID "${TARGET_NAME}")
  if (TARGET ${_UID})
    set (${RESULT_VARIABLE} TRUE)
  else ()
    set (${RESULT_VARIABLE} FALSE)
  endif ()
  unset (_UID)
endmacro ()

# ----------------------------------------------------------------------------
## @brief Get default subdirectory prefix of scripted library modules.
#
# @param [out] PREFIX   Name of variable which is set to the default library
#                       prefix, i.e., subdirectory relative to the library
#                       root directory as used for the @c PREFIX property of
#                       scripted module libraries (see basis_add_script_library())
#                       or relative to the include directory in case of C++.
#                       Note that this prefix includes a trailing slash to
#                       indicate that the prefix is a subdirectory.
# @param [in]  LANGUAGE Programming language (case-insenitive), e.g.,
#                       @c CXX, @c Python, @c Matlab...
macro (basis_library_prefix PREFIX LANGUAGE)
  string (TOUPPER "${LANGUAGE}" _LANGUAGE_U)
  if (PROJECT_NAMESPACE_${_LANGUAGE_U})
    basis_sanitize_for_regex (_RE "${BASIS_NAMESPACE_DELIMITER_${_LANGUAGE_U}}")
    string (REGEX REPLACE "${_RE}" "/" ${PREFIX} "${PROJECT_NAMESPACE_${_LANGUAGE_U}}")
    set (${PREFIX} "${${PREFIX}}/")
    unset (_RE)
  else ()
    message (FATAL_ERROR "basis_library_prefix(): PROJECT_NAMESPACE_${_LANGUAGE_U} not set!"
                         " Make sure that the LANGUAGE argument is supported and in"
                         " uppercase letters only.")
  endif ()
  unset (_LANGUAGE_U)
endmacro ()

# ----------------------------------------------------------------------------
## @brief Get file name of compiled script.
#
# @param [out] CFILE  File path of compiled script file.
# @param [in]  SOURCE Script source file.
# @param [in]  ARGV2  Language of script file. If not specified, the language
#                     is derived from the file name extension and shebang of
#                     the script source file.
function (basis_get_compiled_file CFILE SOURCE)
  if (ARGC GREATER 2)
    set (LANGUAGE "${ARGV2}")
  else ()
    basis_get_source_language (LANGUAGE "${SOURCE}")
  endif ()
  set (${CFILE} "" PARENT_SCOPE)
  if (SOURCE)
    if (LANGUAGE MATCHES "PYTHON")
      set (${CFILE} "${SOURCE}c" PARENT_SCOPE)
    elseif (LANGUAGE MATCHES "JYTHON")
      if (SOURCE MATCHES "(.*)\\.([^.]+)$")
        set (${CFILE} "${CMAKE_MATCH_1}$${CMAKE_MATCH_2}.class" PARENT_SCOPE)
      endif ()
    endif ()
  endif ()
endfunction ()

# ----------------------------------------------------------------------------
## @brief Get file path of Jython file compiled from the given Python module.
#
# Python modules are also compiled using Jython. This macro returns the file
# path of the compiled Jython file in the build tree which corresponds to the
# specified Python module.
#
# @param [out] CFILE  Path of corresponding compiled Jython file.
# @param [in]  MODULE Path of input Python module in build tree.
macro (basis_get_compiled_jython_file_of_python_module CFILE MODULE)
  if (BINARY_PYTHON_LIBRARY_DIR AND BINARY_JYTHON_LIBRARY_DIR)
    file (RELATIVE_PATH _GCJFOPM_REL "${BINARY_PYTHON_LIBRARY_DIR}" "${MODULE}")
  else ()
    set (_GCJFOPM_REL)
  endif ()
  if (NOT _GCJFOPM_REL MATCHES "^$|^\\.\\./")
    basis_get_compiled_file (${CFILE} "${BINARY_JYTHON_LIBRARY_DIR}/${_GCJFOPM_REL}" JYTHON)
  else ()
    basis_get_compiled_file (${CFILE} "${MODULE}" JYTHON)
  endif ()
  unset (_GCJFOPM_REL)
endmacro ()

# ----------------------------------------------------------------------------
## @brief Whether to compile Python modules for Jython interpreter.
#
# This macro returns a boolean value stating whether Python modules shall also
# be compiled for use by Jython interpreter if BASIS_COMPILE_SCRIPTS is ON.
#
# @param [out] FLAG Set to either TRUE or FALSE depending on whether Python
#                   modules shall be compiled using Jython or not.
macro (basis_compile_python_modules_for_jython FLAG)
  if (BASIS_COMPILE_SCRIPTS AND JYTHON_EXECUTABLE)
    set (${FLAG} TRUE)
  else ()
    set (${FLAG} FALSE)
  endif ()
  if (DEFINED USE_JythonInterp AND NOT USE_JythonInterp)
    set (${FLAG} FALSE)
  endif ()
endmacro ()

# ----------------------------------------------------------------------------
## @brief Glob source files.
#
# This function gets a list of source files and globbing expressions, evaluates
# the globbing expression, and replaces these by the found source files such
# that the resulting list of source files contains only absolute paths of
# source files. It is used by basis_add_executable() and basis_add_library()
# to get a list of all source files. The syntax for the glob expressions
# corresponds to the one used by CMake's
# <a href="http://www.cmake.org/cmake/help/v2.8.8/cmake.html#command:file">
# file(GLOB)</a> command. Additionally, if the pattern <tt>**</tt> is found
# in a glob expression, it is replaced by a single <tt>*</tt> and the
# recursive version, i.e., <tt>file(GLOB_RECURSE)</tt>, is used instead.
#
# @param [in]  TARGET_UID UID of build target which builds the globbed source files.
#                         The custom target which re-globs the source files
#                         before each build of this target is named after this
#                         build target with two leading underscores (__).
# @param [out] SOURCES    List of absolute source paths.
# @param [in]  ARGN       Input file paths and/or globbing expressions.
#
# @sa basis_add_executable()
# @sa basis_add_library()
function (basis_add_glob_target TARGET_UID SOURCES)
  # prepare globbing expressions
  # make paths absolute and turn directories into recursive globbing expressions
  set (EXPRESSIONS)
  foreach (EXPRESSION IN LISTS ARGN)
    if (NOT EXPRESSION MATCHES "^\\$<") # preserve generator expressions
      if (NOT IS_ABSOLUTE "${EXPRESSION}")
        # prefer configured/generated files in the build tree, but disallow
        # globbing within the build tree; glob only files in source tree
        if (NOT EXPRESSION MATCHES "[*?]|\\[[0-9]+-[0-9]+\\]"            AND
            EXISTS           "${CMAKE_CURRENT_BINARY_DIR}/${EXPRESSION}" AND
            NOT IS_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${EXPRESSION}")
          set (EXPRESSION "${CMAKE_CURRENT_BINARY_DIR}/${EXPRESSION}")
        else ()
          set (EXPRESSION "${CMAKE_CURRENT_SOURCE_DIR}/${EXPRESSION}")
        endif ()
      endif ()
      if (IS_DIRECTORY "${EXPRESSION}")
        set (EXPRESSION "${EXPRESSION}/**")
      endif ()
    endif ()
    list (APPEND EXPRESSIONS "${EXPRESSION}")
  endforeach ()
  # only if at least one globbing expression is given we need to go through this hassle
  if (EXPRESSIONS MATCHES "[*?]|\\[[0-9]+-[0-9]+\\]")
    set (BUILD_DIR    "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${TARGET_UID}.dir")
    set (SOURCES_FILE "${BUILD_DIR}/sources.txt")
    # get initial list of source files
    execute_process (
      COMMAND "${CMAKE_COMMAND}"
                  "-DEXPRESSIONS:STRING=${EXPRESSIONS}"
                  "-DINIT:BOOLEAN=TRUE"
                  "-DSOURCES_FILE:FILEPATH=${SOURCES_FILE}"
                  -P "${BASIS_MODULE_PATH}/glob.cmake"
      RESULT_VARIABLE RETVAL
    )
    if (NOT RETVAL EQUAL 0 OR NOT EXISTS "${SOURCES_FILE}")
      message (FATAL_ERROR "Target ${TARGET_UID}: Failed to glob source files!")
    endif ()
    # note that including this file here, which is modified whenever a
    # source file is added or removed, triggers a re-configuration of the
    # build system which is required to re-execute this function.
    include ("${SOURCES_FILE}")
    set (${SOURCES} "${INITIAL_SOURCES}" PARENT_SCOPE)
    # add custom target to re-glob files before each build
    set (ERRORMSG "You have either added, removed, or renamed a source file which"
                  " matches one of the globbing expressions specified for the"
                  " list of source files from which to build the ${TARGET_UID} target."
                  " Therefore, the build system must be re-configured. Either try to"
                  " build again which should trigger CMake and re-configure the build"
                  " system or run CMake manually.")
    basis_list_to_string (ERRORMSG ${ERRORMSG})
    add_custom_target (
      __${TARGET_UID}
      COMMAND "${CMAKE_COMMAND}"
                  "-DEXPRESSIONS:STRING=${EXPRESSIONS}"
                  "-DINIT:BOOLEAN=FALSE"
                  "-DSOURCES_FILE:FILEPATH=${SOURCES_FILE}"
                  "-DERRORMSG:STRING=${ERRORMSG}"
                  -P "${BASIS_MODULE_PATH}/glob.cmake"
      COMMENT "Checking if source files for target ${TARGET_UID} were added or removed"
      VERBATIM
    )
  # otherwise, just return the given absolute source file paths
  else ()
    set (${SOURCES} "${EXPRESSIONS}" PARENT_SCOPE)
  endif ()
endfunction ()

# ----------------------------------------------------------------------------
## @brief Detect programming language of given source code files.
#
# This function determines the programming language in which the given source
# code files are written. If no common programming language could be determined,
# "AMBIGUOUS" is returned. If none of the following programming languages
# could be determined, "UNKNOWN" is returned: CXX (i.e., C++), JAVA, MATLAB,
# PYTHON, JYTHON, PERL, BASH, BATCH.
#
# @param [out] LANGUAGE Detected programming language.
# @param [in]  ARGN     List of source code files.
function (basis_get_source_language LANGUAGE)
  set (LANGUAGE_OUT)
  # iterate over source files
  foreach (SOURCE_FILE ${ARGN})
    # skip generator expressions
    if (NOT SOURCE_FILE MATCHES "^\\$<")
      get_filename_component (SOURCE_FILE "${SOURCE_FILE}" ABSOLUTE)

      if (IS_DIRECTORY "${SOURCE_FILE}")

        file (GLOB_RECURSE SOURCE_FILES "${SOURCE_FILE}/*")
        list (APPEND ARGN ${SOURCE_FILES})

      else ()

        # ------------------------------------------------------------------------
        # determine language based on extension for those without shebang
        set (LANG)
        # C++
        if (SOURCE_FILE MATCHES "\\.(c|cc|cpp|cxx|h|hh|hpp|hxx|txx|inl)(\\.in)?$")
          set (LANG "CXX")
        # Java
        elseif (SOURCE_FILE MATCHES "\\.java(\\.in)?$")
          set (LANG "JAVA")
        # MATLAB
        elseif (SOURCE_FILE MATCHES "\\.m(\\.in)?$")
          set (LANG "MATLAB")
        endif ()

        # ------------------------------------------------------------------------
        # determine language from shebang directive
        #
        # Note that some scripting languages may use identical file name extensions.
        # This is in particular the case for Python and Jython. The only way we
        # can distinguish these two is by looking at the shebang directive.
        if (NOT LANG)
          
          if (NOT EXISTS "${SOURCE_FILE}" AND EXISTS "${SOURCE_FILE}.in")
            set (SOURCE_FILE "${SOURCE_FILE}.in")
          endif ()
          if (EXISTS "${SOURCE_FILE}")
            file (STRINGS "${SOURCE_FILE}" FIRST_LINE LIMIT_COUNT 1)
            if (FIRST_LINE MATCHES "^#!")
              if (FIRST_LINE MATCHES "^#! */usr/bin/env +([^ ]+)")
                set (INTERPRETER "${CMAKE_MATCH_1}")
              elseif (FIRST_LINE MATCHES "^#! *([^ ]+)")
                set (INTERPRETER "${CMAKE_MATCH_1}")
                get_filename_component (INTERPRETER "${INTERPRETER}" NAME)
              else ()
                set (INTERPRETER)
              endif ()
              if (INTERPRETER MATCHES "^(python|jython|perl|bash)$")
                string (TOUPPER "${INTERPRETER}" LANG)
              endif ()
            endif ()
          endif ()
        endif ()

        # ------------------------------------------------------------------------
        # determine language from further known extensions
        if (NOT LANG)
          # Python
          if (SOURCE_FILE MATCHES "\\.py(\\.in)?$")
            set (LANG "PYTHON")
          # Perl
          elseif (SOURCE_FILE MATCHES "\\.(pl|pm|t)(\\.in)?$")
            set (LANG "PERL")
          # BASH
          elseif (SOURCE_FILE MATCHES "\\.sh(\\.in)?$")
            set (LANG "BASH")
          # Batch
          elseif (SOURCE_FILE MATCHES "\\.bat(\\.in)?$")
            set (LANG "BATCH")
          # unknown
          else ()
            set (LANGUAGE_OUT "UNKNOWN")
            break ()
          endif ()
        endif ()

        # ------------------------------------------------------------------------
        # detect ambiguity
        if (LANGUAGE_OUT AND NOT "^${LANG}$" STREQUAL "^${LANGUAGE_OUT}$")
          if (LANGUAGE_OUT MATCHES "CXX" AND LANG MATCHES "MATLAB")
            # MATLAB Compiler can handle this...
          elseif (LANGUAGE_OUT MATCHES "MATLAB" AND LANG MATCHES "CXX")
            set (LANG "MATLAB") # language stays MATLAB
          elseif (LANGUAGE_OUT MATCHES "PYTHON" AND LANG MATCHES "JYTHON")
            # Jython can deal with Python scripts/modules
          elseif (LANGUAGE_OUT MATCHES "JYTHON" AND LANG MATCHES "PYTHON")
            set (LANG "JYTHON") # language stays JYTHON
          else ()
            # ambiguity
            set (LANGUAGE_OUT "AMBIGUOUS")
            break ()
          endif ()
        endif ()

        # update current language
        set (LANGUAGE_OUT "${LANG}")
      endif ()
    endif ()
  endforeach ()
  # return
  if (LANGUAGE_OUT)
    set (${LANGUAGE} "${LANGUAGE_OUT}" PARENT_SCOPE)
  else ()
    message (FATAL_ERROR "basis_get_source_language called without arguments!")
  endif ()
endfunction ()

# ----------------------------------------------------------------------------
## @brief Configure .in source files.
#
# This function configures each source file in the given argument list with
# a .in file name suffix and stores the configured file in the build tree
# with the same relative directory as the template source file itself.
# The first argument names the CMake variable of the list of configured
# source files where each list item is the absolute file path of the
# corresponding (configured) source file.
#
# @param [out] LIST_NAME Name of output list.
# @param [in]  ARGN      These arguments are parsed and the following
#                        options recognized. All remaining arguments are
#                        considered to be source file paths.
# @par
# <table border="0">
#   <tr>
#     @tp @b BINARY_DIRECTORY @endtp
#     <td>Explicitly specify directory in build tree where configured
#         source files should be written to.</td>
#   </tr>
#   <tr>
#     @tp @b KEEP_DOT_IN_SUFFIX @endtp
#     <td>By default, after a source file with the .in extension has been
#         configured, the .in suffix is removed from the file name.
#         This can be omitted by giving this option.</td>
#   </tr>
# </table>
#
# @returns Nothing.
function (basis_configure_sources LIST_NAME)
  # parse arguments
  CMAKE_PARSE_ARGUMENTS (ARGN "KEEP_DOT_IN_SUFFIX" "BINARY_DIRECTORY" "" ${ARGN})

  # ensure that specified BINARY_DIRECTORY is inside build tree of project
  if (ARGN_BINARY_DIRECTORY)
    get_filename_component (_binpath "${ARGN_BINARY_DIRECTORY}" ABSOLUTE)
    file (RELATIVE_PATH _relpath "${PROJECT_BINARY_DIR}" "${_binpath}")
    if (_relpath MATCHES "^\\.\\./")
      message (FATAL_ERROR "Specified BINARY_DIRECTORY must be inside the build tree!")
    endif ()
    unset (_binpath)
    unset (_relpath)
  endif ()

  # configure source files
  set (CONFIGURED_SOURCES)
  foreach (SOURCE ${ARGN_UNPARSED_ARGUMENTS})
    # The .in suffix is optional, add it here if a .in file exists for this
    # source file, but only if the source file itself does not name an actually
    # existing source file.
    #
    # If the source file path is relative, prefer possibly already configured
    # sources in build tree such as the test driver source file created by
    # create_test_sourcelist() or a manual use of configure_file().
    #
    # Note: Make path absolute, otherwise EXISTS check will not work!
    if (NOT IS_ABSOLUTE "${SOURCE}")
      if (EXISTS "${CMAKE_CURRENT_BINARY_DIR}/${SOURCE}")
        set (SOURCE "${CMAKE_CURRENT_BINARY_DIR}/${SOURCE}")
      elseif (EXISTS "${CMAKE_CURRENT_BINARY_DIR}/${SOURCE}.in")
        set (SOURCE "${CMAKE_CURRENT_BINARY_DIR}/${SOURCE}.in")
      elseif (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${SOURCE}")
        set (SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/${SOURCE}")
      elseif (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${SOURCE}.in")
        set (SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/${SOURCE}.in")
      endif ()
    else ()
      if (NOT EXISTS "${SOURCE}" AND EXISTS "${SOURCE}.in")
        set (SOURCE "${SOURCE}.in")
      endif ()
    endif ()
    # configure source file if filename ends in .in suffix
    if (SOURCE MATCHES "\\.in$")
      # if binary directory was given explicitly, use it
      if (ARGN_BINARY_DIRECTORY)
        get_filename_component (SOURCE_NAME "${SOURCE}" NAME)
        if (NOT ARGN_KEEP_DOT_IN_SUFFIX)
          string (REGEX REPLACE "\\.in$" "" SOURCE_NAME "${SOURCE_NAME}")
        endif ()
        set (CONFIGURED_SOURCE "${ARGN_BINARY_DIRECTORY}/${SOURCE_NAME}")
      # otherwise,
      else ()
        # if source is in project's source tree use relative binary directory
        if ("^${SOURCE}$" STREQUAL "^${PROJECT_SOURCE_DIR}$")
          basis_get_relative_path (CONFIGURED_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}" "${SOURCE}")
          get_filename_component (CONFIGURED_SOURCE "${CMAKE_CURRENT_BINARY_DIR}/${CONFIGURED_SOURCE}" ABSOLUTE)
          if (NOT ARGN_KEEP_DOT_IN_SUFFIX)
            string (REGEX REPLACE "\\.in$" "" CONFIGURED_SOURCE "${CONFIGURED_SOURCE}")
          endif ()
        # otherwise, use current binary directory
        else ()
          get_filename_component (SOURCE_NAME "${SOURCE}" NAME)
          if (NOT ARGN_KEEP_DOT_IN_SUFFIX)
            string (REGEX REPLACE "\\.in$" "" SOURCE_NAME "${SOURCE_NAME}")
          endif ()
          set (CONFIGURED_SOURCE "${CMAKE_CURRENT_BINARY_DIR}/${SOURCE_NAME}")
        endif ()
      endif ()
      # configure source file
      configure_file ("${SOURCE}" "${CONFIGURED_SOURCE}" @ONLY)
      if (BASIS_DEBUG)
        message ("** Configured source file with .in extension")
      endif ()
    # otherwise, skip configuration of this source file
    else ()
      set (CONFIGURED_SOURCE "${SOURCE}")
      if (BASIS_DEBUG)
        message ("** Skipped configuration of source file")
      endif ()
    endif ()
    if (BASIS_DEBUG)
      message ("**     Source:            ${SOURCE}")
      message ("**     Configured source: ${CONFIGURED_SOURCE}")
    endif ()
    list (APPEND CONFIGURED_SOURCES "${CONFIGURED_SOURCE}")
  endforeach ()
  # return
  set (${LIST_NAME} "${CONFIGURED_SOURCES}" PARENT_SCOPE)
endfunction ()

# ----------------------------------------------------------------------------
## @brief Remove one blank line from top of string
macro (basis_remove_blank_line STRVAR)
  if (${STRVAR} MATCHES "(^|(.*)\n)[ \t]*\n(.*)")
    set (${STRVAR} "${CMAKE_MATCH_1}${CMAKE_MATCH_3}")
  endif ()
endmacro ()

# ----------------------------------------------------------------------------
## @brief Configure and optionally compile script file.
#
# This function is used to configure script files during the build. It is
# called by the build script generated by basis_add_script_target() for each script
# target. It is also used to configure the modules of the packages
# implemented in supported scripting languages which are located in the
# @c PROJECT_LIBRARY_DIR of the source tree.
#
# In case of executable scripts, this function automatically prepends the
# module search paths such that the modules of this software package are found
# (and preferred in case of potential name conflicts with other packages).
# Moreover, it adds (or replaces) the shebang directive on Unix such that the
# configured interpreter version is used. On Windows, it converts the executable
# script into a Windows Command instead which executes the proper interpreter
# with the code section of the input script.
#
# @param [in] INPUT  Input script file.
# @param [in] OUTPUT Configured output script file.
# @param [in] ARGN   Optional arguments:
# @par
# <table border=0>
#   <tr>
#     @tp @b COMPILE @endtp
#     <td>Whether to compile module scripts if suitable, i.e., an intermediate
#         format exists for the specific scripting language. For example,
#         Python modules can be compiled.</td>
#   </tr>
#     @tp @b CONFIGURATION name @endtp
#     <td>Name of build configuration.</td>
#   </tr>
#   <tr>
#     @tp @b COPYONLY @endtp
#     <td>Whether to only copy the script file without replacing CMake variables
#         within the file. This option is passed on to CMake's configure_file()
#         command used to configure the script file. By default, the option
#         \@ONLY is used instead.</td>
#   </tr>
#   <tr>
#     @tp @b EXECUTABLE @endtp
#     <td>Specifies that the given script file is an executable script and not a
#         module script. Otherwise, if this option is not given and the output
#         file name contains a file name extension, the given script file is
#         configured as module script. A script file with an output file name
#         that has no extension, is always considered to be an executable.</td>
#   </tr>
#   <tr>
#     @tp @b DIRECTORY dir @endtp
#     <td>Build tree directory of configured script. By default, the directory
#         of the output script file is used to set the @c __DIR__ variable
#         used to make absolute file paths relative to the configured build
#         tree script file in the script configuration file
#         (see basis_set_script_path()).
#   </tr>
#   <tr>
#     @tp @b DESTINATION dir @endtp
#     <td>Installation directory for configured script. If this option is given,
#         the @c BUILD_INSTALL_SCRIPT variable is set to @c TRUE before including
#         any specified script configuration files (see @p CONFIG_FILE option).
#         Moreover, the @c __DIR__ variable is set to the specified directory.
#         Otherwise, if this option is omitted, the @c BUILD_INSTALL_SCRIPT variable
#         is set to @c FALSE instead and @c __DIR__ is set to the directory of
#         the configured @p OUTPUT file. Note that the @c BUILD_INSTALL_SCRIPT and
#         @c __DIR__ variables are in particular used by basis_set_script_path()
#         to convert the given paths to paths relative to the location of the
#         configured/installed script.</td>
#   </tr>
#   <tr>
#     @tp @b CACHE_FILE file1 [file2...] @endtp
#     <td>List of CMake files with dump of variables which should be included
#         before configuring the script. The cache files can be generated using
#         the basis_dump_variables() function.</td>
#   </tr>
#   <tr>
#     @tp @b CONFIG_FILE file1 [file2...] @endtp
#     <td>List of script configuration files to include before the configuration
#         of the script. See also the documentation of the @p DESTINATION option.</td>
#   </tr>
#   <tr>
#     @tp @b LINK_DEPENDS dep1 [dep2...] @endtp
#     <td>List of "link" dependencies, i.e., modules and script/module libraries
#         required by this script. For executable scripts, the paths to these
#         modules/packages is added to the module search path. If the prefix
#         "relative " is given before a file path, it is made relative to the
#         output/installation directory of the script file. All given input paths
#         must be absolute, however, as the relative location depends on
#         whether the script will be installed, i.e., the @c DESTINATION
#         is specified, or not.</td>
#   </tr>
# </table>
function (basis_configure_script INPUT OUTPUT)
  # rename arguments to avoid conflict with script configuration
  set (_INPUT_FILE  "${INPUT}")
  set (_OUTPUT_FILE "${OUTPUT}")
  # --------------------------------------------------------------------------
  # parse arguments
  CMAKE_PARSE_ARGUMENTS (
    ARGN
      "COMPILE;COPYONLY;EXECUTABLE"
      "DIRECTORY;DESTINATION;LANGUAGE;CONFIGURATION"
      "CACHE_FILE;CONFIG_FILE;LINK_DEPENDS"
    ${ARGN}
  )
  if (ARGN_UNPARSED_ARGUMENTS)
    message (FATAL_ERROR "Unrecognized arguments given: ${ARGN_UNPARSED_ARGUMENTS}")
  endif ()
  if (NOT ARGN_LANGUAGE)
    basis_get_source_language (ARGN_LANGUAGE "${_INPUT_FILE}")
  endif ()
  # --------------------------------------------------------------------------
  # include cache files
  foreach (_F IN LISTS ARGN_CACHE_FILE)
    get_filename_component (_F "${_F}" ABSOLUTE)
    if (NOT EXISTS "${_F}")
      message (FATAL_ERROR "Cache file ${_F} does not exist!")
    endif ()
    include ("${_F}")
  endforeach ()
  # --------------------------------------------------------------------------
  # general variables for use in scripts and those intended for use in script
  # configurations. The __DIR__ is in particular used by basis_set_script_path()
  # to make absolute paths relative to the script file location
  if (ARGN_DESTINATION)
    if (NOT IS_ABSOLUTE "${ARGN_DESTINATION}")
      set (ARGN_DESTINATION "${CMAKE_INSTALL_PREFIX}/${ARGN_DESTINATION}")
    endif ()
    set (BUILD_INSTALL_SCRIPT TRUE)
    set (__DIR__ "${ARGN_DESTINATION}")
  else ()
    set (BUILD_INSTALL_SCRIPT FALSE)
    if (ARGN_DIRECTORY)
      if (NOT IS_ABSOLUTE "${ARGN_DIRECTORY}")
        message (FATAL_ERROR "Build tree script DIRECTORY must be an absolute path")
      endif ()
      set (__DIR__ "${ARGN_DIRECTORY}")
      string (REPLACE "$<${BASIS_GE_CONFIG}>" "${ARGN_CONFIGURATION}" __DIR__ "${__DIR__}")
    else ()
      get_filename_component (__DIR__ "${_OUTPUT_FILE}" PATH)
    endif ()
  endif ()
  get_filename_component (__NAME__ "${_OUTPUT_FILE}" NAME)
  set (__FILE__   "${__DIR__}/${__NAME__}")
  set (__CONFIG__ "${ARGN_CONFIGURATION}")
  # --------------------------------------------------------------------------
  # replace $<CONFIG> (deprecated $<CONFIGURATION>) in paths of link dependencies
  string (REPLACE "$<${BASIS_GE_CONFIG}>" "${ARGN_CONFIGURATION}" ARGN_LINK_DEPENDS "${ARGN_LINK_DEPENDS}")
  # --------------------------------------------------------------------------
  # include script configuration files
  foreach (_F IN LISTS ARGN_CONFIG_FILE)
    get_filename_component (_F "${_F}" ABSOLUTE)
    if (NOT EXISTS "${_F}")
      message (FATAL_ERROR "Script configuration file ${_F} does not exist!")
    endif ()
    include ("${_F}")
  endforeach ()
  # --------------------------------------------------------------------------
  # configure executable script
  if (ARGN_EXECUTABLE)
    # Attention: Every line of code added/removed will introduce a mismatch
    #            between error messages of the interpreter and the original
    #            source file. To not confuse/mislead developers too much,
    #            keep number of lines added/removed at a minimum or at least
    #            try to balance the number of lines added and removed.
    #            Moreover, blank lines can be used to insert code without
    #            changing the number of source code lines.
    file (READ "${_INPUT_FILE}" SCRIPT)
    # (temporarily) remove existing shebang directive
    file (STRINGS "${_INPUT_FILE}" FIRST_LINE LIMIT_COUNT 1)
    if (FIRST_LINE MATCHES "^#!")
      basis_sanitize_for_regex (FIRST_LINE_RE "${FIRST_LINE}")
      string (REGEX REPLACE "^${FIRST_LINE_RE}\n?" "" SCRIPT "${SCRIPT}")
      set (SHEBANG "${FIRST_LINE}")
    endif ()
    # replace CMake variables used in script
    if (NOT ARGN_COPYONLY)
      string (CONFIGURE "${SCRIPT}" SCRIPT @ONLY)
    endif ()
    # add code to set module search path
    if (ARGN_LANGUAGE MATCHES "[JP]YTHON")
      if (ARGN_LINK_DEPENDS)
        set (PYTHON_CODE "import sys; import os.path; __dir__ = os.path.dirname(os.path.realpath(__file__))")
        list (REVERSE ARGN_LINK_DEPENDS)
        foreach (DIR ${ARGN_LINK_DEPENDS})
          if (DIR MATCHES "^relative +(.*)$")
            basis_get_relative_path (DIR "${__DIR__}" "${CMAKE_MATCH_1}")
          endif ()
          if (DIR MATCHES "\\.(py|class)$")
            get_filename_component (DIR "${DIR}" PATH)
          endif ()
          if (IS_ABSOLUTE "${DIR}")
            set (PYTHON_CODE "${PYTHON_CODE}; sys.path.insert(0, os.path.realpath('${DIR}'))")
          else ()
            set (PYTHON_CODE "${PYTHON_CODE}; sys.path.insert(0, os.path.realpath(os.path.join(__dir__, '${DIR}')))")
          endif ()
        endforeach ()
        # insert extra Python code near top, but after any future statement
        # (http://docs.python.org/2/reference/simple_stmts.html#future)
        set (FUTURE_STATEMENTS)
        if (SCRIPT MATCHES "^(.*from[ \t]+__future__[ \t]+import[ \t]+[a-z_]+[^\n]*\n)(.*)$")
          set (FUTURE_STATEMENTS "${CMAKE_MATCH_1}")
          set (SCRIPT            "${CMAKE_MATCH_2}")
        endif ()
        basis_remove_blank_line (SCRIPT) # remove a blank line therefore
        set (SCRIPT "${FUTURE_STATEMENTS}${PYTHON_CODE} # <-- added by BASIS\n${SCRIPT}")
      endif ()
    elseif (ARGN_LANGUAGE MATCHES "PERL")
      if (ARGN_LINK_DEPENDS)
        set (PERL_CODE "use Cwd qw(realpath); use File::Basename;")
        foreach (DIR ${ARGN_LINK_DEPENDS})
          if (DIR MATCHES "^relative +(.*)$")
            basis_get_relative_path (DIR "${__DIR__}" "${CMAKE_MATCH_1}")
          endif ()
          if (DIR MATCHES "\\.pm$")
            get_filename_component (DIR "${DIR}" PATH)
          endif ()
          if (IS_ABSOLUTE "${DIR}")
            set (PERL_CODE "${PERL_CODE} use lib '${DIR}';")
          else ()
            set (PERL_CODE "${PERL_CODE} use lib dirname(realpath(__FILE__)) . '/${DIR}';")
          endif ()
        endforeach ()
        basis_remove_blank_line (SCRIPT) # remove a blank line therefore
        set (SCRIPT "${PERL_CODE} # <-- added by BASIS\n${SCRIPT}")
      endif ()
    elseif (ARGN_LANGUAGE MATCHES "BASH")
      basis_library_prefix (PREFIX BASH)
      # In case of Bash, set BASIS_BASH_UTILITIES which is required to first source the
      # BASIS utilities modules (in particular core.sh). This variable should be set to
      # the utilities.sh module of BASIS by default as part of the BASIS installation
      # (environment variable) and is here set to the project-specific basis.sh module.
      #
      # Note that folks at SBIA may submit a Bash script directly to a batch queuing
      # system such as the Oracle Grid Engine (SGE) instead of writing a separate submit
      # script. To avoid not finding the BASIS utilities in this case only because the
      # Bash file was copied by SGE to a temporary file, consider the <PROJECT>_DIR
      # environment variable as an alternative.
      set (BASH_CODE
# Note: Code formatted such that it can be on single line. Use no comments within!
"__FILE__=\"$(cd -P -- \"$(dirname -- \"$BASH_SOURCE\")\" && pwd -P)/$(basename -- \"$BASH_SOURCE\")\"
if [[ -n \"$SGE_ROOT\" ]] && [[ $__FILE__ =~ $SGE_ROOT/.* ]] && [[ -n \"\${${PROJECT_NAME}_DIR}\" ]] && [[ -f \"\${${PROJECT_NAME}_DIR}/bin/${__NAME__}\" ]]
then __FILE__=\"\${${PROJECT_NAME}_DIR}/bin/${__NAME__}\"
fi
i=0
lnk=\"$__FILE__\"
while [[ -h \"$lnk\" ]] && [[ $i -lt 100 ]]
do dir=`dirname -- \"$lnk\"`
lnk=`readlink -- \"$lnk\"`
lnk=`cd \"$dir\" && cd $(dirname -- \"$lnk\") && pwd`/`basename -- \"$lnk\"`
let i++
done
[[ $i -lt 100 ]] && __FILE__=\"$lnk\"
unset -v i dir lnk
__DIR__=\"$(dirname -- \"$__FILE__\")\"
BASIS_BASH_UTILITIES=\"$__DIR__/${BASH_LIBRARY_DIR}/${PREFIX}basis.sh\""
      )
      string (REPLACE "\n" "; " BASH_CODE "${BASH_CODE}")
      # set BASHPATH which is used by import() function provided by core.sh module of BASIS
      set (BASHPATH)
      foreach (DIR ${ARGN_LINK_DEPENDS})
        if (DIR MATCHES "^relative +(.*)$")
          basis_get_relative_path (DIR "${__DIR__}" "${CMAKE_MATCH_1}")
        endif ()
        if (DIR MATCHES "\\.sh$")
          get_filename_component (DIR "${DIR}" PATH)
        endif ()
        if (IS_ABSOLUTE "${DIR}")
          list (APPEND BASHPATH "${DIR}")
        else ()
          list (APPEND BASHPATH "$__DIR__/${DIR}")
        endif ()
      endforeach ()
      if (BASHPATH)
        list (REMOVE_DUPLICATES BASHPATH)
        list (APPEND BASHPATH "$BASHPATH")
        basis_list_to_delimited_string (BASHPATH ":" NOAUTOQUOTE ${BASHPATH})
        set (BASH_CODE "${BASH_CODE}; BASHPATH=\"${BASHPATH}\"")
      endif ()
      basis_remove_blank_line (SCRIPT) # remove a blank line therefore
      set (SCRIPT "${BASH_CODE} # <-- added by BASIS\n${SCRIPT}")
    endif ()
    # replace shebang directive
    if (PYTHON_EXECUTABLE AND ARGN_LANGUAGE MATCHES "PYTHON")
      if (WIN32)
        set (SHEBANG "@setlocal enableextensions & \"${PYTHON_EXECUTABLE}\" -x \"%~f0\" %* & goto :EOF")
      else ()
        set (SHEBANG "#! ${PYTHON_EXECUTABLE}")
      endif ()
    elseif (JYTHON_EXECUTABLE AND ARGN_LANGUAGE MATCHES "JYTHON")
      if (WIN32)
        set (SHEBANG "@setlocal enableextensions & \"${JYTHON_EXECUTABLE}\" -x \"%~f0\" %* & goto :EOF")
      else ()
        # Attention: It is IMPORTANT to not use "#! <interpreter>" even if the <interpreter>
        #            is given as full path in case of jython. Otherwise, the Jython executable
        #            fails to execute from within a Python script using the os.system(),
        #            subprocess.popen(), subprocess.call() or similar function!
        #            Don't ask me for an explanation, but possibly the used shell otherwise does
        #            not recognize the shebang as being valid. Using /usr/bin/env helps out here,
        #            -schuha
        set (SHEBANG "#! /usr/bin/env ${JYTHON_EXECUTABLE}")
      endif ()
    elseif (PERL_EXECUTABLE AND ARGN_LANGUAGE MATCHES "PERL")
      if (WIN32)
        set (SHEBANG "@goto = \"START_OF_BATCH\" ;\n@goto = ();")
        set (SCRIPT "${SCRIPT}\n\n__END__\n\n:\"START_OF_BATCH\"\n@\"${PERL_EXECUTABLE}\" -w -S \"%~f0\" %*")
      else ()
        set (SHEBANG "#! ${PERL_EXECUTABLE} -w")
      endif ()
    elseif (BASH_EXECUTABLE AND ARGN_LANGUAGE MATCHES "BASH")
      set (SHEBANG "#! ${BASH_EXECUTABLE}")
    endif ()
    # add (modified) shebang directive again
    if (SHEBANG)
      set (SCRIPT "${SHEBANG}\n${SCRIPT}")
    endif ()
    # write configured script
    file (WRITE "${_OUTPUT_FILE}" "${SCRIPT}")
    # make script executable on Unix
    if (UNIX AND NOT ARGN_DESTINATION)
      execute_process (COMMAND /bin/chmod +x "${_OUTPUT_FILE}")
    endif ()
  # --------------------------------------------------------------------------
  # configure module script
  else ()
    # configure module - do not use configure_file() as it will not update the
    #                    file if nothing has changed. the update of the modification
    #                    time is however in particular required for the
    #                    configure_script.cmake build command which uses this
    #                    function to build script targets. otherwise, the custom
    #                    build command is reexecuted only because the output files
    #                    never appear to be more recent than the dependencies
    file (READ "${_INPUT_FILE}" SCRIPT)
    if (NOT ARGN_COPYONLY)
      string (CONFIGURE "${SCRIPT}" SCRIPT @ONLY)
    endif ()
    file (WRITE "${_OUTPUT_FILE}" "${SCRIPT}")
    # compile module if requested
    if (ARGN_COMPILE)
      if (ARGN_LANGUAGE MATCHES "PYTHON" AND PYTHON_EXECUTABLE)
        basis_get_compiled_file (CFILE "${_OUTPUT_FILE}" PYTHON)
        execute_process (COMMAND "${PYTHON_EXECUTABLE}" -E -c "import py_compile; py_compile.compile('${_OUTPUT_FILE}', '${CFILE}')")
        basis_compile_python_modules_for_jython (RV)
        if (RV)
          basis_get_compiled_jython_file_of_python_module (CFILE "${_OUTPUT_FILE}")
          get_filename_component (CDIR "${CFILE}" PATH)
          file (MAKE_DIRECTORY "${CDIR}")
          execute_process (COMMAND "${JYTHON_EXECUTABLE}" -c "import py_compile; py_compile.compile('${_OUTPUT_FILE}', '${CFILE}')")
        endif ()
      elseif (ARGN_LANGUAGE MATCHES "JYTHON" AND JYTHON_EXECUTABLE)
        basis_get_compiled_file (CFILE "${_OUTPUT_FILE}" JYTHON)
        execute_process (COMMAND "${JYTHON_EXECUTABLE}" -c "import py_compile; py_compile.compile('${_OUTPUT_FILE}', '${CFILE}')")
      endif ()
    endif ()
  endif ()
endfunction ()

# ----------------------------------------------------------------------------
## @brief Get type name of target.
#
# @param [out] TYPE       The target's type name or NOTFOUND.
# @param [in]  TARGET_UID The UID of the target.
function (_basis_get_target_type TYPE TARGET_UID)
  get_target_property (IMPORTED ${TARGET_UID} IMPORTED)
  if (IMPORTED)
    get_target_property (TYPE_OUT ${TARGET_UID} TYPE)
  else ()
    get_target_property (TYPE_OUT ${TARGET_UID} BASIS_TYPE)
    if (NOT TYPE_OUT)
      get_target_property (TYPE_OUT ${TARGET_UID} TYPE)
    endif ()
  endif ()
  set ("${TYPE}" "${TYPE_OUT}" PARENT_SCOPE)
endfunction ()

# ----------------------------------------------------------------------------
## @brief Get type name of target.
#
# @param [out] TYPE        The target's type name or NOTFOUND.
# @param [in]  TARGET_NAME The name of the target.
function (basis_get_target_type TYPE TARGET_NAME)
  basis_get_target_uid (TARGET_UID "${TARGET_NAME}")
  if (TARGET ${TARGET_UID})
    _basis_get_target_type(TYPE_OUT ${TARGET_UID})
  else ()
    set (TYPE_OUT "NOTFOUND")
  endif ()
  set ("${TYPE}" "${TYPE_OUT}" PARENT_SCOPE)
endfunction ()

# ----------------------------------------------------------------------------
## @brief Get location of build target output file(s).
#
# This convenience function can be used to get the full path of the output
# file(s) generated by a given build target. It is similar to the read-only
# @c LOCATION property of CMake targets and should be used instead of
# reading this porperty. In case of scripted libraries, this function returns
# the path of the root directory of the library that has to be added to the
# module search path.
#
# @note If the target is a binary built from C++ source files and the CMake
#       generator is an IDE such as Visual Studio or Xcode, the absolute
#       directory of the target location ends with the generator expression
#       "/$<${BASIS_GE_CONFIG}>" which is to be substituted by the respective
#       build configuration.
#
# @param [out] VAR         Path of build target output file.
# @param [in]  TARGET_NAME Name of build target.
# @param [in]  PART        Which file name component of the @c LOCATION
#                          property to return. See get_filename_component().
#                          If POST_INSTALL_RELATIVE is given as argument,
#                          @p VAR is set to the path of the installed file
#                          relative to the installation prefix. Similarly,
#                          POST_INSTALL sets @p VAR to the absolute path
#                          of the installed file post installation.
#
# @returns Path of output file similar to @c LOCATION property of CMake targets.
#
# @sa http://www.cmake.org/cmake/help/cmake-2-8-docs.html#prop_tgt:LOCATION
function (basis_get_target_location VAR TARGET_NAME PART)
  basis_get_target_uid (TARGET_UID "${TARGET_NAME}")
  if (TARGET "${TARGET_UID}")
    _basis_get_target_name (TARGET_NAME "${TARGET_UID}")
    _basis_get_target_type (TYPE        "${TARGET_UID}")
    get_target_property (IMPORTED ${TARGET_UID} IMPORTED)
    # ------------------------------------------------------------------------
    # imported targets
    #
    # Note: This might not be required though as even custom executable
    #       and library targets can be imported using CMake's
    #       add_executable(<NAME> IMPORTED) and add_library(<NAME> <TYPE> IMPORTED)
    #       commands. Such executable can, for example, also be a BASH
    #       script built by basis_add_script().
    if (IMPORTED)
      # 1. Try IMPORTED_LOCATION_<CMAKE_BUILD_TYPE>
      if (CMAKE_BUILD_TYPE)
        string (TOUPPER "${CMAKE_BUILD_TYPE}" U)
      else ()
        set (U "NOCONFIG")
      endif ()
      get_target_property (LOCATION ${TARGET_UID} IMPORTED_LOCATION_${U})
      # 2. Try IMPORTED_LOCATION
      if (NOT LOCATION)
        get_target_property (LOCATION ${TARGET_UID} IMPORTED_LOCATION)
      endif ()
      # 3. Prefer Release over all other configurations
      if (NOT LOCATION)
        get_target_property (LOCATION ${TARGET_UID} IMPORTED_LOCATION_RELEASE)
      endif ()
      # 4. Just use any of the imported configurations
      if (NOT LOCATION)
        get_property (CONFIGS TARGET ${TARGET_UID} PROPERTY IMPORTED_CONFIGURATIONS)
        foreach (C IN LISTS CONFIGS)
          string (TOUPPER "${C}" C)
          get_target_property (LOCATION ${TARGET_UID} IMPORTED_LOCATION_${C})
          if (LOCATION)
            break ()
          endif ()
        endforeach ()
      endif ()
      # make path relative to CMAKE_INSTALL_PREFIX if POST_INSTALL_RELATIVE given
      if (LOCATION AND ARGV2 MATCHES "POST_INSTALL_RELATIVE")
        file (RELATIVE_PATH LOCATION "${CMAKE_INSTALL_PREFIX}" "${LOCATION}")
      endif ()
    # ------------------------------------------------------------------------
    # non-imported targets
    else ()
      # Attention: The order of the matches/if cases matters here!
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # scripts
      if (TYPE MATCHES "^SCRIPT_(EXECUTABLE|MODULE)$")
        if (PART MATCHES "POST_INSTALL")
          get_target_property (DIRECTORY ${TARGET_UID} INSTALL_DIRECTORY)
        else ()
          get_target_property (DIRECTORY ${TARGET_UID} OUTPUT_DIRECTORY)
          if (DIRECTORY AND CMAKE_GENERATOR MATCHES "Visual Studio|Xcode")
            set (DIRECTORY "${DIRECTORY}/$<${BASIS_GE_CONFIG}>")
          endif ()
        endif ()
        get_target_property (FNAME ${TARGET_UID} OUTPUT_NAME)
      elseif (TYPE STREQUAL "^SCRIPT_MODULE$")
        if (PART MATCHES "POST_INSTALL")
          get_target_property (DIRECTORY ${TARGET_UID} INSTALL_DIRECTORY)
        else ()
          get_target_property (COMPILE   ${TARGET_UID} COMPILE)
          get_target_property (DIRECTORY ${TARGET_UID} OUTPUT_DIRECTORY)
          if (DIRECTORY AND (COMPILE OR CMAKE_GENERATOR MATCHES "Visual Studio|Xcode"))
            set (DIRECTORY "${DIRECTORY}/$<${BASIS_GE_CONFIG}>")
          endif ()
        endif ()
        get_target_property (FNAME ${TARGET_UID} OUTPUT_NAME)
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # libraries
      elseif (TYPE MATCHES "LIBRARY|MODULE|MEX")
        if (TYPE MATCHES "STATIC")
          if (PART MATCHES "POST_INSTALL")
            get_target_property (DIRECTORY ${TARGET_UID} ARCHIVE_INSTALL_DIRECTORY)
          else ()
            get_target_property (DIRECTORY ${TARGET_UID} ARCHIVE_OUTPUT_DIRECTORY)
          endif ()
          get_target_property (FNAME ${TARGET_UID} ARCHIVE_OUTPUT_NAME)
        else ()
          if (PART MATCHES "POST_INSTALL")
            get_target_property (DIRECTORY ${TARGET_UID} LIBRARY_INSTALL_DIRECTORY)
          else ()
            get_target_property (DIRECTORY ${TARGET_UID} LIBRARY_OUTPUT_DIRECTORY)
          endif ()
          get_target_property (FNAME ${TARGET_UID} LIBRARY_OUTPUT_NAME)
        endif ()
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # executables
      else ()
        if (PART MATCHES "POST_INSTALL")
          get_target_property (DIRECTORY ${TARGET_UID} RUNTIME_INSTALL_DIRECTORY)
        else ()
          get_target_property (DIRECTORY ${TARGET_UID} RUNTIME_OUTPUT_DIRECTORY)
        endif ()
        get_target_property (FNAME ${TARGET_UID} RUNTIME_OUTPUT_NAME)
      endif ()
      if (DIRECTORY MATCHES "NOTFOUND")
        message (FATAL_ERROR "Failed to get directory of ${TYPE} ${TARGET_UID}!"
                             " Check implementation of basis_get_target_location()"
                             " and make sure that the required *INSTALL_DIRECTORY"
                             " property is set on the target!")
      endif ()
      if (DIRECTORY)
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # get output name of built file (if applicable)
        if (NOT FNAME)
          get_target_property (FNAME ${TARGET_UID} OUTPUT_NAME)
        endif ()
        if (NOT "^${TYPE}$" STREQUAL "^SCRIPT_LIBRARY$")
          get_target_property (PREFIX ${TARGET_UID} PREFIX)
          get_target_property (SUFFIX ${TARGET_UID} SUFFIX)
          if (FNAME)
            set (TARGET_FILE "${FNAME}")
          else ()
            set (TARGET_FILE "${TARGET_NAME}")
          endif ()
          if (PREFIX)
            set (TARGET_FILE "${PREFIX}${TARGET_FILE}")
          endif ()
          if (SUFFIX)
            set (TARGET_FILE "${TARGET_FILE}${SUFFIX}")
          elseif (WIN32 AND "^${TYPE}$" STREQUAL "^EXECUTABLE$")
            set (TARGET_FILE "${TARGET_FILE}.exe")
          endif ()
        endif ()
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # prepend $<CONFIG> "generator expression" for non-custom binaries
        # when built with an IDE such as Visual Studio or Xcode
        if ("^${TYPE}$" STREQUAL "^EXECUTABLE$" OR "^${TYPE}$" STREQUAL "^LIBRARY$")
          if (NOT PART MATCHES "INSTALL")
            if (CMAKE_GENERATOR MATCHES "Visual Studio|Xcode")
              set (DIRECTORY "${DIRECTORY}/$<${BASIS_GE_CONFIG}>")
            endif ()
          endif ()
        endif ()
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # assemble final path
        if (PART MATCHES "POST_INSTALL_RELATIVE")
          if (IS_ABSOLUTE "${DIRECTORY}")
            file (RELATIVE_PATH DIRECTORY "${CMAKE_INSTALL_PREFIX}" "${DIRECTORY}")
            if (NOT DIRECTORY)
              set (DIRECTORY ".")
            endif ()
          endif ()
        elseif (PART MATCHES "POST_INSTALL")
          if (NOT IS_ABSOLUTE "${DIRECTORY}")
            set (DIRECTORY "${CMAKE_INSTALL_PREFIX}/${DIRECTORY}")
          endif ()
        endif ()
        if (TARGET_FILE)
          set (LOCATION "${DIRECTORY}/${TARGET_FILE}")
        else ()
          set (LOCATION "${DIRECTORY}")
        endif ()
      else ()
        set (LOCATION "${DIRECTORY}")
      endif ()
    endif ()
    # get filename component
    if (LOCATION AND PART MATCHES "(^|_)(PATH|NAME|NAME_WE)$")
      get_filename_component (LOCATION "${LOCATION}" "${CMAKE_MATCH_2}")
    endif ()
  else ()
    message (FATAL_ERROR "basis_get_target_location(): Unknown target ${TARGET_UID}")
  endif ()
  # return
  set ("${VAR}" "${LOCATION}" PARENT_SCOPE)
endfunction ()

# ----------------------------------------------------------------------------
## @brief Get link libraries/dependencies of (imported) target.
#
# This function recursively adds the dependencies of the dependencies as well
# and returns them together with the list of the direct link dependencies.
# Moreover, for script targets, if any of the dependencies uses the BASIS
# utilities for the given language (@c BASIS_UTILITIES property), the
# corresponding utilities library is added to the list of dependencies.
# Note that therefore the BASIS utilities targets have to be added already,
# which is only the case during the finalization of script targets.
#
# @param [out] LINK_DEPENDS List of all link dependencies. In case of scripts,
#                           the dependencies are the required modules or
#                           paths to required packages, respectively.
# @param [in]  TARGET_NAME  Name of the target.
function (basis_get_target_link_libraries LINK_DEPENDS TARGET_NAME)
  basis_get_target_uid (TARGET_UID "${TARGET_NAME}")
  if (NOT TARGET "${TARGET_UID}")
    message (FATAL_ERROR "basis_get_target_link_libraries(): Unknown target: ${TARGET_UID}")
  endif ()
  if (BASIS_DEBUG AND BASIS_VERBOSE)
    message ("** basis_get_target_link_libraries():")
    message ("**   TARGET_NAME:     ${TARGET_NAME}")
    message ("**   CURRENT_DEPENDS: ${ARGN}")
  endif ()
  # get type of target
  get_target_property (BASIS_TYPE ${TARGET_UID} BASIS_TYPE)
  # get direct link dependencies of target
  get_target_property (IMPORTED ${TARGET_UID} IMPORTED)
  if (IMPORTED)
    # 1. Try IMPORTED_LINK_INTERFACE_LIBRARIES_<CMAKE_BUILD_TYPE>
    if (CMAKE_BUILD_TYPE)
      string (TOUPPER "${CMAKE_BUILD_TYPE}" U)
    else ()
      set (U "NOCONFIG")
    endif ()
    get_target_property (DEPENDS ${TARGET_UID} "IMPORTED_LINK_INTERFACE_LIBRARIES_${U}")
    # 2. Try IMPORTED_LINK_INTERFACE_LIBRARIES
    if (NOT DEPENDS)
      get_target_property (DEPENDS ${TARGET_UID} "IMPORTED_LINK_INTERFACE_LIBRARIES")
    endif ()
    # 3. Prefer Release over all other configurations
    if (NOT DEPENDS)
      get_target_property (DEPENDS ${TARGET_UID} "IMPORTED_LINK_INTERFACE_LIBRARIES_RELEASE")
    endif ()
    # 4. Just use any of the imported configurations
    if (NOT DEPENDS)
      get_property (CONFIGS TARGET "${TARGET_UID}" PROPERTY IMPORTED_CONFIGURATIONS)
      foreach (C IN LISTS CONFIGS)
        get_target_property (DEPENDS ${TARGET_UID} "IMPORTED_LINK_INTERFACE_LIBRARIES_${C}")
        if (DEPENDS)
          break ()
        endif ()
      endforeach ()
    endif ()
  # otherwise, get LINK_DEPENDS property value
  elseif (BASIS_TYPE MATCHES "^EXECUTABLE$|^(SHARED|STATIC|MODULE)_LIBRARY$")
    get_target_property (DEPENDS ${TARGET_UID} BASIS_LINK_DEPENDS)
  else ()
    get_target_property (DEPENDS ${TARGET_UID} LINK_DEPENDS)
  endif ()
  if (NOT DEPENDS)
    set (DEPENDS)
  endif ()
  # prepend BASIS utilities if used (and added)
  if (BASIS_TYPE MATCHES "SCRIPT")
    set (BASIS_UTILITIES_TARGETS)
    foreach (UID IN ITEMS ${TARGET_UID} ${DEPENDS})
      if (TARGET "${UID}")
        get_target_property (BASIS_UTILITIES ${UID} BASIS_UTILITIES)
        get_target_property (LANGUAGE        ${UID} LANGUAGE)
        if (BASIS_UTILITIES)
          set (BASIS_UTILITIES_TARGET)
          if (LANGUAGE MATCHES "[JP]YTHON")
            basis_get_source_target_name (BASIS_UTILITIES_TARGET "basis.py" NAME)
          elseif (LANGUAGE MATCHES "PERL")
            basis_get_source_target_name (BASIS_UTILITIES_TARGET "Basis.pm" NAME)
          elseif (LANGUAGE MATCHES "BASH")
            basis_get_source_target_name (BASIS_UTILITIES_TARGET "basis.sh" NAME)
          endif ()
          if (BASIS_UTILITIES_TARGET)
            basis_get_target_uid (BASIS_UTILITIES_TARGET ${BASIS_UTILITIES_TARGET})
          endif ()
          if (TARGET ${BASIS_UTILITIES_TARGET})
            list (APPEND BASIS_UTILITIES_TARGETS ${BASIS_UTILITIES_TARGET})
          endif ()
        endif ()
      endif ()
    endforeach ()
    if (BASIS_UTILITIES_TARGETS)
      list (INSERT DEPENDS 0 ${BASIS_UTILITIES_TARGETS})
    endif ()
  endif ()
  # convert target names to UIDs
  set (_DEPENDS)
  foreach (LIB IN LISTS DEPENDS)
    basis_get_target_uid (UID "${LIB}")
    if (TARGET ${UID})
      list (APPEND _DEPENDS "${UID}")
    else ()
      list (APPEND _DEPENDS "${LIB}")
    endif ()
  endforeach ()
  set (DEPENDS "${_DEPENDS}")
  unset (_DEPENDS)
  # recursively add link dependencies of dependencies
  # TODO implement it non-recursively for better performance
  foreach (LIB IN LISTS DEPENDS)
    if (TARGET ${LIB})
      list (FIND ARGN "${LIB}" IDX) # avoid recursive loop
      if (IDX EQUAL -1)
        basis_get_target_link_libraries (LIB_DEPENDS ${LIB} ${ARGN} ${DEPENDS})
        list (APPEND DEPENDS ${LIB_DEPENDS})
      endif ()
    endif ()
  endforeach ()
  # remove duplicate entries
  if (DEPENDS)
    list (REMOVE_DUPLICATES DEPENDS)
  endif ()
  # return
  set (${LINK_DEPENDS} "${DEPENDS}" PARENT_SCOPE)
endfunction ()

# ============================================================================
# generator expressions
# ============================================================================

# ----------------------------------------------------------------------------
## @brief Process generator expressions in arguments.
#
# This command evaluates the $&lt;TARGET_FILE:tgt&gt; and related generator
# expressions also for custom targets such as scripts and MATLAB Compiler
# targets. For other generator expressions whose argument is a target name,
# this function replaces the target name by the target UID, i.e., the actual
# CMake target name such that the expression can be evaluated by CMake.
# The following generator expressions are directly evaluated by this function:
# <table border=0>
#   <tr>
#     @tp <b><tt>$&lt;TARGET_FILE:tgt&gt;</tt></b> @endtp
#     <td>Absolute file path of built target.</td>
#   </tr>
#   <tr>
#     @tp <b><tt>$&lt;TARGET_FILE_POST_INSTALL:tgt&gt;</tt></b> @endtp
#     <td>Absolute path of target file after installation using the
#         current @c CMAKE_INSTALL_PREFIX.</td>
#   </tr>
#   <tr>
#     @tp <b><tt>$&lt;TARGET_FILE_POST_INSTALL_RELATIVE:tgt&gt;</tt></b> @endtp
#     <td>Path of target file after installation relative to @c CMAKE_INSTALL_PREFIX.</td>
#   </tr>
# </table>
# Additionally, the suffix <tt>_NAME</tt> or <tt>_DIR</tt> can be appended
# to the name of each of these generator expressions to get only the basename
# of the target file including the extension or the corresponding directory
# path, respectively.
#
# Generator expressions are in particular supported by basis_add_test().
#
# @param [out] ARGS Name of output list variable.
# @param [in]  ARGN List of arguments to process.
#
# @sa basis_add_test()
# @sa http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:add_test
function (basis_process_generator_expressions ARGS)
  set (ARGS_OUT)
  foreach (ARG IN LISTS ARGN)
    string (REGEX MATCHALL "\\$<.*TARGET.*:.*>" EXPRS "${ARG}")
    foreach (EXPR IN LISTS EXPRS)
      if (EXPR MATCHES "\\$<(.*):(.*)>")
        set (EXPR_NAME   "${CMAKE_MATCH_1}")
        set (TARGET_NAME "${CMAKE_MATCH_2}")
        # TARGET_FILE* expression, including custom targets
        if (EXPR_NAME MATCHES "^TARGET_FILE(.*)")
          if (NOT CMAKE_MATCH_1)
            set (CMAKE_MATCH_1 "ABSOLUTE")
          endif ()
          string (REGEX REPLACE "^_" "" PART "${CMAKE_MATCH_1}")
          basis_get_target_location (ARG "${TARGET_NAME}" ${PART})
        # other generator expression supported by CMake
        # only replace target name, but do not evaluate expression
        else ()
          basis_get_target_uid (TARGET_UID "${CMAKE_MATCH_2}")
          string (REPLACE "${EXPR}" "$<${CMAKE_MATCH_1}:${TARGET_UID}>" ARG "${ARG}")
        endif ()
        if (BASIS_DEBUG AND BASIS_VERBOSE)
          message ("** basis_process_generator_expressions():")
          message ("**   Expression:  ${EXPR}")
          message ("**   Keyword:     ${EXPR_NAME}")
          message ("**   Argument:    ${TARGET_NAME}")
          message ("**   Replaced by: ${ARG}")
        endif ()
      endif ()
    endforeach ()
    list (APPEND ARGS_OUT "${ARG}")
  endforeach ()
  set (${ARGS} "${ARGS_OUT}" PARENT_SCOPE)
endfunction ()


##
#  @brief basis_append_to_each takes an input list and appends a single element to each item in that list and appends it to the output list.
#                For example, this is useful for adding relative paths to the end of a list of paths.
#
#  @param OUTPUT_LIST Name of list that will be filled with appended names.
#  @param INPUT_LIST  Name of list that contains items to have text appended.
#  @param ITEM_TO_APPEND text to append to each item in the input list.
#
function(basis_append_to_each OUTPUT_LIST INPUT_LIST ITEM_TO_APPEND)
  foreach(PATH IN LISTS ${INPUT_LIST})
    list(APPEND ${OUTPUT_LIST} ${PATH}${ITEM_TO_APPEND} )
  endforeach()

  if(${OUTPUT_LIST})
    set(${OUTPUT_LIST} ${${OUTPUT_LIST}} PARENT_SCOPE)
  endif()
endfunction()


## @}
# end of Doxygen group
