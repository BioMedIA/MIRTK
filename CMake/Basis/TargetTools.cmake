# ============================================================================
# Copyright (c) 2011-2012 University of Pennsylvania
# Copyright (c) 2013-2016 Andreas Schuh
# All rights reserved.
#
# See COPYING file for license information or visit
# https://cmake-basis.github.io/download.html#license
# ============================================================================

##############################################################################
# @file  TargetTools.cmake
# @brief Functions and macros to add executable and library targets.
#
# @ingroup CMakeTools
##############################################################################

if (__BASIS_TARGETTOOLS_INCLUDED)
  return ()
else ()
  set (__BASIS_TARGETTOOLS_INCLUDED TRUE)
endif ()


## @addtogroup CMakeUtilities
#  @{


# ============================================================================
# properties
# ============================================================================

# ----------------------------------------------------------------------------
## @brief Set properties on a target.
#
# This function replaces CMake's
# <a href="http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:set_target_properties">
# set_target_properties()</a> command and extends its functionality.
# In particular, it maps the given target names to the corresponding target UIDs.
#
# @note If @c BASIS_USE_TARGET_UIDS is @c OFF and is not required by a project,
#       it is recommended to use set_target_properties() instead (note that
#       set_target_properties is overriden by the ImportTools.cmake module of BASIS).
#       This will break the build configuration scripts when @c BASIS_USE_TARGET_UIDS
#       is set to @c ON later. It should thus only be used if the project will
#       never use the target UID feature of BASIS. A project can possibly define
#       a global macro which either calls set_target_properties or
#       basis_set_target_properties. But be aware of the related CMake bugs
#       which prevent basis_set_target_properties to do the same already.
#       ARGV/ARGN do not preserve empty arguments nor list arguments!
#
# @note Due to a bug in CMake (http://www.cmake.org/Bug/view.php?id=12303),
#       except of the first property given directly after the @c PROPERTIES keyword,
#       only properties listed in @c BASIS_PROPERTIES_ON_TARGETS can be set.
#
# @param [in] ARGN List of arguments. See
#                  <a href="http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:set_target_properties">
#                  set_target_properties()</a>.
#
# @returns Sets the specified properties on the given target.
#
# @sa http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:set_target_properties
#
# @ingroup CMakeAPI
function (basis_set_target_properties)
  # convert target names to UIDs
  set (TARGET_UIDS)
  list (LENGTH ARGN N)
  if (N EQUAL 0)
    message (FATAL_ERROR "basis_set_target_properties(): Missing arguments!")
  endif ()
  list (GET ARGN 0 ARG)
  while (NOT ARG MATCHES "^PROPERTIES$")
    basis_get_target_uid (TARGET_UID "${ARG}")
    list (APPEND TARGET_UIDS "${TARGET_UID}")
    list (REMOVE_AT ARGN 0)
    list (LENGTH ARGN N)
    if (N EQUAL 0)
      break ()
    else ()
      list (GET ARGN 0 ARG)
    endif ()
  endwhile ()
  if (NOT ARG MATCHES "^PROPERTIES$")
    message (FATAL_ERROR "Missing PROPERTIES argument!")
  elseif (NOT TARGET_UIDS)
    message (FATAL_ERROR "No target specified!")
  endif ()
  # remove PROPERTIES keyword
  list (REMOVE_AT ARGN 0)
  math (EXPR N "${N} - 1")
  # set targets properties
  #
  # Note: By iterating over the properties, the empty property values
  #       are correctly passed on to CMake's set_target_properties()
  #       command, while
  #       _set_target_properties(${TARGET_UIDS} PROPERTIES ${ARGN})
  #       (erroneously) discards the empty elements in ARGN.
  if (BASIS_DEBUG)
    message ("** basis_set_target_properties:")
    message ("**   Target(s):  ${TARGET_UIDS}")
    message ("**   Properties: [${ARGN}]")
  endif ()
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
    if (PROPERTY MATCHES "^$") # remember: STREQUAL is buggy and evil!
      message (FATAL_ERROR "Empty property name given!")
    endif ()
    # set target property
    if (COMMAND _set_target_properties) # i.e. ImportTools.cmake included
      _set_target_properties (${TARGET_UIDS} PROPERTIES ${PROPERTY} "${VALUE}")
    else ()
      set_target_properties (${TARGET_UIDS} PROPERTIES ${PROPERTY} "${VALUE}")
    endif ()
  endwhile ()
  # make sure that every property had a corresponding value
  if (NOT N EQUAL 0)
    message (FATAL_ERROR "No value given for target property ${ARGN}")
  endif ()
endfunction ()

# ----------------------------------------------------------------------------
## @brief Get value of property set on target.
#
# This function replaces CMake's
# <a href="http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:set_target_properties">
# get_target_properties()</a> command and extends its functionality.
# In particular, it maps the given @p TARGET_NAME to the corresponding target UID.
#
# @param [out] VAR         Name of output variable.
# @param [in]  TARGET_NAME Name of build target.
# @param [in]  ARGN        Remaining arguments for
#                          <a href="http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:get_target_properties">
#                          get_target_properties()</a>.
#
# @returns Sets @p VAR to the value of the requested property.
#
# @sa http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:get_target_property
#
# @ingroup CMakeAPI
function (basis_get_target_property VAR TARGET_NAME)
  basis_get_target_uid (TARGET_UID "${TARGET_NAME}")
  get_target_property (VALUE "${TARGET_UID}" ${ARGN})
  set (${VAR} "${VALUE}" PARENT_SCOPE)
endfunction ()

# ============================================================================
# definitions
# ============================================================================

# ----------------------------------------------------------------------------
## @brief Add compile definitions.
#
# This function replaces CMake's
# <a href="http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:add_definitions">
# add_definitions()</a> command.
#
# @param [in] ARGN List of arguments for
#                  <a href="http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:add_definitions">
#                  add_definitions()</a>.
#
# @returns Adds the given definitions.
#
# @sa http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:add_definitions
#
# @ingroup CMakeAPI
function (basis_add_definitions)
  add_definitions (${ARGN})
endfunction ()

# ----------------------------------------------------------------------------
## @brief Remove previously added compile definitions.
#
# This function replaces CMake's
# <a href="http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:remove_definitions">
# remove_definitions()</a> command.
#
# @param [in] ARGN List of arguments for
#                  <a href="http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:remove_definitions">
#                  remove_definitions()</a>.
#
# @returns Removes the specified definitions.
#
# @sa http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:remove_definition
#
# @ingroup CMakeAPI
function (basis_remove_definitions)
  remove_definitions (${ARGN})
endfunction ()

# ============================================================================
# directories
# ============================================================================

# ----------------------------------------------------------------------------
## @brief Add directories to search path for include files.
#
# Overwrites CMake's
# <a href="http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:include_directories">
# include_directories()</a> command. This is required because the
# basis_include_directories() function is not used by other projects in their
# package use files. Therefore, this macro is an alias for
# basis_include_directories().
#
# @param [in] ARGN List of arguments for basis_include_directories().
#
# @returns Adds the given paths to the search path for include files.
#
# @sa http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:include_directories
macro (include_directories)
  basis_include_directories (${ARGN})
endmacro ()

# ----------------------------------------------------------------------------
## @brief Add directories to search path for include files.
#
# This function replaces CMake's
# <a href="http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:include_directories">
# include_directories()</a> command. Besides invoking CMake's internal command
# with the given arguments, it updates the @c PROJECT_INCLUDE_DIRECTORIES
# property on the current project (see basis_set_project_property()). This list
# contains a list of all include directories used by a project, regardless of
# the directory in which the basis_include_directories() function was used.
#
# @param ARGN List of arguments for
#             <a href="http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:include_directories">
#             include_directories()</a> command.
#
# @returns Nothing.
#
# @sa http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:include_directories
#
# @ingroup CMakeAPI
function (basis_include_directories)
  # CMake's include_directories ()
  _include_directories (${ARGN})

  # parse arguments
  CMAKE_PARSE_ARGUMENTS (ARGN "AFTER;BEFORE;SYSTEM" "" "" ${ARGN})

  # make relative paths absolute
  set (DIRS)
  foreach (P IN LISTS ARGN_UNPARSED_ARGUMENTS)
    if (NOT P MATCHES "^\\$<") # preserve generator expressions
      get_filename_component (P "${P}" ABSOLUTE)
    endif ()
    list (APPEND DIRS "${P}")
  endforeach ()

  if (DIRS)
    # append directories to "global" list of include directories
    basis_get_project_property (INCLUDE_DIRS PROPERTY PROJECT_INCLUDE_DIRS)
    if (BEFORE)
      list (INSERT INCLUDE_DIRS 0 ${DIRS})
    else ()
      list (APPEND INCLUDE_DIRS ${DIRS})
    endif ()
    if (INCLUDE_DIRS)
      list (REMOVE_DUPLICATES INCLUDE_DIRS)
    endif ()
    if (BASIS_DEBUG)
      message ("** basis_include_directories():")
      if (BEFORE)
        message ("**    Add before:  ${DIRS}")
      else ()
        message ("**    Add after:   ${DIRS}")
      endif ()
      if (BASIS_VERBOSE)
        message ("**    Directories: ${INCLUDE_DIRS}")
      endif ()
    endif ()
    basis_set_project_property (PROPERTY PROJECT_INCLUDE_DIRS ${INCLUDE_DIRS})
  endif ()
endfunction ()

# ----------------------------------------------------------------------------
## @brief Add directories to search path for libraries.
#
# Overwrites CMake's
# <a href="http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:link_directories">
# link_directories()</a> command. This is required because the
# basis_link_directories() function is not used by other projects in their
# package use files. Therefore, this macro is an alias for
# basis_link_directories().
#
# @param [in] ARGN List of arguments for basis_link_directories().
#
# @returns Adds the given paths to the search path for libraries.
#
# @sa http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:link_directories
macro (link_directories)
  basis_link_directories (${ARGN})
endmacro ()

# ----------------------------------------------------------------------------
## @brief Add directories to search path for libraries.
#
# This function replaces CMake's
# <a href="http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:link_directories">
# link_directories()</a> command. Even though this function yet only invokes
# CMake's internal command, it should be used in BASIS projects to enable the
# extension of this command's functionality as part of BASIS if required.
#
# @param [in] ARGN List of arguments for
#                  <a href="http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:link_directories">
#                  link_directories()</a>.
#
# @returns Adds the given paths to the search path for libraries.
#
# @sa http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:link_directories
#
# @ingroup CMakeAPI
function (basis_link_directories)
  # CMake's link_directories() command
  _link_directories (${ARGN})
  # make relative paths absolute
  set (DIRS)
  foreach (P IN LISTS ARGN)
    if (NOT P MATCHES "^\\$<") # preserve generator expressions
      get_filename_component (P "${P}" ABSOLUTE)
    endif ()
    list (APPEND DIRS "${P}")
  endforeach ()
  if (DIRS)
    # append directories to "global" list of link directories
    basis_get_project_property (LINK_DIRS PROPERTY PROJECT_LINK_DIRS)
    list (APPEND LINK_DIRS ${DIRS})
    if (LINK_DIRS)
      list (REMOVE_DUPLICATES LINK_DIRS)
    endif ()
    if (BASIS_DEBUG)
      message ("** basis_link_directories():")
      message ("**    Add:         [${DIRS}]")
      if (BASIS_VERBOSE)
        message ("**    Directories: [${LINK_DIRS}]")
      endif ()
    endif ()
    basis_set_project_property (PROPERTY PROJECT_LINK_DIRS "${LINK_DIRS}")
    # if the directories are added by an external project's <Pkg>Use.cmake
    # file which is part of the same superbuild as this project, add the
    # directories further to the list of directories that may be added to
    # the RPATH. see basis_set_target_install_rpath().
    if (BUNDLE_PROJECT) # set in basis_use_package()
      basis_get_project_property (LINK_DIRS PROPERTY BUNDLE_LINK_DIRS)
      list (APPEND LINK_DIRS ${DIRS})
      if (LINK_DIRS)
        list (REMOVE_DUPLICATES LINK_DIRS)
      endif ()
      basis_set_project_property (PROPERTY BUNDLE_LINK_DIRS "${LINK_DIRS}")
    endif ()
  endif ()
endfunction ()

# ============================================================================
# dependencies
# ============================================================================

# ----------------------------------------------------------------------------
## @brief Add dependencies to build target.
#
# This function replaces CMake's
# <a href="http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:add_dependencies">
# add_dependencies()</a> command and extends its functionality.
# In particular, it maps the given target names to the corresponding target UIDs.
#
# @param [in] ARGN Arguments for
#                  <a href="http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:add_dependencies">
#                  add_dependencies()</a>.
#
# @returns Adds the given dependencies of the specified build target.
#
# @sa http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:add_dependencies
#
# @ingroup CMakeAPI
if (BASIS_USE_TARGET_UIDS)
  function (basis_add_dependencies)
    set (ARGS)
    foreach (ARG ${ARGN})
      basis_get_target_uid (UID "${ARG}")
      if (TARGET "${UID}")
        list (APPEND ARGS "${UID}")
      else ()
        list (APPEND ARGS "${ARG}")
      endif ()
    endforeach ()
    add_dependencies (${ARGS})
  endfunction ()
else ()
  macro (basis_add_dependencies)
    add_dependencies (${ARGV})
  endmacro ()
endif ()

# ----------------------------------------------------------------------------
## @brief Add link dependencies to build target.
#
# This function replaces CMake's
# <a href="http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:target_link_libraries">
# target_link_libraries()</a> command.
#
# The main reason for replacing this function is to treat libraries such as
# MEX-files which are supposed to be compiled into a MATLAB executable added
# by basis_add_executable() special. In this case, these libraries are added
# to the LINK_DEPENDS property of the given MATLAB Compiler target. Similarly,
# executable scripts and modules written in a scripting language may depend
# on other modules.
#
# Another reason is the mapping of build target names to fully-qualified
# build target names as used by BASIS (see basis_get_target_uid()).
#
# Only link dependencies added with this function are considered for the setting
# of the INSTALL_RPATH of executable targets (see basis_set_target_install_rpath()).
#
# Example:
# @code
# basis_add_library (MyMEXFunc MEX myfunc.c)
# basis_add_executable (MyMATLABApp main.m)
# basis_target_link_libraries (MyMATLABApp MyMEXFunc OtherMEXFunc.mexa64)
# @endcode
#
# @param [in] TARGET_NAME Name of the target.
# @param [in] ARGN        Link libraries.
#
# @returns Adds link dependencies to the specified build target.
#          For custom targets, the given libraries are added to the
#          @c LINK_DEPENDS property of these targets, in particular.
#
# @sa http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:target_link_libraries
#
# @ingroup CMakeAPI
function (basis_target_link_libraries TARGET_NAME)
  basis_get_target_uid (TARGET_UID "${TARGET_NAME}")
  if (NOT TARGET "${TARGET_UID}")
    message (FATAL_ERROR "basis_target_link_libraries(): Unknown target ${TARGET_UID}.")
  endif ()
  # get type of named target
  get_target_property (BASIS_TYPE ${TARGET_UID} BASIS_TYPE)
  # substitute non-fully qualified target names
  set (ARGS)
  foreach (ARG ${ARGN})
    if ("^${ARG}$" STREQUAL "^basis$")
      get_target_property (LANGUAGE ${TARGET_UID} LANGUAGE)
      if (NOT LANGUAGE OR "^${LANGUAGE}$" STREQUAL "^UNKNOWN$")
        message (FATAL_ERROR "Target ${TARGET_UID} is of unknown LANGUAGE! Cannot add dependency on \"basis\" utilities.")
      endif ()
      basis_add_utilities_library (BASIS_UTILITIES_TARGET ${LANGUAGE})
      list (APPEND ARGS ${BASIS_UTILITIES_TARGET})
      set_target_properties (${TARGET_UID} PROPERTIES BASIS_UTILITIES TRUE)
    else ()
      basis_get_target_uid (UID "${ARG}")
      if (TARGET "${UID}")
        if ("^${UID}$" STREQUAL "^${TARGET_UID}$")
          message (FATAL_ERROR "Cannot add link library ${UID} as dependency of itself!")
        endif ()
        list (APPEND ARGS "${UID}")
      else ()
        list (APPEND ARGS "${ARG}")
      endif ()
    endif ()
  endforeach ()
  # get current link libraries
  if (BASIS_TYPE MATCHES "^EXECUTABLE$|^(SHARED|STATIC|MODULE)_LIBRARY$")
    get_target_property (DEPENDS ${TARGET_UID} BASIS_LINK_DEPENDS)
  else ()
    get_target_property (DEPENDS ${TARGET_UID} LINK_DEPENDS)
  endif ()
  if (NOT DEPENDS)
    set (DEPENDS)
  endif ()
  # note that MCC does itself a dependency check and in case of scripts
  # the basis_get_target_link_libraries() function is used
  if (BASIS_TYPE MATCHES "MCC|SCRIPT")
    list (APPEND DEPENDS ${ARGS})
  # otherwise
  else ()
    list (APPEND DEPENDS ${ARGS})
    # pull implicit dependencies (e.g., ITK uses this)
    set (DEPENDENCY_ADDED 1)
    while (DEPENDENCY_ADDED)
      set (DEPENDENCY_ADDED 0)
      foreach (LIB IN LISTS DEPENDS)
        foreach (LIB_DEPEND IN LISTS ${LIB}_LIB_DEPENDS)
          if (NOT LIB_DEPEND MATCHES "^$|^general$")
            string (REGEX REPLACE "^-l" "" LIB_DEPEND "${LIB_DEPEND}")
            list (FIND DEPENDS ${LIB_DEPEND} IDX)
            if (IDX EQUAL -1)
              list (APPEND DEPENDS ${LIB_DEPEND})
              set (DEPENDENCY_ADDED 1)
            endif ()
          endif ()
        endforeach ()
      endforeach ()
    endwhile ()
  endif ()
  # update LINK_DEPENDS
  if (BASIS_TYPE MATCHES "^EXECUTABLE$|^(SHARED|STATIC|MODULE)_LIBRARY$")
    set_target_properties (${TARGET_UID} PROPERTIES BASIS_LINK_DEPENDS "${DEPENDS}")
    target_link_libraries (${TARGET_UID} ${ARGS})
  else ()
    # FIXME cannot use LINK_DEPENDS for CMake targets as these depends are
    #       otherwise added as is to the Makefile. therefore, consider renaming
    #       LINK_DEPENDS in general to BASIS_LINK_DEPENDS.
    set_target_properties (${TARGET_UID} PROPERTIES LINK_DEPENDS "${DEPENDS}")
  endif ()
endfunction ()

# ============================================================================
# add targets
# ============================================================================

# ----------------------------------------------------------------------------
## @brief Add custom target.
macro (basis_add_custom_target TARGET_NAME)
  basis_check_target_name ("${TARGET_NAME}")
  basis_make_target_uid (_UID "${TARGET_NAME}")
  add_custom_target (${_UID} ${ARGN})
  unset (_UID)
endmacro ()

# ----------------------------------------------------------------------------
## @brief Determine language of source files.
# @sa basis_add_executable(), basis_add_library()
macro (_basis_target_source_language)
  if (NOT ARGN_LANGUAGE)
    basis_get_source_language (ARGN_LANGUAGE ${SOURCES})
    if (ARGN_LANGUAGE MATCHES "AMBIGUOUS|UNKNOWN")
      set (_FILES)
      foreach (SOURCE IN LISTS SOURCES)
        set (_FILES "${_FILES}\n  ${SOURCE}")
      endforeach ()
      if (ARGN_LANGUAGE MATCHES "AMBIGUOUS")
        message (FATAL_ERROR "Target ${TARGET_UID}: Ambiguous source code files! Try to set LANGUAGE manually and make sure that no unknown option was given:${_FILES}")
      elseif (ARGN_LANGUAGE MATCHES "UNKNOWN")
        message (FATAL_ERROR "Target ${TARGET_UID}: Unknown source code language! Try to set LANGUAGE manually and make sure that no unknown option was given:${_FILES}")
      endif ()
    endif ()
  else ()
    string (TOUPPER "${ARGN_LANGUAGE}" ARGN_LANGUAGE)
  endif ()
endmacro ()

# ----------------------------------------------------------------------------
## @brief Add executable target.
#
# This is the main function to add an executable target to the build system,
# where an executable can be a binary file or a script written in a scripting
# language. In general we refer to any output file which is part of the software
# (i.e., excluding configuration files) and which can be executed
# (e.g., a binary file in the ELF format) or interpreted (e.g., a Python script)
# directly, as executable file. Natively, CMake supports only executables built
# from C/C++ source code files. This function extends CMake's capabilities
# by adding custom build commands for non-natively supported programming
# languages and further standardizes the build of executable targets.
# For example, by default, it is not necessary to specify installation rules
# separately as these are added by this function already (see below).
#
# @par Programming languages
# Besides adding usual executable targets build by the set <tt>C/CXX</tt>
# language compiler, this function inspects the list of source files given and
# detects whether this list contains sources which need to be build using a
# different compiler. In particular, it supports the following languages:
# @n
# <table border="0">
#   <tr>
#     @tp @b CXX @endtp
#     <td>The default behavior, adding an executable target build from C/C++
#         source code. The target is added via CMake's add_executable() command.</td>
#   </tr>
#   <tr>
#     @tp <b>PYTHON</b>|<b>JYTHON</b>|<b>PERL</b>|<b>BASH</b> @endtp
#     <td>Executables written in one of the named scripting languages are built by
#         configuring and/or copying the script files to the build tree and
#         installation tree, respectively. During the build step, certain strings
#         of the form \@VARIABLE\@ are substituted by the values set during the
#         configure step. How these CMake variables are set is specified by a
#         so-called script configuration, which itself is either a CMake script
#         file or a string of CMake code set as value of the @c SCRIPT_DEFINITIONS
#         property of the executable target.</td>
#   </tr>
#   <tr>
#     @tp @b MATLAB @endtp
#     <td>Standalone application built from MATLAB sources using the
#         MATLAB Compiler (mcc). This language option is used when the list
#         of source files contains one or more *.m files. A custom target is
#         added which depends on custom command(s) that build the executable.</td>
#         @n@n
#         Attention: The *.m file with the entry point/main function of the
#                    executable has to be given before any other *.m file.
#   </tr>
# </table>
#
# @par Helper functions
# If the programming language of the input source files is not specified
# explicitly by providing the @p LANGUAGE argument, the extensions of the
# source files and if necessary the first line of script files are inspected
# by the basis_get_source_language() function. Once the programming language is
# known, this function invokes the proper subcommand which adds the respective
# build target. In particular, it calls basis_add_executable_target() for C++
# sources (.cxx), basis_add_mcc_target() for MATLAB scripts (.m), and
# basis_add_script() for all other source files.
#
# @note DO NOT use the mentioned subcommands directly. Always use
#       basis_add_executable() to add an executable target to your project.
#       Only refer to the documentation of the subcommands to learn about the
#       available options of the particular subcommand and considered target
#       properties.
#
# @par Output directories
# The built executable file is output to the @c BINARY_RUNTIME_DIR or
# @c BINARY_LIBEXEC_DIR if the @p LIBEXEC option is given.
# If this function is used within the @c PROJECT_TESTING_DIR, however,
# the built executable is output to the @c TESTING_RUNTIME_DIR or
# @c TESTING_LIBEXEC_DIR instead.
#
# @par Installation
# An install command for the added executable target is added by this function
# as well. The executable will be installed as part of the specified @p COMPONENT
# in the directory @c INSTALL_RUNTIME_DIR or @c INSTALL_LIBEXEC_DIR if the option
# @p LIBEXEC is given. Executable targets are exported by default such that they
# can be imported by other CMake-aware projects by including the CMake
# configuration file of this package (&lt;Package&gt;Config.cmake file).
# No installation rules are added, however, if this function is used within the
# @c PROJECT_TESTING_DIR or if "none" (case-insensitive) is given as
# @p DESTINATION. Test executables are further only exported as part of the
# build tree, but not the installation as they are by default not installed.
#
# @param [in] TARGET_NAME Name of the target. If an existing source file is given
#                         as first argument, it is added to the list of source files
#                         and the build target name is derived from the name of this file.
# @param [in] ARGN        This argument list is parsed and the following
#                         arguments are extracted, all other arguments are passed
#                         on to add_executable() or the respective custom
#                         commands used to add an executable build target.
# @par
# <table border="0">
#   <tr>
#     @tp @b COMPONENT name @endtp
#     <td>Name of component as part of which this executable will be installed
#         if the specified @c DESTINATION is not "none".
#         (default: @c BASIS_RUNTIME_COMPONENT)</td>
#   </tr>
#   <tr>
#     @tp @b DESTINATION dir @endtp
#     <td>Installation directory relative to @c CMAKE_INSTALL_PREFIX.
#         If "none" (case-insensitive) is given as argument, no default
#         installation rules are added for this executable target.
#         (default: @c INSTALL_RUNTIME_DIR or @c INSTALL_LIBEXEC_DIR
#         if the @p LIBEXEC option is given)</td>
#   </tr>
#   <tr>
#     @tp @b LANGUAGE lang @endtp
#     <td>Programming language in which source files are written (case-insensitive).
#         If not specified, the programming language is derived from the file name
#         extensions of the source files and, if applicable, the shebang directive
#         on the first line of the script file. If the programming language could
#         not be detected automatically, check the file name extensions of the
#         source files and whether no unrecognized additional arguments were given
#         or specify the programming language using this option.
#         (default: auto-detected)</td>
#   </tr>
#   <tr>
#     @tp @b LIBEXEC @endtp
#     <td>Specifies that the built executable is an auxiliary executable which
#         is only called by other executables. (default: @c FALSE)</td>
#   </tr>
#   <tr>
#     @tp @b [NO]EXPORT @endtp
#     <td>Whether to export this target. (default: @c TRUE)</td>
#   </tr>
#   <tr>
#     @tp @b NO_BASIS_UTILITIES @endtp
#     <td>Specify that the BASIS utilities are not used by this executable and hence
#         no link dependency on the BASIS utilities shall be added.
#         (default: @c NOT @c BASIS_UTILITIES)</td>
#   </tr>
#   <tr>
#     @tp @b USE_BASIS_UTILITIES @endtp
#     <td>Specify that the BASIS utilities are used and required by this executable
#         and hence a link dependency on the BASIS utilities has to be added.
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
# @returns Adds an executable build target. In case of an executable which is
#          not build from C++ source files, the function basis_finalize_targets()
#          has to be invoked to finalize the addition of the custom build target.
#          This is done by the basis_project_end() macro.
#
# @sa basis_add_executable_target()
# @sa basis_add_script()
# @sa basis_add_mcc_target()
#
# @ingroup CMakeAPI
function (basis_add_executable TARGET_NAME)
  # --------------------------------------------------------------------------
  # parse arguments
  CMAKE_PARSE_ARGUMENTS (
    ARGN
      "EXECUTABLE;LIBEXEC;NO_BASIS_UTILITIES;USE_BASIS_UTILITIES;EXPORT;NOEXPORT;FINAL"
      "COMPONENT;DESTINATION;LANGUAGE"
      ""
    ${ARGN}
  )
  # derive target name from path if existing source path is given as first argument instead
  # and get list of library source files
  get_filename_component (S "${TARGET_NAME}" ABSOLUTE)
  if (IS_DIRECTORY "${S}" AND NOT ARGN_UNPARSED_ARGUMENTS)
    set (SOURCES "${S}")
    basis_get_source_target_name (TARGET_NAME "${TARGET_NAME}" NAME_WE)
  elseif (EXISTS "${S}" AND NOT IS_DIRECTORY "${S}" OR (NOT S MATCHES "\\.in$" AND EXISTS "${S}.in" AND NOT IS_DIRECTORY "${S}.in"))
    set (SOURCES "${S}")
    basis_get_source_target_name (TARGET_NAME "${TARGET_NAME}" NAME_WE)
  else ()
    set (SOURCES)
  endif ()
  if (ARGN_UNPARSED_ARGUMENTS)
    list (APPEND SOURCES ${ARGN_UNPARSED_ARGUMENTS})
  endif ()
  if (NOT SOURCES)
    message (FATAL_ERROR "basis_add_executable called with only one argument which does however not"
                         " appear to be a file name. Note that the filename extension must be"
                         " included if the target name should be derived from the base filename"
                         " of the source file.")
  endif ()
  # --------------------------------------------------------------------------
  # make target UID
  basis_check_target_name ("${TARGET_NAME}")
  basis_make_target_uid (TARGET_UID "${TARGET_NAME}")
  # --------------------------------------------------------------------------
  # process globbing expressions to get complete list of source files
  basis_add_glob_target (${TARGET_UID} SOURCES ${SOURCES})
  # --------------------------------------------------------------------------
  # determine programming language
  _basis_target_source_language ()
  # --------------------------------------------------------------------------
  # prepare arguments for subcommand
  foreach (ARG IN LISTS ARGN_UNPARSED_ARGUMENTS)
    list (REMOVE_ITEM ARGN "${ARG}")
  endforeach ()
  list (APPEND ARGN ${SOURCES})
  # --------------------------------------------------------------------------
  # C++
  if (ARGN_LANGUAGE MATCHES "CXX")
    basis_add_executable_target (.${TARGET_UID} ${ARGN})
  # --------------------------------------------------------------------------
  # MATLAB
  elseif (ARGN_LANGUAGE MATCHES "MATLAB")
    if (ARGN_LIBEXEC)
      list (REMOVE_ITEM ARGN LIBEXEC)
      basis_add_mcc_target (.${TARGET_UID} LIBEXEC ${ARGN})
    else ()
      list (REMOVE_ITEM ARGN EXECUTABLE)
      basis_add_mcc_target (.${TARGET_UID} EXECUTABLE ${ARGN})
    endif ()
  # --------------------------------------------------------------------------
  # others
  else ()
    if (ARGN_LIBEXEC)
      list (REMOVE_ITEM ARGN LIBEXEC)
      basis_add_script (.${TARGET_UID} LIBEXEC ${ARGN})
    else ()
      list (REMOVE_ITEM ARGN EXECUTABLE)
      basis_add_script (.${TARGET_UID} EXECUTABLE ${ARGN})
    endif ()
  endif ()
  # --------------------------------------------------------------------------
  # re-glob source files before each build (if necessary)
  if (TARGET __${TARGET_UID})
    if (TARGET _${TARGET_UID})
      add_dependencies (_${TARGET_UID} __${TARGET_UID})
    endif ()
    add_dependencies (${TARGET_UID} __${TARGET_UID})
  endif ()
endfunction ()

# ----------------------------------------------------------------------------
## @brief Add library target.
#
# This is the main function to add a library target to the build system, where
# a library can be a binary archive, shared library, a MEX-file or module(s)
# written in a scripting language. In general we refer to any output file which
# is part of the software (i.e., excluding configuration files), but cannot be
# executed (e.g., a binary file in the ELF format) or interpreted
# (e.g., a Python module) directly, as library file. Natively, CMake supports only
# libraries built from C/C++ source code files. This function extends CMake's
# capabilities by adding custom build commands for non-natively supported
# programming languages and further standardizes the build of library targets.
# For example, by default, it is not necessary to specify installation rules
# separately as these are added by this function already (see below).
#
# @par Programming languages
# Besides adding usual library targets built from C/C++ source code files,
# this function can also add custom build targets for libraries implemented
# in other programming languages. It therefore tries to detect the programming
# language of the given source code files and delegates the addition of the
# build target to the proper helper functions. It in particular supports the
# following languages:
# @n
# <table border="0">
#   <tr>
#     @tp @b CXX @endtp
#     <td>Source files written in C/C++ are by default built into either
#         @p STATIC, @p SHARED, or @p MODULE libraries. If the @p MEX option
#         is given, however, a MEX-file (a shared library) is build using
#         the MEX script instead of using the default C++ compiler directly.</td>
#   </tr>
#   <tr>
#     @tp <b>PYTHON</b>|<b>JYTHON</b>|<b>PERL</b>|<b>BASH</b> @endtp
#     <td>Modules written in one of the named scripting languages are built similar
#         to executable scripts except that the file name extension is preserved
#         and no executable file permission is set on Unix. These modules are
#         intended for import/inclusion in other modules or executables written
#         in the particular scripting language only.</td>
#   </tr>
#   <tr>
#     @tp @b MATLAB @endtp
#     <td>Libraries of M-files or shared libraries built using the MATLAB Compiler (mcc).
#         This language option is used when the list of source files contains one or
#         more *.m files. A custom target is added which depends on custom command(s)
#         that build the library. If the type of the library is @c SHARED, a shared
#         library is build using the MATLAB Compiler. Otherwise, the M-files are
#         configured and installed such that they can be used in MATLAB.</td>
#   </tr>
# </table>
#
# @par Helper functions
# If the programming language of the input source files is not specified
# explicitly by providing the @p LANGUAGE argument, the extensions of the
# source files are inspected using basis_get_source_language(). Once the
# programming language is known, this function invokes the proper subcommand.
# In particular, it calls basis_add_library_target() for C++ sources (.cxx)
# if the target is not a MEX-file target, basis_add_mex_file() for C++ sources
# if the @p MEX option is given, basis_add_mcc_target() for MATLAB scripts (.m),
# and basis_add_script_library() for all other source files.
#
# @note DO NOT use the mentioned subcommands directly. Always use
#       basis_add_library() to add a library target to your project. Only refer
#       to the documentation of the subcommands to learn about the available
#       options of the particular subcommand and the considered target properties.
#
# @par Output directories
# In case of modules written in a scripting language, the libraries are output to
# the <tt>BINARY_&lt;LANGUAGE&gt;_LIBRARY_DIR</tt> if defined. Otherwise,
# the built libraries are output to the @c BINARY_RUNTIME_DIR, @c BINARY_LIBRARY_DIR,
# and/or @c BINARY_ARCHIVE_DIR. If this command is used within the @c PROJECT_TESTING_DIR,
# however, the files are output to the corresponding directories in the testing tree,
# instead.
#
# @par Installation
# An installation rule for the added library target is added by this function
# if the destination is not "none" (case-insensitive). Runtime libraries are
# installed as part of the @p RUNTIME_COMPONENT to the @p RUNTIME_DESTINATION.
# Library components are installed as part of the @p LIBRARY_COMPONENT to the
# @p LIBRARY_DESTINATION. Library targets are further exported such that they
# can be imported by other CMake-aware projects by including the CMake
# configuration file of this package (&lt;Package&gt;Config.cmake file).
# If this function is used within the @c PROJECT_TESTING_DIR, however, no
# installation rules are added. Test library targets are further only exported
# as part of the build tree.
#
# @par Example
# @code
# basis_add_library (MyLib1 STATIC mylib.cxx)
# basis_add_library (MyLib2 STATIC mylib.cxx COMPONENT dev)
#
# basis_add_library (
#   MyLib3 SHARED mylib.cxx
#   RUNTIME_COMPONENT bin
#   LIBRARY_COMPONENT dev
# )
#
# basis_add_library (MyMex MEX mymex.cxx)
# basis_add_library (PythonModule MyModule.py.in)
# basis_add_library (ShellModule MODULE MyModule.sh.in)
# @endcode
#
# @param [in] TARGET_NAME Name of build target. If an existing file is given as
#                         argument, it is added to the list of source files and
#                         the target name is derived from the name of this file.
# @param [in] ARGN        This argument list is parsed and the following
#                         arguments are extracted. All unparsed arguments are
#                         treated as source files.
# @par
# <table border="0">
#   <tr>
#     @tp <b>STATIC</b>|<b>SHARED</b>|<b>MODULE</b>|<b>MEX</b> @endtp
#     <td>Type of the library. (default: @c SHARED for C++ libraries if
#         @c BUILD_SHARED_LIBS evaluates to true or @c STATIC otherwise,
#         and @c MODULE in all other cases)</td>
#   </tr>
#   <tr>
#     @tp @b COMPONENT name @endtp
#     <td>Name of component as part of which this library will be installed
#         if the @c RUNTIME_DESTINATION or @c LIBRARY_DESTINATION is not "none".
#         Used only if @p RUNTIME_COMPONENT or @p LIBRARY_COMPONENT not specified.
#         (default: see @p RUNTIME_COMPONENT and @p LIBRARY_COMPONENT)</td>
#   </tr>
#   <tr>
#     @tp @b DESTINATION dir @endtp
#     <td>Installation directory for runtime and library component relative
#         to @c CMAKE_INSTALL_PREFIX. Used only if @p RUNTIME_DESTINATION or
#         @p LIBRARY_DESTINATION not specified. If "none" (case-insensitive)
#         is given as argument, no default installation rules are added.
#         (default: see @p RUNTIME_DESTINATION and @p LIBRARY_DESTINATION)</td>
#   </tr>
#   <tr>
#     @tp @b LANGUAGE lang @endtp
#     <td>Programming language in which source files are written (case-insensitive).
#         If not specified, the programming language is derived from the file name
#         extensions of the source files and, if applicable, the shebang directive
#         on the first line of the script file. If the programming language could
#         not be detected automatically, check the file name extensions of the
#         source files and whether no unrecognized additional arguments were given
#         or specify the programming language using this option.
#         (default: auto-detected)</td>
#   </tr>
#   <tr>
#     @tp @b LIBRARY_COMPONENT name @endtp
#     <td>Name of component as part of which import/static library will be intalled
#         if @c LIBRARY_DESTINATION is not "none".
#         (default: @c COMPONENT if specified or @c BASIS_LIBRARY_COMPONENT otherwise)</td>
#   </tr>
#   <tr>
#     @tp @b LIBRARY_DESTINATION dir @endtp
#     <td>Installation directory of the library component relative to
#         @c CMAKE_INSTALL_PREFIX. If "none" (case-insensitive) is given as argument,
#         no installation rule for the library component is added.
#         (default: @c INSTALL_ARCHIVE_DIR)</td>
#   </tr>
#   <tr>
#     @tp @b RUNTIME_COMPONENT name @endtp
#     <td>Name of component as part of which runtime library will be installed
#         if @c RUNTIME_DESTINATION is not "none".
#         (default: @c COMPONENT if specified or @c BASIS_RUNTIME_COMPONENT otherwise)</td>
#   </tr>
#   <tr>
#     @tp @b RUNTIME_DESTINATION dir @endtp
#     <td>Installation directory of the runtime component relative to
#         @c CMAKE_INSTALL_PREFIX. If "none" (case-insensitive) is given as argument,
#         no installation rule for the runtime library is added.
#         (default: @c INSTALL_LIBRARY_DIR on Unix or @c INSTALL_RUNTIME_DIR Windows)</td>
#   </tr>
#   <tr>
#     @tp @b [NO]EXPORT @endtp
#     <td>Whether to export this target. (default: @c TRUE)</td>
#   </tr>
#   <tr>
#     @tp @b NO_BASIS_UTILITIES @endtp
#     <td>Specify that the BASIS utilities are not used by this executable and hence
#         no link dependency on the BASIS utilities shall be added.
#         (default: @c NOT @c BASIS_UTILITIES)</td>
#   </tr>
#   <tr>
#     @tp @b USE_BASIS_UTILITIES @endtp
#     <td>Specify that the BASIS utilities are used and required by this executable
#         and hence a link dependency on the BASIS utilities has to be added.
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
# @returns Adds a library build target. In case of a library not written in C++
#          or MEX-file targets, basis_finalize_targets() has to be invoked
#          to finalize the addition of the build target(s). This is done
#          by the basis_project_end() macro.
#
# @sa basis_add_library_target()
# @sa basis_add_script_library()
# @sa basis_add_mex_file()
# @sa basis_add_mcc_target()
#
# @ingroup CMakeAPI
function (basis_add_library TARGET_NAME)
  # --------------------------------------------------------------------------
  # parse arguments
  CMAKE_PARSE_ARGUMENTS (
    ARGN
      "STATIC;SHARED;MODULE;MEX;USE_BASIS_UTILITIES;NO_BASIS_UTILITIES;EXPORT;NOEXPORT;FINAL"
      "COMPONENT;RUNTIME_COMPONENT;LIBRARY_COMPONENT;DESTINATION;RUNTIME_DESTINATION;LIBRARY_DESTINATION;LANGUAGE"
      ""
    ${ARGN}
  )
  # derive target name from path if existing source path is given as first argument instead
  # and get list of library source files
  get_filename_component (S "${TARGET_NAME}" ABSOLUTE)
  if (IS_DIRECTORY "${S}" AND NOT ARGN_UNPARSED_ARGUMENTS)
    set (SOURCES "${S}")
    basis_get_source_target_name (TARGET_NAME "${TARGET_NAME}" NAME)
  elseif (EXISTS "${S}" AND NOT IS_DIRECTORY "${S}" OR (NOT S MATCHES "\\.in$" AND EXISTS "${S}.in" AND NOT IS_DIRECTORY "${S}.in"))
    set (SOURCES "${S}")
    if (ARGN_MEX)
      basis_get_source_target_name (TARGET_NAME "${TARGET_NAME}" NAME_WE)
    else ()
      _basis_target_source_language ()
      if ("$${ARGN_LANGUAGE}" STREQUAL "$CXX")
        basis_get_source_target_name (TARGET_NAME "${TARGET_NAME}" NAME_WE)
      else ()
        basis_get_source_target_name (TARGET_NAME "${TARGET_NAME}" NAME)
      endif ()
    endif ()
  else ()
    set (SOURCES)
  endif ()
  if (ARGN_UNPARSED_ARGUMENTS)
    list (APPEND SOURCES ${ARGN_UNPARSED_ARGUMENTS})
  endif ()
  # --------------------------------------------------------------------------
  # make target UID
  basis_check_target_name ("${TARGET_NAME}")
  basis_make_target_uid (TARGET_UID "${TARGET_NAME}")
  # --------------------------------------------------------------------------
  # process globbing expressions to get complete list of source files
  basis_add_glob_target (${TARGET_UID} SOURCES ${SOURCES})
  # --------------------------------------------------------------------------
  # determine programming language
  _basis_target_source_language ()
  # --------------------------------------------------------------------------
  # prepare arguments for subcommand
  foreach (ARG IN LISTS ARGN_UNPARSED_ARGUMENTS)
    list (REMOVE_ITEM ARGN "${ARG}")
  endforeach ()
  list (APPEND ARGN ${SOURCES})
  # --------------------------------------------------------------------------
  # C++
  if (ARGN_LANGUAGE MATCHES "CXX")
    # MEX-file
    if (ARGN_MEX)
      if (ARGN_STATIC)
        message (FATAL_ERROR "Target ${TARGET_UID}: Invalid library type! Only modules or shared libraries can be built by the MEX script.")
      endif ()
      list (REMOVE_ITEM ARGN MODULE)
      list (REMOVE_ITEM ARGN SHARED)
      list (REMOVE_ITEM ARGN MEX)
      basis_add_mex_file (.${TARGET_UID} ${ARGN})
    # library
    else ()
      basis_add_library_target (.${TARGET_UID} ${ARGN})
    endif ()
  # --------------------------------------------------------------------------
  # MATLAB
  elseif (ARGN_LANGUAGE MATCHES "MATLAB")
    if (ARGN_STATIC OR ARGN_MEX)
      message (FATAL_ERROR "Target ${TARGET_UID}: Invalid library type! Only shared libraries can be built by the MATLAB Compiler.")
    endif ()
    if (ARGN_SHARED)
      list (REMOVE_ITEM ARGN SHARED)
      basis_add_mcc_target (.${TARGET_UID} SHARED ${ARGN})
    else ()
      list (REMOVE_ITEM ARGN MODULE) # optional
      basis_add_script_library (.${TARGET_UID} ${ARGN})
    endif ()
  # --------------------------------------------------------------------------
  # other
  else ()
    if (ARGN_STATIC OR ARGN_SHARED OR ARGN_MEX)
      message (FATAL_ERROR "Target ${TARGET_UID}: Invalid library type! Only modules can be built from scripts.")
    endif ()
    list (REMOVE_ITEM ARGN MODULE)
    basis_add_script_library (.${TARGET_UID} ${ARGN})
  endif ()
  # --------------------------------------------------------------------------
  # re-glob source files before each build (if necessary)
  if (TARGET __${TARGET_UID})
    if (TARGET _${TARGET_UID})
      add_dependencies (_${TARGET_UID} __${TARGET_UID})
    endif ()
    add_dependencies (${TARGET_UID} __${TARGET_UID})
  endif ()
endfunction ()

# ----------------------------------------------------------------------------
## @brief Add single arbitrary or executable script.
#
# @note This function should not be used directly for executable scripts or
#       module libraries. Use basis_add_executable() or basis_add_library()
#       in such (most) cases instead.
#
# This function can be used to add a single arbitrary script file (i.e., any
# text file which is input to a program), such as a CTest script for example,
# to the build if neither basis_add_executable() nor basis_add_library() are
# appropriate choices. In all other cases, either basis_add_executable() or
# basis_add_library() should be used. Note that the script file is by default
# not considered to be an executable. Instead it is assumed that the program
# which interprets/processes the script must be executed explicitly with this
# script as argument. Only scripts built with the @p EXECUTABLE or @p LIBEXEC
# type option are treated as executable files, where in case of Unix a shebang
# directive implicitly states the program used by the shell to interpret the
# script and on Windows a Windows Command which imitates the behavior of Unix
# shells is generated by BASIS. Do not use these type options, however, but
# only use the default @p MODULE option. The basis_add_executable() function
# should be used instead to add an executable script. The basis_add_script()
# function shall only be used for none-executable arbitrary script files which
# cannot be built by basis_add_executable() or basis_add_library().
#
# If the script name ends in <tt>.in</tt>, the <tt>.in</tt> suffix is removed
# from the output name. Further, in case of executable scripts, the file name
# extension is removed from the output file name. Instead, a shebang directive
# is added on Unix to the built script. In order to enable the convenient
# execution of Python and Perl scripts also on Windows without requiring the
# user to setup a proper associate between the filename extension with the
# corresponding interpreter executable, a few lines of Batch code are added at
# the top and bottom of executable Python and Perl scripts. This Batch code
# invokes the configured interpreter with the script file and the given script
# arguments as command-line arguments. Note that both the original script source
# code and the Batch code are stored within the single file. The file name
# extension of such modified scripts is by default set to <tt>.cmd</tt>, the
# common extension for Windows NT Command Scripts. Scripts in other languages
# are not modified and the extension of the original scripts script file is
# preserved on Windows in this case. In case of non-executable scripts, the
# file name extension is kept in any case.
#
# Certain CMake variables within the source file are replaced during the
# built of the script. See the
# <a href="https://cmake-basis.github.io/scripttargets/>
# Build System Standard</a> for details.
# Note, however, that source files are only configured if the file name
# ends in the <tt>.in</tt> suffix.
#
# A custom CMake build target with the following properties is added by this
# function to the build system. These properties are used by basis_build_script()
# to generate a build script written in CMake code which is executed by a custom
# CMake command. Before the invokation of basis_build_script(), the target
# properties can be modified using basis_set_target_properties().
#
# @note Custom BASIS build targets are finalized by BASIS using basis_project_end(),
#       i.e., the end of the root CMake configuration file of the (sub-)project.
#
# @par Properties on script targets
# <table border=0>
#   <tr>
#     @tp @b BASIS_TYPE @endtp
#     <td>Read-only property with value "SCRIPT_FILE" for arbitrary scripts,
#         "SCRIPT_EXECUTABLE" for executable scripts, and "SCRIPT_LIBEXEC" for
#          auxiliary executable scripts.
#          (default: see @p MODULE, @p EXECUTABLE, @p LIBEXEC options)</td>
#   </tr>
#   <tr>
#     @tp @b BASIS_UTILITIES @endtp
#     <td>Whether the BASIS utilities are used by this script. For the supported
#         scripting languages for which BASIS utilities are implemented, BASIS
#         will in most cases automatically detect whether these utilities are
#         used by a script or not. Otherwise, set this property manually or use
#         either the @p USE_BASIS_UTILITIES or the @p NO_BASIS_UTILITIES option
#         when adding the script target. (default: auto-detected or @c UNKNOWN)</td>
#   </tr>
#   <tr>
#     @tp @b BINARY_DIRECTORY @endtp
#     <td>Build tree directory of this target. (default: @c CMAKE_CURRENT_BINARY_DIR)</td>
#   </tr>
#   <tr>
#     @tp @b COMPILE @endtp
#     <td>Whether to compile the script if the programming language allows such
#         pre-compilation as in case of Python, for example. If @c TRUE, only the
#         compiled file is installed. (default: @c BASIS_COMPILE_SCRIPTS)</td>
#   </tr>
#   <tr>
#     @tp @b SCRIPT_DEFINITIONS @endtp
#     <td>CMake code which is evaluated after the inclusion of the default script
#         configuration files. This code can be used to set the replacement text of the
#         CMake variables ("@VAR@" patterns) used in the source file.
#         See <a href="https://cmake-basis.github.io/standard/scripttargets.html#script-configuration">
#         Build System Standard</a> for details. (default: "")</td>
#   </tr>
#   <tr>
#     @tp @b SCRIPT_DEFINITIONS_FILE @endtp
#     <td>CMake script file with compile definitions, also referred to as script
#         configuration file. The named files are included after the default BASIS
#         script configuration and before the @c SCRIPT_DEFINITIONS code is being
#         evaluated. (default: @c BINARY_CONFIG_DIR/ScriptConfig.cmake)</td>
#   </tr>
#   <tr>
#     @tp @b COMPONENT @endtp
#     <td>Name of component as part of which this script is installed if
#         @c INSTALL_DIRECTORY is not set to "none".
#         (default: see @p COMPONENT argument)</td>
#   </tr>
#   <tr>
#     @tp @b EXPORT @endtp
#     <td>Whether to export this build target in which case an import library
#         target is added to the custom exports file with the path to the
#         built/installed script set as @c IMPORT_LOCATION. (default: @c TRUE)</td>
#   </tr>
#   <tr>
#     @tp @b INSTALL_DIRECTORY @endtp
#     <td>Installation directory of script file configured for use in installation tree
#         relative to @c CMAKE_INSTALL_PREFIX. Set to "none" (case-insensitive) to skip the
#         addition of an installation rule. (default: see @p DESTINATION argument)</td>
#   </tr>
#   <tr>
#     @tp @b LANGUAGE @endtp
#     <td>Read-only property of programming language of script file in uppercase letters.
#         (default: see @p LANGUAGE argument)</td>
#   </tr>
#   <tr>
#     @tp @b LINK_DEPENDS @endtp
#     <td>Paths or target names of script modules and libraries used by this script.
#         In case of an (auxiliary) executable script, the directories of these modules
#         are added to the search path for modules of the given programming language
#         if such search paths are supported by the language and BASIS knows how to set
#         these (as in case of Python/Jython, Perl, and MATLAB, in particular).
#         Moreover, for each listed build target a dependency is added between this
#         script target and the named build targets. Use basis_target_link_libraries()
#         to add additional link dependencies.
#         (default: BASIS utilities module if used or empty list otherwise)</td>
#   </tr>
#   <tr>
#     @tp @b OUTPUT_DIRECTORY @endtp
#     <td>Output directory for built script file configured for use in build tree.
#         (default: @c BINARY_LIBRARY_DIR for arbitrary scripts, @c BINARY_RUNTIME_DIR
#         for executable scripts, and @c BINARY_LIBEXEC_DIR for auxiliary executables)</td>
#   </tr>
#   <tr>
#     @tp @b OUTPUT_NAME @endtp
#     <td>Name of built script file including file name extension (if any).
#         (default: basename of script file for arbitrary scripts, without extension
#         for executable scripts on Unix, and <tt>.cmd</tt> extension on Windows
#         in case of executable Python/Jython or Perl script)</td>
#   </tr>
#   <tr>
#     @tp @b SOURCE_DIRECTORY @endtp
#     <td>Source directory of this target. (default: @c CMAKE_CURRENT_SOURCE_DIR)</td>
#   </tr>
#   <tr>
#     @tp @b SOURCES @endtp
#     <td>Read-only property which lists the source file of this script target.
#         Note that the first element in this list actually names a directory
#         in the build, the one where the build script for this target is located
#         instead of a source file and thus should be ignored. The second entry
#         corresponds to the source file of this script target.</td>
#   </tr>
# </table>
#
# @attention Properties documented as read-only must not be modified.
#
# @note If this function is used within the @c PROJECT_TESTING_DIR, the built
#       executable is output to the @c BINARY_TESTING_DIR directory tree instead.
#       Moreover, no installation rules are added. Test executables are further
#       not exported, regardless of the @c EXPORT property.
#
# @param [in] TARGET_NAME Name of build target. If an existing file is given as
#                         argument, it is added to the list of source files and
#                         the target name is derived from the name of this file.
# @param [in] ARGN        The remaining arguments are parsed and the following arguments
#                         recognized. All unparsed arguments are treated as source files,
#                         where in particular exactly one source file is required if the
#                         @p TARGET_NAME argument does not name an existing source file.
# @par
# <table border=0>
#   <tr>
#     @tp <b>MODULE</b>|<b>EXECUTABLE</b>|<b>LIBEXEC</b> @endtp
#     <td>Type of script to built, i.e., either arbitrary module script which
#         cannot be executed directly, an executable script with proper shebang
#         directive and execute permissions on Unix or Windows Command on Windows,
#         or an auxiliary executable. The type of the script mainly changes the
#         default values of the target properties such as the output and installation
#         directories. To add an (auxiliary) executable script, use
#         basis_add_executable(), however, instead of this function.
#         The @c EXECUTABLE and @c LIBEXEC options are only intended for
#         internal use by BASIS. (default: MODULE)</td>
#   </tr>
#   <tr>
#     @tp @b COMPONENT name @endtp
#     <td>Name of installation component as part of which this script is being
#         installed if the @c INSTALL_DIRECTORY property is not "none".
#         (default: @c BASIS_LIBRARY_COMPONENT for arbitrary scripts or
#         @c BASIS_RUNTIME_COMPONENT for executable scripts)</td>
#   </tr>
#   <tr>
#     @tp @b DESTINATION dir @endtp
#     <td>Installation directory for script file relative to @c CMAKE_INSTALL_PREFIX.
#         If an absolute path is given as argument, it is made relative to the
#         configured installation prefix.
#         (default: @c INSTALL_LIBRARY_DIR for arbitrary scripts,
#         @c INSTALL_RUNTIME_DIR for executable scripts, and @c INSTALL_LIBEXEC_DIR
#         for auxiliary executable scripts)</td>
#   </tr>
#   <tr>
#     @tp @b LANGUAGE lang @endtp
#     <td>Programming language in which script file is written (case-insensitive).
#         If not specified, the programming language is derived from the file name
#         extension of the source file and the shebang directive on the first line
#         of the script if any. If the programming language could not be detected
#         automatically, the @c LANGUAGE property is set to @c UNKNOWN. Note that
#         for arbitrary script targets, the script file will still be built correctly
#         even if the scripting language was not recognized. The automatic detection
#         whether the BASIS utilities are used and required will fail, however.
#         In this case, specify the programming language using this option.
#         (default: auto-detected or @c UNKNOWN)</td>
#   </tr>
#   <tr>
#     @tp @b [NO]EXPORT @endtp
#     <td>Whether to export this target. (default: @c TRUE)</td>
#   </tr>
#   <tr>
#     @tp @b NO_BASIS_UTILITIES @endtp
#     <td>Specify that the BASIS utilities are not used by this script. If the
#         programming language of the script is known and BASIS utilities are
#         available for this language, BASIS will in most cases automatically
#         detect whether these utilities are used by a script or not. Use this
#         option to skip this check because the script does not make use of the
#         BASIS utilities.</td>
#   </tr>
#   <tr>
#     @tp @b USE_BASIS_UTILITIES @endtp
#     <td>Specify that the BASIS utilities are used and thus required by this script.
#         If the programming language of the script is known and BASIS utilities are
#         available for this language, BASIS will in most cases automatically
#         detect whether these utilities are used by a script or not. Use this option
#         to skip this check because it is already known that the script makes use of
#         the BASIS utilities. Note that an error is raised if this option is given,
#         but no BASIS utilities are available for the programming language of this
#         script or if the programming language is unknown, respectively, not detected
#         correctly. In this case, consider the use of the @p LANGUAGE argument.</td>
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
# @returns Adds a custom CMake target with the documented properties. The actual custom
#          command to build the script is added by basis_build_script().
#
# @ingroup CMakeAPI
function (basis_add_script TARGET_NAME)
  # parse arguments
  CMAKE_PARSE_ARGUMENTS (
    ARGN
      "MODULE;EXECUTABLE;LIBEXEC;NO_BASIS_UTILITIES;USE_BASIS_UTILITIES;EXPORT;NOEXPORT;FINAL"
      "COMPONENT;DESTINATION;LANGUAGE"
      ""
    ${ARGN}
  )
  if (NOT ARGN_MODULE AND NOT ARGN_EXECUTABLE AND NOT ARGN_LIBEXEC)
    set (ARGN_MODULE TRUE)
  endif ()
  if (ARGN_MODULE)
    set (TYPE MODULE)
  else ()
    set (TYPE EXECUTABLE)
  endif ()
  string (TOLOWER "${TYPE}" type)
  # derive target name from file name if existing source file given as first argument
  get_filename_component (S "${TARGET_NAME}" ABSOLUTE)
  if (EXISTS "${S}" AND NOT IS_DIRECTORY "${S}" OR (NOT S MATCHES "\\.in$" AND EXISTS "${S}.in" AND NOT IS_DIRECTORY "${S}.in"))
    set (SOURCES "${S}")
    if (ARGN_MODULE)
      basis_get_source_target_name (TARGET_NAME "${TARGET_NAME}" NAME)
    else ()
      basis_get_source_target_name (TARGET_NAME "${TARGET_NAME}" NAME_WE)
    endif ()
    set (SET_OUTPUT_NAME_TO_TARGET_NAMEE FALSE)
  else ()
    set (SOURCES)
    set (SET_OUTPUT_NAME_TO_TARGET_NAME TRUE)
  endif ()
  # check target name
  basis_check_target_name ("${TARGET_NAME}")
  basis_make_target_uid (TARGET_UID "${TARGET_NAME}")
  message (STATUS "Adding ${type} script ${TARGET_UID}...")
  if (ARGN_MODULE AND TYPE MATCHES "EXECUTABLE")
    message (FATAL_ERROR "Target ${TARGET_UID}: MODULE and EXECUTABLE or LIBEXEC options are mutually exclusive!")
  endif ()
  # check/set parsed arguments
  basis_set_flag (ARGN EXPORT ${BASIS_EXPORT_DEFAULT})
  if (ARGN_USE_BASIS_UTILITIES AND ARGN_NO_BASIS_UTILITIES)
    message (FATAL_ERROR "Options USE_BASIS_UTILITIES and NO_BASIS_UTILITIES are mutually exclusive!")
  endif ()
  list (LENGTH ARGN_UNPARSED_ARGUMENTS N)
  if (SOURCES)
    math (EXPR N "${N} + 1")
  endif ()
  if (N GREATER 1)
    if (NOT SOURCES)
      list (REMOVE_AT ARGN_UNPARSED_ARGUMENTS 0)
    endif ()
    message (FATAL_ERROR "Target ${TARGET_UID}: Too many or unrecognized arguments: ${ARGN_UNPARSED_ARGUMENTS}!\n"
                         " Only one script can be built by each script target.")
  elseif (NOT SOURCES)
    set (SOURCES "${ARGN_UNPARSED_ARGUMENTS}")
    get_filename_component (SOURCES "${SOURCES}" ABSOLUTE)
  endif ()
  if (NOT EXISTS "${SOURCES}" AND NOT SOURCES MATCHES "\\.in$" AND EXISTS "${SOURCES}.in")
    set (SOURCES "${SOURCES}.in")
  endif ()
  if (NOT EXISTS "${SOURCES}")
    string (REGEX REPLACE "\\.in$" "" SOURCES "${SOURCES}")
    message (FATAL_ERROR "Target ${TARGET_UID}: Source file ${SOURCES}[.in] does not exist!")
  endif ()
  # add custom target
  add_custom_target (${TARGET_UID} ALL SOURCES ${SOURCES})
  # dump CMake variables for configuration of script
  set (BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${TARGET_UID}")
  basis_dump_variables ("${BUILD_DIR}.dir/cache.cmake.tmp")
  # auto-detect programming language (may be as well UNKNOWN)
  if (ARGN_LANGUAGE)
    string (TOUPPER "${ARGN_LANGUAGE}" ARGN_LANGUAGE)
  else ()
    basis_get_source_language (ARGN_LANGUAGE ${SOURCES})
  endif ()
  # IS_TEST flag
  basis_sanitize_for_regex (RE "${PROJECT_TESTING_DIR}")
  if (CMAKE_CURRENT_SOURCE_DIR MATCHES "^${RE}")
    set (IS_TEST TRUE)
  else ()
    set (IS_TEST FALSE)
  endif ()
  # default directory infix used below
  if (ARGN_MODULE)
    set (TYPE_INFIX "LIBRARY")
  elseif (ARGN_LIBEXEC)
    set (TYPE_INFIX "LIBEXEC")
  else ()
    set (TYPE_INFIX "RUNTIME")
  endif ()
  # output name
  string (REGEX REPLACE "\\.in$" "" SOURCE_NAME "${SOURCES}")
  if (SET_OUTPUT_NAME_TO_TARGET_NAME)
    basis_get_target_name (OUTPUT_NAME ${TARGET_UID})
  else ()
    get_filename_component (OUTPUT_NAME "${SOURCE_NAME}" NAME_WE)
  endif ()
  if (ARGN_MODULE)
    get_filename_component (SUFFIX "${SOURCE_NAME}" EXT)
  else ()
    if (WIN32)
      if (ARGN_LANGUAGE MATCHES "[JP]YTHON|PERL")
        set (SUFFIX ".cmd")
      else ()
        get_filename_component (SUFFIX "${SOURCE_NAME}" EXT)
      endif ()
    else ()
      set (SUFFIX)
    endif ()
  endif ()
  # output directory
  if (IS_TEST)
    set (OUTPUT_DIRECTORY "${TESTING_${TYPE_INFIX}_DIR}")
  else ()
    set (OUTPUT_DIRECTORY "${BINARY_${TYPE_INFIX}_DIR}")
  endif ()
  # installation component
  if (NOT ARGN_COMPONENT)
    if (ARGN_MODULE)
      set (ARGN_COMPONENT "${BASIS_LIBRARY_COMPONENT}")
    else ()
      set (ARGN_COMPONENT "${BASIS_RUNTIME_COMPONENT}")
    endif ()
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
  elseif (IS_TEST)
    set (ARGN_DESTINATION) # do not install
  else ()
    set (ARGN_DESTINATION "${INSTALL_${TYPE_INFIX}_DIR}")
  endif ()
  # script configuration ("compile definitions")
  if (EXISTS "${BINARY_CONFIG_DIR}/ScriptConfig.cmake")
    set (CONFIG_FILE "${BINARY_CONFIG_DIR}/ScriptConfig.cmake")
  else ()
    set (CONFIG_FILE)
  endif ()
  # auto-detect use of BASIS utilities
  set (LINK_DEPENDS)
  if (ARGN_LANGUAGE MATCHES "[JP]YTHON")
    set (UTILITIES_LANGUAGE "PYTHON")
  else ()
    set (UTILITIES_LANGUAGE "${ARGN_LANGUAGE}")
  endif ()
  if (ARGN_USE_BASIS_UTILITIES)
    if (NOT BASIS_UTILITIES_ENABLED MATCHES "${UTILITIES_LANGUAGE}")
      message (FATAL_ERROR "Target ${TARGET_UID} requires the BASIS utilities for ${UTILITIES_LANGUAGE}"
                           " but BASIS was either build without the build of these utilities enabled"
                           " or no utilities for this programming language are implemented. Remove the"
                           " USE_BASIS_UTILITIES option if no BASIS utilities are used by the script"
                           " ${SOURCES} or specify the correct programming language if it was not"
                           " detected correctly.")
    endif ()
    set (USES_BASIS_UTILITIES TRUE)
  elseif (NOT ARGN_NO_BASIS_UTILITIES AND NOT UTILITIES_LANGUAGE MATCHES "UNKNOWN")
    basis_utilities_check (USES_BASIS_UTILITIES ${SOURCES} ${UTILITIES_LANGUAGE})
  else ()
    set (USES_BASIS_UTILITIES FALSE)
  endif ()
  if (USES_BASIS_UTILITIES)
    basis_set_project_property (PROPERTY PROJECT_USES_${UTILITIES_LANGUAGE}_UTILITIES TRUE)
    if (BASIS_DEBUG)
      message ("** Target ${TARGET_UID} uses the BASIS utilities for ${UTILITIES_LANGUAGE}.")
    endif ()
  endif ()
  # set properties of custom build target
  set_target_properties (
    ${TARGET_UID}
    PROPERTIES
      LANGUAGE                 ${ARGN_LANGUAGE}
      BASIS_TYPE               SCRIPT_${TYPE}
      BASIS_UTILITIES          ${USES_BASIS_UTILITIES}
      BUILD_DIRECTORY          "${BUILD_DIR}"
      SOURCE_DIRECTORY         "${CMAKE_CURRENT_SOURCE_DIR}"
      BINARY_DIRECTORY         "${CMAKE_CURRENT_BINARY_DIR}"
      OUTPUT_DIRECTORY         "${OUTPUT_DIRECTORY}"
      INSTALL_DIRECTORY        "${ARGN_DESTINATION}"
      COMPONENT                "${ARGN_COMPONENT}"
      OUTPUT_NAME              "${OUTPUT_NAME}"
      PREFIX                   ""
      SUFFIX                   "${SUFFIX}"
      SCRIPT_DEFINITIONS       ""
      SCRIPT_DEFINITIONS_FILE  "${CONFIG_FILE}"
      LINK_DEPENDS             "${LINK_DEPENDS}"
      EXPORT                   ${EXPORT}
      COMPILE                  ${BASIS_COMPILE_SCRIPTS}
      TEST                     ${IS_TEST}
      LIBEXEC                  ${ARGN_LIBEXEC}
  )
  # finalize target
  if (ARGN_FINAL)
    basis_finalize_targets (${TARGET_UID})
  endif ()
  # add target to list of targets
  basis_set_project_property (APPEND PROPERTY TARGETS "${TARGET_UID}")
  message (STATUS "Adding ${type} script ${TARGET_UID}... - done")
endfunction ()

# ----------------------------------------------------------------------------
## @brief Finalize custom targets by adding the missing build commands.
#
# This function is called by basis_project_end() in order to finalize the
# addition of the custom build targets such as, for example, build targets
# for the build of executable scripts, Python packages, MATLAB Compiler
# executables and shared libraries, and MEX-files. It can, however, also
# be called explicitly either at the end of each CMakeLists.txt file that
# adds new build targets or with the name of the target(s) to finalize in
# some of the CMakeLists.txt files. This is to ensure that the custom target
# which executes the actual build command is added in the same directory as
# the original custom target with the respective target properties.
#
# @param[in] ARGN List of targets to finalize. If none specified, all custom
#                 targets that were added before and are not finalized already
#                 will be finalized using the current binary directory.
#
# @returns Generates the CMake build scripts and adds custom build commands
#          and corresponding targets for the execution of these scripts.
#
# @sa basis_build_script()
# @sa basis_build_script_library()
# @sa basis_build_mcc_target()
# @sa basis_build_mex_file()
function (basis_finalize_targets)
  if (ARGN)
    set (TARGETS)
    foreach (TARGET_NAME ${ARGN})
      basis_get_target_uid (TARGET_UID ${TARGET_NAME})
      list (APPEND TARGETS ${TARGET_UID})
    endforeach ()
  else ()
    basis_get_project_property (TARGETS PROPERTY TARGETS)
    if (NOT TARGETS)
      return()
    endif ()
    # targets of BASIS utilities are finalized separately
    # because some properties still have to be set by
    # basis_configure_utilities, see UtilitiesTools module
    basis_make_target_uid (TARGET_UID_basis_sh basis_sh)
    basis_make_target_uid (TARGET_UID_basis_py basis_py)
    basis_make_target_uid (TARGET_UID_Basis_pm Basis_pm)
    list (REMOVE_ITEM TARGETS
      ${TARGET_UID_basis_sh}
      ${TARGET_UID_basis_py}
      ${TARGET_UID_Basis_pm}
    )
  endif ()
  basis_get_project_property (FINALIZED_TARGETS PROPERTY FINALIZED_TARGETS)
  if (FINALIZED_TARGETS)
    list (REMOVE_ITEM TARGETS ${FINALIZED_TARGETS})
  endif ()
  if (BASIS_DEBUG)
    message (
      "** basis_finalize_targets:"
      "\n     TARGETS:   [${TARGETS}]"
      "\n     FINALIZED: [${FINALIZED_TARGETS}]"
    )
  endif ()
  foreach (TARGET_UID ${TARGETS})
    get_target_property (BASIS_TYPE ${TARGET_UID} BASIS_TYPE)
    if (BASIS_TYPE MATCHES "^EXECUTABLE$|^(SHARED|MODULE)_LIBRARY$")
      if (BASIS_INSTALL_RPATH AND NOT CMAKE_SKIP_RPATH)
        # Only if BASIS is allowed to take care of the INSTALL_RPATH property
        # and the use of this property was not disabled by the project
        basis_set_target_install_rpath (${TARGET_UID})
      endif ()
    elseif (BASIS_TYPE MATCHES "SCRIPT_LIBRARY")
      basis_build_script_library (${TARGET_UID})
    elseif (BASIS_TYPE MATCHES "SCRIPT")
      basis_build_script (${TARGET_UID})
    elseif (BASIS_TYPE MATCHES "MEX")
      basis_build_mex_file (${TARGET_UID})
    elseif (BASIS_TYPE MATCHES "MCC")
      basis_build_mcc_target (${TARGET_UID})
    endif ()
    list (APPEND FINALIZED_TARGETS ${TARGET_UID})
  endforeach ()
  basis_set_project_property (PROPERTY FINALIZED_TARGETS ${FINALIZED_TARGETS})
endfunction ()

# ============================================================================
# internal helpers
# ============================================================================

# ----------------------------------------------------------------------------
## @brief Add executable target.
#
# This BASIS function overwrites CMake's
# <a href="http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:add_executable">
# add_executable()</a> command in order to store information of imported targets
# which is in particular used to generate the source code of the ExecutableTargetInfo
# modules which are part of the BASIS utilities.
#
# @note Use basis_add_executable() instead where possible!
#
# @param [in] TARGET_UID Name of the target.
# @param [in] ARGN       Further arguments of CMake's add_executable().
#
# @sa http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:add_executable
function (add_executable TARGET_UID)
  _add_executable (${TARGET_UID} ${ARGN})
  list (FIND ARGN IMPORTED IDX)
  if (IDX EQUAL 0)
    if (COMMAND basis_add_imported_target)
      basis_add_imported_target ("${TARGET_UID}" EXECUTABLE)
    endif ()
  else ()
    basis_set_project_property (APPEND PROPERTY TARGETS "${TARGET_UID}")
  endif ()
endfunction ()

# ----------------------------------------------------------------------------
## @brief Add library target.
#
# This BASIS function overwrites CMake's
# <a href="http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:add_library">
# add_library()</a> command in order to store information of imported targets.
#
# @note Use basis_add_library() instead where possible!
#
# @param [in] TARGET_UID Name of the target.
# @param [in] ARGN       Further arguments of CMake's add_library().
#
# @sa http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:add_library
function (add_library TARGET_UID)
  _add_library (${TARGET_UID} ${ARGN})
  list (FIND ARGN IMPORTED IDX)
  if (IDX EQUAL 1)
    if (COMMAND basis_add_imported_target)
      basis_add_imported_target ("${TARGET_UID}" "${ARGV1}")
    endif ()
  else ()
    basis_set_project_property (APPEND PROPERTY TARGETS "${TARGET_UID}")
  endif ()
endfunction ()

# ----------------------------------------------------------------------------
## @brief Add executable built from C++ source code.
#
# @note This function should not be used directly. Instead, it is called
#       by basis_add_executable() if the (detected) programming language
#       of the given source code files is @c CXX (i.e., C/C++).
#
# This function adds an executable target for the build of an executable from
# C++ source code files. Refer to the documentation of basis_add_executable()
# for a description of general options for adding an executable target.
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
# @param [in] TARGET_NAME Name of build target.
# @param [in] ARGN        This argument list is parsed and the following
#                         arguments are extracted, all other arguments are
#                         considered to be source code files and simply passed
#                         on to CMake's add_executable() command.
# @par
# <table border=0>
#   <tr>
#     @tp @b COMPONENT name @endtp
#     <td>Name of component as part of which this executable will be installed
#         if the specified @c DESTINATION is not "none".
#         (default: @c BASIS_RUNTIME_COMPONENT)</td>
#   </tr>
#   <tr>
#     @tp @b DESTINATION dir @endtp
#     <td>Installation directory relative to @c CMAKE_INSTALL_PREFIX.
#         If "none" (case-insensitive) is given as argument, no default
#         installation rules are added for this executable target.
#         (default: @c INSTALL_RUNTIME_DIR or @c INSTALL_LIBEXEC_DIR
#         if @p LIBEXEC is given)</td>
#   </tr>
#   <tr>
#     @tp @b LIBEXEC @endtp
#     <td>Specifies that the built executable is an auxiliary executable which
#         is only called by other executables. (default: @c FALSE)</td>
#   </tr>
#   <tr>
#     @tp @b [NO]EXPORT @endtp
#     <td>Whether to export this target. (default: @c TRUE)</td>
#   </tr>
#   <tr>
#     @tp @b NO_BASIS_UTILITIES @endtp
#     <td>Specify that the BASIS utilities are not used by this executable and hence
#         no link dependency on the BASIS utilities library shall be added.
#         (default: @c NOT @c BASIS_UTILITIES)</td>
#   </tr>
#   <tr>
#     @tp @b USE_BASIS_UTILITIES @endtp
#     <td>Specify that the BASIS utilities are used and required by this executable
#         and hence a link dependency on the BASIS utilities library has to be added.
#         (default: @c BASIS_UTILITIES)</td>
#   </tr>
# </table>
#
# @returns Adds executable target using CMake's add_executable() command.
#
# @sa basis_add_executable()
function (basis_add_executable_target TARGET_NAME)
  # check target name
  basis_check_target_name (${TARGET_NAME})
  basis_make_target_uid (TARGET_UID "${TARGET_NAME}")
  message (STATUS "Adding executable ${TARGET_UID}...")
  # parse arguments
  CMAKE_PARSE_ARGUMENTS (
    ARGN
      "USE_BASIS_UTILITIES;NO_BASIS_UTILITIES;EXPORT;NOEXPORT;LIBEXEC"
      "COMPONENT;DESTINATION"
      ""
    ${ARGN}
  )
  set (SOURCES ${ARGN_UNPARSED_ARGUMENTS})
  basis_set_flag (ARGN EXPORT  ${BASIS_EXPORT_DEFAULT})
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
  # IS_TEST flag
  basis_sanitize_for_regex (RE "${PROJECT_TESTING_DIR}")
  if (CMAKE_CURRENT_SOURCE_DIR MATCHES "^${RE}")
    set (IS_TEST TRUE)
  else ()
    set (IS_TEST FALSE)
  endif ()
  # installation component
  if (NOT ARGN_COMPONENT)
    set (ARGN_COMPONENT "${BASIS_RUNTIME_COMPONENT}")
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
  elseif (ARGN_LIBEXEC)
    set (ARGN_DESTINATION "${INSTALL_LIBEXEC_DIR}")
  else ()
    set (ARGN_DESTINATION "${INSTALL_RUNTIME_DIR}")
  endif ()
  # configure (.in) source files
  basis_configure_sources (SOURCES ${SOURCES})
  # add executable target
  add_executable (${TARGET_UID} ${SOURCES})
  basis_make_target_uid (HEADERS_TARGET headers)
  if (TARGET "${HEADERS_TARGET}")
    add_dependencies (${TARGET_UID} ${HEADERS_TARGET})
  endif ()
  basis_get_target_name (OUTPUT_NAME ${TARGET_UID})
  set_target_properties (${TARGET_UID} PROPERTIES BASIS_TYPE "EXECUTABLE" LANGUAGE "CXX" OUTPUT_NAME "${OUTPUT_NAME}")
  if (ARGN_LIBEXEC)
    set_target_properties (${TARGET_UID} PROPERTIES LIBEXEC 1 COMPILE_DEFINITIONS LIBEXEC SCRIPT_DEFINITIONS LIBEXEC)
  else ()
    set_target_properties (${TARGET_UID} PROPERTIES LIBEXEC 0)
  endif ()
  set_target_properties (${TARGET_UID} PROPERTIES TEST ${IS_TEST})
  # output directory
  if (IS_TEST)
    if (ARGN_LIBEXEC)
      set_target_properties (${TARGET_UID} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${TESTING_LIBEXEC_DIR}")
    else ()
      set_target_properties (${TARGET_UID} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${TESTING_RUNTIME_DIR}")
    endif ()
  elseif (ARGN_LIBEXEC)
    set_target_properties (${TARGET_UID} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${BINARY_LIBEXEC_DIR}")
  else ()
    set_target_properties (${TARGET_UID} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${BINARY_RUNTIME_DIR}")
  endif ()
  # installation directory
  set_target_properties (${TARGET_UID} PROPERTIES RUNTIME_INSTALL_DIRECTORY "${ARGN_DESTINATION}")
  # link to BASIS utilities
  if (USES_BASIS_UTILITIES)
    basis_target_link_libraries (.${TARGET_UID} basis)
  else ()
    set_target_properties (${TARGET_UID} PROPERTIES BASIS_UTILITIES FALSE)
  endif ()
  # export
  set (EXPORT_OPT)
  if (EXPORT)
    basis_add_export_target (EXPORT_OPT ${TARGET_UID} "${IS_TEST}" ${ARGN_DESTINATION})
  endif ()
  # installation
  if (ARGN_DESTINATION)
    if (IS_TEST)
      # TODO install (selected?) tests
    else ()
      install (
        TARGETS ${TARGET_UID} ${EXPORT_OPT}
        DESTINATION "${ARGN_DESTINATION}"
        COMPONENT   "${ARGN_COMPONENT}"
      )
    endif ()
  endif ()
  # done
  message (STATUS "Adding executable ${TARGET_UID}... - done")
endfunction ()

# ----------------------------------------------------------------------------
## @brief Add library built from C++ source code.
#
# @note This function should not be used directly. Instead, it is called
#       by basis_add_library() if the (detected) programming language
#       of the given source code files is @c CXX (i.e., C/C++) and the
#       option @c MEX is not given.
#
# This function adds a library target which builds a library from C++ source
# code files. Refer to the documentation of basis_add_library() for a
# description of the general options for adding a library target.
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
# @param [in] TARGET_NAME Name of build target.
# @param [in] ARGN        This argument list is parsed and the following
#                         arguments are extracted. All other arguments are
#                         considered to be source code files and simply
#                         passed on to CMake's add_library() command.
# @par
# <table border=0>
#   <tr>
#     @tp <b>STATIC</b>|<b>SHARED</b>|<b>MODULE</b> @endtp
#     <td>Type of the library. (default: @c SHARED if @c BUILD_SHARED_LIBS
#         evaluates to true or @c STATIC otherwise)</td>
#   </tr>
#   <tr>
#     @tp @b COMPONENT name @endtp
#     <td>Name of component as part of which this library will be installed
#         if either the @c RUNTIME_INSTALL_DIRECTORY or
#         @c LIBRARY_INSTALL_DIRECTORY property is not "none". Used only if
#         either @p RUNTIME_COMPONENT or @p LIBRARY_COMPONENT not specified.
#         (default: see @p RUNTIME_COMPONENT and @p LIBRARY_COMPONENT)</td>
#   </tr>
#   <tr>
#     @tp @b DESTINATION dir @endtp
#     <td>Installation directory for runtime and library component relative
#         to @c CMAKE_INSTALL_PREFIX. Used only if either @p RUNTIME_DESTINATION
#         or @p LIBRARY_DESTINATION not specified. If "none" (case-insensitive)
#         is given as argument, no default installation rules are added.
#         (default: see @p RUNTIME_DESTINATION and @p LIBRARY_DESTINATION)</td>
#   </tr>
#   <tr>
#     @tp @b LIBRARY_COMPONENT name @endtp
#     <td>Name of component as part of which import/static library will be intalled
#         if @c LIBRARY_INSTALL_DIRECTORY property is not "none".
#         (default: @c COMPONENT if specified or @c BASIS_LIBRARY_COMPONENT otherwise)</td>
#   </tr>
#   <tr>
#     @tp @b LIBRARY_DESTINATION dir @endtp
#     <td>Installation directory of the library component relative to
#         @c CMAKE_INSTALL_PREFIX. If "none" (case-insensitive) is given as argument,
#         no installation rule for the library component is added.
#         (default: @c INSTALL_ARCHIVE_DIR)</td>
#   </tr>
#   <tr>
#     @tp @b RUNTIME_COMPONENT name @endtp
#     <td>Name of component as part of which runtime library will be installed
#         if @c RUNTIME_INSTALL_DIRECTORY property is not "none".
#         (default: @c COMPONENT if specified or @c BASIS_RUNTIME_COMPONENT otherwise)</td>
#   </tr>
#   <tr>
#     @tp @b RUNTIME_DESTINATION dir @endtp
#     <td>Installation directory of the runtime component relative to
#         @c CMAKE_INSTALL_PREFIX. If "none" (case-insensitive) is given as argument,
#         no installation rule for the runtime library is added.
#         (default: @c INSTALL_LIBRARY_DIR on Unix or @c INSTALL_RUNTIME_DIR Windows)</td>
#   </tr>
#   <tr>
#     @tp @b [NO]EXPORT @endtp
#     <td>Whether to export this target. (default: @c TRUE)</td>
#   </tr>
#   <tr>
#     @tp @b NO_BASIS_UTILITIES @endtp
#     <td>Specify that the BASIS utilities are not used by this executable and hence
#         no link dependency on the BASIS utilities library shall be added.
#         (default: @c NOT BASIS_UTILITIES)</td>
#   </tr>
#   <tr>
#     @tp @b USE_BASIS_UTILITIES @endtp
#     <td>Specify that the BASIS utilities are used and required by this executable
#         and hence a link dependency on the BASIS utilities library shall be added.
#         (default: @c BASIS_UTILITIES)</td>
#   </tr>
# </table>
#
# @returns Adds library target using CMake's add_library() command.
#
# @sa basis_add_library()
function (basis_add_library_target TARGET_NAME)
  # On UNIX-based systems setting the VERSION property only creates
  # annoying files with the version string as suffix.
  # Moreover, MEX-files may NEVER have a suffix after the MEX extension!
  # Otherwise, the MATLAB Compiler when using the symbolic link
  # without this suffix will create code that fails on runtime
  # with an .auth file missing error.
  #
  # Thus, do NOT set VERSION and SOVERSION properties on library targets!

  # check target name
  basis_check_target_name (${TARGET_NAME})
  basis_make_target_uid (TARGET_UID "${TARGET_NAME}")
  # parse arguments
  CMAKE_PARSE_ARGUMENTS (
    ARGN
      "STATIC;SHARED;MODULE;USE_BASIS_UTILITIES;NO_BASIS_UTILITIES;EXPORT;NOEXPORT"
      "COMPONENT;RUNTIME_COMPONENT;LIBRARY_COMPONENT;DESTINATION;RUNTIME_DESTINATION;LIBRARY_DESTINATION"
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
  # IS_TEST flag
  basis_sanitize_for_regex (RE "${PROJECT_TESTING_DIR}")
  if (CMAKE_CURRENT_SOURCE_DIR MATCHES "^${RE}")
    set (IS_TEST TRUE)
  else ()
    set (IS_TEST FALSE)
  endif ()
  # library type
  if (NOT ARGN_SHARED AND NOT ARGN_STATIC AND NOT ARGN_MODULE)
    if (BUILD_SHARED_LIBS)
      set (ARGN_SHARED TRUE)
    else ()
      set (ARGN_STATIC TRUE)
    endif ()
  endif ()
  set (TYPE)
  if (ARGN_STATIC)
    if (TYPE)
      message (FATAL_ERROR "More than one library type specified for target ${TARGET_UID}!")
    endif ()
    set (TYPE "STATIC")
  endif ()
  if (ARGN_SHARED)
    if (TYPE)
      message (FATAL_ERROR "More than one library type specified for target ${TARGET_UID}!")
    endif ()
    set (TYPE "SHARED")
  endif ()
  if (ARGN_MODULE)
    if (TYPE)
      message (FATAL_ERROR "More than one library type specified for target ${TARGET_UID}!")
    endif ()
    set (TYPE "MODULE")
  endif ()
  string (TOLOWER "${TYPE}" type)
  # status message
  message (STATUS "Adding ${type} library ${TARGET_UID}...")
  # installation component
  if (ARGN_COMPONENT)
    if (NOT ARGN_RUNTIME_COMPONENT)
      set (ARGN_RUNTIME_COMPONENT "${ARGN_COMPONENT}")
    endif ()
    if (NOT ARGN_LIBRARY_COMPONENT)
      set (ARGN_LIBRARY_COMPONENT "${ARGN_COMPONENT}")
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
    if (NOT ARGN_STATIC AND NOT ARGN_RUNTIME_DESTINATION)
      set (ARGN_RUNTIME_DESTINATION "${ARGN_DESTINATION}")
    endif ()
    if (NOT ARGN_LIBRARY_DESTINATION)
      set (ARGN_LIBRARY_DESTINATION "${ARGN_DESTINATION}")
    endif ()
  endif ()
  if (NOT ARGN_RUNTIME_DESTINATION)
    set (ARGN_RUNTIME_DESTINATION "${INSTALL_RUNTIME_DIR}")
  endif ()
  if (NOT ARGN_LIBRARY_DESTINATION)
    set (ARGN_LIBRARY_DESTINATION "${INSTALL_LIBRARY_DIR}")
  endif ()
  if (ARGN_STATIC OR ARGN_RUNTIME_DESTINATION MATCHES "^[nN][oO][nN][eE]$")
    set (ARGN_RUNTIME_DESTINATION)
  endif ()
  if (ARGN_LIBRARY_DESTINATION MATCHES "^[nN][oO][nN][eE]$")
    set (ARGN_LIBRARY_DESTINATION)
  endif ()
  # configure (.in) source files
  basis_configure_sources (SOURCES ${SOURCES})
  # add library target
  add_library (${TARGET_UID} ${TYPE} ${SOURCES})
  basis_make_target_uid (HEADERS_TARGET headers)
  if (TARGET ${HEADERS_TARGET})
    add_dependencies (${TARGET_UID} ${HEADERS_TARGET})
  endif ()
  basis_get_target_name (OUTPUT_NAME ${TARGET_UID})
  set_target_properties (${TARGET_UID} PROPERTIES BASIS_TYPE "${TYPE}_LIBRARY" LANGUAGE "CXX" OUTPUT_NAME "${OUTPUT_NAME}")
  # output directory
  if (IS_TEST)
    set_target_properties (
      ${TARGET_UID}
      PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${TESTING_RUNTIME_DIR}"
        LIBRARY_OUTPUT_DIRECTORY "${TESTING_LIBRARY_DIR}"
        ARCHIVE_OUTPUT_DIRECTORY "${TESTING_ARCHIVE_DIR}"
    )
  else ()
    set_target_properties (
      ${TARGET_UID}
      PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${BINARY_RUNTIME_DIR}"
        LIBRARY_OUTPUT_DIRECTORY "${BINARY_LIBRARY_DIR}"
        ARCHIVE_OUTPUT_DIRECTORY "${BINARY_ARCHIVE_DIR}"
    )
  endif ()
  # installation directory
  # these properties are used by basis_get_target_location() in particular
  set_target_properties (
    ${TARGET_UID}
    PROPERTIES
      RUNTIME_INSTALL_DIRECTORY "${ARGN_RUNTIME_DESTINATION}"
      LIBRARY_INSTALL_DIRECTORY "${ARGN_LIBRARY_DESTINATION}"
      ARCHIVE_INSTALL_DIRECTORY "${ARGN_LIBRARY_DESTINATION}"
  )
  # link to BASIS utilities
  if (USES_BASIS_UTILITIES)
    basis_target_link_libraries (.${TARGET_UID} basis)
  else ()
    set_target_properties (${TARGET_UID} PROPERTIES BASIS_UTILITIES FALSE)
  endif ()
  # export
  set (EXPORT_OPT)
  if (EXPORT)
    basis_add_export_target (EXPORT_OPT ${TARGET_UID} "${IS_TEST}" ${ARGN_RUNTIME_DESTINATION} ${ARGN_LIBRARY_DESTINATION})
  endif ()
  # installation
  set (DESTINATION_OPTS)
  if (IS_TEST)
    # TODO At the moment, no tests are installed. Once there is a way to
    #      install selected tests, the shared libraries they depend on
    #      need to be installed as well.
  else ()
    if (ARGN_RUNTIME_DESTINATION)
      list (APPEND DESTINATION_OPTS
        RUNTIME
          DESTINATION "${ARGN_RUNTIME_DESTINATION}"
          COMPONENT   "${ARGN_RUNTIME_COMPONENT}"
      )
    endif ()
    if (ARGN_LIBRARY_DESTINATION)
      list (APPEND DESTINATION_OPTS
        LIBRARY
          DESTINATION "${ARGN_LIBRARY_DESTINATION}"
          COMPONENT   "${ARGN_LIBRARY_COMPONENT}"
        ARCHIVE
          DESTINATION "${ARGN_LIBRARY_DESTINATION}"
          COMPONENT   "${ARGN_LIBRARY_COMPONENT}"
      )
    endif ()
  endif ()
  if (DESTINATION_OPTS)
    install (TARGETS ${TARGET_UID} ${EXPORT_OPT} ${DESTINATION_OPTS})
  endif ()
  # done
  message (STATUS "Adding ${type} library ${TARGET_UID}... - done")
endfunction ()

# ----------------------------------------------------------------------------
## @brief Add script library target.
#
# @note This function should not be used directly. Instead, it is called
#       by basis_add_library() if the (detected) programming language
#       of the given source code files is neither @c CXX (i.e., C/C++) nor
#       @c MATLAB.
#
# This function adds a build target for libraries which are a collection of
# one or more modules written in a scripting language. The relative paths
# of the modules relative to the library's @p SOURCE_DIRECTORY property are
# preserved. This is important for the most widely used scripting languages
# such as Python, Perl, or MATLAB, where the file path relative to the
# package root directory defines the package namespace.
#
# A custom CMake build target with the following properties is added by this
# function to the build system. These properties are used by
# basis_build_script_library() to generate a build script written in CMake
# code which is executed by a custom CMake command. Before the invokation of
# basis_build_script_library(), the target properties can be modified using
# basis_set_target_properties().
#
# @note Custom BASIS build targets are finalized by BASIS using basis_project_end(),
#       i.e., the end of the root CMake configuration file of the (sub-)project.
#
# @par Properties on script library targets
# <table border=0>
#   <tr>
#     @tp @b BASIS_TYPE @endtp
#     <td>Read-only property with value "SCRIPT_LIBRARY" for script library targets.</td>
#   </tr>
#   <tr>
#     @tp @b BASIS_UTILITIES @endtp
#     <td>Whether the BASIS utilities are used by any module of this library.
#         For the supported scripting languages for which BASIS utilities are
#         implemented, BASIS will in most cases automatically detect whether
#         these utilities are used by a module or not. Otherwise, set this
#         property manually or use either the @p USE_BASIS_UTILITIES or the
#         @p NO_BASIS_UTILITIES option when adding the library target.
#         (default: auto-detected or @c UNKNOWN)</td>
#   </tr>
#   <tr>
#     @tp @b BINARY_DIRECTORY @endtp
#     <td>Build tree directory of this target. (default: @c CMAKE_CURRENT_BINARY_DIR)</td>
#   </tr>
#   <tr>
#     @tp @b COMPILE @endtp
#     <td>Whether to compile the library, respectively, it's modules if the
#         programming language allows such pre-compilation as in case of Python,
#         for example. If @c TRUE, only the compiled files are installed.
#         (default: @c BASIS_COMPILE_SCRIPTS)</td>
#   </tr>
#   <tr>
#     @tp @b SCRIPT_DEFINITIONS @endtp
#     <td>CMake code which is evaluated after the inclusion of the default script
#         configuration files. This code can be used to set the replacement text of the
#         CMake variables ("@VAR@" patterns) used in the source files.
#         See <a href="https://cmake-basis.github.io/standard/scripttargets.html#script-configuration">
#         Build System Standard</a> for details. (default: "")</td>
#   </tr>
#   <tr>
#     @tp @b SCRIPT_DEFINITIONS_FILE @endtp
#     <td>CMake script file with compile definitions, also referred to as script
#         configuration file. The named files are included after the default BASIS
#         script configuration and before the @c SCRIPT_DEFINITIONS code is being
#         evaluated. (default: @c BINARY_CONFIG_DIR/ScriptConfig.cmake)</td>
#   </tr>
#   <tr>
#     @tp @b EXPORT @endtp
#     <td>Whether to export this build target in which case an import library
#         target is added to the custom exports file with the path to the
#         built/installed modules set as @c IMPORT_LOCATION. (default: @c TRUE)</td>
#   </tr>
#   <tr>
#     @tp @b LANGUAGE @endtp
#     <td>Read-only property of programming language of modules in uppercase letters.
#         (default: see @p LANGUAGE argument)</td>
#   </tr>
#   <tr>
#     @tp @b LIBRARY_COMPONENT @endtp
#     <td>Name of component as part of which this library is installed if
#         @c LIBRARY_INSTALL_DIRECTORY is not set to "none".
#         (default: see @p COMPONENT argument)</td>
#   </tr>
#   <tr>
#     @tp @b LIBRARY_INSTALL_DIRECTORY @endtp
#     <td>Installation directory of library configured for use in installation tree
#         relative to @c CMAKE_INSTALL_PREFIX. Set to "none" (case-insensitive) to skip the
#         addition of an installation rule.
#         (default: <tt>INSTALL_&lt;LANGUAGE&gt;_LIBRARY_DIR</tt> if defined or
#         @c INSTALL_LIBRARY_DIR otherwise)</td>
#   </tr>
#   <tr>
#     @tp @b LIBRARY_OUTPUT_DIRECTORY @endtp
#     <td>Output directory of library configured for use within the build tree.
#         (default: <tt>BINARY_&lt;LANGUAGE&gt;_LIBRARY_DIR</tt> if defined or
#         @c BINARY_LIBRARY_DIR otherwise)</td>
#   </tr>
#   <tr>
#     @tp @b LINK_DEPENDS @endtp
#     <td>Paths or target names of script modules and libraries used by this script.
#         For each listed build target, a dependency is added between this
#         library target and the named build targets. Use basis_target_link_libraries()
#         to add additional link dependencies. Further note that if this library is
#         a link dependency of an executable script added by basis_add_executable()
#         (i.e., basis_add_script() actually), the link dependencies of this library
#         are inherited by the executable script.
#         (default: BASIS utilities module if used or empty list otherwise)</td>
#   </tr>
#   <tr>
#     @tp @b PREFIX @endtp
#     <td>Common module prefix. The given directory path is appended to both
#         @c LIBRAR_OUTPUT_DIRECTORY and @c LIBRARY_INSTALL_DIRECTORY and can,
#         for example, be used to install modules of a Python package as part of
#         another Python package, where @c LIBRARY_OUTPUT_DIRECTORY or
#         @c LIBRARY_INSTALL_DIRECTORY, respectively, is the directory of the
#         main package which is added to the @c PYTHONPATH. Possibly missing
#         __init__.py files in case of Python are generated by the _initpy target
#         which is automatically added by BASIS in that case and further added to
#         the dependencies of this library target.
#         (default: @c PROJECT_NAMESPACE_PYTHON if @c LANGUAGE is @c PYTHON with
#         periods (.) replaced by slashes (/), @c PROJECT_NAMESPACE_PERL if
#         @c LANGUAGE is @c PERL with <tt>::</tt> replaced by slashes (/),
#         and "" otherwise)</td>
#   </tr>
#   <tr>
#     @tp @b SOURCE_DIRECTORY @endtp
#     <td>Source directory of this target. This directory is in particular
#         used to convert the paths of the given source files to relative paths.
#         The built modules within the build and installation tree will have the
#         same relative path (relative to the @c LIBRARY_OUTPUT_DIRECTORY or
#         @c LIBRARY_INSTALL_DIRECTORY, respectively).
#         (default: @c CMAKE_CURRENT_SOURCE_DIR)</td>
#   </tr>
#   <tr>
#     @tp @b SOURCES @endtp
#     <td>Read-only property which lists the source files of this library.
#         Note that the first element in this list actually names a directory
#         in the build, the one where the build script for this target is located
#         instead of a source file and thus should be ignored.</td>
#   </tr>
# </table>
#
# @attention Properties documented as read-only must not be modified.
#
# @param [in] TARGET_NAME Name of build target.
# @param [in] ARGN        The remaining arguments are parsed and the following
#                         arguments extracted. All unparsed arguments are treated
#                         as the module files of the script library.
# @par
# <table border=0>
#   <tr>
#     @tp @b COMPONENT name @endtp
#     <td>Name of installation component as part of which this library is being
#         installed if the @c LIBRARY_INSTALL_DIRECTORY property is not "none".
#         (default: @c BASIS_LIBRARY_COMPONENT)</td>
#   </tr>
#   <tr>
#     @tp @b DESTINATION dir @endtp
#     <td>Installation directory for library relative to @c CMAKE_INSTALL_PREFIX.
#         If an absolute path is given as argument, it is made relative to the
#         configured installation prefix. (default: @c INSTALL_LIBRARY_DIR)</td>
#   </tr>
#   <tr>
#     @tp @b LANGUAGE lang @endtp
#     <td>Programming language in which modules are written (case-insensitive).
#         If not specified, the programming language is derived from the file name
#         extensions of the source files and the shebang directive on the first line
#         of each module if any. If the programming language could not be detected
#         automatically, the @c LANGUAGE property is set to @c UNKNOWN. Note that
#         for script library targets, the library may still be built correctly
#         even if the scripting language was not recognized. The automatic detection
#         whether the BASIS utilities are used and required will fail, however.
#         In this case, specify the programming language using this option.
#         (default: auto-detected or @c UNKNOWN)</td>
#   </tr>
#   <tr>
#     @tp @b [NO]EXPORT @endtp
#     <td>Whether to export this target. (default: @c TRUE)</td>
#   </tr>
#   <tr>
#     @tp @b NO_BASIS_UTILITIES @endtp
#     <td>Specify that the BASIS utilities are not used by this library. If the
#         programming language of the modules is known and BASIS utilities are
#         available for this language, BASIS will in most cases automatically
#         detect whether these utilities are used by any module of this library.
#         Use this option to skip this check in the case that no module makes
#         use of the BASIS utilities.</td>
#   </tr>
#   <tr>
#     @tp @b USE_BASIS_UTILITIES @endtp
#     <td>Specify that the BASIS utilities are used and thus required by this library.
#         If the programming language of the modules is known and BASIS utilities are
#         available for this language, BASIS will in most cases automatically
#         detect whether these utilities are used by any module of this library.
#         Use this option to skip this check when it is already known that no module
#         makes use of the BASIS utilities. Note that an error is raised if this option
#         is given, but no BASIS utilities are available for the programming language
#         of this script or if the programming language is unknown, respectively, not
#         detected correctly. In this case, consider the use of the @p LANGUAGE argument.</td>
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
# @returns Adds a custom CMake target with the documented properties. The actual custom
#          command to build the library is added by basis_build_script_library().
#
# @sa basis_add_library()
function (basis_add_script_library TARGET_NAME)
  # check target name
  basis_check_target_name ("${TARGET_NAME}")
  basis_make_target_uid (TARGET_UID "${TARGET_NAME}")
  message (STATUS "Adding script library ${TARGET_UID}...")
  # dump CMake variables for configuration of script
  set (BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${TARGET_UID}")
  basis_dump_variables ("${BUILD_DIR}.dir/cache.cmake.tmp")
  # parse arguments
  CMAKE_PARSE_ARGUMENTS (
    ARGN
      "NO_BASIS_UTILITIES;USE_BASIS_UTILITIES;EXPORT;NOEXPORT;FINAL"
      "COMPONENT;DESTINATION;LANGUAGE"
      ""
    ${ARGN}
  )
  basis_set_flag (ARGN EXPORT ${BASIS_EXPORT_DEFAULT})
  set (SOURCES "${ARGN_UNPARSED_ARGUMENTS}")
  # IS_TEST flag
  basis_sanitize_for_regex (RE "${PROJECT_TESTING_DIR}")
  if (CMAKE_CURRENT_SOURCE_DIR MATCHES "^${RE}")
    set (IS_TEST TRUE)
  else ()
    set (IS_TEST FALSE)
  endif ()
  # check source files
  set (_SOURCES)
  foreach (S IN LISTS SOURCES)
    get_filename_component (S "${S}" ABSOLUTE)
    if (NOT EXISTS "${S}" AND NOT S MATCHES "\\.in$" AND EXISTS "${S}.in" AND NOT IS_DIRECTORY "${S}.in")
      set (S "${S}.in")
    elseif (IS_DIRECTORY "${S}")
      message (FATAL_ERROR "Target ${TARGET_UID}: Directory ${S} given where file name expected!")
    endif ()
    if (NOT EXISTS "${S}")
      string (REGEX REPLACE "\\.in$" "" S "${S}")
      message (FATAL_ERROR "Target ${TARGET_UID}: Source file ${S}[.in] does not exist!")
    endif ()
    list (APPEND _SOURCES "${S}")
  endforeach ()
  if (NOT _SOURCES)
    message (FATAL_ERROR "Target ${TARGET_UID}: No source files specified!")
  endif ()
  set (SOURCES "${_SOURCES}")
  unset (_SOURCES)
  # auto-detect programming language (may be as well UNKNOWN)
  string (TOUPPER "${ARGN_LANGUAGE}" ARGN_LANGUAGE)
  if (NOT ARGN_LANGUAGE)
    basis_get_source_language (ARGN_LANGUAGE ${SOURCES})
    if (ARGN_LANGUAGE MATCHES "AMBIGUOUS|UNKNOWN")
      message (FATAL_ERROR "Target ${TARGET_UID}: Failed to determine programming"
                           " language of modules! Make sure that all modules are"
                           " written in the same language and that the used programming"
                           " language is supported by BASIS, i.e., either Python (Jython),"
                           " Perl, Bash, or MATLAB. Otherwise, try to specify the language"
                           " explicitly using the LANGUAGE option.")
    endif ()
  endif ()
  # output directory
  if (IS_TEST)
    if (DEFINED TESTING_${ARGN_LANGUAGE}_LIBRARY_DIR)
      set (OUTPUT_DIRECTORY "${TESTING_${ARGN_LANGUAGE}_LIBRARY_DIR}")
    else ()
      set (OUTPUT_DIRECTORY "${TESTING_LIBRARY_DIR}")
    endif ()
  else ()
    if (DEFINED BINARY_${ARGN_LANGUAGE}_LIBRARY_DIR)
      set (OUTPUT_DIRECTORY "${BINARY_${ARGN_LANGUAGE}_LIBRARY_DIR}")
    else ()
      set (OUTPUT_DIRECTORY "${BINARY_LIBRARY_DIR}")
    endif ()
  endif ()
  # installation component
  if (NOT ARGN_COMPONENT)
    set (ARGN_COMPONENT "${BASIS_LIBRARY_COMPONENT}")
  endif ()
  if (NOT ARGN_COMPONENT)
    set (ARGN_COMPONENT "Unspecified")
  endif ()
  # installation directory
  if (IS_TEST)
    if (ARGN_DESTINATION)
      message (WARNING "Target ${TARGET_UID} is a library used for testing only."
                       " Installation to the specified directory will be skipped.")
      set (ARGN_DESTINATION)
    endif ()
  else ()
    if (ARGN_DESTINATION)
      if (IS_ABSOLUTE "${ARGN_DESTINATION}")
        file (RELATIVE_PATH ARGN_DESTINATION "${CMAKE_INSTALL_PREFIX}" "${ARGN_DESTINATION}")
      endif ()
    else ()
      if (DEFINED INSTALL_${ARGN_LANGUAGE}_LIBRARY_DIR)
        set (ARGN_DESTINATION "${INSTALL_${ARGN_LANGUAGE}_LIBRARY_DIR}")
      else ()
        set (ARGN_DESTINATION "${INSTALL_LIBRARY_DIR}")
      endif ()
    endif ()
  endif ()
  # common module prefix
  if (ARGN_LANGUAGE MATCHES "^([JP]YTHON|PERL|MATLAB|BASH)$")
    basis_library_prefix (PREFIX ${ARGN_LANGUAGE})
  else ()
    set (PREFIX)
  endif ()
  # script configuration ("compile definitions")
  if (EXISTS "${BINARY_CONFIG_DIR}/ScriptConfig.cmake")
    set (CONFIG_FILE "${BINARY_CONFIG_DIR}/ScriptConfig.cmake")
  else ()
    set (CONFIG_FILE)
  endif ()
  # auto-detect use of BASIS utilities
  if (ARGN_LANGUAGE MATCHES "[JP]YTHON")
    set (UTILITIES_LANGUAGE "PYTHON")
  else ()
    set (UTILITIES_LANGUAGE "${ARGN_LANGUAGE}")
  endif ()
  if (ARGN_USE_BASIS_UTILITIES)
    if (NOT BASIS_UTILITIES_ENABLED MATCHES "${UTILITIES_LANGUAGE}")
      message (FATAL_ERROR "Target ${TARGET_UID} requires the BASIS utilities for ${UTILITIES_LANGUAGE}"
                           " but BASIS was either build without the build of these utilities enabled"
                           " or no utilities for this programming language are implemented. Remove the"
                           " USE_BASIS_UTILITIES option if no BASIS utilities are used by the modules"
                           " of the library or specify the correct programming language if it was not"
                           " detected correctly.")
    endif ()
    set (USES_BASIS_UTILITIES TRUE)
  elseif (BASIS_UTILITIES AND NOT ARGN_NO_BASIS_UTILITIES AND NOT UTILITIES_LANGUAGE MATCHES "UNKNOWN")
    set (USES_BASIS_UTILITIES FALSE)
    foreach (M IN LISTS SOURCES)
      basis_utilities_check (USES_BASIS_UTILITIES "${M}" ${UTILITIES_LANGUAGE})
      if (USES_BASIS_UTILITIES)
        break ()
      endif ()
    endforeach ()
  else ()
    set (USES_BASIS_UTILITIES FALSE)
  endif ()
  # add custom target
  add_custom_target (${TARGET_UID} ALL SOURCES ${SOURCES})
  set_target_properties (
    ${TARGET_UID}
    PROPERTIES
      LANGUAGE                  "${ARGN_LANGUAGE}"
      BASIS_TYPE                "SCRIPT_LIBRARY"
      BASIS_UTILITIES           "${USES_BASIS_UTILITIES}"
      BUILD_DIRECTORY           "${BUILD_DIR}"
      SOURCE_DIRECTORY          "${CMAKE_CURRENT_SOURCE_DIR}"
      BINARY_DIRECTORY          "${CMAKE_CURRENT_BINARY_DIR}"
      LIBRARY_OUTPUT_DIRECTORY  "${OUTPUT_DIRECTORY}"
      LIBRARY_INSTALL_DIRECTORY "${ARGN_DESTINATION}"
      LIBRARY_COMPONENT         "${BASIS_LIBRARY_COMPONENT}"
      PREFIX                    "${PREFIX}"
      SCRIPT_DEFINITIONS        ""
      SCRIPT_DEFINITIONS_FILE   "${CONFIG_FILE}"
      LINK_DEPENDS              ""
      EXPORT                    "${EXPORT}"
      COMPILE                   "${BASIS_COMPILE_SCRIPTS}"
      TEST                      "${IS_TEST}"
  )
  # link to BASIS utilities
  if (USES_BASIS_UTILITIES)
    basis_target_link_libraries (.${TARGET_UID} basis)
    if (BASIS_DEBUG)
      message ("** Target ${TARGET_UID} uses the BASIS utilities for ${UTILITIES_LANGUAGE}.")
    endif ()
  endif ()
  # finalize target
  if (ARGN_FINAL)
    basis_finalize_targets (${TARGET_UID})
  endif ()
  # add target to list of targets
  basis_set_project_property (APPEND PROPERTY TARGETS "${TARGET_UID}")
  message (STATUS "Adding script library ${TARGET_UID}... - done")
endfunction ()

# ============================================================================
# custom build commands
# ============================================================================

# ----------------------------------------------------------------------------
## @brief Set INSTALL_RPATH property of executable or shared library target.
#
# This function sets the @c INSTALL_RPATH property of a specified executable or
# shared library target using the @c LINK_DEPENDS obtained using the
# basis_get_target_link_libraries() function. It determines the installation
# location of each dependency using the basis_get_target_location() function.
#
# @returns Sets the @c INSTALL_RPATH property of the specified target.
#
# @sa basis_get_target_link_libraries()
function (basis_set_target_install_rpath TARGET_NAME)
  basis_get_target_uid (TARGET_UID "${TARGET_NAME}")
  if (BASIS_VERBOSE)
    message (STATUS "Setting INSTALL_RPATH property of ${TARGET_UID}...")
  endif ()
  if (NOT TARGET "${TARGET_UID}")
    message (FATAL_ERROR "Unknown target: ${TARGET_UID}")
  endif ()
  if (BASIS_DEBUG)
    message ("** basis_set_target_install_rpath():")
    message ("**    TARGET_NAME:   ${TARGET_UID}")
  endif ()
  if (CMAKE_HOST_APPLE)
    set (ORIGIN "@loader_path")
  else ()
    set (ORIGIN "\$ORIGIN")
  endif ()
  # always prefer libraries located within the same directory
  set (INSTALL_RPATH "${ORIGIN}/.")
  # common default RPATH (rarely used)
  if (CMAKE_INSTALL_RPATH)
    set (INSTALL_RPATH "${INSTALL_RPATH};${CMAKE_INSTALL_RPATH}")
  endif ()
  # get location of target used to make paths relative to this $ORIGIN
  basis_get_target_location (TARGET_LOCATION ${TARGET_UID} POST_INSTALL_PATH)
  # directories of external projects belonging to same bundle which
  # were added using [basis_]link_directories() command
  basis_get_project_property (BUNDLE_LINK_DIRS PROPERTY BUNDLE_LINK_DIRS)
  foreach (LINK_DIR ${BUNDLE_LINK_DIRS})
    list (FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${LINK_DIR}" IDX)
    if (IDX EQUAL -1)
      if (BASIS_DEBUG AND BASIS_VERBOSE)
        message ("**    BUNDLE_LINK_DIR: ${LINK_DIR}")
      endif ()
      basis_get_relative_path (RPATH "${TARGET_LOCATION}" "${LINK_DIR}")
      list (APPEND INSTALL_RPATH "${ORIGIN}/${RPATH}")
    endif ()
  endforeach ()
  # directories of link libraries
  #
  # only the libraries of this project and targets imported
  # from other projects which are part of the same bundle
  basis_get_target_link_libraries (LINK_DEPENDS ${TARGET_UID})
  if (BASIS_DEBUG AND BASIS_VERBOSE)
    message ("**    LINK_DEPENDS: [${LINK_DEPENDS}]")
  endif ()
  foreach (LINK_DEPEND ${LINK_DEPENDS})
    set (DEPEND_LOCATION)
    if (TARGET "${LINK_DEPEND}")
      basis_get_target_type (LINK_TYPE ${LINK_DEPEND})
      if ("^${LINK_TYPE}$" STREQUAL "^SHARED_LIBRARY$")
        basis_get_target_property (BUNDLED  ${LINK_DEPEND} BUNDLED)
        basis_get_target_property (IMPORTED ${LINK_DEPEND} IMPORTED)
        if (NOT IMPORTED OR BUNDLED)
          basis_get_target_location (DEPEND_LOCATION ${LINK_DEPEND} POST_INSTALL_PATH)
          if (BASIS_DEBUG AND BASIS_VERBOSE)
            message ("**    LOCATION(${LINK_DEPEND}): ${DEPEND_LOCATION}")
          endif ()
        endif ()
      endif ()
    elseif (IS_ABSOLUTE "${LINK_DEPEND}")
      if (IS_DIRECTORY "${LINK_DEPEND}")
        set (DEPEND_LOCATION "${LINK_DEPEND}")
      else ()
        get_filename_component (DEPEND_LOCATION "${LINK_DEPEND}" PATH)
      endif ()
      list (FIND BUNDLE_LINK_DIRS "${DEPEND_LOCATION}" IDX)
      if (IDX EQUAL -1)
        set (DEPEND_LOCATION)
      endif ()
    endif ()
    if (DEPEND_LOCATION)
      list (FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${DEPEND_LOCATION}" IDX)
      if (IDX EQUAL -1)
        basis_get_relative_path (RPATH "${TARGET_LOCATION}" "${DEPEND_LOCATION}")
        list (APPEND INSTALL_RPATH "${ORIGIN}/${RPATH}")
      endif ()
    endif ()
  endforeach ()
  # remove duplicates
  if (INSTALL_RPATH)
    list (REMOVE_DUPLICATES INSTALL_RPATH)
  endif ()
  # set INSTALL_RPATH property
  set_target_properties (${TARGET_UID} PROPERTIES INSTALL_RPATH "${INSTALL_RPATH}")
  if (BASIS_DEBUG)
    message ("**    INSTALL_RPATH: [${INSTALL_RPATH}]")
  endif ()
  if (BASIS_VERBOSE)
    message (STATUS "Setting INSTALL_RPATH property of ${TARGET_UID}... - done")
  endif ()
endfunction ()

# ----------------------------------------------------------------------------
## @brief Add custom command for build of single script.
#
# This function is called by basis_finalize_targets() which in turn is called
# by basis_project_end(), i.e., the end of the root CMake configuration file
# of the (sub-)project.
#
# @param [in] TARGET_UID Name/UID of custom target added by basis_add_script().
#
# @sa basis_add_script()
function (basis_build_script TARGET_UID)
  # does this target exist ?
  basis_get_target_uid (TARGET_UID "${TARGET_UID}")
  if (NOT TARGET "${TARGET_UID}")
    message (FATAL_ERROR "Unknown build target: ${TARGET_UID}")
  endif ()
  if (BASIS_VERBOSE)
    message (STATUS "Adding build command for target ${TARGET_UID}...")
  endif ()
  # get target properties
  basis_get_target_link_libraries (LINK_DEPENDS ${TARGET_UID}) # paths of script modules/packages
                                                               # including BASIS utilities if used
  set (
    PROPERTIES
      LANGUAGE                 # programming language of script
      BASIS_TYPE               # must match "^SCRIPT_(EXECUTABLE|LIBEXEC|MODULE)$"
      BUILD_DIRECTORY          # CMakeFiles build directory
      SOURCE_DIRECTORY         # CMake source directory
      BINARY_DIRECTORY         # CMake binary directory
      OUTPUT_DIRECTORY         # output directory for built script
      INSTALL_DIRECTORY        # installation directory for built script
      COMPONENT                # installation component
      OUTPUT_NAME              # name of built script including extension (if any)
      PREFIX                   # name prefix
      SUFFIX                   # name suffix (e.g., extension for executable script)
      SCRIPT_DEFINITIONS       # CMake code to set variables used to configure script
      SCRIPT_DEFINITIONS_FILE  # script configuration file
      EXPORT                   # whether this target shall be exported
      COMPILE                  # whether to compile script if applicable
      SOURCES                  # path of script source file
  )
  get_target_property (IS_TEST ${TARGET_UID} TEST) # whether this script is used for testing only
  foreach (PROPERTY ${PROPERTIES})
    get_target_property (${PROPERTY} ${TARGET_UID} ${PROPERTY})
  endforeach ()
  set (EXECUTABLE FALSE)
  set (LIBEXEC    FALSE)
  set (MODULE     FALSE)
  if (BASIS_TYPE MATCHES "^SCRIPT_(EXECUTABLE|LIBEXEC|MODULE)$")
    set (${CMAKE_MATCH_1} TRUE)
    if (LIBEXEC)
      set (EXECUTABLE TRUE)
    endif ()
  else ()
    message (FATAL_ERROR "Target ${TARGET_UID}: Unexpected BASIS_TYPE: ${BASIS_TYPE}")
  endif ()
  if (NOT BINARY_DIRECTORY)
    message (FATAL_ERROR "Target ${TARGET_UID}: Missing BINARY_DIRECTORY property!")
  endif ()
  file (RELATIVE_PATH _relpath "${CMAKE_BINARY_DIR}" "${BINARY_DIRECTORY}")
  if (_relpath MATCHES "^\\.\\./")
    message (FATAL_ERROR "Target ${TARGET_UID}: BINARY_DIRECTORY must be inside of build tree!")
  endif ()
  unset (_relpath)
  if (INSTALL_DIRECTORY AND NOT COMPONENT)
    set (COMPONENT "Unspecified")
  endif ()
  list (GET SOURCES 0 BUILD_DIR) # CMake <3.1 stores path to internal build directory here
  if (BUILD_DIR MATCHES "CMakeFiles")
    list (REMOVE_AT SOURCES 0)
  endif ()
  list (LENGTH SOURCES L)
  if (NOT L EQUAL 1)
    message (FATAL_ERROR "Target ${TARGET_UID}: Expected one element in SOURCES list!"
                         " Have you modified this (read-only) property or is your"
                         " (newer) CMake version not compatible with BASIS?")
  endif ()
  set (SOURCE_FILE "${SOURCES}")
  set (BUILD_DIR   "${BUILD_DIRECTORY}.dir")
  # output directory and name
  if (NOT OUTPUT_NAME)
    basis_get_target_name (OUTPUT_NAME ${TARGET_UID})
  endif ()
  if (PREFIX)
    set (OUTPUT_NAME "${PREFIX}${OUTPUT_NAME}")
  endif ()
  if (SUFFIX)
    set (OUTPUT_NAME "${OUTPUT_NAME}${SUFFIX}")
  endif ()
  if (CMAKE_GENERATOR MATCHES "Visual Studio|Xcode")
    set (OUTPUT_FILE "${BUILD_DIR}/build/${OUTPUT_NAME}")
    set (OUTPUT_DIR  "${OUTPUT_DIRECTORY}/$<${BASIS_GE_CONFIG}>")
  elseif (MODULE AND COMPILE)
    set (OUTPUT_FILE "${BUILD_DIR}/build/${OUTPUT_NAME}")
    set (OUTPUT_DIR  "${OUTPUT_DIRECTORY}")
  else ()
    set (OUTPUT_FILE "${OUTPUT_DIRECTORY}/${OUTPUT_NAME}")
    unset (OUTPUT_DIR)
  endif ()
  # arguments of build script
  if (INSTALL_DIRECTORY)
    if (MODULE)
      get_filename_component (SOURCE_NAME "${SOURCE_FILE}" NAME)
      set (INSTALL_FILE "${BUILD_DIR}/install/${SOURCE_NAME}")
      string (REGEX REPLACE "\\.in$" "" INSTALL_FILE "${INSTALL_FILE}")
    else ()
      set (INSTALL_FILE "${BUILD_DIR}/install/${OUTPUT_NAME}")
    endif ()
    set (DESTINATION "${INSTALL_DIRECTORY}")
    if (NOT IS_ABSOLUTE "${DESTINATION}")
      set (DESTINATION "${CMAKE_INSTALL_PREFIX}/${DESTINATION}")
    endif ()
  else ()
    set (INSTALL_FILE)
    set (DESTINATION)
  endif ()
  set (CONFIG_FILES)
  if (EXISTS "${BASIS_SCRIPT_CONFIG_FILE}")
    list (APPEND CONFIG_FILES "${BASIS_SCRIPT_CONFIG_FILE}")
  endif ()
  if (SCRIPT_DEFINITIONS_FILE)
    list (APPEND CONFIG_FILES "${SCRIPT_DEFINITIONS_FILE}")
  endif ()
  if (SCRIPT_DEFINITIONS)
    set (SCRIPT_CONFIG_FILE "${BUILD_DIR}/ScriptConfig.cmake")
    file (WRITE "${SCRIPT_CONFIG_FILE}.tmp" "# DO NOT edit. Automatically generated by BASIS.\n${SCRIPT_DEFINITIONS}\n")
    execute_process (COMMAND "${CMAKE_COMMAND}" -E copy_if_different "${SCRIPT_CONFIG_FILE}.tmp" "${SCRIPT_CONFIG_FILE}")
    file (REMOVE "${SCRIPT_CONFIG_FILE}.tmp")
    list (APPEND CONFIG_FILES "${SCRIPT_CONFIG_FILE}")
  endif ()
  set (CACHE_FILE "${BUILD_DIR}/cache.cmake")
  execute_process (COMMAND "${CMAKE_COMMAND}" -E copy_if_different "${CACHE_FILE}.tmp" "${CACHE_FILE}")
  file (REMOVE "${CACHE_FILE}.tmp")
  set (OPTIONS)
  foreach (FLAG IN ITEMS COMPILE EXECUTABLE)
    if (${FLAG})
      list (APPEND OPTIONS ${FLAG})
    endif ()
  endforeach ()
  # link dependencies - module search paths
  set (BUILD_LINK_DEPENDS)
  set (INSTALL_LINK_DEPENDS)
  foreach (LINK_DEPEND IN LISTS LINK_DEPENDS)
    basis_get_target_uid (UID "${LINK_DEPEND}")
    if (TARGET "${UID}")
      get_target_property (IMPORTED ${UID} IMPORTED)
      get_target_property (BUNDLED  ${UID} BUNDLED)
      if (IMPORTED AND NOT BUNDLED)
        basis_get_target_location (LOCATION "${UID}" ABSOLUTE)
        if (LOCATION)
          list (APPEND BUILD_LINK_DEPENDS   "${LOCATION}")
          list (APPEND INSTALL_LINK_DEPENDS "${LOCATION}")
        else ()
          message (WARNING "Could not determine build tree location of file corresponding to target ${UID}")
        endif ()
      else ()
        basis_get_target_location (LOCATION "${UID}" ABSOLUTE)
        if (LOCATION)
          list (APPEND BUILD_LINK_DEPENDS "${LOCATION}")
        else ()
          message (WARNING "Could not determine build tree location of file corresponding to target ${UID}")
        endif ()
        basis_get_target_property (LINK_DEPEND_IS_TEST "${UID}" TEST)
        if (NOT LINK_DEPEND_IS_TEST)
          basis_get_target_location (LOCATION "${UID}" POST_INSTALL)
          if (LOCATION)
            list (APPEND INSTALL_LINK_DEPENDS "relative ${LOCATION}")
          else ()
            message (WARNING "Could not determine installation location of file corresponding to target ${UID}")
          endif ()
        endif ()
      endif ()
    else ()
      list (APPEND BUILD_LINK_DEPENDS   "${LINK_DEPEND}")
      list (APPEND INSTALL_LINK_DEPENDS "${LINK_DEPEND}")
    endif ()
  endforeach ()
  # prepend own module search paths - if dependencies among own modules
  #                                   not specified or to ensure that
  #                                   these are preferred
  if (LANGUAGE MATCHES "JYTHON")
    if (BINARY_PYTHON_LIBRARY_DIR)
      list (INSERT BUILD_LINK_DEPENDS 0 "${BINARY_PYTHON_LIBRARY_DIR}")
    endif ()
    if (INSTALL_PYTHON_LIBRARY_DIR)
      list (INSERT INSTALL_LINK_DEPENDS 0 "relative ${CMAKE_INSTALL_PREFIX}/${INSTALL_PYTHON_LIBRARY_DIR}")
    endif ()
  endif ()
  if (BINARY_${LANGUAGE}_LIBRARY_DIR)
    list (INSERT BUILD_LINK_DEPENDS 0 "${BINARY_${LANGUAGE}_LIBRARY_DIR}")
  endif ()
  if (INSTALL_${LANGUAGE}_SITE_DIR AND NOT BASIS_INSTALL_SITE_PACKAGES)
    if (IS_ABSOLUTE "${INSTALL_${LANGUAGE}_SITE_DIR}")
      list (INSERT INSTALL_LINK_DEPENDS 0 "${INSTALL_${LANGUAGE}_SITE_DIR}")
    else ()
      list (INSERT INSTALL_LINK_DEPENDS 0 "relative ${CMAKE_INSTALL_PREFIX}/${INSTALL_${LANGUAGE}_SITE_DIR}")
    endif ()
  endif ()
  if (INSTALL_${LANGUAGE}_LIBRARY_DIR)
    list (INSERT INSTALL_LINK_DEPENDS 0 "relative ${CMAKE_INSTALL_PREFIX}/${INSTALL_${LANGUAGE}_LIBRARY_DIR}")
  endif ()
  list (INSERT BUILD_LINK_DEPENDS   0 "${BINARY_LIBRARY_DIR}")
  list (INSERT INSTALL_LINK_DEPENDS 0 "relative ${CMAKE_INSTALL_PREFIX}/${INSTALL_LIBRARY_DIR}")
  if (IS_TEST)
    list (INSERT BUILD_LINK_DEPENDS 0 "${TESTING_LIBRARY_DIR}")
  endif ()
  if (BUILD_LINK_DEPENDS)
    list (REMOVE_DUPLICATES BUILD_LINK_DEPENDS)
  endif ()
  if (INSTALL_LINK_DEPENDS)
    list (REMOVE_DUPLICATES INSTALL_LINK_DEPENDS)
  endif ()
  # remove default site-packages directories
  if (LANGUAGE MATCHES "[JP]YTHON|PERL")
    list (REMOVE_ITEM BUILD_LINK_DEPENDS   "${${LANGUAGE}_SITELIB}")
    list (REMOVE_ITEM INSTALL_LINK_DEPENDS "${${LANGUAGE}_SITELIB}")
  endif ()
  # configure build script
  set (BUILD_SCRIPT "${BUILD_DIR}/build.cmake")
  configure_file ("${BASIS_MODULE_PATH}/configure_script.cmake.in" "${BUILD_SCRIPT}" @ONLY)
  # list of all output files
  set (OUTPUT_FILES "${OUTPUT_FILE}")
  if (INSTALL_FILE)
    list (APPEND OUTPUT_FILES "${INSTALL_FILE}")
  endif ()
  if (MODULE AND COMPILE)
    basis_get_compiled_file (OUTPUT_CFILE  "${OUTPUT_FILE}"  ${LANGUAGE})
    basis_get_compiled_file (INSTALL_CFILE "${INSTALL_FILE}" ${LANGUAGE})
    if (OUTPUT_CFILE)
      list (APPEND OUTPUT_FILES "${OUTPUT_CFILE}")
    endif ()
    if (INSTALL_CFILE)
      list (APPEND OUTPUT_FILES "${INSTALL_CFILE}")
    endif ()
    if (LANGUAGE MATCHES "PYTHON")
      basis_compile_python_modules_for_jython (RV)
      if (RV)
        basis_get_compiled_jython_file_of_python_module (JYTHON_OUTPUT_CFILE "${OUTPUT_FILE}")
        basis_get_compiled_file (JYTHON_INSTALL_CFILE "${INSTALL_FILE}" JYTHON)
        if (JYTHON_OUTPUT_CFILE)
          list (APPEND OUTPUT_FILES "${JYTHON_OUTPUT_CFILE}")
        endif ()
        if (JYTHON_INSTALL_CFILE)
          list (APPEND OUTPUT_FILES "${JYTHON_INSTALL_CFILE}")
        endif ()
      endif ()
    endif ()
  endif ()
  # add build command for script
  if (OUTPUT_CFILE)
    file (RELATIVE_PATH REL "${CMAKE_BINARY_DIR}" "${OUTPUT_CFILE}")
  else ()
    file (RELATIVE_PATH REL "${CMAKE_BINARY_DIR}" "${OUTPUT_FILE}")
  endif ()
  if (LANGUAGE MATCHES "UNKNOWN")
    set (COMMENT "Building script ${REL}...")
  elseif (MODULE)
    set (COMMENT "Building ${LANGUAGE} module ${REL}...")
  else ()
    set (COMMENT "Building ${LANGUAGE} executable ${REL}...")
  endif ()
  add_custom_command (
    OUTPUT          ${OUTPUT_FILES}
    COMMAND         "${CMAKE_COMMAND}" -D "CONFIGURATION:STRING=$<${BASIS_GE_CONFIG}>" -P "${BUILD_SCRIPT}"
    MAIN_DEPENDENCY "${SOURCE_FILE}"
    DEPENDS         "${BUILD_SCRIPT}" "${BASIS_MODULE_PATH}/CommonTools.cmake" # basis_configure_script() definition
    COMMENT         "${COMMENT}"
    VERBATIM
  )
  # add custom target
  add_custom_target (_${TARGET_UID} DEPENDS ${OUTPUT_FILES})
  foreach (T IN LISTS LINK_DEPENDS)
    if (TARGET ${T})
      add_dependencies (_${TARGET_UID} ${T})
    endif ()
  endforeach ()
  add_dependencies (${TARGET_UID} _${TARGET_UID})
  # cleanup on "make clean" - including compiled files regardless of COMPILE flag
  set_property (DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES ${OUTPUT_FILES})
  foreach (OUTPUT_FILE IN LISTS OUTPUT_FILES)
    basis_get_compiled_file (CFILE "${OUTPUT_FILE}" ${LANGUAGE})
    if (CFILE)
      list (FIND OUTPUT_FILES "${CFILE}" IDX)
      if (IDX EQUAL -1)
        set_property (DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES "${CFILE}")
      endif ()
    endif ()
  endforeach ()
  # copy configured (and compiled) script to (configuration specific) output directory
  if (OUTPUT_DIR)
    if (OUTPUT_CFILE)
      get_filename_component (OUTPUT_CNAME "${OUTPUT_CFILE}" NAME)
      add_custom_command (
        TARGET _${TARGET_UID} POST_BUILD
        COMMAND "${CMAKE_COMMAND}" -E copy "${OUTPUT_CFILE}" "${OUTPUT_DIR}/${OUTPUT_CNAME}"
      )
      set_property (DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES "${OUTPUT_DIR}/${OUTPUT_CNAME}")
      if (JYTHON_OUTPUT_CFILE)
        get_filename_component (OUTPUT_CNAME "${JYTHON_OUTPUT_CFILE}" NAME)
        add_custom_command (
          TARGET _${TARGET_UID} POST_BUILD
          COMMAND "${CMAKE_COMMAND}" -E copy "${JYTHON_OUTPUT_CFILE}" "${OUTPUT_DIR}/${OUTPUT_CNAME}"
        )
        set_property (DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES "${OUTPUT_DIR}/${OUTPUT_CNAME}")
      endif ()
    else ()
      add_custom_command (
        TARGET _${TARGET_UID} POST_BUILD
        COMMAND "${CMAKE_COMMAND}" -E copy "${OUTPUT_FILE}" "${OUTPUT_DIR}/${OUTPUT_NAME}"
      )
      set_property (DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES "${OUTPUT_DIR}/${OUTPUT_NAME}")
    endif ()
  endif ()
  # export target
  if (EXPORT)
    basis_add_custom_export_target (${TARGET_UID} "${IS_TEST}")
  endif ()
  # install script
  if (INSTALL_DIRECTORY)
    if (INSTALL_CFILE)
      if (MODULE AND LANGUAGE MATCHES "PYTHON")
        basis_compile_python_modules_for_jython (RV)
        if (RV)
          basis_get_compiled_script (INSTALL_JYTHON_CFILE "${INSTALL_FILE}" JYTHON)
          basis_sanitize_for_regex (LIBRE  "${INSTALL_PYTHON_LIBRARY_DIR}")
          basis_sanitize_for_regex (SITERE "${INSTALL_PYTHON_SITE_DIR}")
          if (INSTALL_CFILE MATCHES "^${LIBRE}/" AND INSTALL_JYTHON_LIBRARY_DIR)
            string (REGEX REPLACE "^${LIBRE}/" "${INSTALL_JYTHON_LIBRARY_DIR}/" INSTALL_DIRECTORY_JYTHON "${INSTALL_DIRECTORY}")
          elseif (INSTALL_CFILE MATCHES "^${SITERE}/" AND INSTALL_JYTHON_SITE_DIR)
            string (REGEX REPLACE "^${SITERE}/" "${INSTALL_JYTHON_SITE_DIR}/" INSTALL_DIRECTORY_JYTHON "${INSTALL_DIRECTORY}")
          else ()
            set (INSTALL_DIRECTORY_JYTHON "${INSTALL_DIRECTORY}")
          endif ()
          install (
            FILES       "${INSTALL_JYTHON_CFILE}"
            DESTINATION "${INSTALL_DIRECTORY_JYTHON}"
            COMPONENT   "${COMPONENT}"
          )
        endif ()
      endif ()
      set (INSTALL_FILE "${INSTALL_CFILE}")
    elseif (NOT INSTALL_FILE)
      set (INSTALL_FILE "${OUTPUT_FILE}")
    endif ()
    if (MODULE)
      set (INSTALLTYPE FILES)
    else ()
      set (INSTALLTYPE PROGRAMS)
    endif ()
    install (
      ${INSTALLTYPE} "${INSTALL_FILE}"
      DESTINATION    "${INSTALL_DIRECTORY}"
      COMPONENT      "${COMPONENT}"
      RENAME         "${OUTPUT_NAME}"
    )
  endif ()
  # done
  if (BASIS_VERBOSE)
    message (STATUS "Adding build command for target ${TARGET_UID}... - done")
  endif ()
endfunction ()

# ----------------------------------------------------------------------------
## @brief Add custom command for build of script library.
#
# This function is called by basis_finalize_targets() which in turn is called
# by basis_project_end(), i.e., the end of the root CMake configuration file
# of the (sub-)project.
#
# @param [in] TARGET_UID Name/UID of custom target added by basis_add_script_library().
#
# @sa basis_add_script_library()
function (basis_build_script_library TARGET_UID)
  # does this target exist ?
  basis_get_target_uid (TARGET_UID "${TARGET_UID}")
  if (NOT TARGET "${TARGET_UID}")
    message (FATAL_ERROR "Unknown target: ${TARGET_UID}")
  endif ()
  if (BASIS_VERBOSE)
    message (STATUS "Adding build command for target ${TARGET_UID}...")
  endif ()
  # get target properties
  basis_get_target_link_libraries (LINK_DEPENDS ${TARGET_UID}) # paths of script modules/packages
                                                               # including BASIS utilities if used
  set (
    PROPERTIES
      LANGUAGE                   # programming language of modules
      BASIS_TYPE                 # must be "SCRIPT_LIBRARY"
      BASIS_UTILITIES            # whether this target requires the BASIS utilities
      BUILD_DIRECTORY            # CMakeFiles build directory
      SOURCE_DIRECTORY           # CMake source directory
      BINARY_DIRECTORY           # CMake binary directory
      LIBRARY_OUTPUT_DIRECTORY   # output directory for built modules
      LIBRARY_INSTALL_DIRECTORY  # installation directory for built modules
      LIBRARY_COMPONENT          # installation component
      PREFIX                     # common prefix for modules
      SCRIPT_DEFINITIONS         # CMake code to set variables used to configure modules
      SCRIPT_DEFINITIONS_FILE    # script configuration file
      LINK_DEPENDS               # paths of script modules/packages used by the modules of this library
      EXPORT                     # whether to export this target
      COMPILE                    # whether to compile the modules/library if applicable
      SOURCES                    # source files of module scripts
  )
  get_target_property (IS_TEST ${TARGET_UID} TEST) # whether this script is used for testing only
  foreach (PROPERTY ${PROPERTIES})
    get_target_property (${PROPERTY} ${TARGET_UID} ${PROPERTY})
  endforeach ()
  if (NOT BASIS_TYPE MATCHES "^SCRIPT_LIBRARY$")
    message (FATAL_ERROR "Target ${TARGET_UID}: Unexpected BASIS_TYPE: ${BASIS_TYPE}")
  endif ()
  if (NOT SOURCE_DIRECTORY)
    message (FATAL_ERROR "Missing SOURCE_DIRECTORY property!")
  endif ()
  if (NOT LIBRARY_OUTPUT_DIRECTORY)
    message (FATAL_ERROR "Missing LIBRARY_OUTPUT_DIRECTORY property!")
  endif ()
  if (NOT IS_ABSOLUTE "${LIBRARY_OUTPUT_DIRECTORY}")
    set (LIBRARY_OUTPUT_DIRECTORY "${TOPLEVEL_PROJECT_BINARY_DIR}/${LIBRARY_OUTPUT_DIRECTORY}")
  endif ()
  file (RELATIVE_PATH _relpath "${CMAKE_BINARY_DIR}" "${LIBRARY_OUTPUT_DIRECTORY}")
  if (_relpath MATCHES "^\\.\\./")
    message (FATAL_ERROR "Output directory LIBRARY_OUTPUT_DIRECTORY is outside the build tree!")
  endif ()
  unset(_relpath)
  if (NOT LIBRARY_COMPONENT)
    set (LIBRARY_COMPONENT "Unspecified")
  endif ()
  list (GET SOURCES 0 BUILD_DIR) # CMake <3.1 stores path to internal build directory here
  if (BUILD_DIR MATCHES "CMakeFiles")
    list (REMOVE_AT SOURCES 0)
  endif ()
  if (NOT SOURCES)
    message (FATAL_ERROR "Target ${TARGET_UID}: Expected at least one element in SOURCES list!"
                         " Have you incorrectly modified this (read-only) property or"
                         " is your (newer) CMake version not compatible with BASIS?")
  endif ()
  set (BUILD_DIR "${BUILD_DIRECTORY}.dir")
  # common arguments of build script
  set (CONFIG_FILES)
  if (EXISTS "${BASIS_SCRIPT_CONFIG_FILE}")
    list (APPEND CONFIG_FILES "${BASIS_SCRIPT_CONFIG_FILE}")
  endif ()
  if (SCRIPT_DEFINITIONS_FILE)
    list (APPEND CONFIG_FILES ${SCRIPT_DEFINITIONS_FILE})
  endif ()
  if (SCRIPT_DEFINITIONS)
    set (SCRIPT_CONFIG_FILE "${BUILD_DIR}/ScriptConfig.cmake")
    file (WRITE "${SCRIPT_CONFIG_FILE}.tmp" "# DO NOT edit. Automatically generated by BASIS.\n${SCRIPT_DEFINITIONS}\n")
    execute_process (COMMAND "${CMAKE_COMMAND}" -E copy_if_different "${SCRIPT_CONFIG_FILE}.tmp" "${SCRIPT_CONFIG_FILE}")
    file (REMOVE "${SCRIPT_CONFIG_FILE}.tmp")
    list (APPEND CONFIG_FILES "${SCRIPT_CONFIG_FILE}")
  endif ()
  set (CACHE_FILE "${BUILD_DIR}/cache.cmake")
  execute_process (COMMAND "${CMAKE_COMMAND}" -E copy_if_different "${CACHE_FILE}.tmp" "${CACHE_FILE}")
  file (REMOVE "${CACHE_FILE}.tmp")
  set (OPTIONS) # no additional options
  if (COMPILE)
    list (APPEND OPTIONS COMPILE)
  endif ()
  # add build command for each module
  set (OUTPUT_FILES)                                    # list of all output files
  set (FILES_TO_COPY)                                   # relative paths of build tree modules
  set (FILES_TO_INSTALL)                                # list of output files for installation
  set (BINARY_INSTALL_DIRECTORY "${BUILD_DIR}/install") # common base directory for files to install 
  if (COMPILE OR CMAKE_GENERATOR MATCHES "Visual Studio|Xcode") # post-build copy to <libdir>/<config>/
    set (BINARY_OUTPUT_DIRECTORY "${BUILD_DIR}/build")  # common base directory for build tree files
    set (POST_BUILD_COPY TRUE)
  else ()
    set (BINARY_OUTPUT_DIRECTORY "${LIBRARY_OUTPUT_DIRECTORY}")
    set (POST_BUILD_COPY FALSE)
  endif ()
  foreach (SOURCE_FILE IN LISTS SOURCES)
    file (RELATIVE_PATH S "${SOURCE_DIRECTORY}" "${SOURCE_FILE}")
    string (REGEX REPLACE "\\.in$" "" S "${S}")
    basis_get_source_target_name (BUILD_SCRIPT_NAME "build_${S}")
    # arguments of build script
    if (PREFIX)
      set (S "${PREFIX}${S}")
    endif ()
    set (OUTPUT_FILE "${BINARY_OUTPUT_DIRECTORY}/${S}")
    get_filename_component (OUTPUT_DIR "${LIBRARY_OUTPUT_DIRECTORY}/${S}" PATH)
    if (LIBRARY_INSTALL_DIRECTORY)
      set (INSTALL_FILE "${BINARY_INSTALL_DIRECTORY}/${S}")
      set (DESTINATION  "${LIBRARY_INSTALL_DIRECTORY}/${S}")
      if (NOT IS_ABSOLUTE "${DESTINATION}")
        set (DESTINATION "${CMAKE_INSTALL_PREFIX}/${DESTINATION}")
      endif ()
      get_filename_component (DESTINATION "${DESTINATION}" PATH)
    else ()
      set (INSTALL_FILE)
      set (DESTINATION)
    endif ()
    # configure build script
    set (BUILD_SCRIPT "${BUILD_DIR}/${BUILD_SCRIPT_NAME}.cmake")
    configure_file ("${BASIS_MODULE_PATH}/configure_script.cmake.in" "${BUILD_SCRIPT}" @ONLY)
    # output files of this command
    set (_OUTPUT_FILES "${OUTPUT_FILE}")
    if (INSTALL_FILE)
      list (APPEND _OUTPUT_FILES "${INSTALL_FILE}")
    endif ()
    if (COMPILE)
      basis_get_compiled_file (OUTPUT_CFILE  "${OUTPUT_FILE}"  ${LANGUAGE})
      basis_get_compiled_file (INSTALL_CFILE "${INSTALL_FILE}" ${LANGUAGE})
      if (OUTPUT_CFILE)
        list (APPEND _OUTPUT_FILES "${OUTPUT_CFILE}")
      endif ()
      if (INSTALL_CFILE)
        list (APPEND _OUTPUT_FILES "${INSTALL_CFILE}")
      endif ()
      if (LANGUAGE MATCHES "PYTHON")
        basis_compile_python_modules_for_jython (RV)
        if (RV)
          basis_get_compiled_jython_file_of_python_module (JYTHON_OUTPUT_CFILE "${OUTPUT_FILE}")
          basis_get_compiled_file (JYTHON_INSTALL_CFILE "${INSTALL_FILE}" JYTHON)
          if (JYTHON_OUTPUT_CFILE)
            list (APPEND _OUTPUT_FILES "${JYTHON_OUTPUT_CFILE}")
          endif ()
          if (JYTHON_INSTALL_CFILE)
            list (APPEND _OUTPUT_FILES "${JYTHON_INSTALL_CFILE}")
          endif ()
        endif ()
      endif ()
    endif ()
    if (POST_BUILD_COPY)
      if (OUTPUT_CFILE)
        file (RELATIVE_PATH _file "${BINARY_OUTPUT_DIRECTORY}" "${OUTPUT_CFILE}")
        list (APPEND FILES_TO_COPY "${_file}")
        if (JYTHON_OUTPUT_CFILE)
          file (RELATIVE_PATH _file "${BINARY_OUTPUT_DIRECTORY}" "${JYTHON_OUTPUT_CFILE}")
          list (APPEND FILES_TO_COPY "${_file}")
        endif ()
      else ()
        list (APPEND FILES_TO_COPY "${S}")
      endif ()
    endif ()
    if (INSTALL_CFILE)
      list (APPEND FILES_TO_INSTALL "${INSTALL_CFILE}")
    elseif (INSTALL_FILE)
      list (APPEND FILES_TO_INSTALL "${INSTALL_FILE}")
    endif ()
    # add build command
    if (OUTPUT_CFILE)
      file (RELATIVE_PATH REL "${CMAKE_BINARY_DIR}" "${OUTPUT_CFILE}")
    else ()
      file (RELATIVE_PATH REL "${CMAKE_BINARY_DIR}" "${OUTPUT_FILE}")
    endif ()
    set (COMMENT "Building ${LANGUAGE} module ${REL}...")
    add_custom_command (
      OUTPUT          ${_OUTPUT_FILES}
      COMMAND         "${CMAKE_COMMAND}" -D "CONFIGURATION=$<${BASIS_GE_CONFIG}>" -P "${BUILD_SCRIPT}"
      MAIN_DEPENDENCY "${SOURCE_FILE}"
      DEPENDS         "${BUILD_SCRIPT}" "${BASIS_MODULE_PATH}/CommonTools.cmake" # basis_configure_script() definition
      COMMENT         "${COMMENT}"
      VERBATIM
    )
    # add output files of command to list of all output files
    list (APPEND OUTPUT_FILES ${_OUTPUT_FILES})
  endforeach ()
  # add custom target to build modules
  add_custom_target (_${TARGET_UID} DEPENDS ${OUTPUT_FILES})
  foreach (T IN LISTS LINK_DEPENDS)
    if (TARGET ${T})
      add_dependencies (_${TARGET_UID} ${T})
    endif ()
  endforeach ()
  add_dependencies (${TARGET_UID} _${TARGET_UID})
  # copy configured modules to output directory
  foreach (OUTPUT_FILE IN LISTS FILES_TO_COPY)
    add_custom_command (
      TARGET _${TARGET_UID} POST_BUILD
      COMMAND "${CMAKE_COMMAND}" -E copy "${BINARY_OUTPUT_DIRECTORY}/${OUTPUT_FILE}" "${LIBRARY_OUTPUT_DIRECTORY}/${OUTPUT_FILE}"
    )
  endforeach ()
  # cleanup on "make clean" - including compiled files regardless of COMPILE flag
  set_property (DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES ${OUTPUT_FILES})
  foreach (OUTPUT_FILE IN LISTS OUTPUT_FILES)
    basis_get_compiled_file (CFILE "${OUTPUT_FILE}" ${LANGUAGE})
    if (CFILE)
      list (FIND OUTPUT_FILES "${CFILE}" IDX)
      if (IDX EQUAL -1)
        set_property (DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES "${CFILE}")
      endif ()
    endif ()
  endforeach ()
  # export target
  if (EXPORT)
    basis_add_custom_export_target (${TARGET_UID} "${IS_TEST}")
  endif ()
  # add installation rule
  foreach (INSTALL_FILE IN LISTS FILES_TO_INSTALL)
    get_filename_component (D "${INSTALL_FILE}" PATH)
    file (RELATIVE_PATH D "${BINARY_INSTALL_DIRECTORY}" "${D}")
    install (
      FILES       "${INSTALL_FILE}"
      DESTINATION "${LIBRARY_INSTALL_DIRECTORY}/${D}"
      COMPONENT   "${LIBRARY_COMPONENT}"
    )
    if (LANGUAGE MATCHES "PYTHON" AND INSTALL_FILE MATCHES "\\.pyc$")
      basis_compile_python_modules_for_jython (RV)
      if (RV)
        basis_sanitize_for_regex (LIBRE  "${INSTALL_PYTHON_LIBRARY_DIR}")
        basis_sanitize_for_regex (SITERE "${INSTALL_PYTHON_SITE_DIR}")
        if (LIBRARY_INSTALL_DIRECTORY MATCHES "^${LIBRE}/*$")
          set (JYTHON_INSTALL_DIRECTORY "${INSTALL_JYTHON_LIBRARY_DIR}")
        elseif (LIBRARY_INSTALL_DIRECTORY MATCHES "^${SITERE}/*$")
          set (JYTHON_INSTALL_DIRECTORY "${INSTALL_JYTHON_SITE_DIR}")
        else ()
          set (JYTHON_INSTALL_DIRECTORY "${LIBRARY_INSTALL_DIRECTORY}")
        endif ()
        string (REGEX REPLACE "c$" "" SOURCE_FILE "${INSTALL_FILE}")
        basis_get_compiled_file (INSTALL_JYTHON_CFILE "${SOURCE_FILE}" JYTHON)
        install (
          FILES       "${INSTALL_JYTHON_CFILE}"
          DESTINATION "${JYTHON_INSTALL_DIRECTORY}/${D}"
          COMPONENT   "${LIBRARY_COMPONENT}"
        )
      endif ()
    endif ()
  endforeach ()
  # done
  if (BASIS_VERBOSE)
    message (STATUS "Adding build command for target ${TARGET_UID}... - done")
  endif ()
endfunction ()

# ----------------------------------------------------------------------------
# @brief Add target to build/install __init__.py files.
function (basis_add_init_py_target)
  if (NOT BASIS_PYTHON_TEMPLATES_DIR)
    message (WARNING "BASIS_PYTHON_TEMPLATES_DIR not set, skipping basis_add_init_py_target")
    return ()
  endif ()
  # constants
  set (BUILD_DIR "${PROJECT_BINARY_DIR}/CMakeFiles/_initpy.dir")
  basis_sanitize_for_regex (BINARY_PYTHON_LIBRARY_DIR_RE  "${BINARY_PYTHON_LIBRARY_DIR}")
  basis_sanitize_for_regex (TESTING_PYTHON_LIBRARY_DIR_RE "${TESTING_PYTHON_LIBRARY_DIR}")
  basis_sanitize_for_regex (INSTALL_PYTHON_LIBRARY_DIR_RE "${INSTALL_PYTHON_LIBRARY_DIR}")
  basis_sanitize_for_regex (BINARY_JYTHON_LIBRARY_DIR_RE  "${BINARY_JYTHON_LIBRARY_DIR}")
  basis_sanitize_for_regex (TESTING_JYTHON_LIBRARY_DIR_RE "${TESTING_JYTHON_LIBRARY_DIR}")
  basis_sanitize_for_regex (INSTALL_JYTHON_LIBRARY_DIR_RE "${INSTALL_JYTHON_LIBRARY_DIR}")
  # collect directories requiring a __init__.py file
  set (DEPENDENTS)          # targets which generate Python/Jython modules and depend on _initpy
  set (PYTHON_DIRS)         # Python library directories requiring a __init__.py file
  set (JYTHON_DIRS)         # Jython library directories requiring a __init__.py file
  set (EXCLUDE)             # exclude these directories
  set (INSTALL_EXCLUDE)     # exclude these directories upon installation
  set (COMPONENTS)          # installation components
  basis_get_project_property (TARGETS PROPERTY TARGETS)
  foreach (TARGET_UID IN LISTS TARGETS)
    get_target_property (BASIS_TYPE ${TARGET_UID} BASIS_TYPE)
    get_target_property (LANGUAGE   ${TARGET_UID} LANGUAGE)
    if (BASIS_TYPE MATCHES "MODULE|LIBRARY" AND LANGUAGE MATCHES "^[JP]YTHON$")
      # get path of built Python modules
      basis_get_target_location (LOCATION         ${TARGET_UID} ABSOLUTE)
      basis_get_target_location (INSTALL_LOCATION ${TARGET_UID} POST_INSTALL_RELATIVE)
      if (BASIS_TYPE MATCHES "^SCRIPT_LIBRARY$")
        get_target_property (PREFIX           ${TARGET_UID} PREFIX)
        get_target_property (SOURCES          ${TARGET_UID} SOURCES)
        get_target_property (SOURCE_DIRECTORY ${TARGET_UID} SOURCE_DIRECTORY)
        string (REGEX REPLACE "/+$" "" PREFIX "${PREFIX}")
        if (PREFIX)
          set (LOCATION         "${LOCATION}/${PREFIX}")
          set (INSTALL_LOCATION "${INSTALL_LOCATION}/${PREFIX}")
        endif ()
        list (GET SOURCES 0 BINARY_DIR) # CMake <3.1 stores path to internal build directory here
        if (BINARY_DIR MATCHES "CMakeFiles")
          list (REMOVE_AT SOURCES 0)
        endif ()
        foreach (SOURCE IN LISTS SOURCES)
          file (RELATIVE_PATH SOURCE "${SOURCE_DIRECTORY}" "${SOURCE}")
          list (APPEND _LOCATION         "${LOCATION}/${SOURCE}")
          list (APPEND _INSTALL_LOCATION "${INSTALL_LOCATION}/${SOURCE}")
        endforeach ()
        set (LOCATION         "${_LOCATION}")
        set (INSTALL_LOCATION "${_INSTALL_LOCATION}")
      endif ()
      # get component (used by installation rule)
      get_target_property (COMPONENT ${TARGET_UID} "LIBRARY_COMPONENT")
      list (FIND COMPONENTS "${COMPONENT}" IDX)
      if (IDX EQUAL -1)
        list (APPEND COMPONENTS "${COMPONENT}")
        set (INSTALL_${LANGUAGE}_DIRS_${COMPONENT}) # list of directories for which to install
                                                    # __init__.py for this component
      endif ()
      # directories for which to build a __init__.py file
      foreach (L IN LISTS LOCATION)
        basis_get_filename_component (DIR "${L}" PATH)
        if (L MATCHES "/__init__\\.py$")
          list (APPEND EXCLUDE "${DIR}")
        else ()
          list (APPEND DEPENDENTS ${TARGET_UID}) # depends on _initpy
          if (BINARY_${LANGUAGE}_LIBRARY_DIR_RE AND DIR MATCHES "^${BINARY_${LANGUAGE}_LIBRARY_DIR_RE}/.+")
            while (NOT "${DIR}" MATCHES "^${BINARY_${LANGUAGE}_LIBRARY_DIR_RE}$")
              list (APPEND ${LANGUAGE}_DIRS "${DIR}")
              get_filename_component (DIR "${DIR}" PATH)
            endwhile ()
          elseif (TESTING_${LANGUAGE}_LIBRARY_DIR_RE AND DIR MATCHES "^${TESTING_${LANGUAGE}_LIBRARY_DIR_RE}/.+")
            while (NOT "${DIR}" MATCHES "^${TESTING_${LANGUAGE}_LIBRARY_DIR_RE}$")
              list (APPEND ${LANGUAGE}_DIRS "${DIR}")
              get_filename_component (DIR "${DIR}" PATH)
            endwhile ()
          endif ()
        endif ()
      endforeach ()
      # directories for which to install a __init__.py file
      foreach (L IN LISTS INSTALL_LOCATION)
        basis_get_filename_component (DIR "${L}" PATH)
        if (L MATCHES "/__init__\\.py$")
          list (APPEND INSTALL_EXCLUDE "${DIR}")
        else ()
          list (APPEND DEPENDENTS ${TARGET_UID}) # depends on _initpy
          if (INSTALL_${LANGUAGE}_LIBRARY_DIR_RE AND DIR MATCHES "^${INSTALL_${LANGUAGE}_LIBRARY_DIR_RE}/.+")
            while (NOT "${DIR}" MATCHES "^${INSTALL_${LANGUAGE}_LIBRARY_DIR_RE}$")
              list (APPEND INSTALL_${LANGUAGE}_DIRS_${COMPONENT} "${DIR}")
              if (BASIS_COMPILE_SCRIPTS AND LANGUAGE MATCHES "PYTHON")
                basis_compile_python_modules_for_jython (RV)
                if (RV)
                  if (INSTALL_JYTHON_LIBRARY_DIR)
                    file (RELATIVE_PATH REL "${CMAKE_INSTALL_PREFIX}/${INSTALL_PYTHON_LIBRARY_DIR}" "${CMAKE_INSTALL_PREFIX}/${DIR}")
                  else ()
                    set (REL)
                  endif ()
                  if (NOT REL MATCHES "^$|^\\.\\./")
                    list (APPEND INSTALL_JYTHON_DIRS_${COMPONENT} "${INSTALL_JYTHON_LIBRARY_DIR}/${REL}")
                  endif ()
                endif ()
              endif ()
              get_filename_component (DIR "${DIR}" PATH)
            endwhile ()
          endif ()
        endif ()
      endforeach ()
    endif ()
  endforeach ()
  if (DEPENDENTS)
    list (REMOVE_DUPLICATES DEPENDENTS)
  endif ()
  # return if nothing to do
  if (NOT PYTHON_DIRS AND NOT JYTHON_DIRS)
    return ()
  endif ()
  if (PYTHON_DIRS)
    list (REMOVE_DUPLICATES PYTHON_DIRS)
  endif ()
  if (JYTHON_DIRS)
    list (REMOVE_DUPLICATES JYTHON_DIRS)
  endif ()
  if (EXCLUDE)
    list (REMOVE_DUPLICATES EXCLUDE)
  endif ()
  if (INSTALL_EXCLUDE)
    list (REMOVE_DUPLICATES INSTALL_EXCLUDE)
  endif ()
  # generate build script
  set (PYTHON_COMPILED_FILES "${CFILE}")
  set (JYTHON_COMPILED_FILES "${CFILE}")
  set (C "configure_file (\"${BASIS_PYTHON_TEMPLATES_DIR}/__init__.py.in\" \"${BUILD_DIR}/__init__.py\" @ONLY)\n")
  foreach (LANGUAGE IN ITEMS PYTHON JYTHON) # Python *must* come first. See JYTHON_COMPILED_FILES list.
    set (${LANGUAGE}_OUTPUT_FILES)
    set (${LANGUAGE}_INSTALL_FILE "${BUILD_DIR}/__init__.py")
    set (C "${C}\nset (${LANGUAGE}_EXECUTABLE \"${${LANGUAGE}_EXECUTABLE}\")\n\n")
    foreach (DIR IN LISTS ${LANGUAGE}_DIRS)
      list (FIND EXCLUDE "${DIR}" IDX)
      if (IDX EQUAL -1)
        set (C "${C}configure_file (\"${BASIS_PYTHON_TEMPLATES_DIR}/__init__.py.in\" \"${DIR}/__init__.py\" @ONLY)\n")
        list (APPEND ${LANGUAGE}_OUTPUT_FILES "${DIR}/__init__.py")
      endif ()
    endforeach ()
    if (BASIS_COMPILE_SCRIPTS AND ${LANGUAGE}_EXECUTABLE)
      set (C "${C}\n")
      string (TOLOWER "${LANGUAGE}" language)
      basis_get_compiled_file (CFILE "${BUILD_DIR}/${language}/__init__.py" ${LANGUAGE})
      set (C "${C}file (MAKE_DIRECTORY \"${BUILD_DIR}/${language}\")\n")
      set (C "${C}execute_process (COMMAND \"${${LANGUAGE}_EXECUTABLE}\" -c \"import py_compile; py_compile.compile('${BUILD_DIR}/__init__.py', '${CFILE}')\")\n")
      set (${LANGUAGE}_INSTALL_FILE "${CFILE}")
      set (C "${C}\n")
      foreach (SFILE IN LISTS ${LANGUAGE}_OUTPUT_FILES)
        basis_get_compiled_file (CFILE "${SFILE}" ${LANGUAGE})
        set (C "${C}execute_process (COMMAND \"${${LANGUAGE}_EXECUTABLE}\" -c \"import py_compile; py_compile.compile('${SFILE}', '${CFILE}')\")\n")
        list (APPEND ${LANGUAGE}_COMPILED_FILES "${CFILE}")
        if (LANGUAGE MATCHES "PYTHON")
          basis_compile_python_modules_for_jython (RV)
          if (RV)
            if (BINARY_PYTHON_LIBRARY_DIR AND BINARY_JYTHON_LIBRARY_DIR)
              file (RELATIVE_PATH REL "${BINARY_PYTHON_LIBRARY_DIR}" "${SFILE}")
            else ()
              set (REL)
            endif ()
            if (NOT REL MATCHES "^$|^\\.\\./")
              basis_get_compiled_file (CFILE "${BINARY_JYTHON_LIBRARY_DIR}/${REL}" JYTHON)
            else ()
              basis_get_compiled_file (CFILE "${SFILE}" JYTHON)
            endif ()
            get_filename_component (CDIR "${CFILE}" PATH)
            set (C "${C}file (MAKE_DIRECTORY \"${CDIR}\")\n")
            set (C "${C}execute_process (COMMAND \"${JYTHON_EXECUTABLE}\" -c \"import py_compile; py_compile.compile('${SFILE}', '${CFILE}')\")\n")
            list (APPEND JYTHON_COMPILED_FILES "${CFILE}")
          endif ()
        endif ()
      endforeach ()
      list (APPEND ${LANGUAGE}_OUTPUT_FILES ${${LANGUAGE}_COMPILED_FILES})
    endif ()
  endforeach ()
  # write/update build script
  set (BUILD_SCRIPT "${BUILD_DIR}/build.cmake")
  if (EXISTS "${BUILD_SCRIPT}")
    file (WRITE "${BUILD_SCRIPT}.tmp" "${C}")
    execute_process (COMMAND "${CMAKE_COMMAND}" -E copy_if_different "${BUILD_SCRIPT}.tmp" "${BUILD_SCRIPT}")
    file (REMOVE "${BUILD_SCRIPT}.tmp")
  else ()
    file (WRITE "${BUILD_SCRIPT}" "${C}")
  endif ()
  # add custom build command
  add_custom_command (
    OUTPUT          "${BUILD_DIR}/__init__.py" ${PYTHON_OUTPUT_FILES} ${JYTHON_OUTPUT_FILES}
    COMMAND         "${CMAKE_COMMAND}" -P "${BUILD_SCRIPT}"
    MAIN_DEPENDENCY "${BASIS_PYTHON_TEMPLATES_DIR}/__init__.py.in"
    COMMENT         "Building PYTHON modules */__init__.py..."
  )
  # add custom target which triggers execution of build script
  add_custom_target (_initpy ALL DEPENDS "${BUILD_DIR}/__init__.py" ${PYTHON_OUTPUT_FILES} ${JYTHON_OUTPUT_FILES})
  if (BASIS_DEBUG)
    message ("** basis_add_init_py_target():")
  endif ()
  foreach (DEPENDENT IN LISTS DEPENDENTS)
    if (BASIS_DEBUG)
      message ("**    Adding dependency on _initpy target to ${DEPENDENT}")
    endif ()
    add_dependencies (${DEPENDENT} _initpy)
  endforeach ()
  # cleanup on "make clean" - including compiled modules regardless of COMPILE flag
  set_property (DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES "${BUILD_DIR}/__init__.py" ${PYTHON_OUTPUT_FILES} ${JYTHON_OUTPUT_FILES})
  # add install rules
  foreach (LANGUAGE IN ITEMS PYTHON JYTHON)
    foreach (COMPONENT IN LISTS COMPONENTS)
      if (INSTALL_${LANGUAGE}_DIRS_${COMPONENT})
        list (REMOVE_DUPLICATES INSTALL_${LANGUAGE}_DIRS_${COMPONENT})
      endif ()
      foreach (DIR IN LISTS INSTALL_${LANGUAGE}_DIRS_${COMPONENT})
        list (FIND INSTALL_EXCLUDE "${DIR}" IDX)
        if (IDX EQUAL -1)
          if (BASIS_DEBUG AND BASIS_VERBOSE)
            message("**    Copy ${${LANGUAGE}_INSTALL_FILE} to ${DIR} upon installation of ${COMPONENT} component")
          endif ()
          install (
            FILES       "${${LANGUAGE}_INSTALL_FILE}"
            DESTINATION "${DIR}"
            COMPONENT   "${COMPONENT}"
          )
        endif ()
      endforeach ()
    endforeach ()
  endforeach ()
endfunction ()


## @}
# end of Doxygen group
