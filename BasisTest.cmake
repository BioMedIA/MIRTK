# ============================================================================
# Copyright (c) 2011-2012 University of Pennsylvania
# Copyright (c) 2013-2016 Andreas Schuh
# All rights reserved.
#
# See COPYING file for license information or visit
# https://cmake-basis.github.io/download.html#license
# ============================================================================

##############################################################################
# @file  BasisTest.cmake
# @brief CTest configuration. Include this module instead of CTest.
#
# @note This module is included by basis_project_initialize().
#
# @ingroup CMakeAPI
##############################################################################

# ============================================================================
# configuration
# ============================================================================

## @brief Request build of tests.
option (BUILD_TESTING "Request build of tests." OFF)

# include CTest module which enables testing, but prevent it from generating
# any configuration file or adding targets yet as we want to adjust the
# default CTest settings--in particular the site name--before
set (RUN_FROM_DART 1)
include (CTest)

# mark timeout option as advanced
mark_as_advanced (DART_TESTING_TIMEOUT)

set (RUN_FROM_CTEST_OR_DART 1)
include (CTestTargets)
set (RUN_FROM_CTEST_OR_DART)

# custom CTest configuration
if (EXISTS "${PROJECT_CONFIG_DIR}/CTestCustom.cmake.in")
  configure_file (
    "${PROJECT_CONFIG_DIR}/CTestCustom.cmake.in"
    "${PROJECT_BINARY_DIR}/CTestCustom.cmake"
    @ONLY
  )
elseif (EXISTS "${PROJECT_CONFIG_DIR}/CTestCustom.cmake")
  configure_file (
    "${PROJECT_CONFIG_DIR}/CTestCustom.cmake"
    "${PROJECT_BINARY_DIR}/CTestCustom.cmake"
    COPYONLY
  )
else ()
  basis_get_relative_path (CONFIG_DIR "${PROJECT_SOURCE_DIR}" "${PROJECT_CONFIG_DIR}")
  if (EXISTS "${BASIS_TEMPLATE_DIR}/${CONFIG_DIR}/CTestCustom.cmake.in")
    # to avoid a Doxygen warning, we need to replace certain patterns used by
    # the basisproject tool to replace them with project related information
    #
    # Note: Do this only on the first pass, otherwise configure_file() will
    #       retrigger CMake every time b/c the modification timestamp of the
    #       source file is newer than the previously configured file.
    #       The use of the "cmake -E copy_if_different" command could be used
    #       here instead, but as we do not expect the CTestCustom.cmake.in
    #       file of the BASIS project template to change, it is sufficient
    #       to only check if we copied the template file already.
    if (NOT EXISTS "${PROJECT_BINARY_DIR}/CTestCustom.cmake.in")
      file (READ "${BASIS_TEMPLATE_DIR}/${CONFIG_DIR}/CTestCustom.cmake.in" _TEMPLATE)
      string (REGEX REPLACE "<year>" "" _TEMPLATE "${_TEMPLATE}")
      file (WRITE "${PROJECT_BINARY_DIR}/CTestCustom.cmake.in" "${_TEMPLATE}")
      unset (_TEMPLATE)
    endif ()
    # configure the modified template file
    configure_file (
      "${PROJECT_BINARY_DIR}/CTestCustom.cmake.in"
      "${PROJECT_BINARY_DIR}/CTestCustom.cmake"
      @ONLY
    )
  endif ()
endif ()

# ============================================================================
# constants
# ============================================================================

## @brief Names of recognized properties on tests.
#
# Unfortunately, the @c ARGV and @c ARGN arguments of a CMake function()
# or macro() does not preserve values which themselves are lists. Therefore,
# it is not possible to distinguish between property names and their values
# in the arguments passed to basis_set_tests_properties().
# To overcome this problem, this list specifies all the possible property names.
# Everything else is considered to be a property value except the first argument
# follwing right after the @c PROPERTIES keyword. Alternatively,
# basis_set_property() can be used as here no disambiguity exists.
#
# @note Placeholders such as &lt;CONFIG&gt; are allowed. These are treated
#       as the regular expression "[^ ]+". See basis_list_to_regex().
#
# @sa http://www.cmake.org/cmake/help/cmake-2-8-docs.html#section_PropertiesonTests
set (BASIS_PROPERTIES_ON_TESTS
  ATTACHED_FILES
  ATTACHED_FILES_ON_FAIL
  COST
  DEPENDS
  ENVIRONMENT
  FAIL_REGULAR_EXPRESSION
  LABELS
  MEASUREMENT
  PASS_REGULAR_EXPRESSION
  PROCESSORS
  REQUIRED_FILES
  RESOURCE_LOCK
  RUN_SERIAL
  TIMEOUT
  WILL_FAIL
  WORKING_DIRECTORY
)

# convert list into regular expression
basis_list_to_regex (BASIS_PROPERTIES_ON_TESTS_RE ${BASIS_PROPERTIES_ON_TESTS})

# ============================================================================
# utilities
# ============================================================================

## @addtogroup CMakeUtilities
#  @{


# ----------------------------------------------------------------------------
## @brief Disable testing if project does not implement any tests.
#
# This function checks if there are test/ subdirectories in the project and
# disables and hides the BUILD_TESTING option if none are found.
function (basis_disable_testing_if_no_tests)
  if (IS_DIRECTORY "${PROJECT_TESTING_DIR}")
    set (DISABLE_TESTING FALSE)
  else ()
    set (DISABLE_TESTING TRUE)
  endif ()
  if (DISABLE_TESTING)
    basis_get_relative_path (TESTING_DIR "${PROJECT_SOURCE_DIR}" "${PROJECT_TESTING_DIR}")
    foreach (M IN LISTS PROJECT_MODULES_ENABLED)
      if (IS_DIRECTORY "${MODULE_${M}_SOURCE_DIR}/${TESTING_DIR}")
        set (DISABLE_TESTING FALSE)
        break ()
      endif ()
    endforeach ()
  endif ()
  if (DISABLE_TESTING)
    set_property (CACHE BUILD_TESTING PROPERTY VALUE OFF)
    set_property (CACHE BUILD_TESTING PROPERTY TYPE  INTERNAL)
  else ()
    set_property (CACHE BUILD_TESTING PROPERTY TYPE BOOL)
  endif ()
endfunction ()

# ----------------------------------------------------------------------------
## @brief Set a property of the tests.
#
# This function replaces CMake's
# <a href="http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:set_tests_properties">
# set_tests_properties()</a> command.
#
# @note Due to a bug in CMake (http://www.cmake.org/Bug/view.php?id=12303),
#       except of the first property given directly after the @c PROPERTIES keyword,
#       only properties listed in @c BASIS_PROPERTIES_ON_TESTS can be set.
#
# @param [in] ARGN List of arguments for
#                  <a href="http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:set_tests_properties">
#                  set_tests_properties()</a>.
#
# @returns Sets the given properties of the specified test.
#
# @sa http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:set_tests_properties
#
# @ingroup CMakeAPI
function (basis_set_tests_properties)
  # convert test names to UIDs
  set (TEST_UIDS)
  list (GET ARGN 0 ARG)
  while (ARG AND NOT ARG MATCHES "^PROPERTIES$")
    basis_get_test_uid (TEST_UID "${ARG}")
    list (APPEND TEST_UIDS "${TEST_UID}")
    list (REMOVE_AT ARGN 0)
    list (GET ARGN 0 ARG)
  endwhile ()
  if (NOT ARG MATCHES "^PROPERTIES$")
    message (FATAL_ERROR "Missing PROPERTIES argument!")
  elseif (NOT TEST_UIDS)
    message (FATAL_ERROR "No test specified!")
  endif ()
  # remove PROPERTIES keyword
  list (REMOVE_AT ARGN 0)
  # set tests properties
  #
  # Note: By iterating over the properties, the empty property values
  #       are correctly passed on to CMake's set_tests_properties()
  #       command, while
  #       set_tests_properties(${TEST_UIDS} PROPERTIES ${ARGN})
  #       (erroneously) discards the empty elements in ARGN.
  if (BASIS_DEBUG)
    message ("** basis_set_tests_properties:")
    message ("**   Test(s)  :  ${TEST_UIDS}")
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
      if (ARG MATCHES "${BASIS_PROPERTIES_ON_TESTS_RE}")
        break ()
      endif ()
      list (APPEND VALUE "${ARG}")
      list (REMOVE_AT ARGN 0)
      list (LENGTH ARGN N)
    endwhile ()
    # check property name
    if ("${PROPERTY}" STREQUAL "")
      message (FATAL_ERROR "Empty property name given!")
    endif ()
    # map test names to test UIDs
    if (PROPERTY MATCHES "^DEPENDS$")
      set (UIDS)
      foreach (V IN LISTS VALUE)
        basis_get_test_uid (UID "${V}")
        list (APPEND UIDS "${UID}")
      endforeach ()
      set (VALUE ${UIDS})
    endif ()
    # set target property
    if (BASIS_DEBUG)
      message ("**   -> ${PROPERTY} = [${VALUE}]")
    endif ()
    set_tests_properties (${TEST_UIDS} PROPERTIES ${PROPERTY} "${VALUE}")
  endwhile ()
  # make sure that every property had a corresponding value
  if (NOT N EQUAL 0)
    message (FATAL_ERROR "No value given for test property ${ARGN}")
  endif ()
endfunction ()

# ----------------------------------------------------------------------------
## @brief Get a property of the test.
#
# This function replaces CMake's
# <a href="http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:get_test_property">
# get_test_property()</a> command.
#
# @param [out] VAR       Property value.
# @param [in]  TEST_NAME Name of test.
# @param [in]  ARGN      Remaining arguments of
#                        <a href="http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:get_test_property">
#                        get_test_property()</a>.
#
# @returns Sets @p VAR to the value of the requested property.
#
# @sa http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:get_test_property
#
# @ingroup CMakeAPI
function (basis_get_test_property VAR TEST_NAME)
  basis_get_test_uid (TEST_UID "${TEST_NAME}")
  get_test_property (VALUE "${TEST_UID}" ${ARGN})
  set (${VAR} "${VALUE}" PARENT_SCOPE)
endfunction ()

# ----------------------------------------------------------------------------
## @brief Create and add a test driver executable.
#
# @param [in] TESTDRIVER_NAME Name of the test driver.
# @param [in] ARGN            List of source files implementing tests.
#
# @ingroup CMakeAPI
function (basis_add_test_driver TESTDRIVER_NAME)
  # parse arguments
  CMAKE_PARSE_ARGUMENTS (
    ARGN
      ""
      "EXTRA_INCLUDE;FUNCTION"
      ""
    ${ARGN}
  )
  if (ARGN_EXTRA_INCLUDE OR ARGN_FUNCTION)
    message (FATAL_ERROR "Invalid argument EXTRA_INCLUDE or FUNCTION! "
                         "Use create_test_sourcelist() if you want to create "
                         "a custom test driver.")
  endif ()
  # choose test driver implementation depending on which packages are available
  set (TESTDRIVER_INCLUDE      "basis/testdriver.h")
  set (TESTDRIVER_LINK_DEPENDS basis.basis)
  if (ITK_FOUND)
    basis_include_directories (BEFORE ${ITK_INCLUDE_DIRS})
    list (APPEND TESTDRIVER_LINK_DEPENDS ${ITK_LIBRARIES})
  endif ()
  if (WIN32)
    list (APPEND TESTDRIVER_LINK_DEPENDS Ws2_32)
  endif ()
  # create test driver source code
  set (TESTDRIVER_SOURCE "${TESTDRIVER_NAME}")
  if (NOT TESTDRIVER_SOURCE MATCHES "\\.cxx")
    set (TESTDRIVER_SOURCE "${TESTDRIVER_SOURCE}.cxx")
  endif ()
  set (CMAKE_TESTDRIVER_BEFORE_TESTMAIN "    #include <basis/testdriver-before-test.inc>")
  set (CMAKE_TESTDRIVER_AFTER_TESTMAIN  "    #include <basis/testdriver-after-test.inc>")
  create_test_sourcelist (
    TESTDRIVER_SOURCES
      ${TESTDRIVER_SOURCE} ${ARGN_UNPARSED_ARGUMENTS}
      EXTRA_INCLUDE "${TESTDRIVER_INCLUDE}"
      FUNCTION      testdriversetup
  )
  # add test driver executable
  basis_add_executable (${TESTDRIVER_NAME} ${TESTDRIVER_SOURCES})
  basis_target_link_libraries (${TESTDRIVER_NAME} ${TESTDRIVER_LINK_DEPENDS})
  if (ITK_FOUND)
    basis_set_target_properties (
      ${TESTDRIVER_NAME}
      PROPERTIES
        COMPILE_DEFINITIONS
          TESTDRIVER_NAME=\"${TESTDRIVER_NAME}\"
          ITK_VERSION=\"${ITK_VERSION_MAJOR}.${ITK_VERSION_MINOR}.${ITK_VERSION_PATCH}\"
          ITK_VERSION_MAJOR=${ITK_VERSION_MAJOR}
    )
  else ()
    basis_set_target_properties (
      ${TESTDRIVER_NAME}
      PROPERTIES
        COMPILE_DEFINITIONS
          TESTDRIVER_NAME=\"${TESTDRIVER_NAME}\"
    )
  endif ()
endfunction ()

# ----------------------------------------------------------------------------
## @brief Add test.
#
# This command is used similar to CMake's
# <a href="http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:add_test">
# add_test()</a> command. It adds a test to the
# <a href="http://www.cmake.org/cmake/help/ctest-2-8-docs.html">CTest</a>-based
# testing system. Unlike CMake's
# <a href="http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:add_test">
# add_test()</a>, this command can, for convenience, implicitly add
# the necessary executable build target to the build system. Therefore,
# instead of the name of the executable command, specify the sources of the
# test implementation. An executable build target is then added by this
# function using basis_add_executable(), and the built executable is used
# as test command. If the @p UNITTEST option is given, the necessary unit
# testing libraries which are part of the BASIS installation are added as
# link dependencies as well as the default implementation of the main()
# function if none of the specified source files has the suffix
# <tt>-main</tt> or <tt>_main</tt> in the file name.
#
# Generator expressions as supported by CMake's
# <a href="http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:add_test">
# add_test()</a> command are also supported by basis_add_test() as arguments of
# the test command. For the argument specifying the test command itself, however,
# only the generator expression $<TARGET_FILE:tgt> is allowed. Alternatively,
# for this special argument, the name of the executable target can be supplied
# directly without the use of the $<TARGET_FILE:tgt> generator expression.
# See documentation of basis_process_generator_expressions() for details on
# the supported generator expressions.
#
# Example:
# @code
# basis_add_test (COMMAND $<TARGET_FILE:basis.testdriver> $<TARGET_FILE:myexe> ...)
# basis_add_test (COMMAND basis.testdriver $<TARGET_FILE:myexe> ...)
# @endcode
#
# @param [in] TEST_NAME Name of the test. If a source file is given
#                       as first argument, the test name is derived
#                       from the name of this source file and the source
#                       file is added to the list of sources which implement
#                       the test command.
# @param [in] ARGN      The following parameters are parsed:
# @par
# <table border="0">
#   <tr>
#     @tp @b COMMAND cmd [arg1 [arg2 ...]] @endtp
#     <td>The command to execute and optionally its arguments. The command
#         can be the name of an executable target (including imported targets),
#         or the name or path of an executable. Alternatively,
#         a test can be build from sources and the build executable
#         used as command. In this case, specify the sources using the
#         @p SOURCES argument. The command name @c cmd if given is used as
#         output name of the built executable. If you do not want to
#         specify the name of the output executable explicitly, but have
#         it derived from the @p TEST_NAME, do not specify the @p COMMAND
#         option and use the @p ARGS option instead to only specify the
#         arguments of the test command.</td>
#   </tr>
#   <tr>
#     @tp @b ARGS arg1 [arg2 ...] @endtp
#     <td>Arguments of the test command. If this option is given, the
#         specified arguments are appended to the arguments specified
#         already as part of the @p COMMAND option, if any.</td>
#   </tr>
#   <tr>
#     @tp @b WORKING_DIRECTORY dir @endtp
#     <td>The working directory of the test command. The generator expression
#         $<TARGET_FILE_DIR:tgt> can be used to specify a working directory
#         which corresponds to the output directory of a given target file.
#         Default: @c TESTING_OUTPUT_DIR / @p TEST_NAME. </td>
#   </tr>
#   <tr>
#     @tp @b CONFIGURATIONS @endtp
#     <td>If a CONFIGURATIONS option is given then the test will be executed
#         only when testing under one of the named configurations.</td>
#   </tr>
#   <tr>
#     @tp @b SOURCES file1 [file2 ...] @endtp
#     <td>The source files of the test. Use the @p UNITTEST option to specify
#         that the sources are an implementation of a unit test. In this case,
#         the default implementation of the main() function is added to the
#         build of the test executable. However, if this list contains a
#         file with the suffix <tt>-main</tt> or <tt>_main</tt> in the name,
#         the default implementation of the main() function is not used.
#         See the documentation of the @p UNITTEST option for further details.</td>
#   </tr>
#   <tr>
#     @tp @b LINK_DEPENDS file1|target1 [file2|target2 ...] @endtp
#     <td>Link dependencies of test executable build from sources.</td>
#   </tr>
#   <tr>
#     @tp @b NO_DEFAULT_MAIN @endtp
#     <td>Force that the implementation of the default main() function
#         is not added to unit tests even if neither of the given source
#         files has the suffix <tt>-main</tt> or <tt>_main</tt> in the
#         file name.</td>
#   </tr>
#   <tr>
#     @tp @b UNITTEST @endtp
#     <td>Specifies that the test is a unit test. In this case, the test
#         implementation is linked to the default unit testing framework
#         for the used programming language which is part of the BASIS
#         installation.</td>
#   </tr>
#   <tr>
#     @tp @b WITH_EXT @endtp
#     <td>Do not strip extension if test name is derived from source file name.</td>
#   </tr>
#   <tr>
#     @tp @b ARGN @endtp
#     <td>All other arguments are passed on to basis_add_executable() if
#         an executable target for the test is added.</td>
#   </tr>
# </table>
#
# @returns Adds build target for test executable if test source files
#          are given and/or adds a CTest test which executes the given
#          test command.
#
# @sa basis_process_generator_expressions()
# @sa http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:add_test
#
# @todo Make use of ExternalData module to fetch remote test data.
#
# @ingroup CMakeAPI
function (basis_add_test TEST_NAME)
  basis_sanitize_for_regex (R "${PROJECT_TESTING_DIR}")
  if (NOT CMAKE_CURRENT_SOURCE_DIR MATCHES "^${R}")
    message (FATAL_ERROR "Tests can only be added inside ${PROJECT_TESTING_DIR}!")
  endif ()

  # --------------------------------------------------------------------------
  # parse arguments
  CMAKE_PARSE_ARGUMENTS (
    ARGN
    "UNITTEST;NO_DEFAULT_MAIN;WITH_EXT;CLEAN_WORKING_DIRECTORY_BEFORE_TEST;CLEAN_WORKING_DIRECTORY_AFTER_TEST"
    "WORKING_DIRECTORY"
    "CONFIGURATIONS;SOURCES;LINK_DEPENDS;COMMAND;ARGS"
    ${ARGN}
  )

  if (BASIS_DEBUG AND BASIS_VERBOSE)
    message ("** basis_add_test():")
    message ("**   Test:      ${TEST_NAME}")
    message ("**   Command:   ${ARGN_COMMAND}")
    message ("**   Arguments: ${ARGN_ARGS}")
  endif ()

  list (LENGTH ARGN_COMMAND N)
  if (N GREATER 1)
    list (GET ARGN_COMMAND 0 CMD)
    list (REMOVE_AT ARGN_COMMAND 0)
    list (APPEND ARGN_ARGS ${ARGN_COMMAND})
    set (ARGN_COMMAND "${CMD}")
  endif ()

  # --------------------------------------------------------------------------
  # test name
  if (NOT ARGN_COMMAND AND NOT ARGN_SOURCES)
    get_filename_component (ARGN_SOURCES "${TEST_NAME}" ABSOLUTE)
    if (ARGN_WITH_EXT)
      basis_get_source_target_name (TEST_NAME "${TEST_NAME}" NAME)
      list (APPEND ARGN_UNPARSED_ARGUMENTS WITH_EXT) # pass on to basis_add_executable()
    else ()
      basis_get_source_target_name (TEST_NAME "${TEST_NAME}" NAME_WE)
    endif ()
  endif ()

  basis_check_test_name ("${TEST_NAME}")
  basis_make_test_uid (TEST_UID "${TEST_NAME}")

  if (BASIS_DEBUG)
    message ("** basis_add_test():")
    message ("**   Name:      ${TEST_NAME}")
    message ("**   Command:   ${ARGN_COMMAND}")
    message ("**   Arguments: ${ARGN_ARGS}")
    message ("**   Sources:   ${ARGN_SOURCES}")
  endif ()

  # --------------------------------------------------------------------------
  # build test executable
  set (LANGUAGE)
  if (ARGN_SOURCES)
    if (ARGN_UNITTEST)
      basis_get_source_language (LANGUAGE ${ARGN_SOURCES})
      if (NOT ARGN_NO_DEFAULT_MAIN)
        foreach (SOURCE ${ARGN_SOURCES})
          if (SOURCE MATCHES "-main\\.|_main\\.")
            set (ARGN_NO_DEFAULT_MAIN 1)
            break ()
          endif ()
        endforeach ()
      endif ()
      if (LANGUAGE MATCHES "CXX")
        if (NOT ARGN_NO_DEFAULT_MAIN)
          if (BASIS_TEST_MAIN_LIBRARY)
            list (APPEND ARGN_LINK_DEPENDS "${BASIS_TEST_MAIN_LIBRARY}")
          else ()
            message (FATAL_ERROR "Test ${TEST_NAME} added as UNITTEST without NO_DEFAULT_MAIN option,"
                                 " but BASIS_TEST_MAIN_LIBRARY not set. This library is part of the"
                                 " CMake BASIS installation. When you are using the CMake BASIS Modules"
                                 " without CMake BASIS installation, add GTest to TEST_DEPENDS and"
                                 " set BASIS_TEST_MAIN_LIBRARY to the value of GTEST_MAIN_LIBRARY"
                                 " before using the basis_add_test command.")
          endif ()
        endif ()
        if (BASIS_TEST_LIBRARY)
          list (APPEND ARGN_LINK_DEPENDS "${BASIS_TEST_LIBRARY}")
        else ()
         message (FATAL_ERROR "Test ${TEST_NAME} added as UNITTEST, but BASIS_TEST_LIBRARY not set."
                              " This library is part of the CMake BASIS installation. When you are"
                              " using the CMake BASIS Modules without CMake BASIS installation,"
                              " add GTest to TEST_DEPENDS and set BASIS_TEST_LIBRARY to the value"
                              " of GTEST_LIBRARY before using the basis_add_test command.")
        endif ()
      endif ()
    endif ()

    basis_add_executable (${TEST_NAME} ${ARGN_SOURCES} ${ARGN_UNPARSED_ARGUMENTS})
    if (ARGN_LINK_DEPENDS)
      basis_target_link_libraries (${TEST_NAME} ${ARGN_LINK_DEPENDS})
    endif ()

    if (ARGN_COMMAND)
      basis_set_target_properties (${TEST_NAME} PROPERTIES OUTPUT_NAME ${CMD})
    endif ()
    set (ARGN_COMMAND ${TEST_NAME})

    if (BASIS_DEBUG)
      message ("** Added test executable ${TEST_UID} (${TEST_NAME})")
      message ("**   Sources:           ${ARGN_SOURCES}")
      message ("**   Link dependencies: ${ARGN_LINK_DEPENDS}")
      message ("**   Options:           ${ARGN_UNPARSED_ARGUMENTS}")
    endif ()
  endif ()

  # --------------------------------------------------------------------------
  # add test
  message (STATUS "Adding test ${TEST_UID}...")

  if (NOT ARGN_WORKING_DIRECTORY)
    set (ARGN_WORKING_DIRECTORY "${TESTING_OUTPUT_DIR}/${TEST_NAME}")
  endif ()
  if (ARGN_WORKING_DIRECTORY MATCHES "^\\$<(.*):(.*)>$")
    if (NOT "${CMAKE_MATCH_1}" STREQUAL "TARGET_FILE_DIR")
      message (FATAL_ERROR "Invalid generator expression used for working directory."
                           " Only $<TARGET_FILE_DIR:tgt> can be used as argument"
                           " of the WORKING_DIRECTORY option.")
    endif ()
    basis_get_target_location (ARGN_WORKING_DIRECTORY "${CMAKE_MATCH_2}" PATH)
  else ()
    if (NOT IS_DIRECTORY "${ARGN_WORKING_DIRECTORY}")
      file (MAKE_DIRECTORY "${ARGN_WORKING_DIRECTORY}")
    endif ()
  endif ()
  set (OPTS "WORKING_DIRECTORY" "${ARGN_WORKING_DIRECTORY}")
  if (ARGN_CONFIGURATIONS)
    list (APPEND OPTS "CONFIGURATIONS")
    foreach (CONFIG ${ARGN_CONFIGURATIONS})
      list (APPEND OPTS "${CONFIG}")
    endforeach ()
  endif ()

  if (ARGN_COMMAND MATCHES "^\\$<(.*):(.*)>$")
    if (NOT "${CMAKE_MATCH_1}" STREQUAL "TARGET_FILE")
      message (FATAL_ERROR "Invalid generator expression used for test command."
                           " Only $<TARGET_FILE:tgt> can be used as first"
                           " argument of the COMMAND option.")
    endif ()
    basis_get_target_location (ARGN_COMMAND "${CMAKE_MATCH_2}" ABSOLUTE)
  else ()
    basis_get_target_uid (COMMAND_UID "${ARGN_COMMAND}")
    if (TARGET "${COMMAND_UID}")
      basis_get_target_location (ARGN_COMMAND "${COMMAND_UID}" ABSOLUTE)
    endif ()
  endif ()

  basis_process_generator_expressions (ARGN_ARGS ${ARGN_ARGS})

  if (BASIS_DEBUG)
    message ("** Add test ${TEST_UID}")
    message ("**   Command:    ${ARGN_COMMAND}")
    message ("**   Arguments:  ${ARGN_ARGS}")
    message ("**   Working in: ${ARGN_WORKING_DIRECTORY}")
  endif ()

  add_test (NAME ${TEST_UID} COMMAND ${ARGN_COMMAND} ${ARGN_ARGS} ${OPTS})

  # especially in case of C++ unit tests, if the linkage of the tests is done
  # incorrectly, no tests are actually run and the unit test passes
  # therefore, add this fail regular expression to identify such issues
  set_tests_properties (${TEST_UID} PROPERTIES FAIL_REGULAR_EXPRESSION "(\\[ *PASSED *\\]|Ran) 0 tests|No tests were found!!!")

  message (STATUS "Adding test ${TEST_UID}... - done")
endfunction ()

# ----------------------------------------------------------------------------
## @brief Add tests of default options for given executable.
#
# @par Default options:
# <table border="0">
#   <tr>
#     @tp <b>--helpshort</b> @endtp
#     <td>Short help. The output has to match the regular expression
#         "[Uu]sage:\n\s*\<executable name\>", where \<executable name\>
#         is the name of the tested executable.</td>
#   </tr>
#   <tr>
#     @tp <b>--help, -h</b> @endtp
#     <td>Help screen. Simply tests if the option is accepted.</td>
#   </tr>
#   <tr>
#     @tp <b>--version</b> @endtp
#     <td>Version information. Output has to include the project version string.</td>
#   </tr>
#   <tr>
#     @tp <b>--verbose, -v</b> @endtp
#     <td>Increase verbosity of output messages. Simply tests if the option is accepted.</td>
#   </tr>
# </table>
#
# @param [in] TARGET_NAME Name of executable or script target.
#
# @returns Adds tests for the default options of the specified executable.
function (basis_add_tests_of_default_options TARGET_NAME)
  basis_get_target_uid (TARGET_UID "${TARGET_NAME}")

  message (STATUS "Adding tests of default options for ${TARGET_UID}...")

  if (NOT TARGET "${TARGET_UID}")
    message (FATAL_ERROR "Unknown target ${TARGET_UID}.")
  endif ()

  # get executable name
  get_target_property (PREFIX      ${TARGET_UID} "PREFIX")
  get_target_property (OUTPUT_NAME ${TARGET_UID} "OUTPUT_NAME")
  get_target_property (SUFFIX      ${TARGET_UID} "SUFFIX")

  if (NOT OUTPUT_NAME)
    set (EXEC_NAME "${TARGET_UID}")
  endif ()
  if (PREFIX)
    set (EXEC_NAME "${PREFIX}${EXEC_NAME}")
  endif ()
  if (SUFFIX)
    set (EXEC_NAME "${EXEC_NAME}${SUFFIX}")
  endif ()

  # get absolute path to executable
  get_target_property (EXEC_DIR ${TARGET_UID} "RUNTIME_OUTPUT_DIRECTORY")

  # executable command
  set (EXEC_CMD "${EXEC_DIR}/${EXEC_NAME}")

  # test option: -v
  basis_add_test (${EXEC}VersionS "${EXEC_CMD}" "-v")

  set_tests_properties (
    ${EXEC}VersionS
    PROPERTIES
      PASS_REGULAR_EXPRESSION "${EXEC} ${PROJECT_VERSION}"
  )

  # test option: --version
  basis_add_test (${EXEC}VersionL "${EXEC_CMD}" "--version")

  set_tests_properties (
    ${EXEC}VersionL
    PROPERTIES
      PASS_REGULAR_EXPRESSION "${EXEC} ${PROJECT_VERSION}"
  )

  # test option: -h
  basis_add_test (${EXEC}HelpS "${EXEC_CMD}" "-h")

  # test option: --help
  basis_add_test (${EXEC}HelpL "${EXEC_CMD}" "--help")

  # test option: --helpshort
  basis_add_test (${EXEC}UsageL "${EXEC_CMD}" "--helpshort")

  set_tests_properties (
    ${EXEC}UsageL
    PROPERTIES
      PASS_REGULAR_EXPRESSION "[Uu]sage:(\n)( )*${EXEC_NAME}"
  )

  message (STATUS "Adding tests of default options for ${EXEC}... - done")
endfunction ()


## @}
# end of Doxygen group
