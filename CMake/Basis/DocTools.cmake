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
# @file  DocTools.cmake
# @brief Tools related to gnerating or adding software documentation.
#
# @ingroup CMakeTools
##############################################################################

if (__BASIS_DOCTOOLS_INCLUDED)
  return ()
else ()
  set (__BASIS_DOCTOOLS_INCLUDED TRUE)
endif ()


# ============================================================================
# adding / generating documentation
# ============================================================================

# ----------------------------------------------------------------------------
## @brief Add documentation target.
#
# This function is used to add a software documentation files to the project
# which are either just copied to the installation or generated from input
# files such as in particular source code files and documentation files
# marked up using one of the supported lightweight markup languages.
#
# The supported generators are:
# <table border="0">
#   <tr>
#     @tp @b None @endtp
#     <td>This generator simply installs the given file or all files within
#         the specified directory.</td>
#   </tr>
#   <tr>
#     @tp @b Doxygen @endtp
#     <td>Used to generate API documentation from in-source code comments and
#         other related files marked up using Doxygen comments. See
#         basis_add_doxygen_doc() for more details.</td>
#   </tr>
#   <tr>
#     @tp @b Sphinx @endtp
#     <td>Used to generate documentation such as a web site from reStructuredText.
#         See basis_add_sphinx_doc() for more details.</td>
#   </tr>
# </table>
#
# @param [in] TARGET_NAME Name of the documentation target or file.
# @param [in] ARGN        Documentation generator as "GENERATOR generator" option
#                         and additional arguments for the particular generator.
#                         The case of the generator name is ignored, i.e.,
#                         @c Doxygen, @c DOXYGEN, @c doxYgen are all valid arguments
#                         which select the @c Doxygen generator. The default generator
#                         is the @c None generator.</td>
#
# @returns Adds a custom target @p TARGET_NAME for the generation of the
#          documentation.
#
# @sa basis_install_doc()
# @sa basis_add_doxygen_doc()
# @sa basis_add_sphinx_doc()
#
# @ingroup CMakeAPI
function (basis_add_doc TARGET_NAME)
  CMAKE_PARSE_ARGUMENTS (ARGN "" "GENERATOR" "" ${ARGN})
  if (NOT ARGN_GENERATOR)
    set (ARGN_GENERATOR "NONE")
  else ()
    string (TOUPPER "${ARGN_GENERATOR}" ARGN_GENERATOR)
  endif ()
  if (ARGN_GENERATOR MATCHES "NONE")
    basis_install_doc (${TARGET_NAME} ${ARGN_UNPARSED_ARGUMENTS})
  elseif (ARGN_GENERATOR MATCHES "DOXYGEN")
    basis_add_doxygen_doc (${TARGET_NAME} ${ARGN_UNPARSED_ARGUMENTS})
  elseif (ARGN_GENERATOR MATCHES "SPHINX")
    basis_add_sphinx_doc (${TARGET_NAME} ${ARGN_UNPARSED_ARGUMENTS})
  else ()
    message (FATAL_ERROR "Unknown documentation generator: ${ARGN_GENERATOR}.")
  endif ()
endfunction ()

# ----------------------------------------------------------------------------
## @brief Install documentation file(s).
#
# This function either adds an installation rule for a single documentation
# file or a directory containing multiple documentation files.
#
# Example:
# @code
# basis_install_doc ("User Manual.pdf" OUTPUT_NAME "BASIS User Manual.pdf")
# basis_install_doc (DeveloperManual.docx COMPONENT dev)
# basis_install_doc (SourceManual.html    COMPONENT src)
# @endcode
#
# @param [in] SOURCE Documentation file or directory to install.
# @param [in] ARGN   List of optional arguments. Valid arguments are:
# @par
# <table border="0">
#   <tr>
#     @tp @b COMPONENT component @endtp
#     <td>Name of the component this documentation belongs to.
#         Defaults to @c BASIS_RUNTIME_COMPONENT.</td>
#   </tr>
#   <tr>
#     @tp @b DESTINATION dir @endtp
#     <td>Installation directory prefix. Defaults to @c INSTALL_DOC_DIR.</td>
#   </tr>
#   <tr>
#     @tp @b OUTPUT_NAME name @endtp
#     <td>Name of file or directory after installation.</td>
#   </tr>
# </table>
#
# @sa basis_add_doc()
function (basis_install_doc SOURCE)
  CMAKE_PARSE_ARGUMENTS (ARGN "" "COMPONENT;DESTINATION;OUTPUT_NAME" "" ${ARGN})

  if (NOT IS_ABSOLUTE "${SOURCE}")
    get_filename_component (SOURCE "${SOURCE}" ABSOLUTE)
  endif ()
  if (NOT ARGN_DESTINATION)
    set (ARGN_DESTINATION "${INSTALL_DOC_DIR}")
  endif ()
  if (NOT ARGN_COMPONENT)
    set (ARGN_COMPONENT "${BASIS_RUNTIME_COMPONENT}")
  endif ()
  if (NOT ARGN_COMPONENT)
    set (ARGN_COMPONENT "Unspecified")
  endif ()
  if (NOT ARGN_OUTPUT_NAME)
    basis_get_filename_component (ARGN_OUTPUT_NAME "${SOURCE}" NAME)
  endif ()

  basis_get_relative_path (
    RELPATH
      "${CMAKE_SOURCE_DIR}"
      "${CMAKE_CURRENT_SOURCE_DIR}/${ARGN_OUTPUT_NAME}"
  )

  message (STATUS "Adding documentation ${RELPATH}...")

  if (IS_DIRECTORY "${SOURCE}")
    basis_install_directory (
      "${SOURCE}" "${ARGN_DESTINATION}/${ARGN_OUTPUT_NAME}"
      COMPONENT "${ARGN_COMPONENT}"
    )
  else ()
    install (
      FILES       "${SOURCE}"
      DESTINATION "${ARGN_DESTINATION}"
      COMPONENT   "${ARGN_COMPONENT}"
      RENAME      "${ARGN_OUTPUT_NAME}"
    )
  endif ()

  message (STATUS "Adding documentation ${RELPATH}... - done")
endfunction ()

# ----------------------------------------------------------------------------
## @brief Add documentation to be generated by Doxygen.
#
# This function adds a build target to generate documentation from in-source
# code comments and other related project pages using
# <a href="http://www.stack.nl/~dimitri/doxygen/index.html">Doxygen</a>.
#
# @param [in] TARGET_NAME Name of the documentation target.
# @param [in] ARGN        List of arguments. The valid arguments are:
# @par
# <table border="0">
#   <tr>
#     @tp @b EXCLUDE_FROM_DOC @endtp
#     <td>By default, the specified target is build as part of the global
#         @c doc target. If this option is given, however, the added
#         documentation will not be build as part of this target.</td>
#   </tr>
#   <tr>
#     @tp @b COMPONENT component @endtp
#     <td>Name of the component this documentation belongs to.
#         Defaults to @c BASIS_LIBRARY_COMPONENT.</td>
#   </tr>
#   <tr>
#     @tp @b DESTINATION dir @endtp
#     <td>Installation directory prefix. Defaults to
#         @c BASIS_INSTALL_&ltTARGET&gt;_DIR in case of HTML output if set.
#         Otherwise, the generated HTML files are not installed.</td>
#   </tr>
#   <tr>
#     @tp @b DOXYFILE file @endtp
#     <td>Name of the template Doxyfile.</td>
#   </tr>
#   <tr>
#     @tp @b PROJECT_NAME name @endtp
#     <td>Value for Doxygen's @c PROJECT_NAME tag which is used to
#         specify the project name.@n
#         Default: @c PROJECT_NAME.</td>
#   </tr>
#   <tr>
#     @tp @b PROJECT_NUMBER version @endtp
#     <td>Value for Doxygen's @c PROJECT_NUMBER tag which is used
#         to specify the project version number.@n
#         Default: @c PROJECT_RELEASE.</td>
#   </tr>
#   <tr>
#     @tp @b PROJECT_WEBSITE url @endtp
#     <td>Used for links to project website.@n
#         Default: @c PROJECT_PACKAGE_WEBSITE </td>
#   </tr>
#   <tr>
#     @tp @b INPUT path1 [path2 ...] @endtp
#     <td>Value for Doxygen's @c INPUT tag which is used to specify input
#         directories/files. Any given input path is added to the default
#         input paths.@n
#         Default: @c PROJECT_CODE_DIRS, @c BINARY_CODE_DIR,
#                  @c PROJECT_INCLUDE_DIRS, @c BINARY_INCLUDE_DIR.</td>
#   </tr>
#   <tr>
#     @tp @b EXCLUDE_BASIS_MODULES @endtp
#     <td>Do not add project CMake files used and generated by BASIS to @b INPUT.</td>
#   </tr>
#   <tr>
#     @tp @b EXCLUDE_BASIS_UTILITIES @endtp
#     <td>Do not add documentation (.dox) files for used CMake BASIS Utilities to @b INPUT.</td>
#   </tr>
#   <tr>
#     @tp @b INPUT_FILTER filter @endtp
#     <td>
#       Value for Doxygen's @c INPUT_FILTER tag which can be used to
#       specify a default filter for all input files.@n
#       Default: @c doxyfilter of BASIS.
#     </td>
#   </tr>
#   <tr>
#     @tp @b FILTER_PATTERNS pattern1 [pattern2...] @endtp
#     <td>Value for Doxygen's @c FILTER_PATTERNS tag which can be used to
#         specify filters on a per file pattern basis.@n
#         Default: None.</td>
#   </tr>
#   <tr>
#     @tp @b INCLUDE_PATH path1 [path2...] @endtp
#     <td>Doxygen's @c INCLUDE_PATH tag can be used to specify one or more
#         directories that contain include files that are not input files
#         but should be processed by the preprocessor. Any given directories
#         are appended to the default include path considered.
#         Default: Directories added by basis_include_directories().</td>
#   </tr>
#   <tr>
#     @tp @b EXCLUDE_PATTERNS pattern1 [pattern2 ...] @endtp
#     <td>Additional patterns used for Doxygen's @c EXCLUDE_PATTERNS tag
#         which can be used to specify files and/or directories that
#         should be excluded from the INPUT source files.@n
#         Default: No exclude patterns.</td>
#   </tr>
#   <tr>
#     @tp @b PREDEFINED name1|name1=value1 [name2|name2=value2...] @endtp
#     <td>Add preprocessor definitions to be expanded by Doxygen.</td>
#   </tr>
#   <tr>
#     @tp @b OUTPUT fmt @endtp
#     <td>Specify output formats in which to generate the documentation.
#         Currently, only @c html and @c xml are supported.</td>
#   </tr>
#   <tr>
#     @tp @b OUTPUT_DIRECTORY dir @endtp
#     <td>Value for Doxygen's @c OUTPUT_DIRECTORY tag which can be used to
#         specify the output directory. The output files are written to
#         subdirectories named "html", "latex", "rtf", and "man".@n
#         Default: <tt>CMAKE_CURRENT_BINARY_DIR/TARGET_NAME</tt>.</td>
#   </tr>
#   <tr>
#     @tp @b COLS_IN_ALPHA_INDEX n @endtp
#     <td>Number of columns in alphabetical index. Default: 3.</td>
#   </tr>
#   <tr>
#     @tp @b IGNORE_PREFIX prefix1 [prefix2...] @endtp
#     <td>In case all classes in a project start with a common prefix, all 
#         classes will be put under the same header in the alphabetical index. 
#         The IGNORE_PREFIX tag can be used to specify one or more prefixes that 
#         should be ignored while generating the index headers.</td>
#   </tr>
#   <tr>
#     @tp @b PROVIDER_NAME name @endtp
#     <td>Value for provider name, such as a company name,  
#         that will be used for pages in the doxygen output.@n
#         Default: @c PROJECT_PROVIDER_NAME.</td>
#   </tr>
#   <tr>
#     @tp @b PROVIDER_WEBSITE url @endtp
#     <td>Value for provider website, such as a company website,  
#         that will be used for pages in the doxygen output.@n
#         Default: @c PROJECT_PROVIDER_WEBSITE.</td>
#   </tr>
#   <tr>
#     @tp @b PROVIDER_LOGO image_file @endtp
#     <td>Value for provider logo file, such as a company logo,  
#         that will be used for pages in the doxygen output.@n
#         Default: @c PROJECT_PROVIDER_LOGO.</td>
#   </tr>
#   <tr>
#     @tp @b DIVISION_NAME name @endtp
#     <td>Value for division name, such as a company division name,  
#         that will be used for pages in the doxygen output.@n
#         Default: @c PROJECT_DIVISION_NAME.</td>
#   </tr>
#   <tr>
#     @tp @b DIVISION_WEBSITE url @endtp
#     <td>Value for division website, such as a company division website,  
#         that will be used for pages in the doxygen output.@n
#         Default: @c PROJECT_DIVISION_WEBSITE.</td>
#   </tr>
#   <tr>
#     @tp @b DIVISION_LOGO image_file @endtp
#     <td>Value for division logo file, such as a company division logo,  
#         that will be used for pages in the doxygen output.@n
#         Default: @c PROJECT_DIVISION_LOGO.</td>
#   </tr>
#   <tr>
#     @tp @b ENABLED_SECTIONS section1 [section2 ...] @endtp
#     <td>The ENABLED_SECTIONS tag can be used to enable conditional  
#         documentation sections, marked by "\if sectionname ... \endif".</td>
#   </tr>
#   <tr>
#     @tp @b HTML_HEADER html_file @endtp
#     <td>The HTML_HEADER tag can be used to specify a personal HTML header for 
#         each generated HTML page. If none specified, the
#         @c "PROJECT_SOURCE_DIR/doc/doxygen_header.html(.in)?" file is used if present.
#         Otherwise, a default header is used. Specify the value @c Doxygen to use the
#         standard header generated by Doxygen instead.</td>
#   </tr>
#   <tr>
#     @tp @b HTML_FOOTER html_file @endtp
#     <td>The HTML_FOOTER tag can be used to specify a personal HTML footer for 
#         each generated HTML page. If none specified, the
#         @c "PROJECT_SOURCE_DIR/doc/doxygen_footer.html(.in)?" file is used if present.
#         Otherwise, a default footer is used. Specify the value @c Doxygen to use the
#         standard footer generated by Doxygen instead.</td>
#   </tr>
#   <tr>
#     @tp @b HTML_EXTRA_STYLESHEET css_file @endtp
#     <td>The HTML_EXTRA_STYLESHEET tag can be used to specify a user-defined cascading 
#         style sheet that is used by each HTML page. It can be used to 
#         fine-tune the look of the HTML output. If none specified, the
#         @c "PROJECT_SOURCE_DIR/doc/doxygen_extra.css(.in)?" file is used if present.</td>
#   </tr>
#   <tr>
#     @tp @b HTML_EXTRA_FILES file1 [file2...] @endtp
#     <td>The HTML_EXTRA_FILES tag can be used to specify additional files needed
#         for the HTML output of the API documentation.</td>
#   </tr>
#   <tr>
#     @tp @b DISABLE_PROJECT_NAME_DISPLAY@endtp
#     <td>The DISABLE_PROJECT_NAME_DISPLAY option causes Doxygen's 
#         @c PROJECT_NAME text not to be displayed in the header.
#         Use this if the project name is already part of the logo 
#         so it won't be there twice in the logo image and title text.</td>
#   </tr>
# </table>
# @n
# See <a href="http://www.stack.nl/~dimitri/doxygen/manual/config.html">here</a> for a
# documentation of the Doxygen tags.
# @n@n
# Example:
# @code
# basis_add_doxygen_doc (
#   apidoc
#   DOXYFILE        "Doxyfile.in"
#   PROJECT_NAME    "${PROJECT_NAME}"
#   PROJECT_VERSION "${PROJECT_VERSION}"
#   COMPONENT       dev
# )
# @endcode
#
# @sa basis_add_doc()
function (basis_add_doxygen_doc TARGET_NAME)
  # check target name
  basis_check_target_name ("${TARGET_NAME}")
  basis_make_target_uid (TARGET_UID "${TARGET_NAME}")
  string (TOLOWER "${TARGET_NAME}" TARGET_NAME_L)
  string (TOUPPER "${TARGET_NAME}" TARGET_NAME_U)
  # verbose output
  message (STATUS "Adding documentation ${TARGET_UID}...")
  # find Doxygen
  find_package (Doxygen QUIET)
  if (NOT DOXYGEN_EXECUTABLE)
    if (BUILD_DOCUMENTATION)
      message (FATAL_ERROR "Doxygen not found! Either install Doxygen and/or set DOXYGEN_EXECUTABLE or disable BUILD_DOCUMENTATION.")
    endif ()
    message (STATUS "Doxygen not found. Generation of ${TARGET_UID} documentation disabled.")
    message (STATUS "Adding documentation ${TARGET_UID}... - skipped")
    return ()
  endif ()
  # parse arguments
  set (VALUEARGS
    PROJECT_NAME
    PROJECT_NUMBER
    PROJECT_WEBSITE
    PROVIDER_NAME
    PROVIDER_WEBSITE
    DIVISION_NAME
    DIVISION_WEBSITE
    COMPONENT
    DESTINATION
    HTML_DESTINATION
    MAN_DESTINATION
    OUTPUT_DIRECTORY
    COLS_IN_ALPHA_INDEX
    MAN_SECTION
  )
  set (OPTIONAL_FILE_OPTIONS
    HTML_FOOTER
    HTML_HEADER
    HTML_EXTRA_STYLESHEET
    PROJECT_LOGO
    PROVIDER_LOGO
    DIVISION_LOGO
    DOXYFILE
    TAGFILE
  )
  CMAKE_PARSE_ARGUMENTS (
    DOXYGEN
      "EXCLUDE_FROM_DOC;DISABLE_PROJECT_NAME_DISPLAY;EXCLUDE_BASIS_MODULES;EXCLUDE_BASIS_UTILITIES"
      "${VALUEARGS};${OPTIONAL_FILE_OPTIONS}"
      "INPUT;OUTPUT;INPUT_FILTER;FILTER_PATTERNS;EXCLUDE_PATTERNS;INCLUDE_PATH;IGNORE_PREFIX;ENABLED_SECTIONS;PREDEFINED;HTML_EXTRA_FILES"
      ${ARGN}
  )
  unset (VALUEARGS)
  # handle special arguments
  set (DOXYGEN_HTML_HEADER_IS_DEFAULT FALSE)
  if (DOXYGEN_HTML_HEADER MATCHES "^(Doxygen|doxygen|DOXYGEN|none|None|NONE)$")
    set (DOXYGEN_HTML_HEADER)
  elseif (NOT DOXYGEN_HTML_HEADER OR DOXYGEN_HTML_HEADER MATCHES "^(Default|default|DEFAULT)$")
    set (DOXYGEN_HTML_HEADER "${BASIS_MODULE_PATH}/doxygen_header.html.in")
    set (DOXYGEN_HTML_HEADER_IS_DEFAULT TRUE)
  endif ()
  if (DOXYGEN_HTML_FOOTER MATCHES "^(Doxygen|doxygen|DOXYGEN|none|None|NONE)$")
    set (DOXYGEN_HTML_FOOTER)
  elseif (NOT DOXYGEN_HTML_FOOTER OR DOXYGEN_HTML_FOOTER MATCHES "^(Default|default|DEFAULT)$")
    set (DOXYGEN_HTML_FOOTER "${BASIS_MODULE_PATH}/doxygen_footer.html.in")
  endif ()
  # make file paths absolute and check if files exist
  foreach (opt IN LISTS OPTIONAL_FILE_OPTIONS)
    if (DOXYGEN_${opt})
      get_filename_component (DOXYGEN_${opt} "${DOXYGEN_${opt}}" ABSOLUTE)
      if (NOT EXISTS "${DOXYGEN_${opt}}")
        message (FATAL_ERROR "File ${DOXYGEN_${opt}} does not exist. Check value of the ${opt} option and make sure the file is present.")
      endif ()
    endif ()
  endforeach ()
  set (TMP_DOXYGEN_HTML_EXTRA_FILES)
  foreach (path IN LISTS DOXYGEN_HTML_EXTRA_FILES)
    get_filename_component (abspath "${path}" ABSOLUTE)
    if (NOT EXISTS "${path}")
      message (FATAL_ERROR "File ${path} does not exist. Check value of the HTML_EXTRA_FILES option and make sure the file is present.")
    endif ()
    list (APPEND TMP_DOXYGEN_HTML_EXTRA_FILES "${path}")
  endforeach ()
  set (DOXYGEN_HTML_EXTRA_FILES "${TMP_DOXYGEN_HTML_EXTRA_FILES}")
  unset (TMP_DOXYGEN_HTML_EXTRA_FILES)
  # default component
  if (NOT DOXYGEN_COMPONENT)
    set (DOXYGEN_COMPONENT "${BASIS_LIBRARY_COMPONENT}")
  endif ()
  if (NOT DOXYGEN_COMPONENT)
    set (DOXYGEN_COMPONENT "Unspecified")
  endif ()
  # configuration file
  if (NOT DOXYGEN_DOXYFILE)
    set (DOXYGEN_DOXYFILE "${BASIS_DOXYGEN_DOXYFILE}")
  endif ()
  if (NOT EXISTS "${DOXYGEN_DOXYFILE}")
    message (FATAL_ERROR "Missing option DOXYGEN_FILE or Doxyfile ${DOXYGEN_DOXYFILE} does not exist.")
  endif ()
  # default project attributes and logos
  if (NOT DOXYGEN_PROJECT_NAME)
    set (DOXYGEN_PROJECT_NAME "${PROJECT_NAME}")
  endif ()
  if (NOT DOXYGEN_PROJECT_NUMBER)
    set (DOXYGEN_PROJECT_NUMBER "${PROJECT_RELEASE}")
  endif ()
  if (NOT DOXYGEN_PROJECT_LOGO)
    set (DOXYGEN_PROJECT_LOGO "${PROJECT_PACKAGE_LOGO}")
  elseif (DOXYGEN_PROJECT_LOGO MATCHES "^None|none|NONE$")
    set (DOXYGEN_PROJECT_LOGO)
  endif ()
  if (NOT DOXYGEN_PROJECT_WEBSITE)
    set (DOXYGEN_PROJECT_WEBSITE "${PROJECT_PACKAGE_WEBSITE}")
  endif ()
  if (NOT DOXYGEN_PROVIDER_NAME)
    set (DOXYGEN_PROVIDER_NAME "${PROJECT_PROVIDER_NAME}")
  endif ()
  if (NOT DOXYGEN_PROVIDER_WEBSITE)
    set (DOXYGEN_PROVIDER_WEBSITE "${PROJECT_PROVIDER_WEBSITE}")
  endif ()
  if (NOT DOXYGEN_PROVIDER_LOGO)
    set (DOXYGEN_PROVIDER_LOGO "${PROJECT_PROVIDER_LOGO}")
  endif ()
  if (NOT DOXYGEN_DIVISION_NAME)
    set (DOXYGEN_DIVISION_NAME "${PROJECT_DIVISION_NAME}")
  endif ()
  if (NOT DOXYGEN_DIVISION_WEBSITE)
    set (DOXYGEN_DIVISION_WEBSITE "${PROJECT_DIVISION_WEBSITE}")
  endif ()
  if (NOT DOXYGEN_DIVISION_LOGO)
    set (DOXYGEN_DIVISION_LOGO "${PROJECT_DIVISION_LOGO}")
  endif ()
  # set visibility property of project logos
  if (DOXYGEN_PROJECT_LOGO)
    set (DOXYGEN_PROJECT_LOGO_DISPLAY "block")
  else ()
    set (DOXYGEN_PROJECT_LOGO_DISPLAY "none")
  endif ()
  if (DOXYGEN_PROVIDER_LOGO)
    set (DOXYGEN_PROVIDER_LOGO_DISPLAY "inline")
  else ()
    set (DOXYGEN_PROVIDER_LOGO_DISPLAY "block")
  endif ()
  if (DOXYGEN_DIVISION_LOGO)
    set (DOXYGEN_DIVISION_LOGO_DISPLAY "inline")
  else ()
    set (DOXYGEN_DIVISION_LOGO_DISPLAY "block")
  endif ()
  # allow the user to disable the text header if desired
  if(DOXYGEN_DISABLE_PROJECT_NAME_DISPLAY)
    set (DOXYGEN_PROJECT_NAME_DISPLAY "none")
  else()
    set (DOXYGEN_PROJECT_NAME_DISPLAY "inline")
  endif()
  # standard input files
  if (NOT EXCLUDE_BASIS_MODULES)
    list (APPEND DOXYGEN_INPUT "${PROJECT_SOURCE_DIR}/BasisProject.cmake")
    if (EXISTS "${PROJECT_CONFIG_DIR}/Depends.cmake")
      list (APPEND DOXYGEN_INPUT "${PROJECT_CONFIG_DIR}/Depends.cmake")
    endif ()
    if (EXISTS "${BINARY_CONFIG_DIR}/Directories.cmake")
      list (APPEND DOXYGEN_INPUT "${BINARY_CONFIG_DIR}/Directories.cmake")
    endif ()
    if (EXISTS "${BINARY_CONFIG_DIR}/BasisSettings.cmake")
      list (APPEND DOXYGEN_INPUT "${BINARY_CONFIG_DIR}/BasisSettings.cmake")
    endif ()
    if (EXISTS "${BINARY_CONFIG_DIR}/ProjectSettings.cmake")
      list (APPEND DOXYGEN_INPUT "${BINARY_CONFIG_DIR}/ProjectSettings.cmake")
    endif ()
    if (EXISTS "${BINARY_CONFIG_DIR}/Settings.cmake")
      list (APPEND DOXYGEN_INPUT "${BINARY_CONFIG_DIR}/Settings.cmake")
    elseif (EXISTS "${PROJECT_CONFIG_DIR}/Settings.cmake")
      list (APPEND DOXYGEN_INPUT "${PROJECT_CONFIG_DIR}/Settings.cmake")
    endif ()
    if (EXISTS "${BASIS_SCRIPT_CONFIG_FILE}")
      list (APPEND DOXYGEN_INPUT "${BASIS_SCRIPT_CONFIG_FILE}")
    endif ()
    if (EXISTS "${BINARY_CONFIG_DIR}/ScriptConfig.cmake")
      list (APPEND DOXYGEN_INPUT "${BINARY_CONFIG_DIR}/ScriptConfig.cmake")
    endif ()
    if (EXISTS "${PROJECT_CONFIG_DIR}/ConfigSettings.cmake")
      list (APPEND DOXYGEN_INPUT "${PROJECT_CONFIG_DIR}/ConfigSettings.cmake")
    endif ()
    if (EXISTS "${PROJECT_SOURCE_DIR}/CTestConfig.cmake")
      list (APPEND DOXYGEN_INPUT "${PROJECT_SOURCE_DIR}/CTestConfig.cmake")
    endif ()
    if (EXISTS "${PROJECT_BINARY_DIR}/CTestCustom.cmake")
      list (APPEND DOXYGEN_INPUT "${PROJECT_BINARY_DIR}/CTestCustom.cmake")
    endif ()
    # package configuration files - only exist *after* this function executed
    list (APPEND DOXYGEN_INPUT "${BINARY_LIBCONF_DIR}/${PROJECT_PACKAGE_CONFIG_PREFIX}Config.cmake")
    list (APPEND DOXYGEN_INPUT "${BINARY_LIBCONF_DIR}/${PROJECT_PACKAGE_CONFIG_PREFIX}ConfigVersion.cmake")
    list (APPEND DOXYGEN_INPUT "${BINARY_LIBCONF_DIR}/${PROJECT_PACKAGE_CONFIG_PREFIX}Use.cmake")
    # add .dox files with definition of BASIS Modules groups
    if (BASIS_DIR)
      list (APPEND DOXYGEN_INPUT "${BASIS_MODULE_PATH}/Modules.dox")
    endif ()
  endif ()
  # input directories
  foreach (_DIR IN LISTS BINARY_INCLUDE_DIR PROJECT_INCLUDE_DIRS BINARY_CODE_DIR PROJECT_CODE_DIRS)
    if (IS_DIRECTORY ${_DIR})
      list (APPEND DOXYGEN_INPUT "${_DIR}")
    endif ()
  endforeach ()
  foreach (M IN LISTS PROJECT_MODULES_ENABLED)
    foreach (_DIR IN LISTS ${M}_INCLUDE_DIRS ${M}_CODE_DIRS)
      if (IS_DIRECTORY ${_DIR})
        list (APPEND DOXYGEN_INPUT "${_DIR}")
      endif ()
    endforeach ()
  endforeach ()
  # in case of scripts, have Doxygen process the configured versions for the
  # installation which are further located in proper subdirectories instead
  # of the original source files
  basis_get_project_property (TARGETS)
  foreach (T IN LISTS TARGETS)
    get_target_property (BASIS_TYPE ${T} BASIS_TYPE)
    get_target_property (IS_TEST    ${T} TEST)
    if (NOT IS_TEST AND BASIS_TYPE MATCHES "SCRIPT")
      get_target_property (SOURCES ${T} SOURCES)
      if (SOURCES)
        list (GET SOURCES 0 BUILD_DIR) # CMake <3.1 stores path to internal build directory here
        if (BUILD_DIR MATCHES "CMakeFiles")
          list (REMOVE_AT SOURCES 0)
        endif ()
        get_target_property (BUILD_DIR ${T} BUILD_DIRECTORY)
        list (APPEND DOXYGEN_INPUT "${BUILD_DIR}.dir/install")
        foreach (S IN LISTS SOURCES)
          list (APPEND DOXYGEN_EXCLUDE_PATTERNS "${S}")
          list (APPEND DOXYGEN_EXCLUDE_PATTERNS "${BUILD_DIR}.dir/build")
        endforeach ()
      endif ()
    endif ()
  endforeach ()
  # add .dox files as input
  file (GLOB_RECURSE DOX_FILES "${PROJECT_DOC_DIR}/*.dox")
  list (SORT DOX_FILES) # alphabetic order
  list (APPEND DOXYGEN_INPUT ${DOX_FILES})
  # add .dox files of used BASIS utilities
  if (BASIS_DIR AND NOT EXCLUDE_BASIS_UTILITIES)
    list (APPEND DOXYGEN_INPUT "${BASIS_MODULE_PATH}/Utilities.dox")
    list (APPEND DOXYGEN_INPUT "${BASIS_MODULE_PATH}/CxxUtilities.dox")
    foreach (L IN ITEMS Cxx Java Python Perl Bash Matlab)
      string (TOUPPER "${L}" U)
      basis_get_project_property (USES_${U}_UTILITIES PROPERTY PROJECT_USES_${U}_UTILITIES)
      if (USES_${U}_UTILITIES)
        list (FIND DOXYGEN_INPUT "${BASIS_MODULE_PATH}/Utilities.dox" IDX)
        if (IDX EQUAL -1)
          list (APPEND DOXYGEN_INPUT "${BASIS_MODULE_PATH}/Utilities.dox")
        endif ()
        list (APPEND DOXYGEN_INPUT "${BASIS_MODULE_PATH}/${L}Utilities.dox")
      endif ()
    endforeach ()
  endif ()
  # include path - Disabled as this increases the runtime of Doxygen but
  #                generally the source of third-party packages are not
  #                really referenced. Only the source files of this
  #                project have to be considered. This code is kept as it
  #                might be used again at a later point once it is figured
  #                how Doxygen can be only rerun if necessary.
  #basis_get_project_property (INCLUDE_DIRS PROPERTY PROJECT_INCLUDE_DIRS)
  #foreach (D IN LISTS INCLUDE_DIRS)
  #  list (FIND DOXYGEN_INPUT "${D}" IDX)
  #  if (IDX EQUAL -1)
  #    list (APPEND DOXYGEN_INCLUDE_PATH "${D}")
  #  endif ()
  #endforeach ()
  #basis_list_to_delimited_string (
  #  DOXYGEN_INCLUDE_PATH "\"\nINCLUDE_PATH          += \"" ${DOXYGEN_INCLUDE_PATH}
  #)
  #set (DOXYGEN_INCLUDE_PATH "\"${DOXYGEN_INCLUDE_PATH}\"")
  # make string from DOXYGEN_INPUT - after include path was set
  basis_list_to_delimited_string (
    DOXYGEN_INPUT "\"\nINPUT                 += \"" ${DOXYGEN_INPUT}
  )
  set (DOXYGEN_INPUT "\"${DOXYGEN_INPUT}\"")
  # preprocessor definitions
  basis_list_to_delimited_string (
    DOXYGEN_PREDEFINED "\"\nPREDEFINED            += \"" ${DOXYGEN_PREDEFINED}
  )
  set (DOXYGEN_PREDEFINED "\"${DOXYGEN_PREDEFINED}\"")
  # outputs
  if (NOT DOXYGEN_OUTPUT_DIRECTORY)
    set (DOXYGEN_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${TARGET_NAME_L}")
  endif ()
  if (DOXYGEN_TAGFILE MATCHES "^(None|NONE|none)$")
    set (DOXYGEN_TAGFILE)
  else ()
    set (DOXYGEN_TAGFILE "${DOXYGEN_OUTPUT_DIRECTORY}/Doxytags.${TARGET_NAME_L}")
  endif ()
  if (NOT DOXYGEN_OUTPUT)
    set (DOXYGEN_OUTPUT html)
  endif ()
  foreach (F IN ITEMS HTML XML RTF LATEX MAN)
    set (DOXYGEN_GENERATE_${F} NO)
  endforeach ()
  foreach (f IN LISTS DOXYGEN_OUTPUT)
    if (NOT f MATCHES "^(html|xml)$")
      message (FATAL_ERROR "Invalid/Unsupported Doxygen output format: ${f}")
    endif ()
    string (TOUPPER "${f}" F)
    set (DOXYGEN_GENERATE_${F} YES)  # enable generation of this output
    set (DOXYGEN_${F}_OUTPUT "${f}") # relative output directory
  endforeach ()
  # input filters
  if (NOT DOXYGEN_INPUT_FILTER)
    basis_get_target_uid (DOXYFILTER "basis.doxyfilter")
    if (TARGET "${DOXYFILTER}")
      basis_get_target_location (DOXYGEN_INPUT_FILTER ${DOXYFILTER} ABSOLUTE)
    else ()
      basis_get_target_uid (DOXYFILTER "doxyfilter")
      if (TARGET "${DOXYFILTER}")
        basis_get_target_location (DOXYGEN_INPUT_FILTER ${DOXYFILTER} ABSOLUTE)
      endif ()
    endif ()
  else ()
    set (DOXYFILTER)
  endif ()
  if (DOXYGEN_INPUT_FILTER)
    if (WIN32)
      # Doxygen on Windows (XP, 32-bit) (at least up to version 1.8.0) seems
      # to have a problem of not calling filters which have a space character
      # in their file path correctly. The doxyfilter.bat Batch program is used
      # as a wrapper for the actual filter which is part of the BASIS build.
      # As this file is in the working directory of Doxygen, it can be
      # referenced relative to this working directory, i.e., without file paths.
      # The Batch program itself then calls the actual Doxygen filter with proper
      # quotes to ensure that spaces in the file path are handled correctly.
      # The file extension .bat shall distinguish this wrapper script from the actual
      # doxyfilter.cmd which is generated by BASIS on Windows.
      configure_file ("${BASIS_MODULE_PATH}/doxyfilter.bat.in" "${DOXYGEN_OUTPUT_DIRECTORY}/doxyfilter.bat" @ONLY)
      set (DOXYGEN_INPUT_FILTER "doxyfilter.bat")
    endif ()
  endif ()
  basis_list_to_delimited_string (
    DOXYGEN_FILTER_PATTERNS "\"\nFILTER_PATTERNS       += \"" ${DOXYGEN_FILTER_PATTERNS}
  )
  if (DOXYGEN_FILTER_PATTERNS)
    set (DOXYGEN_FILTER_PATTERNS "\"${DOXYGEN_FILTER_PATTERNS}\"")
  endif ()
  # exclude patterns
  list (APPEND DOXYGEN_EXCLUDE_PATTERNS "cmake_install.cmake")
  list (APPEND DOXYGEN_EXCLUDE_PATTERNS "CTestTestfile.cmake")
  basis_list_to_delimited_string (
    DOXYGEN_EXCLUDE_PATTERNS "\"\nEXCLUDE_PATTERNS      += \"" ${DOXYGEN_EXCLUDE_PATTERNS}
  )
  set (DOXYGEN_EXCLUDE_PATTERNS "\"${DOXYGEN_EXCLUDE_PATTERNS}\"")
  # section for man pages
  if (NOT DOXYGEN_MAN_SECTION)
    set (DOXYGEN_MAN_SECTION 3)
  endif ()
  # other settings
  if (NOT DOXYGEN_COLS_IN_ALPHA_INDEX OR DOXYGEN_COLS_IN_ALPHA_INDEX MATCHES "[^0-9]")
    set (DOXYGEN_COLS_IN_ALPHA_INDEX 1)
  endif ()
  basis_list_to_delimited_string (DOXYGEN_IGNORE_PREFIX " " ${DOXYGEN_IGNORE_PREFIX})
  # click & jump in emacs and Visual Studio
  if (CMAKE_BUILD_TOOL MATCHES "(msdev|devenv)")
    set (DOXYGEN_WARN_FORMAT "\"$file($line) : $text \"")
  else ()
    set (DOXYGEN_WARN_FORMAT "\"$file:$line: $text \"")
  endif ()
  # installation directories
  set (BASIS_INSTALL_${TARGET_NAME_U}_DIR "" CACHE PATH "Installation directory for Doxygen ${TARGET_NAME} target.")
  mark_as_advanced (BASIS_INSTALL_${TARGET_NAME_U}_DIR)
  foreach (f IN LISTS DOXYGEN_OUTPUT)
    string (TOUPPER "${f}" F)
    if (BASIS_INSTALL_${TARGET_NAME_U}_DIR)
      set (DOXYGEN_${F}_DESTINATION "${BASIS_INSTALL_${TARGET_NAME_U}_DIR}") # user setting
    endif ()
    if (NOT DOXYGEN_${F}_DESTINATION)
      if (DOXYGEN_DESTINATION)
        set (DOXYGEN_${F}_DESTINATION "${DOXYGEN_DESTINATION}") # common destination
      elseif (f MATCHES "man")
        if (INSTALL_MAN_DIR)
          set (DOXYGEN_MAN_DESTINATION "${INSTALL_MAN_DIR}/man${DOXYGEN_MAN_SECTION}") # default for manual pages
        endif ()
      elseif (NOT f MATCHES "html") # do not install excludes by default
        set (DOXYGEN_${F}_DESTINATION "${INSTALL_DOC_DIR}") # default destination
      endif ()
    endif ()
  endforeach ()
  # determine tool to generate pdf documentation, see USE_PDFLATEX in Doxyfile.in
  find_package (LATEX QUIET)
  if (PDFLATEX_COMPILER)
    set (DOXYGEN_USE_PDFLATEX YES)
  else ()
    set (DOXYGEN_USE_PDFLATEX NO)
    message (STATUS "pdflatex not found. For higher quality PDFs make sure pdflatex is installed and found on the PATH.")
  endif ()
  # use default custom HTML files if available and none explicitly specified
  if (NOT DOXYGEN_HTML_EXTRA_STYLESHEET)
    if (EXISTS "${PROJECT_DOCRES_DIR}/doxygen_extra.css.in")
      set (DOXYGEN_HTML_EXTRA_STYLESHEET "${PROJECT_DOCRES_DIR}/doxygen_extra.css.in")
    elseif (EXISTS "${PROJECT_DOCRES_DIR}/doxygen_extra.css")
      set (DOXYGEN_HTML_EXTRA_STYLESHEET "${PROJECT_DOCRES_DIR}/doxygen_extra.css")
    elseif (EXISTS "${PROJECT_DOC_DIR}/doxygen_extra.css.in")
      set (DOXYGEN_HTML_EXTRA_STYLESHEET "${PROJECT_DOC_DIR}/doxygen_extra.css.in")
    elseif (EXISTS "${PROJECT_DOC_DIR}/doxygen_extra.css")
      set (DOXYGEN_HTML_EXTRA_STYLESHEET "${PROJECT_DOC_DIR}/doxygen_extra.css")
    elseif (EXISTS "${PROJECT_DOCRES_DIR}/doxygen.css.in")
      set (DOXYGEN_HTML_EXTRA_STYLESHEET "${PROJECT_DOCRES_DIR}/doxygen.css.in")
    elseif (EXISTS "${PROJECT_DOCRES_DIR}/doxygen.css")
      set (DOXYGEN_HTML_EXTRA_STYLESHEET "${PROJECT_DOCRES_DIR}/doxygen.css")
    elseif (EXISTS "${PROJECT_DOC_DIR}/doxygen.css.in")
      set (DOXYGEN_HTML_EXTRA_STYLESHEET "${PROJECT_DOC_DIR}/doxygen.css.in")
    elseif (EXISTS "${PROJECT_DOC_DIR}/doxygen.css")
      set (DOXYGEN_HTML_EXTRA_STYLESHEET "${PROJECT_DOC_DIR}/doxygen.css")
    endif ()
  endif ()
  if (DOXYGEN_HTML_HEADER_IS_DEFAULT AND NOT DOXYGEN_HTML_EXTRA_STYLESHEET)
    set (DOXYGEN_HTML_EXTRA_STYLESHEET "${BASIS_MODULE_PATH}/doxygen_extra.css.in")
  endif ()
  if (NOT DOXYGEN_HTML_HEADER)
    if (EXISTS "${PROJECT_DOCRES_DIR}/doxygen_header.html.in")
      set (DOXYGEN_HTML_HEADER "${PROJECT_DOCRES_DIR}/doxygen_header.html.in")
    elseif (EXISTS "${PROJECT_DOCRES_DIR}/doxygen_header.html")
      set (DOXYGEN_HTML_HEADER "${PROJECT_DOCRES_DIR}/doxygen_header.html")
    elseif (EXISTS "${PROJECT_DOC_DIR}/doxygen_header.html.in")
      set (DOXYGEN_HTML_HEADER "${PROJECT_DOC_DIR}/doxygen_header.html.in")
    elseif (EXISTS "${PROJECT_DOC_DIR}/doxygen_header.html")
      set (DOXYGEN_HTML_HEADER "${PROJECT_DOC_DIR}/doxygen_header.html")
    endif ()
  endif ()
  if (NOT DOXYGEN_HTML_FOOTER)
    if (EXISTS "${PROJECT_DOCRES_DIR}/doxygen_footer.html.in")
      set (DOXYGEN_HTML_FOOTER "${PROJECT_DOCRES_DIR}/doxygen_footer.html.in")
    elseif (EXISTS "${PROJECT_DOCRES_DIR}/doxygen_footer.html")
      set (DOXYGEN_HTML_FOOTER "${PROJECT_DOCRES_DIR}/doxygen_footer.html")
    elseif (EXISTS "${PROJECT_DOC_DIR}/doxygen_footer.html.in")
      set (DOXYGEN_HTML_FOOTER "${PROJECT_DOC_DIR}/doxygen_footer.html.in")
    elseif (EXISTS "${PROJECT_DOC_DIR}/doxygen_footer.html")
      set (DOXYGEN_HTML_FOOTER "${PROJECT_DOC_DIR}/doxygen_footer.html")
    endif ()
  endif ()
  # configure/copy custom HTML ressource files
  if (DOXYGEN_GENERATE_HTML)
    set (DOXYGEN_HTML_EXTRA_STYLESHEET_NAME "doxygen_extra.css")
    set (OUTPUT_HTML_EXTRA_STYLESHEET "${DOXYGEN_OUTPUT_DIRECTORY}/${DOXYGEN_HTML_EXTRA_STYLESHEET_NAME}")
    set (OUTPUT_HTML_HEADER           "${DOXYGEN_OUTPUT_DIRECTORY}/doxygen_header.html")
    set (OUTPUT_HTML_FOOTER           "${DOXYGEN_OUTPUT_DIRECTORY}/doxygen_footer.html")
    foreach (res IN ITEMS EXTRA_STYLESHEET HEADER FOOTER)
      if (EXISTS "${DOXYGEN_HTML_${res}}")
        if (DOXYGEN_HTML_${res} MATCHES "\\.in$")
          configure_file (${DOXYGEN_HTML_${res}} "${OUTPUT_HTML_${res}}" @ONLY)
        elseif (DOXYGEN_HTML_${res})
          configure_file (${DOXYGEN_HTML_${res}} "${OUTPUT_HTML_${res}}" COPYONLY)
        else ()
          set (OUTPUT_HTML_${res})
        endif ()
      endif ()
    endforeach ()
    set  (DOXYGEN_HTML_EXTRA_STYLESHEET "${OUTPUT_HTML_EXTRA_STYLESHEET}")
    set  (DOXYGEN_HTML_HEADER           "${OUTPUT_HTML_HEADER}")
    set  (DOXYGEN_HTML_FOOTER           "${OUTPUT_HTML_FOOTER}")
  else ()
    set (DOXYGEN_HTML_EXTRA_STYLESHEET)
    set (DOXYGEN_HTML_FOOTER)
    set (DOXYGEN_HTML_HEADER)
  endif ()
  if (DOXYGEN_PROVIDER_LOGO)
    list (APPEND DOXYGEN_HTML_EXTRA_FILES "${DOXYGEN_PROVIDER_LOGO}")
  endif()
  if (DOXYGEN_DIVISION_LOGO)
    list (APPEND DOXYGEN_HTML_EXTRA_FILES "${DOXYGEN_DIVISION_LOGO}")
  endif()
  if (DOXYGEN_HTML_EXTRA_FILES)
    basis_list_to_delimited_string (
      DOXYGEN_HTML_EXTRA_FILES "\"\nHTML_EXTRA_FILES      += \"" ${DOXYGEN_HTML_EXTRA_FILES}
    )
  endif ()
  # list of enabled Doxygen comment sections
  basis_join ("${DOXYGEN_ENABLED_SECTIONS}" " " DOXYGEN_ENABLED_SECTIONS)
  # configure Doxygen configuration file
  set (DOXYFILE "${DOXYGEN_OUTPUT_DIRECTORY}/Doxyfile.${TARGET_NAME_L}")
  configure_file ("${DOXYGEN_DOXYFILE}" "${DOXYFILE}" @ONLY)
  if (CMAKE_GENERATOR MATCHES "Visual Studio|Xcode")
    file (READ "${DOXYFILE}" DOXYFILE_CONTENT)
    foreach (CONFIG IN LISTS CMAKE_CONFIGURATION_TYPES)
      string (REPLACE "$<${BASIS_GE_CONFIG}>" "${CONFIG}" DOXYFILE_CONTENT_CONFIG "${DOXYFILE_CONTENT}")
      file (WRITE "${DOXYFILE}.${CONFIG}" "${DOXYFILE_CONTENT_CONFIG}")
    endforeach ()
    unset (DOXYFILE_CONTENT)
    unset (DOXYFILE_CONTENT_CONFIG)
    set (DOXYFILE "${DOXYFILE}.$<${BASIS_GE_CONFIG}>")
  endif ()
  # add build target
  set (OPTALL)
  if (BUILD_DOCUMENTATION AND BASIS_ALL_DOC)
    set (OPTALL "ALL")
  endif ()
  file (MAKE_DIRECTORY "${DOXYGEN_OUTPUT_DIRECTORY}")
  add_custom_target (
    ${TARGET_UID} ${OPTALL} "${DOXYGEN_EXECUTABLE}" "${DOXYFILE}"
    WORKING_DIRECTORY "${DOXYGEN_OUTPUT_DIRECTORY}"
    COMMENT "Building documentation ${TARGET_UID}..."
  )
  # memorize certain settings which might be useful to know by other functions
  # in particular, in case of the use of the XML output by other documentation
  # build tools such as Sphinx, the function that wants to make use of this
  # output can check if the Doxygen target has been configured properly and
  # further requires to know the location of the XML output
  set_target_properties (
    ${TARGET_UID}
    PROPERTIES
      BASIS_TYPE       Doxygen
      OUTPUT_DIRECTORY "${DOXYGEN_OUTPUT_DIRECTORY}"
      DOXYFILE         "${DOXYGEN_DOXYFILE}"
      TAGFILE          "${DOXYGEN_TAGFILE}"
      OUTPUT           "${DOXYGEN_OUTPUT}"
  )
  foreach (f IN LISTS DOXYGEN_OUTPUT)
    string (TOUPPER "${f}" F)
    set_target_properties (
      ${TARGET_UID}
      PROPERTIES
        ${F}_INSTALL_DIRECTORY "${DOXYGEN_${F}_DESTINATION}"
        ${F}_OUTPUT_DIRECTORY  "${DOXYGEN_OUTPUT_DIRECTORY}/${DOXYGEN_${F}_OUTPUT}"
    )
    set_property (
      DIRECTORY
      APPEND PROPERTY
        ADDITIONAL_MAKE_CLEAN_FILES
          "${DOXYGEN_OUTPUT_DIRECTORY}/${DOXYGEN_${F}_OUTPUT}"
    )
  endforeach ()
  if (DOXYGEN_TAGFILE)
    set_property (
      DIRECTORY
      APPEND PROPERTY
        ADDITIONAL_MAKE_CLEAN_FILES
          "${DOXYGEN_TAGFILE}"
    )
  endif ()
  # The Doxygen filter, if a build target of this project, has to be build
  # before the documentation can be generated.
  if (TARGET "${DOXYFILTER}")
    add_dependencies (${TARGET_UID} ${DOXYFILTER})
  endif ()
  # The public header files shall be configured/copied before.
  if (TARGET headers)
    add_dependencies (${TARGET_UID} headers)
  endif ()
  # The documentation shall be build after all other executable and library
  # targets have been build. For example, a .py.in script file shall first
  # be "build", i.e., configured before the documentation is being generated
  # from the configured .py file.
  basis_get_project_property (TARGETS PROPERTY TARGETS)
  foreach (_UID ${TARGETS})
    get_target_property (BASIS_TYPE ${_UID} "BASIS_TYPE")
    if (BASIS_TYPE MATCHES "SCRIPT|EXECUTABLE|LIBRARY")
      add_dependencies (${TARGET_UID} ${_UID})
    endif ()
  endforeach ()
  # add general "doc" target
  if (NOT DOXYGEN_EXCLUDE_FROM_DOC)
    if (NOT TARGET doc)
      add_custom_target (doc)
    endif ()
    add_dependencies (doc ${TARGET_UID})
  endif ()
  # install documentation
  install (
    CODE
      "
      set (HTML_DESTINATION \"${DOXYGEN_HTML_DESTINATION}\")
      set (MAN_DESTINATION  \"${DOXYGEN_MAN_DESTINATION}\")

      function (install_doxydoc FMT)
        string (TOUPPER \"\${FMT}\" FMT_U)
        set (CMAKE_INSTALL_PREFIX \"\${\${FMT_U}_DESTINATION}\")
        if (NOT CMAKE_INSTALL_PREFIX)
          return ()
        elseif (NOT IS_ABSOLUTE \"\${CMAKE_INSTALL_PREFIX}\")
          set (CMAKE_INSTALL_PREFIX \"${CMAKE_INSTALL_PREFIX}/\${CMAKE_INSTALL_PREFIX}\")
        endif ()
        set (EXT)
        set (DIR \"\${FMT}\")
        if (FMT MATCHES \".pdf\")
          set (EXT \".pdf\")
          set (DIR \"latex\")
        elseif (FMT MATCHES \".rtf\")
          set (EXT \".rtf\")
        elseif (FMT MATCHES \"man\")
          set (EXT \".?\")
        endif ()
        file (
          GLOB_RECURSE
            FILES
          RELATIVE \"${DOXYGEN_OUTPUT_DIRECTORY}/\${DIR}\"
            \"${DOXYGEN_OUTPUT_DIRECTORY}/\${DIR}/*\${EXT}\"
        )
        foreach (F IN LISTS FILES)
          execute_process (
            COMMAND \"${CMAKE_COMMAND}\" -E compare_files
                \"${DOXYGEN_OUTPUT_DIRECTORY}/\${DIR}/\${F}\"
                \"\${CMAKE_INSTALL_PREFIX}/\${F}\"
            RESULT_VARIABLE RC
            OUTPUT_QUIET
            ERROR_QUIET
          )
          if (RC EQUAL 0)
            message (STATUS \"Up-to-date: \${CMAKE_INSTALL_PREFIX}/\${F}\")
          else ()
            message (STATUS \"Installing: \${CMAKE_INSTALL_PREFIX}/\${F}\")
            execute_process (
              COMMAND \"${CMAKE_COMMAND}\" -E copy_if_different
                  \"${DOXYGEN_OUTPUT_DIRECTORY}/\${DIR}/\${F}\"
                  \"\${CMAKE_INSTALL_PREFIX}/\${F}\"
              RESULT_VARIABLE RC
              OUTPUT_QUIET
              ERROR_QUIET
            )
            if (RC EQUAL 0)
              list (APPEND CMAKE_INSTALL_MANIFEST_FILES \"\${CMAKE_INSTALL_PREFIX}/\${F}\")
            else ()
              message (STATUS \"Failed to install \${CMAKE_INSTALL_PREFIX}/\${F}\")
            endif ()
          endif ()
        endforeach ()
        if (FMT MATCHES \"html\" AND EXISTS \"${DOXYGEN_TAGFILE}\")
          get_filename_component (DOXYGEN_TAGFILE_NAME \"${DOXYGEN_TAGFILE}\" NAME)
          execute_process (
            COMMAND \"${CMAKE_COMMAND}\" -E copy_if_different
              \"${DOXYGEN_TAGFILE}\"
              \"\${CMAKE_INSTALL_PREFIX}/\${DOXYGEN_TAGFILE_NAME}\"
          )
          list (APPEND CMAKE_INSTALL_MANIFEST_FILES \"\${CMAKE_INSTALL_PREFIX}/\${DOXYGEN_TAGFILE_NAME}\")
        endif ()
      endfunction ()

      foreach (FMT IN ITEMS html pdf rtf man)
        install_doxydoc (\${FMT})
      endforeach ()
      "
  )
  # done
  message (STATUS "Adding documentation ${TARGET_UID}... - done")
endfunction ()

# ----------------------------------------------------------------------------
## @brief Add documentation target to be generated by Sphinx (sphinx-build).
#
# This function adds a build target to generate documentation from
# <a href="http://docutils.sourceforge.net/rst.html">reStructuredText</a>
# (.rst files) using <a href="http://sphinx.pocoo.org/">Sphinx</a>.
#
# @param [in] TARGET_NAME Name of the documentation target.
# @param [in] ARGN        List of arguments. The valid arguments are:
# @par
# <table border="0">
#   <tr>
#     @tp @b EXCLUDE_FROM_DOC @endtp
#     <td>By default, the specified target is build as part of the global
#         @c doc target. If this option is given, however, the added
#         documentation will not be build as part of this target.</td>
#   </tr>
#   <tr>
#     @tp @b BUILDER(S) builder... @endtp
#     <td>Sphinx builders to use. For each named builder, a build target
#         named &lt;TARGET_NAME&gt;_&lt;builder&gt; is added.</td>
#   </tr>
#   <tr>
#     @tp @b DEFAULT_BUILDER builder @endtp
#     <td>Default Sphinx builder to associated with the @c TARGET_NAME
#         build target. Defaults to the first builder named by @c BUILDERS.</td>
#   </tr>
#   <tr>
#     @tp @b AUTHOR(S) name @endtp
#     <td>Names of authors who wrote this documentation.
#         (default: @c PROJECT_AUTHORS)</td>
#   </tr>
#   <tr>
#     @tp @b COPYRIGHT text @endtp
#     <td>Copyright statement for generated files. (default: @c PROJECT_COPYRIGHT)</td>
#   </tr>
#   <tr>
#     @tp @b COMPONENT component @endtp
#     <td>Name of the component this documentation belongs to.
#         Defaults to @c BASIS_RUNTIME_COMPONENT.</td>
#   </tr>
#   <tr>
#     @tp @b DESTINATION dir @endtp
#     <td>Installation directory prefix. Used whenever there is no specific
#         destination specified for a particular Sphinx builder. Defaults to
#         @c BASIS_INSTALL_&ltTARGET&gt;_DIR in case of HTML output if set.
#         Otherwise, the generated HTML files are not installed.</td>
#   </tr>
#   <tr>
#     @tp @b &lt;BUILDER&gt;_DESTINATION dir @endtp
#     <td>Installation directory for files generated by the specific builder.<td>
#   </tr>
#   <tr>
#     @tp @b EXTENSIONS ext... @endtp
#     <td>Names of Sphinx extensions to enable.</td>
#   </tr>
#   <tr>
#     @tp @b BREATHE target... @endtp
#     <td>Adds a project for the breathe extension which allows the
#         inclusion of in-source code documentation extracted by Doxygen.
#         For this to work, the specified Doxygen target has to be
#         configured with the XML output enabled.</td>
#   </tr>
#   <tr>
#     @tp @b DOXYLINK target... @endtp
#     <td>Adds a role for the doxylink Sphinx extension which allows to cross-reference
#         generated HTML API documentation generated by Doxygen.</td>
#   </tr>
#   <tr>
#     @tp @b DOXYLINK_URL url @endtp
#     <td>URL to Doxygen documentation. Use DOXYLINK_PREFIX and/or DOXYLINK_SUFFIX
#         instead if you use multiple Doxygen targets, where the target name is
#         part of the URL.</td>
#   </tr>
#   <tr>
#     @tp @b DOXYLINK_PREFIX url @endtp
#     <td>Prefix to use for links to Doxygen generated documentation pages
#         as generated by the doxylink Sphinx extension. If this prefix does
#         not start with a protocol such as http:// or https://, it is prefixed
#         to the default path determined by this function relative to the build
#         or installed Doxygen documentation.</td>
#   </tr>
#   <tr>
#     @tp @b DOXYLINK_SUFFIX suffix @endtp
#     <td>Suffix for links to Doxygen generated documentation pages as generated
#         by the doxylink Sphinx extension.</td>
#   </tr>
#   <tr>
#     @tp @b DOXYDOC target... @endtp
#     <td>Alias for both @c BREATHE and @c DOXYLINK options.</td>
#   </tr>
#   <tr>
#     @tp @b CONFIG_FILE file @endtp
#     <td>Sphinx configuration file. Defaults to @c BASIS_SPHINX_CONFIG.</td>
#   </tr>
#   <tr>
#     @tp @b SOURCE_DIRECTORY @endtp
#     <td>Root directory of Sphinx source files.
#         Defaults to the current source directory or, if a subdirectory
#         named @c TARGET_NAME in lowercase only exists, to this subdirectory.</td>
#   </tr>
#   <tr>
#     @tp @b OUTPUT_NAME @endtp
#     <td>Output name for generated documentation such as PDF document or MAN page.
#         Defaults to @c PROJECT_NAME.</td>
#   </tr>
#   <tr>
#     @tp @b OUTPUT_DIRECTORY @endtp
#     <td>Root output directory for generated files. Defaults to the binary
#         directory corresponding to the set @c SOURCE_DIRECTORY.</td>
#   </tr>
#   <tr>
#     @tp @b TAG tag @endtp
#     <td>Tag argument of <tt>sphinx-build</tt>.</td>
#   </tr>
#   <tr>
#     @tp @b TEMPLATES_PATH @endtp
#     <td>Path to template files. Defaults to <tt>SOURCE_DIRECTORY/templates/</tt>.</td>
#   </tr>
#   <tr>
#     @tp @b MASTER_DOC name @endtp
#     <td>Name of master document. Defaults to <tt>index</tt>.</td>
#   </tr>
#   <tr>
#     @tp @b EXCLUDE_PATTERN pattern @endtp
#     <td>A glob-style pattern that should be excluded when looking for source files.
#         Specify this option more than once to specify multiple exclude patterns.
#         They are matched against the source file names relative to the source directory,
#         using slashes as directory separators on all platforms.</td>
#   </tr>
#   <tr>
#     @tp @b HTML_TITLE title @endtp
#     <td>Title of HTML web site.</td>
#   </tr>
#   <tr>
#     @tp @b HTML_THEME theme @endtp
#     <td>Name of HTML theme. Defaults to the @c sbia theme included with BASIS.</td>
#   </tr>
#   <tr>
#     @tp @b HTML_THEME_PATH dir @endtp
#     <td>Directory of HTML theme. Defaults to @c BASIS_SPHINX_HTML_THEME_PATH.</td>
#   </tr>
#   <tr>
#     @tp @b HTML_LOGO file @endtp
#     <td>Logo to display in sidebar of HTML pages.</td>
#   </tr>
#   <tr>
#     @tp @b HTML_FAVICON file @endtp
#     <td>Favorite square icon often displayed by browsers in the tab bar.
#         Should be a @c .ico file.</td>
#   </tr>
#   <tr>
#     @tp @b HTML_STATIC_PATH dir @endtp
#     <td>Directory for static files of HTML pages. Defaults to <tt>SOURCE_DIRECTORY/static/</tt>.</td>
#   </tr>
#   <tr>
#     @tp @b HTML_STYLE css @endtp
#     <td>The style sheet to use for HTML pages. A file of that name must exist either in Sphinx'
#         default static/ path or the specified @c HTML_STATIC_PATH. Default is the stylesheet
#         given by the selected theme. If you only want to add or override a few things compared
#         to the theme’s stylesheet, use CSS \@import to import the theme’s stylesheet.</td>
#   </tr>
#   <tr>
#     @tp @b HTML_SIDEBARS name... @endtp
#     <td>Names of HTML template files for sidebar(s). Defaults to none if not specified.
#         Valid default templates are @c localtoc, @c globaltoc, @c searchbox, @c relations,
#         @c sourcelink. See <a href="http://sphinx.pocoo.org/config.html#confval-html_sidebars">
#         Shinx documentation of html_sidebars option</a>. Custom templates can be used as
#         well by copying the template <tt>.html</tt> file to the @c TEMPLATES_PATH directory.</td>
#   </tr>
#   <tr>
#     @tp @b NO_HTML_DOMAIN_INDICES @endtp
#     <td>Set Sphinx configuration option @c html_domain_indices to @c False. (Default: @c True)</td>
#   </tr>
#   <tr>
#     @tp @b NO_HTML_MODINDEX @endtp
#     <td>Set Sphinx configuration option @c html_use_modindex to @c False. (Default: @c True)</td>
#   </tr>
#   <tr>
#     @tp @b NO_HTML_INDEX @endtp
#     <td>Set Sphinx configuration option @c html_use_index to @c False. (Default: @c True)</td>
#   </tr>
#   <tr>
#     @tp @b LATEX_MASTER_DOC name @endtp
#     <td>Name of master document for LaTeX builder. Defaults to <tt>MASTER_DOC</tt>.</td>
#   </tr>
#   <tr>
#     @tp @b LATEX_TITLE title @endtp
#     <td>Title for LaTeX/PDF output. Defaults to title of <tt>index.rst</tt>.</td>
#   </tr>
#   <tr>
#     @tp @b LATEX_LOGO file @endtp
#     <td>Logo to display above title in generated LaTeX/PDF output.</td>
#   </tr>
#   <tr>
#     @tp @b LATEX_DOCUMENT_CLASS howto|manual @endtp
#     <td>Document class to use by @c latex builder.</td>
#   </tr>
#   <tr>
#     @tp @b LATEX_SHOW_URLS @endtp
#     <td>See Sphinx documentation of the
#         <a href="http://sphinx.pocoo.org/config.html#confval-latex_show_urls">latex_show_urls</a> option.</td>
#   </tr>
#   <tr>
#     @tp @b LATEX_SHOW_PAGEREFS @endtp
#     <td>See Sphinx documentation of the
#         <a href="">latex_show_pagerefs</a> option.</td>
#   </tr>
#   <tr>
#     @tp @b MAN_SECTION num @endtp
#     <td>Section number for manual pages generated by @c man builder.</td>
#   </tr>
# </table>
#
# @sa basis_add_doc()
function (basis_add_sphinx_doc TARGET_NAME)
  # check target name
  basis_check_target_name ("${TARGET_NAME}")
  basis_make_target_uid (TARGET_UID "${TARGET_NAME}")
  string (TOLOWER "${TARGET_NAME}" TARGET_NAME_L)
  string (TOUPPER "${TARGET_NAME}" TARGET_NAME_U)
  # verbose output
  message (STATUS "Adding documentation ${TARGET_UID}...")
  # parse arguments
  set (ONE_ARG_OPTIONS
    COMPONENT
    DEFAUL_BUILDER
    DESTINATION HTML_DESTINATION MAN_DESTINATION TEXINFO_DESTINATION
    CONFIG_FILE
    SOURCE_DIRECTORY OUTPUT_DIRECTORY OUTPUT_NAME TAG
    COPYRIGHT MASTER_DOC
    HTML_TITLE HTML_THEME HTML_LOGO HTML_FAVICON HTML_THEME_PATH HTML_STYLE
    LATEX_MASTER_DOC LATEX_TITLE LATEX_LOGO LATEX_DOCUMENT_CLASS LATEX_SHOW_URLS LATEX_SHOW_PAGEREFS
    MAN_SECTION
    DOXYLINK_URL DOXYLINK_PREFIX DOXYLINK_SUFFIX
  )
  # note that additional multiple value arguments are parsed later on below
  # this is necessary b/c all unparsed arguments are considered to be options
  # of the used HTML theme
  CMAKE_PARSE_ARGUMENTS (SPHINX
    "EXCLUDE_FROM_DOC;NO_HTML_DOMAIN_INDICES;NO_HTML_MODINDEX;NO_HTML_INDEX"
    "${ONE_ARG_OPTIONS}"
    ""
    ${ARGN}
  )
  # source directory
  if (NOT SPHINX_SOURCE_DIRECTORY)
    if (IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${TARGET_NAME}")
      set (SPHINX_SOURCE_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${TARGET_NAME}")
    else ()
      set (SPHINX_SOURCE_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
    endif ()
  elseif (NOT IS_ABSOLUTE "${SPHINX_SOURCE_DIRECTORY}")
    get_filename_component (SPHINX_SOURCE_DIRECTORY "${SPHINX_SOURCE_DIRECTORY}" ABSOLUTE)
  endif ()
  # component
  if (NOT SPHINX_COMPONENT)
    set (SPHINX_COMPONENT "${BASIS_RUNTIME_COMPONENT}")
  endif ()
  if (NOT SPHINX_COMPONENT)
    set (SPHINX_COMPONENT "Unspecified")
  endif ()
  # find Sphinx
  find_package (Sphinx COMPONENTS build QUIET)
  if (NOT Sphinx-build_EXECUTABLE)
    if (BUILD_DOCUMENTATION)
      message (FATAL_ERROR "Command sphinx-build not found! Either install Sphinx and/or set Sphinx-build_EXECUTABLE or disable BUILD_DOCUMENTATION.")
    endif ()
    message (STATUS "Command sphinx-build not found. Generation of ${TARGET_UID} documentation disabled.")
    message (STATUS "Adding documentation ${TARGET_UID}... - skipped")
    return ()
  endif ()
  if (DEFINED Sphinx_VERSION_MAJOR AND Sphinx_VERSION_MAJOR LESS 1)
    if (BUILD_DOCUMENTATION)
      message (FATAL_ERROR "Found sphinx-build is too old (v${Sphinx_VERSION_STRING})! Please install a more recent version and/or set Sphinx-build_EXECUTABLE to a newer version or disable BUILD_DOCUMENTATION.")
    endif ()
    message (STATUS "Command sphinx-build is too old (v${Sphinx_VERSION_STRING}). Generation of ${TARGET_UID} documentation disabled.")
    message (STATUS "Adding documentation ${TARGET_UID}... - skipped")
    return ()
  endif ()
  # parse remaining arguments
  set (SPHINX_HTML_THEME_OPTIONS)
  set (SPHINX_BUILDERS)
  set (SPHINX_AUTHORS)
  set (SPHINX_EXTENSIONS)
  set (SPHINX_BREATHE_TARGETS)
  set (SPHINX_DOXYLINK_TARGETS)
  set (SPHINX_HTML_SIDEBARS)
  set (SPHINX_TEMPLATES_PATH)
  set (SPHINX_HTML_STATIC_PATH)
  set (SPHINX_EXCLUDE_PATTERNS)
  set (SPHINX_DEPENDS)
  set (OPTION_NAME)
  set (OPTION_VALUE)
  set (OPTION_PATTERN "(authors?|builders?|extensions|breathe|doxylink|doxydoc|html_sidebars|templates_path|html_static_path|exclude_pattern)")
  foreach (ARG IN LISTS SPHINX_UNPARSED_ARGUMENTS)
    if (NOT OPTION_NAME OR ARG MATCHES "^[A-Z_]+$")
      # SPHINX_HTML_THEME_OPTIONS
      if (OPTION_NAME AND NOT OPTION_NAME MATCHES "^${OPTION_PATTERN}$")
        if (NOT OPTION_VALUE)
          message (FATAL_ERROR "Option ${OPTION_NAME} is missing an argument!")
        endif ()
        list (LENGTH OPTION_VALUE NUM)
        if (NUM GREATER 1)
          basis_list_to_delimited_string (OPTION_VALUE ", " NOAUTOQUOTE ${OPTION_VALUE})
          set (OPTION_VALUE "[${OPTION_VALUE}]")
        endif ()
        list (APPEND SPHINX_HTML_THEME_OPTIONS "'${OPTION_NAME}': ${OPTION_VALUE}")
      endif ()
      # name of next option
      set (OPTION_NAME "${ARG}")
      set (OPTION_VALUE)
      string (TOLOWER "${OPTION_NAME}" OPTION_NAME)
    # BUILDER option
    elseif (OPTION_NAME MATCHES "^builders?$")
      if (ARG MATCHES "html dirhtml singlehtml pdf latex man text texinfo linkcheck")
        message (FATAL_ERROR "Invalid/Unsupported Sphinx builder: ${ARG}")
      endif ()
      list (APPEND SPHINX_BUILDERS "${ARG}")
    # AUTHORS option
    elseif (OPTION_NAME MATCHES "^authors?$")
      list (APPEND SPHINX_AUTHORS "'${ARG}'")
    # EXTENSIONS option
    elseif (OPTION_NAME MATCHES "^extensions$")
      # built-in extension
      if (ARG MATCHES "^(autodoc|autosummary|doctest|intersphinx|pngmath|jsmath|mathjax|graphvis|inheritance_graph|ifconfig|coverage|todo|extlinks|viewcode)$")
        set (ARG "sphinx.ext.${CMAKE_MATCH_0}")
      # map originial name of extensions included with BASIS
      elseif (BASIS_SPHINX_EXTENSIONS_PATH AND ARG MATCHES "^sphinxcontrib.(doxylink)$")
        set (ARG "${CMAKE_MATCH_1}")
      endif ()
      list (APPEND SPHINX_EXTENSIONS "'${ARG}'")
    # DOXYDOC
    elseif (OPTION_NAME MATCHES "^doxydoc$")
      list (APPEND SPHINX_BREATHE_TARGETS  "${ARG}")
      list (APPEND SPHINX_DOXYLINK_TARGETS "${ARG}")
    # BREATHE
    elseif (OPTION_NAME MATCHES "^breathe$")
      list (APPEND SPHINX_BREATHE_TARGETS "${ARG}")
    # DOXYLINK
    elseif (OPTION_NAME MATCHES "^doxylink$")
      list (APPEND SPHINX_DOXYLINK_TARGETS "${ARG}")
    # HTML_SIDEBARS
    elseif (OPTION_NAME MATCHES "^html_sidebars$")
      if (NOT ARG MATCHES "\\.html?$")
        set (ARG "${ARG}.html")
      endif ()
      list (APPEND SPHINX_HTML_SIDEBARS "'${ARG}'")
    # TEMPLATES_PATH
    elseif (OPTION_NAME MATCHES "^templates_path$")
      if (NOT IS_ABSOLUTE "${ARG}")
        set (ARG "${SPHINX_SOURCE_DIRECTORY}/${ARG}")
      endif ()
      list (APPEND SPHINX_TEMPLATES_PATH "'${ARG}'")
    # HTML_STATIC_PATH
    elseif (OPTION_NAME MATCHES "^html_static_path$")
      if (NOT IS_ABSOLUTE "${ARG}")
        set (ARG "${SPHINX_SOURCE_DIRECTORY}/${ARG}")
      endif ()
      list (APPEND SPHINX_HTML_STATIC_PATH "'${ARG}'")
    # EXCLUDE_PATTERN
    elseif (OPTION_NAME MATCHES "^exclude_pattern$")
      list (APPEND SPHINX_EXCLUDE_PATTERNS "'${ARG}'")
    # value of theme option
    else ()
      if (ARG MATCHES "^(TRUE|FALSE)$")
        string (TOLOWER "${ARG}" "${ARG}")
      endif ()
      if (NOT ARG MATCHES "^\\[.*\\]$|^{.*}$")
        set (ARG "'${ARG}'")
      endif ()
      list (APPEND OPTION_VALUE "${ARG}")
    endif ()
  endforeach ()
  # append parsed option setting to SPHINX_HTML_THEME_OPTIONS
  if (OPTION_NAME AND NOT OPTION_NAME MATCHES "^${OPTION_PATTERN}$")
    if (NOT OPTION_VALUE)
      message (FATAL_ERROR "Option ${OPTION_NAME} is missing an argument!")
    endif ()
    list (LENGTH OPTION_VALUE NUM)
    if (NUM GREATER 1)
      basis_list_to_delimited_string (OPTION_VALUE ", " NOAUTOQUOTE ${OPTION_VALUE})
      set (OPTION_VALUE "[${OPTION_VALUE}]")
    endif ()
    list (APPEND SPHINX_HTML_THEME_OPTIONS "'${OPTION_NAME}': ${OPTION_VALUE}")
  endif ()
  # authors
  if (NOT SPHINX_AUTHORS)
    foreach (AUTHOR IN LISTS PROJECT_AUTHORS)
      list (APPEND SPHINX_AUTHORS "'${AUTHOR}'")
    endforeach ()
  endif ()
  if (NOT SPHINX_COPYRIGHT)
    set (SPHINX_COPYRIGHT "${PROJECT_COPYRIGHT}")
  endif ()
  # default builders
  if (NOT SPHINX_BUILDERS)
    set (SPHINX_BUILDERS html dirhtml singlehtml man pdf texinfo text linkcheck)
  endif ()
  if (SPHINX_DEFAULT_BUILDER)
    list (FIND SPHINX_BUILDERS "${SPHINX_DEFAULT_BUILDER}" IDX)
    if (IDX EQUAL -1)
      list (INSERT SPHINX_BUILDERS 0 "${SPHINX_DEFAULT_BUILDER}")
    endif ()
  else ()
    list (GET SPHINX_BUILDERS 0 SPHINX_DEFAULT_BUILDER)
  endif ()
  # output directories
  if (NOT SPHINX_OUTPUT_NAME)
    set (SPHINX_OUTPUT_NAME "${PROJECT_NAME}")
  endif ()
  if (NOT SPHINX_OUTPUT_DIRECTORY)
    if (IS_ABSOLUTE "${SPHINX_OUTPUT_NAME}")
      get_filename_component (SPHINX_OUTPUT_DIRECTORY "${SPHINX_OUTPUT_NAME}" PATH)
    else ()
      basis_get_relative_path (SPHINX_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}" "${SPHINX_SOURCE_DIRECTORY}")
    endif ()
  endif ()
  if (NOT IS_ABSOLUTE "${SPHINX_OUTPUT_DIRECTORY}")
    set (SPHINX_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${SPHINX_OUTPUT_DIRECTORY}")
  endif ()
  foreach (b IN LISTS SPHINX_BUILDERS)
    string (TOUPPER "${b}" B)
    if (SPHINX_${B}_OUTPUT_DIRECTORY)
      if (NOT IS_ABSOLUTE "${SPHINX_${B}_OUTPUT_DIRECTORY}")
        set (SPHINX_${B}_OUTPUT_DIRECTORY "${SPHINX_OUTPUT_DIRECTORY}/${SPHINX_${B}_OUTPUT_DIRECTORY}")
      endif ()
    else ()
      set (SPHINX_${B}_OUTPUT_DIRECTORY "${SPHINX_OUTPUT_DIRECTORY}/${b}")
    endif ()
  endforeach ()
  if (IS_ABSOLUTE "${SPHINX_OUTPUT_NAME}")
    basis_get_relative_path (SPHINX_OUTPUT_NAME "${SPHINX_OUTPUT_DIRECTORY}" NAME_WE)
  endif ()
  # configuration directory
  basis_get_relative_path (SPHINX_CONFIG_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}" "${SPHINX_SOURCE_DIRECTORY}")
  set (SPHINX_CONFIG_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${SPHINX_CONFIG_DIRECTORY}")
  # build configuration
  if (NOT SPHINX_MASTER_DOC)
    set (SPHINX_MASTER_DOC "index")
  endif ()
  if (NOT SPHINX_LATEX_MASTER_DOC)
    set (SPHINX_LATEX_MASTER_DOC "${SPHINX_MASTER_DOC}")
  endif ()
  if (NOT SPHINX_TEMPLATES_PATH AND EXISTS "${SPHINX_SOURCE_DIRECTORY}/templates")
    set (SPHINX_TEMPLATES_PATH "'${SPHINX_SOURCE_DIRECTORY}/templates'")
  endif ()
  if (NOT SPHINX_HTML_STATIC_PATH AND EXISTS "${SPHINX_SOURCE_DIRECTORY}/static")
    set (SPHINX_HTML_STATIC_PATH "'${SPHINX_SOURCE_DIRECTORY}/static'")
  endif ()
  if (NOT SPHINX_HTML_THEME)
    set (SPHINX_HTML_THEME "${BASIS_SPHINX_HTML_THEME}")
  endif ()
  if (NOT SPHINX_UNPARSED_ARGUMENTS AND SPHINX_HTML_THEME STREQUAL BASIS_SPHINX_HTML_THEME)
    set (SPHINX_UNPARSED_ARGUMENTS ${BASIS_SPHINX_HTML_THEME_OPTIONS})
  endif () 
  if (NOT SPHINX_LATEX_DOCUMENT_CLASS)
    set (SPHINX_LATEX_DOCUMENT_CLASS "howto")
  endif ()
  if (NOT SPHINX_MAN_SECTION)
    set (SPHINX_MAN_SECTION 1)
  endif ()
  # installation directories
  set (BASIS_INSTALL_${TARGET_NAME_U}_DIR "" CACHE PATH "Installation directory for documentation ${TARGET_NAME} target.")
  mark_as_advanced (BASIS_INSTALL_${TARGET_NAME_U}_DIR)
  foreach (b IN LISTS SPHINX_BUILDERS)
    string (TOUPPER "${b}" B)
    if (BASIS_INSTALL_${TARGET_NAME_U}_DIR)
      set (SPHINX_${B}_DESTINATION "${BASIS_INSTALL_${TARGET_NAME_U}_DIR}") # user setting
    endif ()
    if (NOT SPHINX_${B}_DESTINATION)
      if (SPHINX_DESTINATION)                           
        set (SPHINX_${B}_DESTINATION "${DESTINATION}") # common destination
      elseif (b MATCHES "text")
        set (SPHINX_${B}_DESTINATION "${INSTALL_DOC_DIR}/${TARGET_NAME_L}")
      elseif (b MATCHES "man")
        if (INSTALL_MAN_DIR)
          set (SPHINX_${B}_DESTINATION "${INSTALL_MAN_DIR}/man${SPHINX_MAN_SECTION}") # default for manual pages
        endif ()
      elseif (b MATCHES "texinfo")
        if (INSTALL_TEXINFO_DIR)
          set (SPHINX_${B}_DESTINATION "${INSTALL_TEXINFO_DIR}") # default for Texinfo files
        endif ()
      elseif (NOT b MATCHES "html") # do not install excludes by default
        set (SPHINX_${B}_DESTINATION "${INSTALL_DOC_DIR}") # default location
      endif ()
    endif ()
  endforeach ()
  if (SPHINX_HTML_DESTINATION)
    foreach (b IN LISTS SPHINX_BUILDERS)
      if (b MATCHES "(dir|single)html")
        string (TOUPPER "${b}" B)
        if (NOT SPHINX_${B}_DESTINATION)
          set (SPHINX_${B}_DESTINATION "${SPHINX_HTML_DESTINATION}")
        endif ()
      endif ()
    endforeach ()
  endif ()
  # enable required extension
  if (SPHINX_DOXYLINK_TARGETS)
    if (BASIS_SPHINX_EXTENSIONS_PATH)
      if (NOT SPHINX_EXTENSIONS MATCHES "(^|;)?doxylink(;|$)?")
        list (APPEND SPHINX_EXTENSIONS "'doxylink'")
      endif ()
    else ()
      if (NOT SPHINX_EXTENSIONS MATCHES "(^|;)?sphinxcontrib.doxylink(;|$)?")
        list (APPEND SPHINX_EXTENSIONS "'sphinxcontrib.doxylink'")
      endif ()
    endif ()
  endif ()
  if (SPHINX_BREATHE_TARGETS AND NOT SPHINX_EXTENSIONS MATCHES "(^|;)?breathe(;|$)?")
    list (APPEND SPHINX_EXTENSIONS "'breathe'")
  endif ()
  # doxylink configuration
  foreach (TARGET IN LISTS SPHINX_DOXYLINK_TARGETS)
    basis_get_target_uid (UID "${TARGET}")
    get_target_property (TYPE ${UID} BASIS_TYPE)
    if (NOT TYPE MATCHES "Doxygen")
      message (FATAL_ERROR "Invalid argument for DOXYLINK: Target ${UID} either unknown or it is not a Doxygen target!")
    endif ()
    get_target_property (DOXYGEN_OUTPUT ${UID} OUTPUT)
    if (NOT DOXYGEN_OUTPUT MATCHES "html")
      message (FATAL_ERROR "Doxygen target ${UID} was not configured to generate HTML output! This output is required by the doxylink Sphinx extension.")
    endif ()
    get_target_property (DOXYGEN_TAGFILE        ${UID} TAGFILE)
    get_target_property (DOXYGEN_HTML_DIRECTORY ${UID} HTML_INSTALL_DIRECTORY)
    set (DOXYLINK_PATH)
    if (SPHINX_DOXYLINK_URL)
      set (DOXYLINK_PATH "${SPHINX_DOXYLINK_URL}")
    elseif (SPHINX_DOXYLINK_PREFIX MATCHES "^[a-z]://")
      set (DOXYLINK_PATH "${SPHINX_DOXYLINK_PREFIX}${TARGET}${SPHINX_DOXYLINK_SUFFIX}")
    else ()
      set (DOXYLINK_BASE_DIR)
      if (DOXYGEN_HTML_DIRECTORY)
        if (SPHINX_DOXYLINK_BASE_DIR)
          set (DOXYLINK_BASE_DIR "${SPHINX_DOXYLINK_BASE_DIR}")
        elseif (SPHINX_HTML_INSTALL_DIRECTORY)
          set (DOXYLINK_BASE_DIR "${SPHINX_HTML_INSTALL_DIRECTORY}")
        endif ()
      else ()
        get_target_property (DOXYGEN_HTML_DIRECTORY ${UID} HTML_OUTPUT_DIRECTORY)
        if (SPHINX_DOXYLINK_BASE_DIR)
          set (DOXYLINK_BASE_DIR "${SPHINX_DOXYLINK_BASE_DIR}")
        else ()
          set (DOXYLINK_BASE_DIR "${SPHINX_HTML_OUTPUT_DIRECTORY}")
        endif ()
      endif ()
      if (DOXYLINK_BASE_DIR)
        basis_get_relative_path (DOXYLINK_PATH "${DOXYLINK_BASE_DIR}" "${DOXYGEN_HTML_DIRECTORY}")
      else ()
        set (DOXYLINK_PATH "${TARGET}") # safe fall back
      endif ()
      set (DOXYLINK_PATH "${SPHINX_DOXYLINK_PREFIX}${DOXYLINK_PATH}${SPHINX_DOXYLINK_SUFFIX}")
    endif ()
    list (APPEND SPHINX_DOXYLINK "'${TARGET}': ('${DOXYGEN_TAGFILE}', '${DOXYLINK_PATH}')")
    list (APPEND SPHINX_DEPENDS ${UID})
  endforeach ()
  # breathe configuration
  set (SPHINX_BREATHE_PROJECTS)
  set (SPHINX_BREATHE_DEFAULT_PROJECT)
  foreach (TARGET IN LISTS SPHINX_BREATHE_TARGETS)
    basis_get_target_uid (UID "${TARGET}")
    get_target_property (TYPE ${UID} BASIS_TYPE)
    if (NOT TYPE MATCHES "Doxygen")
      message (FATAL_ERROR "Invalid argument for BREATHE_PROJECTS: Target ${UID} either unknown or it is not a Doxygen target!")
    endif ()
    get_target_property (DOXYGEN_OUTPUT ${UID} OUTPUT)
    if (NOT DOXYGEN_OUTPUT MATCHES "xml")
      message (FATAL_ERROR "Doxygen target ${UID} was not configured to generate XML output! This output is required by the Sphinx extension breathe.")
    endif ()
    get_target_property (DOXYGEN_OUTPUT_DIRECTORY ${UID} XML_OUTPUT_DIRECTORY)
    list (APPEND SPHINX_BREATHE_PROJECTS "'${TARGET}': '${DOXYGEN_OUTPUT_DIRECTORY}'")
    if (NOT SPHINX_BREATHE_DEFAULT_PROJECT)
      set (SPHINX_BREATHE_DEFAULT_PROJECT "${TARGET}")
    endif ()
    list (APPEND SPHINX_DEPENDS ${UID})
  endforeach ()
  # HTML output options
  if (SPHINX_NO_HTML_DOMAIN_INDICES)
    set (SPHINX_USE_DOMAIN_INDICES False)
  else ()
    set (SPHINX_USE_DOMAIN_INDICES True)
  endif ()
  if (SPHINX_NO_HTML_MODINDEX)
    set (SPHINX_USE_MODINDEX False)
  else ()
    set (SPHINX_USE_MODINDEX True)
  endif ()
  if (SPHINX_NO_HTML_INDEX)
    set (SPHINX_USE_INDEX False)
  else ()
    set (SPHINX_USE_INDEX True)
  endif ()
  # LaTeX output options
  if (NOT SPHINX_LATEX_SHOW_URLS)
    set (SPHINX_LATEX_SHOW_URLS "no")
  endif ()
  if (SPHINX_LATEX_SHOW_PAGEREFS)
    set (SPHINX_LATEX_SHOW_PAGEREFS "True")
  else ()
    set (SPHINX_LATEX_SHOW_PAGEREFS "False")
  endif ()
  # turn html_logo, html_favicon, and latex_logo into absolute file path
  foreach (L IN ITEMS HTML_LOGO HTML_FAVICON LATEX_LOGO)
    if (SPHINX_${L} AND NOT IS_ABSOLUTE "${SPHINX_${L}}")
      if (EXISTS "${SPHINX_SOURCE_DIRECTORY}/${SPHINX_${L}}")
        set (SPHINX_${L} "${SPHINX_SOURCE_DIRECTORY}/${SPHINX_${L}}")
      elseif (L MATCHES "^HTML|^LATEX")
        foreach (D IN LISTS SPHINX_${CMAKE_MATCH_0}_STATIC_PATH)
          string (REGEX REPLACE "^'|'$" "" D "${D}")
          if (EXISTS "${D}/${SPHINX_${L}}")
            set (SPHINX_${L} "${D}/${SPHINX_${L}}")
            break ()
          endif ()
        endforeach ()
      endif ()
    endif ()
  endforeach ()
  # turn CMake lists into Python lists
  basis_list_to_delimited_string (SPHINX_EXTENSIONS         ", " NOAUTOQUOTE ${SPHINX_EXTENSIONS})
  basis_list_to_delimited_string (SPHINX_HTML_THEME_OPTIONS ", " NOAUTOQUOTE ${SPHINX_HTML_THEME_OPTIONS})
  basis_list_to_delimited_string (SPHINX_AUTHORS            ", " NOAUTOQUOTE ${SPHINX_AUTHORS})
  basis_list_to_delimited_string (SPHINX_DOXYLINK           ", " NOAUTOQUOTE ${SPHINX_DOXYLINK})
  basis_list_to_delimited_string (SPHINX_BREATHE_PROJECTS   ", " NOAUTOQUOTE ${SPHINX_BREATHE_PROJECTS})
  basis_list_to_delimited_string (SPHINX_HTML_SIDEBARS      ", " NOAUTOQUOTE ${SPHINX_HTML_SIDEBARS})
  basis_list_to_delimited_string (SPHINX_TEMPLATES_PATH     ", " NOAUTOQUOTE ${SPHINX_TEMPLATES_PATH})
  basis_list_to_delimited_string (SPHINX_HTML_STATIC_PATH   ", " NOAUTOQUOTE ${SPHINX_HTML_STATIC_PATH})
  basis_list_to_delimited_string (SPHINX_EXCLUDE_PATTERNS   ", " NOAUTOQUOTE ${SPHINX_EXCLUDE_PATTERNS})
  # configuration file
  if (NOT SPHINX_CONFIG_FILE)
    set (SPHINX_CONFIG_FILE "${BASIS_SPHINX_CONFIG}")
  endif ()
  get_filename_component (SPHINX_CONFIG_FILE "${SPHINX_CONFIG_FILE}" ABSOLUTE)
  if (EXISTS "${SPHINX_CONFIG_FILE}")
    configure_file ("${SPHINX_CONFIG_FILE}" "${SPHINX_CONFIG_DIRECTORY}/conf.py" @ONLY)
  elseif (EXISTS "${SPHINX_CONFIG_FILE}.in")
    configure_file ("${SPHINX_CONFIG_FILE}.in" "${SPHINX_CONFIG_DIRECTORY}/conf.py" @ONLY)
  else ()
    message (FATAL_ERROR "Missing Sphinx configuration file ${SPHINX_CONFIG_FILE}!")
  endif ()
  # add target to build documentation
  set (OPTIONS -a -N -n)
  if (NOT BASIS_VERBOSE)
    list (APPEND OPTIONS "-q")
  endif ()
  foreach (TAG IN LISTS SPHINX_TAG)
    list (APPEND OPTIONS "-t" "${TAG}")
  endforeach ()
  add_custom_target (${TARGET_UID}_all) # target to run all builders
  foreach (BUILDER IN LISTS SPHINX_BUILDERS)
    set (SPHINX_BUILDER "${BUILDER}")
    set (SPHINX_POST_COMMAND)
    if (BUILDER MATCHES "pdf|texinfo")
      if (UNIX)
        set (SPHINX_MAKE_COMMAND "make")
        if (BUILDER MATCHES "pdf")
          set (SPHINX_BUILDER "latex")
          configure_file ("${BASIS_MODULE_PATH}/sphinx_make.sh.in" "${SPHINX_CONFIG_DIRECTORY}/make.sh" @ONLY)
          set (SPHINX_MAKE_COMMAND "${SPHINX_CONFIG_DIRECTORY}/make.sh")
          execute_process (COMMAND chmod +x "${SPHINX_CONFIG_DIRECTORY}/make.sh")
        endif ()
        set (SPHINX_POST_COMMAND COMMAND "${SPHINX_MAKE_COMMAND}" -C "${SPHINX_OUTPUT_DIRECTORY}/${SPHINX_BUILDER}")
      elseif (BUILDER MATCHES "pdf")
        message (WARNING "Target ${TARGET_UID} requires the execution of pdflatex which is currently"
                         " only executed after sphinx-build on Unix platforms. On Windows, pdflatex"
                         " must be executed manually in the directory\n\n\t${SPHINX_OUTPUT_DIRECTORY}/${SPHINX_BUILDER}")
      endif ()
    endif ()
    if (Sphinx_PYTHON_EXECUTABLE)
      add_custom_target (
        ${TARGET_UID}_${BUILDER}
            "${Sphinx_PYTHON_EXECUTABLE}" ${Sphinx_PYTHON_OPTIONS}
                "${Sphinx-build_EXECUTABLE}" ${OPTIONS}
                  -b ${SPHINX_BUILDER}
                  -c "${SPHINX_CONFIG_DIRECTORY}"
                  -d "${SPHINX_CONFIG_DIRECTORY}/doctrees"
                  "${SPHINX_SOURCE_DIRECTORY}"
                  "${SPHINX_OUTPUT_DIRECTORY}/${SPHINX_BUILDER}"
            ${SPHINX_POST_COMMAND}
            ${OPTDEPENDS}
        WORKING_DIRECTORY "${SPHINX_CONFIG_DIRECTORY}"
        COMMENT "Building documentation ${TARGET_UID} (${BUILDER})..."
      )
    elseif (UNIX)
      add_custom_target (
        ${TARGET_UID}_${BUILDER}
            "${Sphinx-build_EXECUTABLE}" ${OPTIONS}
                  -b ${SPHINX_BUILDER}
                  -c "${SPHINX_CONFIG_DIRECTORY}"
                  -d "${SPHINX_CONFIG_DIRECTORY}/doctrees"
                  "${SPHINX_SOURCE_DIRECTORY}"
                  "${SPHINX_OUTPUT_DIRECTORY}/${SPHINX_BUILDER}"
            ${SPHINX_POST_COMMAND}
            ${OPTDEPENDS}
        WORKING_DIRECTORY "${SPHINX_CONFIG_DIRECTORY}"
        COMMENT "Building documentation ${TARGET_UID} (${BUILDER})..."
      )
    else ()
      message (FATAL_ERROR "Python executable required to run Sphinx could not be determined by FindSphinx.cmake!"
                           " Set the advanced PYTHON_EXECUTABLE variable and try again.")
    endif ()
    if (SPHINX_DEPENDS)
      add_dependencies (${TARGET_UID}_${BUILDER} ${SPHINX_DEPENDS})
    endif ()
    add_dependencies (${TARGET_UID}_all ${TARGET_UID}_${BUILDER})
    # cleanup on "make clean"
    set_property (
      DIRECTORY
      APPEND PROPERTY
        ADDITIONAL_MAKE_CLEAN_FILES
          "${SPHINX_OUTPUT_DIRECTORY}"
    )
  endforeach ()
  # add general target which depends on default builder only
  if (BUILD_DOCUMENTATION AND BASIS_ALL_DOC)
    add_custom_target (${TARGET_UID} ALL)
  else ()
    add_custom_target (${TARGET_UID})
  endif ()
  add_dependencies (${TARGET_UID} ${TARGET_UID}_${SPHINX_DEFAULT_BUILDER})
  # add general "doc" target
  if (NOT SPHINX_EXCLUDE_FROM_DOC)
    if (NOT TARGET doc)
      add_custom_target (doc)
    endif ()
    add_dependencies (doc ${TARGET_UID}_${SPHINX_DEFAULT_BUILDER})
  endif ()
  # memorize important target properties
  set_target_properties (
    ${TARGET_UID}
    PROPERTIES
      BASIS_TYPE       Sphinx
      BUILDERS         "${SPHINX_BUILDERS}"
      SOURCE_DIRECTORY "${SPHINX_SOURCE_DIRECTORY}"
      OUTPUT_DIRECTORY "${SPHINX_OUTPUT_DIRECTORY}"
      CONFIG_DIRECTORY "${SPHINX_CONFIG_DIRECTORY}"
  )
  foreach (b IN LISTS SPHINX_BUILDERS)
    string (TOUPPER ${b} B)
    set_target_properties (${TARGET_UID} PROPERTIES ${B}_INSTALL_DIRECTORY "${SPHINX_${B}_DESTINATION}")
  endforeach ()
  # cleanup on "make clean"
  set_property (
    DIRECTORY
    APPEND PROPERTY
      ADDITIONAL_MAKE_CLEAN_FILES
        "${SPHINX_CONFIG_DIRECTORY}/doctrees"
  )
  # install documentation
  install (
    CODE
      "
      set (HTML_DESTINATION    \"${SPHINX_HTML_DESTINATION}\")
      set (PDF_DESTINATION     \"${SPHINX_PDF_DESTINATION}\")
      set (LATEX_DESTINATION   \"${SPHINX_LATEX_DESTINATION}\")
      set (MAN_DESTINATION     \"${SPHINX_MAN_DESTINATION}\")
      set (TEXINFO_DESTINATION \"${SPHINX_TEXINFO_DESTINATION}\")
      set (TEXT_DESTINATION    \"${SPHINX_TEXT_DESTINATION}\")

      function (install_sphinx_doc BUILDER)
        if (BUILDER MATCHES \"pdf\")
          set (SPHINX_BUILDER \"latex\")
        else ()
          set (SPHINX_BUILDER \"\${BUILDER}\")
        endif ()
        string (TOUPPER \"\${BUILDER}\" BUILDER_U)
        set (CMAKE_INSTALL_PREFIX \"\${\${BUILDER_U}_DESTINATION}\")
        if (NOT CMAKE_INSTALL_PREFIX)
          return ()
        elseif (NOT IS_ABSOLUTE \"\${CMAKE_INSTALL_PREFIX}\")
          set (CMAKE_INSTALL_PREFIX \"${CMAKE_INSTALL_PREFIX}/\${CMAKE_INSTALL_PREFIX}\")
        endif ()
        set (EXT)
        if (BUILDER MATCHES \"pdf\")
          set (EXT \".pdf\")
        elseif (BUILDER MATCHES \"man\")
          set (EXT \".?\")
        elseif (BUILDER MATCHES \"texinfo\")
          set (EXT \".info\")
        endif ()
        file (
          GLOB_RECURSE
            FILES
          RELATIVE \"${SPHINX_OUTPUT_DIRECTORY}/\${SPHINX_BUILDER}\"
            \"${SPHINX_OUTPUT_DIRECTORY}/\${SPHINX_BUILDER}/*\${EXT}\"
        )
        foreach (F IN LISTS FILES)
          if (NOT F MATCHES \"\\\\.buildinfo\")
            set (RC 1)
            if (NOT BUILDER MATCHES \"texinfo\")
              execute_process (
                COMMAND \"${CMAKE_COMMAND}\" -E compare_files
                    \"${SPHINX_OUTPUT_DIRECTORY}/\${SPHINX_BUILDER}/\${F}\"
                    \"\${CMAKE_INSTALL_PREFIX}/\${F}\"
                RESULT_VARIABLE RC
                OUTPUT_QUIET
                ERROR_QUIET
              )
            endif ()
            if (RC EQUAL 0)
              message (STATUS \"Up-to-date: \${CMAKE_INSTALL_PREFIX}/\${F}\")
            else ()
              message (STATUS \"Installing: \${CMAKE_INSTALL_PREFIX}/\${F}\")
              if (BUILDER MATCHES \"texinfo\")
                if (EXISTS \"\${CMAKE_INSTALL_PREFIX}/dir\")
                  execute_process (
                    COMMAND install-info
                        \"${SPHINX_OUTPUT_DIRECTORY}/\${SPHINX_BUILDER}/\${F}\"
                        \"\${CMAKE_INSTALL_PREFIX}/dir\"
                    RESULT_VARIABLE RC
                    OUTPUT_QUIET
                    ERROR_QUIET
                  )
                else ()
                  execute_process (
                    COMMAND \"${CMAKE_COMMAND}\" -E copy_if_different
                        \"${SPHINX_OUTPUT_DIRECTORY}/\${SPHINX_BUILDER}/\${F}\"
                        \"\${CMAKE_INSTALL_PREFIX}/dir\"
                    RESULT_VARIABLE RC
                    OUTPUT_QUIET
                    ERROR_QUIET
                  )
                endif ()
              else ()
                execute_process (
                  COMMAND \"${CMAKE_COMMAND}\" -E copy_if_different
                      \"${SPHINX_OUTPUT_DIRECTORY}/\${SPHINX_BUILDER}/\${F}\"
                      \"\${CMAKE_INSTALL_PREFIX}/\${F}\"
                  RESULT_VARIABLE RC
                  OUTPUT_QUIET
                  ERROR_QUIET
                )
              endif ()
              if (RC EQUAL 0)
                # also remember .info files for deinstallation via install-info --delete
                list (APPEND CMAKE_INSTALL_MANIFEST_FILES \"\${CMAKE_INSTALL_PREFIX}/\${F}\")
              else ()
                message (STATUS \"Failed to install \${CMAKE_INSTALL_PREFIX}/\${F}\")
              endif ()
            endif ()
          endif ()
        endforeach ()
      endfunction ()

      set (BUILDERS \"${SPHINX_BUILDERS}\")
      set (HTML_INSTALLED FALSE)
      foreach (BUILDER IN LISTS BUILDERS)
        if ((BUILDER MATCHES \"html\" AND NOT HTML_INSTALLED) OR
              (BUILDER MATCHES \"texinfo|man\" AND UNIX) OR
              NOT BUILDER MATCHES \"html|texinfo|man|latex|linkcheck\")
          install_sphinx_doc (\${BUILDER})
          if (BUILDER MATCHES \"html\")
            set (HTML_INSTALLED TRUE)
          endif ()
        endif ()
      endforeach ()
      "
  )
  # done
  message (STATUS "Adding documentation ${TARGET_UID}... - done")
endfunction ()

# ============================================================================
# change log
# ============================================================================

# ----------------------------------------------------------------------------
## @brief Add target for generation of ChangeLog file.
#
# The ChangeLog is either generated from the Subversion or Git log depending
# on which revision control system is used by the project. Moreover, the
# project's source directory must be either a Subversion working copy or
# the root of a Git repository, respectively. In case of Subversion, if the
# command-line tool svn2cl(.sh) is installed, it is used to output a nicer
# formatted change log.
function (basis_add_changelog)
  _basis_make_target_uid (TARGET_UID changelog)

  option (BUILD_CHANGELOG "Request build and/or installation of the ChangeLog." OFF)
  mark_as_advanced (BUILD_CHANGELOG)
  set (CHANGELOG_FILE "${PROJECT_BINARY_DIR}/ChangeLog")

  message (STATUS "Adding ChangeLog...")

  if (NOT PROJECT_IS_MODULE AND BUILD_CHANGELOG)
    set (_ALL "ALL")
  else ()
    set (_ALL)
  endif ()

  set (DISABLE_BUILD_CHANGELOG FALSE)

  # --------------------------------------------------------------------------
  # generate ChangeLog from Subversion history
  if (EXISTS "${PROJECT_SOURCE_DIR}/.svn")
    find_package (Subversion QUIET)
    if (Subversion_FOUND)

      if (_ALL)
        message ("Generation of ChangeLog enabled as part of ALL."
                 " Be aware that the ChangeLog generation from the Subversion"
                 " commit history can take several minutes and may require the"
                 " input of your Subversion repository credentials during the"
                 " build. If you would like to build the ChangeLog separate"
                 " from the rest of the software package, disable the option"
                 " BUILD_CHANGELOG. You can then build the changelog target"
                 " separate from ALL.")
      endif ()

      # using svn2cl command
      find_program (
        SVN2CL_EXECUTABLE
          NAMES svn2cl svn2cl.sh
          DOC   "The command line tool svn2cl."
      )
      mark_as_advanced (SVN2CL_EXECUTABLE)
      if (SVN2CL_EXECUTABLE)
        add_custom_target (
          ${TARGET_UID} ${_ALL}
          COMMAND "${SVN2CL_EXECUTABLE}"
              "--output=${CHANGELOG_FILE}"
              "--linelen=79"
              "--reparagraph"
              "--group-by-day"
              "--include-actions"
              "--separate-daylogs"
              "${PROJECT_SOURCE_DIR}"
          COMMAND "${CMAKE_COMMAND}"
              "-DCHANGELOG_FILE:FILE=${CHANGELOG_FILE}" -DINPUTFORMAT=SVN2CL
              -P "${BASIS_MODULE_PATH}/PostprocessChangeLog.cmake"
          WORKING_DIRECTORY "${PROJECT_BINARY_DIR}"
          COMMENT "Generating ChangeLog from Subversion log (using svn2cl)..."
        )
      # otherwise, use svn log output directly
      else ()
        add_custom_target (
          ${TARGET_UID} ${_ALL}
          COMMAND "${CMAKE_COMMAND}"
              "-DCOMMAND=${Subversion_SVN_EXECUTABLE};log"
              "-DWORKING_DIRECTORY=${PROJECT_SOURCE_DIR}"
              "-DOUTPUT_FILE=${CHANGELOG_FILE}"
              -P "${BASIS_SCRIPT_EXECUTE_PROCESS}"
          COMMAND "${CMAKE_COMMAND}"
              "-DCHANGELOG_FILE:FILE=${CHANGELOG_FILE}" -DINPUTFORMAT=SVN
              -P "${BASIS_MODULE_PATH}/PostprocessChangeLog.cmake"
          COMMENT "Generating ChangeLog from Subversion log..."
          VERBATIM
        )
      endif ()

    else ()
      if (BASIS_VERBOSE)
        message (STATUS "Project is SVN working copy but Subversion executable was not found."
                        " Generation of ChangeLog disabled.")
      endif ()
      set (DISABLE_BUILD_CHANGELOG TRUE)
    endif ()

  # --------------------------------------------------------------------------
  # generate ChangeLog from Git log
  elseif (EXISTS "${PROJECT_SOURCE_DIR}/.git")
    find_package (Git QUIET)
    if (GIT_FOUND)

      add_custom_target (
        ${TARGET_UID} ${_ALL}
        COMMAND "${CMAKE_COMMAND}"
            "-DCOMMAND=${GIT_EXECUTABLE};log;--date-order;--date=short;--pretty=format:%ad\ \ %an%n%n%w(79,8,10)* %s%n%n%b%n"
            "-DWORKING_DIRECTORY=${PROJECT_SOURCE_DIR}"
            "-DOUTPUT_FILE=${CHANGELOG_FILE}"
            -P "${BASIS_SCRIPT_EXECUTE_PROCESS}"
        COMMAND "${CMAKE_COMMAND}"
            "-DCHANGELOG_FILE=${CHANGELOG_FILE}" -DINPUTFORMAT=GIT
            -P "${BASIS_MODULE_PATH}/PostprocessChangeLog.cmake"
        COMMENT "Generating ChangeLog from Git log..."
        VERBATIM
      )

    else ()
      if (BASIS_VERBOSE)
        message (STATUS "Project is Git repository but Git executable was not found."
                        " Generation of ChangeLog disabled.")
      endif ()
      set (DISABLE_BUILD_CHANGELOG TRUE)
    endif ()

  # --------------------------------------------------------------------------
  # neither SVN nor Git repository
  else ()
    if (BASIS_VERBOSE)
      message (STATUS "Project source directory ${PROJECT_SOURCE_DIR} is neither"
                      " SVN working copy nor Git repository. Generation of ChangeLog disabled.")
    endif ()
    set (DISABLE_BUILD_CHANGELOG TRUE)
  endif ()

  # --------------------------------------------------------------------------
  # disable changelog target
  if (DISABLE_BUILD_CHANGELOG)
    set (BUILD_CHANGELOG OFF CACHE INTERNAL "" FORCE)
    message (STATUS "Adding ChangeLog... - skipped")
    return ()
  endif ()

  # --------------------------------------------------------------------------
  # cleanup on "make clean"
  set_property (DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES "${CHANGELOG_FILE}")

  # --------------------------------------------------------------------------
  # install ChangeLog
  get_filename_component (CHANGELOG_NAME "${CHANGELOG_FILE}" NAME)
  if (PROJECT_IS_MODULE)
    set (CHANGELOG_NAME "${CHANGELOG_NAME}-${PROJECT_NAME}")
  endif ()

  install (
    FILES       "${CHANGELOG_FILE}"
    DESTINATION "${TOPLEVEL_INSTALL_DOC_DIR}"
    COMPONENT   "${BASIS_RUNTIME_COMPONENT}"
    RENAME      "${CHANGELOG_NAME}"
    OPTIONAL
  )

  message (STATUS "Adding ChangeLog... - done")
endfunction ()
