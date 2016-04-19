# ============================================================================
# Copyright (c) 2011-2012 University of Pennsylvania
# Copyright (c) 2013-2016 Andreas Schuh
# All rights reserved.
#
# See COPYING file for license information or visit
# https://cmake-basis.github.io/download.html#license
# ============================================================================

##############################################################################
# @file  BasisSettings.cmake
# @brief Default project-independent settings.
#
# This module defines global CMake constants and variables which are used
# by the BASIS CMake functions and macros. Hence, these values can be used
# to configure the behavior of these functions to some extent without the
# need to modify the functions themselves.
#
# @note As this file also sets the CMake policies to be used, it has to
#       be included using the @c NO_POLICY_SCOPE in order for these policies
#       to take effect also in the including file and its subdirectories.
#
# @attention Be careful when caching any of the variables. Usually, this
#            file is included in the root CMake configuration file of the
#            project which may also be a module of another project and hence
#            may overwrite this project's settings.
#
# @attention Keep in mind that this file is included before any other
#            BASIS module. Further, project-specific information such as
#            the project name are not defined yet.
#
# @ingroup BasisSettings
##############################################################################

if (__BASIS_SETTINGS_INCLUDED)
  return ()
else ()
  set (__BASIS_SETTINGS_INCLUDED TRUE)
endif ()

## @addtogroup BasisSettings
# @{


# ============================================================================
# CMake version and policies
# ============================================================================

cmake_minimum_required (VERSION 2.8.12 FATAL_ERROR)

# Add policies introduced with CMake versions newer than the one specified
# above. These policies would otherwise trigger a policy not set warning by
# newer CMake versions.

if (POLICY CMP0016)
  cmake_policy (SET CMP0016 NEW)
endif ()

if (POLICY CMP0017)
  cmake_policy (SET CMP0017 NEW)
endif ()

if (POLICY CMP0042)
  cmake_policy (SET CMP0042 NEW)
endif ()

# ============================================================================
# required modules
# ============================================================================

include ("${CMAKE_CURRENT_LIST_DIR}/CommonTools.cmake")

# ============================================================================
# generator expressions
# ============================================================================

## @brief Name of build configuration ("$<CONFIG>") generator expression
if (CMAKE_MAJOR_VERSION LESS 3)
  set (BASIS_GE_CONFIG "CONFIGURATION")
else ()
  set (BASIS_GE_CONFIG "CONFIG")
endif ()

# ============================================================================
# meta-data lists
# ============================================================================

## @brief Names of project meta-data switches.
set (
  BASIS_METADATA_LIST_SWITCH
    EXCLUDE_FROM_ALL # exclude module from BUILD_ALL_MODULES
)

## @brief Names of project meta-data with only one argument.
#  @see basis_project() in ProjectTools.cmake
set (
  BASIS_METADATA_LIST_SINGLE
    AUTHOR
    NAME
    SUBPROJECT
    PACKAGE        # short alias for PACKAGE_NAME
    WEBSITE        # short alias for PACKAGE_WEBSITE
    VENDOR         # short alias for PACKAGE_VENDOR
    PROVIDER       # short alias for PACKAGE_VENDOR (not PROVIDER_NAME, see basis_project!)
    PACKAGE_NAME
    PACKAGE_VENDOR # package vendor ID used for installation path
    PACKAGE_WEBSITE
    PACKAGE_LOGO
    PROVIDER_NAME
    PROVIDER_WEBSITE
    PROVIDER_LOGO
    DIVISION_NAME
    DIVISION_WEBSITE
    DIVISION_LOGO
    COPYRIGHT
    LICENSE
    CONTACT
    VERSION
    SOVERSION
    TEMPLATE       # used by basisproject tool
    INCLUDE_DIR    # alias for INCLUDE_DIRS
    CODE_DIR       # alias for CODE_DIRS
    TOOLS_DIR      # alias for TOOLS_DIRS
    MODULES_DIR    # single directory containing multiple modules, see also MODULE_DIRS
    CONFIG_DIR     # directory containing the CMake/BASIS configuration
    DATA_DIR       # directory containing the auxiliary program data
    DOC_DIR        # directory containing the documentation
    DOCRES_DIR     # directory containing the ressource files such as a project logo
    EXAMPLE_DIR    # directory containing some example files
    LIBRARY_DIR    # directory containing script libraries such as Perl or Python modules
    TESTING_DIR    # directory containing the source code and data of the software tests
)

## @brief Names of project meta-data with multiple arguments.
#  @see basis_project() in ProjectTools.cmake
set (
  BASIS_METADATA_LIST_MULTI
    AUTHORS
    DESCRIPTION
    LANGUAGES
    DEFAULT_MODULES  # list of modules to be enabled by default
    EXTERNAL_MODULES # list of external project modules
    DEPENDS
    OPTIONAL_DEPENDS
    TOOLS_DEPENDS
    OPTIONAL_TOOLS_DEPENDS
    TEST_DEPENDS
    OPTIONAL_TEST_DEPENDS
    INCLUDE_DIRS   # list of directories containing public header files
    CODE_DIRS      # list of directories containing source code files, see also CODE_DIR
    TOOLS_DIRS     # list of directories containing source code files of applications, see also TOOLS_DIR
    MODULE_DIRS    # list of separate module directories, see also MODULES_DIR
    OTHER_DIRS     # list of additional (generic) project subdirectories
    SUBDIRS        # deprecated option name for OTHER_DIRS
)

## @brief Names of project meta-data.
#  @see basis_project() in ProjectTools.cmake
set (
  BASIS_METADATA_LIST
    ${BASIS_METADATA_LIST_SWITCH}
    ${BASIS_METADATA_LIST_SINGLE}
    ${BASIS_METADATA_LIST_MULTI}
)

## @brief Names of additional meta-data for Slicer modules with only one argument.
set (
  BASIS_SLICER_METADATA_LIST_SINGLE
     HOMEPAGE
     ICONURL
     STATUS
     SCREENSHOTURLS
)

## @brief Names of additional meta-data for Slicer modules with multiple arguments.
set (
  BASIS_SLICER_METADATA_LIST_MULTI
     ACKNOWLEDGEMENTS
     CATEGORY
     CONTRIBUTORS
     LICENSE_SHORT_DESCRIPTION
)

## @brief Names of additional meta-data for Slicer modules.
set (
  BASIS_SLICER_METADATA_LIST
    ${BASIS_SLICER_METADATA_LIST_SINGLE}
    ${BASIS_SLICER_METADATA_LIST_MULTI}
)

# ============================================================================
# constants and global settings
# ============================================================================

## @brief List of name patterns used for special purpose targets.
#
# Contains a list of target name patterns that are used by the BASIS functions
# for special purposes and are hence not to be used for project targets.
set (
  BASIS_RESERVED_TARGET_NAMES
    "all"
    "bundle"
    "bundle_source"
    "changelog"
    "clean"
    "depend"
    "doc"
    "headers"
    "headers_check"
    "package"
    "package_source"
    "scripts"
    "test"
    "uninstall"
)

## @brief Names of recognized properties on targets.
#
# Unfortunately, the @c ARGV and @c ARGN arguments of a CMake function()
# or macro() does not preserve values which themselves are lists. Therefore,
# it is not possible to distinguish between property names and their values
# in the arguments passed to set_target_properties() or
# basis_set_target_properties(). To overcome this problem, this list specifies
# all the possible property names. Everything else is considered to be a
# property value except the first argument follwing right after the
# @c PROPERTIES keyword. Alternatively, basis_set_property() can be used
# as here no disambiguity exists.
#
# @note Placeholders such as &lt;CONFIG&gt; are allowed. These are treated
#       as the regular expression "[^ ]+". See basis_list_to_regex().
#
# @sa https://cmake.org/cmake/help/v3.4/manual/cmake-properties.7.html#properties-on-targets
set (BASIS_PROPERTIES_ON_TARGETS
  # CMake
  <CONFIG>_OUTPUT_NAME
  <CONFIG>_POSTFIX
  ALIASED_TARGET
  ANDROID_ANT_ADDITIONAL_OPTIONS
  ANDROID_API
  ANDROID_API_MIN
  ANDROID_ARCH
  ANDROID_ASSETS_DIRECTORIES
  ANDROID_GUI
  ANDROID_JAR_DEPENDENCIES
  ANDROID_JAR_DIRECTORIES
  ANDROID_JAVA_SOURCE_DIR
  ANDROID_NATIVE_LIB_DEPENDENCIES
  ANDROID_NATIVE_LIB_DIRECTORIES
  ANDROID_PROCESS_MAX
  ANDROID_PROGUARD
  ANDROID_PROGUARD_CONFIG_PATH
  ANDROID_SECURE_PROPS_PATH
  ANDROID_SKIP_ANT_STEP
  ANDROID_STL_TYPE
  ARCHIVE_OUTPUT_DIRECTORY
  ARCHIVE_OUTPUT_DIRECTORY_<CONFIG>
  ARCHIVE_OUTPUT_NAME
  ARCHIVE_OUTPUT_NAME_<CONFIG>
  AUTOGEN_TARGET_DEPENDS
  AUTOMOC
  AUTOUIC
  AUTOUIC_OPTIONS
  AUTORCC
  AUTORCC_OPTIONS
  BINARY_DIR
  BUILD_WITH_INSTALL_RPATH
  BUNDLE
  BUNDLE_EXTENSION
  C_EXTENSIONS
  C_STANDARD
  C_STANDARD_REQUIRED
  COMPATIBLE_INTERFACE_BOOL
  COMPATIBLE_INTERFACE_NUMBER_MAX
  COMPATIBLE_INTERFACE_NUMBER_MIN
  COMPATIBLE_INTERFACE_STRING
  COMPILE_DEFINITIONS
  COMPILE_DEFINITIONS_<CONFIG>
  COMPILE_FEATURES
  COMPILE_FLAGS
  COMPILE_OPTIONS
  COMPILE_PDB_NAME
  COMPILE_PDB_NAME_<CONFIG>
  <CONFIG>_OUTPUT_NAME
  <CONFIG>_PREFIX
  CROSSCOMPILING_EMULATOR
  CXX_EXTENSIONS
  CXX_STANDARD
  CXX_STANDARD_REQUIRED
  DEBUG_POSTFIX
  DEFINE_SYMBOL
  EchoString
  ENABLE_EXPORTS
  EXCLUDE_FROM_ALL
  EXCLUDE_FROM_DEFAULT_BUILD
  EXCLUDE_FROM_DEFAULT_BUILD_<CONFIG>
  EXPORT_NAME
  FOLDER
  Fortran_FORMAT
  Fortran_MODULE_DIRECTORY
  FRAMEWORK
  FRAMEWORK_VERSION
  GENERATOR_FILE_NAME
  GNUtoMS
  HAS_CXX
  IMPLICIT_DEPENDS_INCLUDE_TRANSFORM
  IMPORTED
  IMPORTED_CONFIGURATIONS
  IMPORTED_IMPLIB
  IMPORTED_IMPLIB_<CONFIG>
  IMPORTED_LINK_DEPENDENT_LIBRARIES
  IMPORTED_LINK_DEPENDENT_LIBRARIES_<CONFIG>
  IMPORTED_LINK_INTERFACE_LANGUAGES
  IMPORTED_LINK_INTERFACE_LANGUAGES_<CONFIG>
  IMPORTED_LINK_INTERFACE_LIBRARIES
  IMPORTED_LINK_INTERFACE_LIBRARIES_<CONFIG>
  IMPORTED_LINK_INTERFACE_MULTIPLICITY
  IMPORTED_LINK_INTERFACE_MULTIPLICITY_<CONFIG>
  IMPORTED_LOCATION
  IMPORTED_LOCATION_<CONFIG>
  IMPORTED_NO_SONAME
  IMPORTED_NO_SONAME_<CONFIG>
  IMPORTED_SONAME
  IMPORTED_SONAME_<CONFIG>
  IMPORT_PREFIX
  IMPORT_SUFFIX
  INCLUDE_DIRECTORIES
  INSTALL_NAME_DIR
  INSTALL_RPATH
  INSTALL_RPATH_USE_LINK_PATH
  INTERFACE_AUTOUIC_OPTIONS
  INTERFACE_COMPILE_DEFINITIONS
  INTERFACE_COMPILE_FEATURES
  INTERFACE_COMPILE_OPTIONS
  INTERFACE_INCLUDE_DIRECTORIES
  INTERFACE_LINK_LIBRARIES
  INTERFACE_POSITION_INDEPENDENT_CODE
  INTERFACE_SOURCES
  INTERFACE_SYSTEM_INCLUDE_DIRECTORIES
  INTERPROCEDURAL_OPTIMIZATION
  INTERPROCEDURAL_OPTIMIZATION_<CONFIG>
  JOB_POOL_COMPILE
  JOB_POOL_LINK
  LABELS
  <LANG>_COMPILER_LAUNCHER
  <LANG>_INCLUDE_WHAT_YOU_USE
  <LANG>_VISIBILITY_PRESET
  LIBRARY_OUTPUT_DIRECTORY
  LIBRARY_OUTPUT_DIRECTORY_<CONFIG>
  LIBRARY_OUTPUT_NAME
  LIBRARY_OUTPUT_NAME_<CONFIG>
  LINKER_LANGUAGE
  LINK_DEPENDS
  LINK_DEPENDS_NO_SHARED
  LINK_FLAGS
  LINK_FLAGS_<CONFIG>
  LINK_INTERFACE_LIBRARIES
  LINK_INTERFACE_LIBRARIES_<CONFIG>
  LINK_INTERFACE_MULTIPLICITY
  LINK_INTERFACE_MULTIPLICITY_<CONFIG>
  LINK_LIBRARIES
  LINK_SEARCH_END_STATIC
  LINK_SEARCH_START_STATIC
  LOCATION
  LOCATION_<CONFIG>
  MACOSX_BUNDLE
  MACOSX_BUNDLE_INFO_PLIST
  MACOSX_FRAMEWORK_INFO_PLIST
  MACOSX_RPATH
  MAP_IMPORTED_CONFIG_<CONFIG>
  NAME
  NO_SONAME
  NO_SYSTEM_FROM_IMPORTED
  OSX_ARCHITECTURES
  OSX_ARCHITECTURES_<CONFIG>
  OUTPUT_NAME
  OUTPUT_NAME_<CONFIG>
  PDB_NAME
  PDB_NAME_<CONFIG>
  PDB_OUTPUT_DIRECTORY
  PDB_OUTPUT_DIRECTORY_<CONFIG>
  POSITION_INDEPENDENT_CODE
  POST_INSTALL_SCRIPT
  PREFIX
  PRE_INSTALL_SCRIPT
  PRIVATE_HEADER
  PROJECT_LABEL
  PUBLIC_HEADER
  RESOURCE
  RULE_LAUNCH_COMPILE
  RULE_LAUNCH_CUSTOM
  RULE_LAUNCH_LINK
  RUNTIME_OUTPUT_DIRECTORY
  RUNTIME_OUTPUT_DIRECTORY_<CONFIG>
  RUNTIME_OUTPUT_NAME
  RUNTIME_OUTPUT_NAME_<CONFIG>
  SKIP_BUILD_RPATH
  SOURCE_DIR
  SOURCES
  SOVERSION
  STATIC_LIBRARY_FLAGS
  STATIC_LIBRARY_FLAGS_<CONFIG>
  SUFFIX
  TYPE
  VERSION
  VISIBILITY_INLINES_HIDDEN
  VS_DESKTOP_EXTENSIONS_VERSION
  VS_DOTNET_REFERENCES
  VS_DOTNET_TARGET_FRAMEWORK_VERSION
  VS_GLOBAL_KEYWORD
  VS_GLOBAL_PROJECT_TYPES
  VS_GLOBAL_ROOTNAMESPACE
  VS_GLOBAL_<variable>
  VS_IOT_EXTENSIONS_VERSION
  VS_IOT_STARTUP_TASK
  VS_KEYWORD
  VS_MOBILE_EXTENSIONS_VERSION
  VS_SCC_AUXPATH
  VS_SCC_LOCALPATH
  VS_SCC_PROJECTNAME
  VS_SCC_PROVIDER
  VS_WINDOWS_TARGET_PLATFORM_MIN_VERSION
  VS_WINRT_COMPONENT
  VS_WINRT_EXTENSIONS
  VS_WINRT_REFERENCES
  WIN32_EXECUTABLE
  WINDOWS_EXPORT_ALL_SYMBOLS
  XCODE_ATTRIBUTE_<an-attribute>
  XCTEST
  # BASIS
  BASIS_INCLUDE_DIRECTORIES    # include directories
  BASIS_LINK_DIRECTORIES       # link directories
  BASIS_LINK_DEPENDS           # link libraries
  BASIS_TYPE                   # BASIS type of target
  BASIS_UTILITIES              # whether BASIS utilities are used by this target
  BUNDLED                      # whether target belongs to same bundle/superbuild
  SCRIPT_DEFINITIONS           # script configuration code
  SCRIPT_DEFINITIONS_FILE      # script configuration file
  LANGUAGE                     # language of source files
  COMPILE                      # enable/disable compilation of script
  EXPORT                       # whether to export target
  LIBEXEC                      # whether the target is an auxiliary executable
  TEST                         # whether the target is a test
  MFILE                        # documentation file of MEX-file
  COMPONENT                    # package component of build target
  LIBRARY_COMPONENT            # package component of the library component
  RUNTIME_COMPONENT            # package component of the runtime component
  ARCHIVE_INSTALL_DIRECTORY    # installation directory of library
  LIBRARY_HEADER_DIRECTORY     # output directory of generated header files
  LIBRARY_INSTALL_DIRECTORY    # installation directory of library
  RUNTIME_INSTALL_DIRECTORY    # installation directory of runtime
  OUTPUT_DIRECTORY             # output directory for generated files
  INSTALL_DIRECTORY            # installation directory for generated files
  HTML_OUTPUT_DIRECTORY        # Doxygen/Sphinx HTML output directory
  HTML_INSTALL_DIRECTORY       # Doxygen/Sphinx HTML installation directory
  DIRHTML_OUTPUT_DIRECTORY     # Sphinx HTML output directory
  DIRHTML_INSTALL_DIRECTORY    # Sphinx HTML installation directory
  SINGLEHTML_OUTPUT_DIRECTORY  # Sphinx HTML output directory
  SINGLEHTML_INSTALL_DIRECTORY # Sphinx HTML installation directory
  LINKCHECK_OUTPUT_DIRECTORY   # Sphinx linkcheck output directory
  LINKCHECK_INSTALL_DIRECTORY  # Sphinx linkcheck installation directory
  XML_OUTPUT_DIRECTORY         # Doxygen XML output directory
  XML_INSTALL_DIRECTORY        # Doxygen XML installation directory
  MAN_OUTPUT_DIRECTORY         # Doxygen/Sphinx MAN output directory
  MAN_INSTALL_DIRECTORY        # Doxygen/Sphinx MAN installation directory
  TEXT_OUTPUT_DIRECTORY        # Sphinx text output directory
  TEXT_INSTALL_DIRECTORY       # Sphinx text installation directory
  TEXINFO_OUTPUT_DIRECTORY     # Sphinx Texinfo output directory
  TEXINFO_INSTALL_DIRECTORY    # Sphinx Texinfo installation directory
  LATEX_OUTPUT_DIRECTORY       # Doxygen/Sphinx LaTeX output directory
  LATEX_INSTALL_DIRECTORY      # Doxygen/Sphinx LaTeX installation directory
  PDF_OUTPUT_DIRECTORY         # Doxygen/Sphinx PDF output directory
  PDF_INSTALL_DIRECTORY        # Doxygen/Sphinx PDF installation directory
  RTF_OUTPUT_DIRECTORY         # Doxygen RTF output directory
  RTF_INSTALL_DIRECTORY        # Doxygen RTF installation directory
  DOXYFILE                     # Doxygen configuration file
  OUTPUT                       # Doxygen output formats
  TAGFILE                      # Doxygen tag file
  BUILD_DIRECTORY              # CMakeFiles build directory of target
  CONFIG_DIRECTORY             # Sphinx configuration directory
  BINARY_DIRECTORY             # CMake build tree directory
  SOURCE_DIRECTORY             # CMake or Sphinx source directory
  BUILDERS                     # Sphinx builders
)

# convert list of property names into regular expression
basis_list_to_regex (BASIS_PROPERTIES_ON_TARGETS_RE ${BASIS_PROPERTIES_ON_TARGETS})

## @brief Whether BASIS shall use target UIDs.
#
# If this option is OFF, target UIDs are idential to the target names
# given as arguments to the "basis_add_*" functions.
#
# The target UIDs ensure that no name conflict between the targets
# of this project and those of an external library which are imported
# occurs. Another reason for using these target UIDs is to avoid
# target name conflicts between modules or subprojects which may
# be developed by different teams.
#
# The downside of using target UIDs is, however, a slower configuration
# of the build system because every target name must be mapped to its
# target UID and possibly vice versa. Moreover, the use of target UIDs
# is less intuitive for those new to BASIS but experienced with CMake.
basis_set_if_not_set (BASIS_USE_TARGET_UIDS OFF)

## @brief Whether BASIS shall use fully-qualified target UIDs.
#
# If this option is OFF, the namespace of the top-level BASIS project is
# not prepended to the actual CMake build target names.
#
# For example, instead of the fully-qualified target UID
# "@PROJECT_NAME_L@.target", the CMake target name will simply
# be "target". Only when the target is referenced from another project,
# the fully-qualified target UID is usually required.
basis_set_if_not_set (BASIS_USE_FULLY_QUALIFIED_UIDS OFF)

## @brief Default component used for library targets when no component is specified.
#
# The default component a library target and its auxiliary files
# are associated with if no component was specified, explicitly.
set (BASIS_LIBRARY_COMPONENT "Development")

## @brief Default component used for executables when no component is specified.
#
# The default component an executable target and its auxiliary files
# are associated with if no component was specified, explicitly.
set (BASIS_RUNTIME_COMPONENT "Runtime")

## @brief Enable the automatic detection of the use of the BASIS utilities.
#
# If @c TRUE, the basis_add_executable() and basis_add_library() commands will try to
# automatically detect whether a given executable or library makes use of the
# BASIS utilities. If so, it configures the utilities for this project and adds
# a respective library build target as well as a link dependency on this target.
# This was the default until BASIS v3.1. Since this version, it is recommended
# to either use the @c USE_BASIS_UTILITIES option of basis_add_executable() and
# basis_add_library() or to add a link dependency on "basis" (recommended):
#
# @code
# basis_add_executable(foo.cxx)
# basis_target_link_libraries(foo basis)
# @endcode
set (BASIS_UTILITIES FALSE)

## @brief Whether to always build the BASIS C++ utilities even if not required by any target
option (BUILD_BASIS_UTILITIES_FOR_CXX    "Force the build of the BASIS C++ Utilities even if not used by this project" OFF)
## @brief Whether to always build the BASIS Python utilities even if not required by any target
option (BUILD_BASIS_UTILITIES_FOR_PYTHON "Force the build of the BASIS Python Utilities even if not used by this project" OFF)
## @brief Whether to always build the BASIS Perl utilities even if not required by any target
option (BUILD_BASIS_UTILITIES_FOR_PERL   "Force the build of the BASIS Perl Utilities even if not used by this project" OFF)
## @brief Whether to always build the BASIS Bash utilities even if not required by any target
option (BUILD_BASIS_UTILITIES_FOR_BASH   "Force the build of the BASIS Bash Utilities even if not used by this project" OFF)

mark_as_advanced (BUILD_BASIS_UTILITIES_FOR_CXX
                  BUILD_BASIS_UTILITIES_FOR_PYTHON
                  BUILD_BASIS_UTILITIES_FOR_PERL
                  BUILD_BASIS_UTILITIES_FOR_BASH)

## @brief Whether to export targets by default.
#
# This global variable specifies the default for the export of build targets if the
# @c EXPORT or @c NO_EXPORT options of the basis_add_executable and basis_add_library
# commands are not given.
set (BASIS_EXPORT_DEFAULT TRUE)

## @brief Suffix used for target exports file "<Package><ExportSuffix>.cmake"
set (BASIS_EXPORT_SUFFIX "Exports")

## @brief Whether to create "<Package><ExportSuffix>.cmake" file so other projects can import the exported targets.
#
# @sa GenerateConfig.cmake, ExportTools.cmake, http://www.cmake.org/cmake/help/v2.8.12/cmake.html#command:export
set (BASIS_EXPORT_ENABLED ON)

basis_is_cached (BASIS_DEPRECATED_CREATE_EXPORTS_FILE_OPTION BASIS_CREATE_EXPORTS_FILE)
if (BASIS_DEPRECATED_CREATE_EXPORTS_FILE_OPTION)
  set (BASIS_EXPORT_ENABLED "${BASIS_CREATE_EXPORTS_FILE}")
endif ()

## @brief Disable use of the revision information obtained from the revision
#         control software such as Subversion or Git.
#
# If this option is @c TRUE, the revision information is not included in the
# @c PROJECT_RELEASE information.
option (BASIS_REVISION_INFO "Enable use of the revision information of the revision control software." ON)
mark_as_advanced (BASIS_REVISION_INFO)

## @brief Enable compilation of scripts if supported by the language.
#
# In particular, Python modules are compiled if this option is enabled and
# only the compiled modules are installed.
#
# @sa basis_add_script_target()
option (BASIS_COMPILE_SCRIPTS "Enable compilation of scripts if supported by the language." OFF)
mark_as_advanced (BASIS_COMPILE_SCRIPTS)


## @brief Enable the installation of scripted modules in site specific default directories.
#
# If this option is @c ON, scripted modules such as Python and Perl modules, in particular,
# are installed in the default installation directories for site packages of the respective
# interpreter. This means that these modules may be installed outside the installation
# root directory as specified by the @c CMAKE_INSTALL_PREFIX. When this option is set to
# @c OFF, all modules are installed in subdirectories of the @c CMAKE_INSTALL_PREFIX.
# These directories may have to be added to the search path for modules manually as they
# might not be in the default search path of the respective interpreter.
#
# The installation directories for public modules which are intended for external use
# can further be set using the -D option of CMake or be modified by editing the respective
# advanced CMake cache variables named <tt>INSTALL_&lt;LANG%gt;_SITE_DIR</tt>.
#
# @note Even though it is more convenient for the user of a package to have the modules
#       being installed in the default directories, this option is set to @c OFF by default.
#       The reasons are that it is in first place expected that the installation will copy
#       files only to directories within the @c CMAKE_INSTALL_PREFIX. Moreover, it is not
#       guaranteed that the user has write permissions for the default site packages directories.
#       Last but not least, when installing public modules located in the @c PROJECT_LIBRARY_DIR
#       source directory, BASIS does not set a default module @c PREFIX which reduces the risk
#       of overwriting existing modules of other packages. If the developer of a BASIS package
#       was not thorough enough and did not follow the guidelines, setting this option to @c ON
#       has the potential risk of overwriting other packages' modules. Therefore,
#       modules are only installed in system default locations if explicitly requested.
option (BASIS_INSTALL_SITE_PACKAGES "Enable the installation of scripted modules in site specific default directories." OFF)
mark_as_advanced (BASIS_INSTALL_SITE_PACKAGES)

## @brief Script used to execute a process in CMake script mode.
#
# In order to be able to assign a timeout to the execution of a custom command
# and to add some error message parsing, this script is used by some build
# rules to actually perform the build step. See for example, the build of
# executables using the MATLAB Compiler.
set (BASIS_SCRIPT_EXECUTE_PROCESS "${BASIS_MODULE_PATH}/ExecuteProcess.cmake")

## @brief File used by default as <tt>--authors</tt> file to <tt>svn2cl</tt>.
#
# This file lists all Subversion users at SBIA and is used by default for
# the mapping of Subversion user names to real names during the generation
# of changelogs.
set (BASIS_SVN_USERS_FILE "${BASIS_MODULE_PATH}/SubversionUsers.txt")

## @brief Force installation of public header files of BASIS C++ utilities.
#
# If this variable is set to FALSE, each header file in the @c PROJECT_INCLUDE_DIR
# is scanned for an include statement which includes one of the public header
# files of the BASIS C++ utilities. If such include statement was found in
# a public header file of the project, the public header files of the BASIS
# C++ utilities are also installed as the project's public header files depend
# on them. You can set this variable to TRUE in the Settings.cmake file of your
# project to force the installation of the public header files of the
# project-specific BASIS C++ utilities.
#
# @sa basis_install_public_headers()
basis_set_if_not_set (BASIS_INSTALL_PUBLIC_HEADERS_OF_CXX_UTILITIES FALSE)

## @brief Whether BASIS should configure any public header file with the .in file name suffix.
#
# If a project does not contain any such public header file (typically one named config.h.in),
# this option can be set to @c FALSE in the @c "config/Settings.cmake" file of the project.
# For better performance, if only one header file needs to be configured, this can be done
# manually by adding a corresponding configure_file() call to the root CMakeLists.txt file
# right after basis_project_begin(). The configured files should be written to the
# @c BINARY_INCLUDE_DIR which is located in the build tree of the project.
set (BASIS_CONFIGURE_PUBLIC_HEADERS FALSE)

## @brief Whether basis_project_begin() should support the configuration of Slicer modules.
#
# This option must be set to @c TRUE in @c "config/Settings.cmake" of a project
# which either itself or one of its modules is a 3D Slicer Extension.
#
# @sa http://www.slicer.org
set (BASIS_SUPPORT_SLICER_MODULES FALSE)

## @brief Enable/disable registration of installed package in CMake registry.
option (BASIS_REGISTER "Request registration of installed package in CMake package registry." ON)
mark_as_advanced (BASIS_REGISTER)

## @brief EXPERIMENTAL - Build project modules as separate external projects.
#
# This may improve performance of the initial configure step but comes with the caveats
# inherent to the superbuild approach as implemented by the ExternalProject module.
option (BASIS_SUPERBUILD_MODULES "EXPERIMENTAL - Build project modules as part of a superbuild. May improve configure speed." OFF)
mark_as_advanced (BASIS_SUPERBUILD_MODULES)

## @brief Enable/disable import of targets using the BASIS ImportTools.
#
# Issues/limitations of CMake's "function" command complicate the definition
# of a custom set_target_properties function which can be used to collect
# "global" information about targets imported from the CMake package
# configuration of project dependencies. The workaround in the custom
# set_target_properties function defined in the ImportTools.cmake is
# extremely inefficient and slows down the configuration step a lot
# (cf. https://github.com/cmake-basis/BASIS/issues/494).
# 
# The only need for collecting this information for all (executable)
# targets imported from dependencies is for generating the executable
# target info table for the BASIS Utilities (cf. UtilitiesTools.cmake).
# Hence, when these are not used, the ImportTools.cmake are not needed.
# Further, when a project does not consist of modules, the imported
# targets are available in the scope of the project.
#
# A project has to set BASIS_IMPORT_TARGETS to TRUE in its root CMakeLists.txt
# file before basis_project_begin() or basis_project_impl(), respectively,
# which in turn include the BASISUse.cmake file. This file includes the
# BASIS ImportTools module when BASIS_IMPORT_TARGETS is true.
basis_set_if_not_set (BASIS_IMPORT_TARGETS FALSE)

# ============================================================================
# programming language specific settings
# ============================================================================

## @brief List of programming languages explicitly supported by BASIS.
#
# @todo Add full support for Java.
set (BASIS_LANGUAGES CMake CXX Python Jython Perl Matlab Bash)

string (TOLOWER "${BASIS_LANGUAGES}" BASIS_LANGUAGES_L)
string (TOUPPER "${BASIS_LANGUAGES}" BASIS_LANGUAGES_U)

# ----------------------------------------------------------------------------
# namespace delimiters
# ----------------------------------------------------------------------------

## @brief Namespace delimiter used in CMake.
set (BASIS_NAMESPACE_DELIMITER_CMAKE .)
## @brief Namespace delimiter used in C++.
set (BASIS_NAMESPACE_DELIMITER_CXX .)
## @brief Namespace delimiter used in Python.
set (BASIS_NAMESPACE_DELIMITER_PYTHON .)
## @brief Namespace delimiter used in Jython.
set (BASIS_NAMESPACE_DELIMITER_JYTHON .)
## @brief Namespace delimiter used in Perl.
set (BASIS_NAMESPACE_DELIMITER_PERL ::)
## @brief Namespace delimiter used in MATLAB.
set (BASIS_NAMESPACE_DELIMITER_MATLAB .)
## @brief Namespace delimiter used in Bash.
#
# @note Bash itself has no concept of namespaces. This namespace delimiter is
#       used by the import() function of the BASIS Utilities for Bash.
#
# @sa BasisBashUtilities
set (BASIS_NAMESPACE_DELIMITER_BASH .)

# ----------------------------------------------------------------------------
# public libraries of script modules
# ----------------------------------------------------------------------------

## @brief Name of library target which builds Python modules in @c PROJECT_LIBRARY_DIR.
#
# This variable is used by basis_configure_script_libraries() which is called
# by basis_project_begin() to add a library target of the given name for the
# build of the Python modules found in the @c PROJECT_LIBRARY_DIR.
#
# @note The given target name is argument to the basis_add_library() command.
#       Overwrite default value in Settings.cmake file of project if desired.
set (BASIS_PYTHON_LIBRARY_TARGET "pythonlib")

## @brief Name of library target which builds Jython modules in @c PROJECT_LIBRARY_DIR.
#
# This variable is used by basis_configure_script_libraries() which is called
# by basis_project_begin() to add a library target of the given name for the
# build of the Jython modules found in the @c PROJECT_LIBRARY_DIR.
#
# @note The given target name is argument to the basis_add_library() command.
#       Overwrite default value in Settings.cmake file of project if desired.
set (BASIS_JYTHON_LIBRARY_TARGET "jythonlib")

## @brief Name of library target which builds Perl modules in @c PROJECT_LIBRARY_DIR.
#
# This variable is used by basis_configure_script_libraries() which is called
# by basis_project_begin() to add a library target of the given name for the
# build of the Perl modules found in the @c PROJECT_LIBRARY_DIR.
#
# @note The given target name is argument to the basis_add_library() command.
#       Overwrite default value in Settings.cmake file of project if desired.
set (BASIS_PERL_LIBRARY_TARGET "perllib")

## @brief Name of library target which builds MATLAB modules in @c PROJECT_LIBRARY_DIR.
#
# This variable is used by basis_configure_script_libraries() which is called
# by basis_project_begin() to add a library target of the given name for the
# build of the MATLAB modules found in the @c PROJECT_LIBRARY_DIR.
#
# @note The given target name is argument to the basis_add_library() command.
#       Overwrite default value in Settings.cmake file of project if desired.
set (BASIS_MATLAB_LIBRARY_TARGET "matlablib")

## @brief Name of library target which builds Bash modules in @c PROJECT_LIBRARY_DIR.
#
# This variable is used by basis_configure_script_libraries() which is called
# by basis_project_begin() to add a library target of the given name for the
# build of the Bash modules found in the @c PROJECT_LIBRARY_DIR.
#
# @note The given target name is argument to the basis_add_library() command.
#       Overwrite default value in Settings.cmake file of project if desired.
set (BASIS_BASH_LIBRARY_TARGET "bashlib")

# ============================================================================
# documentation
# ============================================================================

## @brief Advanced non-cached variable to request build of documentation targets as part of ALL target.
basis_set_if_not_set (BASIS_ALL_DOC OFF)

## @brief Default Doxygen configuration.
set (BASIS_DOXYGEN_DOXYFILE "${CMAKE_CURRENT_LIST_DIR}/Doxyfile.in")

## @brief Default Sphinx configuration.
set (BASIS_SPHINX_CONFIG "${CMAKE_CURRENT_LIST_DIR}/sphinx_conf.py.in")

## @brief Default Sphinx theme.
set (BASIS_SPHINX_HTML_THEME "default")

## @brief Default Sphinx theme options.
set (BASIS_SPHINX_HTML_THEME_OPTIONS
  PROJECT_LOGO   None
  SHOW_SBIA_LOGO false
  SHOW_PENN_LOGO false
  #SHOW_RELBAR2   false
  #ROOT_RELLINKS  "[('home', 'index')]"
)

# ============================================================================
# common options
# ============================================================================

## @brief Request verbose messages from BASIS functions.
option (BASIS_VERBOSE "Request BASIS functions to output verbose messages." OFF)
mark_as_advanced (BASIS_VERBOSE)

## @brief Request debugging messages from BASIS functions.
option (BASIS_DEBUG "Request BASIS functions to help debugging." OFF)
mark_as_advanced (BASIS_DEBUG)

## @brief Request configuration of software build only, skipping steps related to packaging and installation.
option (BASIS_BUILD_ONLY "Request configuration of software build only, skipping steps related to packaging and installation." OFF)
mark_as_advanced (BASIS_BUILD_ONLY)

# ============================================================================
# build configuration
# ============================================================================

if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CXX_FLAGS AND NOT CMAKE_C_FLAGS)
  set (
    CMAKE_BUILD_TYPE "Release"
    CACHE STRING "Choose the type of build, options are: None (CMAKE_C_FLAGS and CMAKE_CXX_FLAGS used) Debug Release RelWithDebInfo MinSizeRel."
    FORCE
  )
endif ()

# the following Mac OS specific variables are yet not further used
# hide them from the normal user, as they are usually not required (yet)
if (DEFINED CMAKE_OSX_ARCHITECTURES)
  mark_as_advanced (CMAKE_OSX_ARCHITECTURES)
endif ()
if (DEFINED CMAKE_OSX_DEPLOYMENT_TARGET)
  mark_as_advanced (CMAKE_OSX_DEPLOYMENT_TARGET)
endif ()
if (DEFINED CMAKE_OSX_SYSROOT)
  mark_as_advanced (CMAKE_OSX_SYSROOT)
endif ()

## @brief Whether to have BASIS set the RPATH of binaries rather than CMake
#
# @sa http://www.cmake.org/Wiki/CMake_RPATH_handling for details on how CMake
#     itself handles the RPATH setting of executables and shared libraries.
option (BASIS_INSTALL_RPATH "Whether to have BASIS set the RPATH of binaries rather than CMake" ON)
mark_as_advanced (BASIS_INSTALL_RPATH)

# use INSTALL_RPATH set by BASIS instead of CMake
if (BASIS_INSTALL_RPATH)
  set (CMAKE_SKIP_RPATH                  FALSE) # use RPATH for installed project own binaries
  set (CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE) # do not add directories outside project to RPATH
endif ()
set (CMAKE_SKIP_BUILD_RPATH              FALSE) # use RPATH for project own binaries
set (CMAKE_BUILD_WITH_INSTALL_RPATH      FALSE) # use different RPATH for build tree executables


## @}
# end of Doxygen group
