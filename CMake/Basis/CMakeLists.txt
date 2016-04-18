# ============================================================================
# Copyright (c) 2011-2012 University of Pennsylvania
# Copyright (c) 2013-2016 Andreas Schuh
# All rights reserved.
#
# See COPYING file for license information or visit
# https://cmake-basis.github.io/download.html#license
# ============================================================================

##############################################################################
# @file  CMakeLists.txt
# @brief Build configuration of CMake BASIS Modules.
##############################################################################

if (NOT INSTALL_BASIS_MODULES_DIR)
  return ()
endif ()

# ----------------------------------------------------------------------------
# install files
set (CMAKE_MODULES
  # main modules
  "BasisPack.cmake"               # packaging, wraps CPack
  "BasisTools.cmake"              # definition of macros and functions
  "BasisTest.cmake"               # testing, wraps CTest
  "BasisSettings.cmake"           # default project-independent settings
  # submodules included by BasisTools.cmake or on demand, respectively
  "CommonTools.cmake"             # common functions/macros
  "DirectoriesSettings.cmake"     # project directory structure
  "DocTools.cmake"                # documentation generation
  "ExportTools.cmake"             # exporting targets
  "ExternalData.cmake"            # storing data externally
  "ExternalData_config.cmake.in"  # config of ExternalData.cmake module
  "ImportTools.cmake"             # importing targets
  "InterpTools.cmake"             # script interpreter tools
  "InstallationTools.cmake"       # software installation
  "MatlabTools.cmake"             # MATLAB support
  "ProjectTools.cmake"            # main project/module helpers
  "RevisionTools.cmake"           # revision control system support
  "SlicerTools.cmake"             # support build of Slicer Extensions
  "SuperBuildTools.cmake"         # super-build related functions
  "TargetTools.cmake"             # addition of build target
  "TopologicalSort.cmake"         # topological sorting of build dependencies
  "UtilitiesTools.cmake"          # generation/configuration of utilities
)

# ----------------------------------------------------------------------------
# auxiliary files
set (CMAKE_MODULES_AUXILIARY_FILES
  # API documentation
  "Doxyfile.in"                         # default Doxygen configuration template
  "Modules.dox"                         # documentation of BASIS modules groups
  "Utilities.dox"                       # documentation of BASIS Utilities group
  "CxxUtilities.dox"                    # documentation of C++ Utilities group
  "JavaUtilities.dox"                   # documentation of Java Utilities group
  "PythonUtilities.dox"                 # documentation of Python Utilities group
  "PerlUtilities.dox"                   # documentation of Perl Utilities group
  "BashUtilities.dox"                   # documentation of BASH Utilities group
  "MatlabUtilities.dox"                 # documentation of MATLAB Utilities group
  # package configuration
  "Config.cmake.in"                     # default <Pkg>Config.cmake template
  "ModuleConfig.cmake.in"               # default <Pkg>Config.cmake template for modules
  "BasisConfigSettings.cmake"           # settings for default <Pkg>Config.cmake
  "ConfigVersion.cmake.in"              # default <Pkg>ConfigVersion.cmake template
  "ConfigUse.cmake.in"                  # default <Pkg>Use.cmake template
  "ModuleConfigUse.cmake.in"            # default <Pkg>Use.cmake template for modules
  "GenerateConfig.cmake"                # generation/configuration of <Pkg>Config.cmake et al.
  # documentation generation
  "doxygen_header.html.in"              # default HTML header of Doxygen documentation
  "doxygen_footer.html.in"              # default HTML footer of Doxygen documentation
  "doxygen_extra.css.in"                # default HTML extra stylesheet of Doxygen documentation
  "sphinx_conf.py.in"                   # default Sphinx configuration
  "sphinx_make.sh.in"                   # custom "make" script for PDF generation
  # other
  "buildtimestamp.cmd"                  # get build timestamp on Windows
  "CheckPublicHeaders.cmake"            # check if header files were added/removed
  "cmake_uninstall.cmake.in"            # uninstall script
  "configure_script.cmake.in"           # used to configure ("build") script files
  "ConfigureIncludeFiles.cmake"         # used to configure public header files
  "Directories.cmake.in"                # documents project directory structure
  "doxyfilter.bat.in"                   # wraps Doxygen filter command on Windows
  "ExecuteProcess.cmake"                # wraps CMake's execute_process() function
  "generate_matlab_executable.cmake.in" # build command for alternative build of MCC targets
  "GenerateExportHeader.cmake"          # GenerateExportHeader module of CMake 3.6
  "exportheader.cmake.in"               # template header file used by GenerateExportHeader
  "get_python_lib.py"                   # auxiliary script to determine site-packages directory
  "glob.cmake"                          # used to glob for source files
  "ScriptConfig.cmake.in"               # default script configuration file
  "PostprocessChangeLog.cmake"          # used to postprocess generated ChangeLog
  "ProjectSettings.cmake.in"            # default project-dependent settings
)

basis_install (
  FILES       ${CMAKE_MODULES} ${CMAKE_MODULES_AUXILIARY_FILES}
  DESTINATION "${INSTALL_BASIS_MODULES_DIR}"
)

# Attention: DO NOT add the uninstall subdirectory here! Instead, the
#            CMakeLists.txt file is installed in a separate directory.
#            This directory is added to the build configuration of each
#            BASIS project (including BASIS) at the end of basis_project_end().
basis_install (
  FILES       "uninstall/CMakeLists.txt"
  DESTINATION "${INSTALL_BASIS_MODULES_DIR}/uninstall"
)

# ----------------------------------------------------------------------------
# make sure that we did not forget to install a file... very easy to forget!
file (GLOB MODULES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" *)
foreach (F IN LISTS MODULES)
  if (NOT F MATCHES "^(CMakeLists\\.txt|COPYING\\.txt|README\\.md|\\..*|.*\\.bak|.*\\.autosave|.*\\.pyc|.*~|uninstall)$")
    list (FIND CMAKE_MODULES "${F}" IDX)
    if (IDX EQUAL -1)
      list (FIND CMAKE_MODULES_AUXILIARY_FILES "${F}" IDX)
    endif ()
    if (IDX EQUAL -1)
      message (WARNING "Seems like ${CMAKE_CURRENT_SOURCE_DIR}/${F} will not be installed!")
    endif ()
  endif ()
endforeach ()
