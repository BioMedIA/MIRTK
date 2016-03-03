# ============================================================================
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
# ============================================================================

##############################################################################
# @file  CMakeLists.txt
# @brief Build configuration of CMake modules and utilities.
#
# This CMake configuration file adds installation rules for the CMake BASIS
# Modules and additional MIRTK CMake modules such that MIRTK Packages can be
# built separately outside the MIRTK top-level project. Note that the
# CMake/Basis directory must contain all CMake BASIS Modules files even if
# not all are used. This is to enable the use of the "git subtree pull" and
# "git subtree split" commands to update these files after an upstream change.
# The CMake/Config directory contains the CMake BASIS configuration files for
# the MIRTK top-level project.
#
# CMake BASIS Modules: https://github.com/schuhschuh/cmake-basis-modules.git
##############################################################################

# ----------------------------------------------------------------------------
# CMake BASIS Modules
set (
  BASIS_MODULES
    # main modules
    "BasisPack.cmake"               # packaging, wraps CPack
    "BasisTools.cmake"              # definition of macros and functions
    "BasisTest.cmake"               # testing, wraps CTest
    "BasisSettings.cmake"           # default project-independent settings
    # submodules included by BasisTools.cmake
    "CommonTools.cmake"             # common functions/macros
    "DirectoriesSettings.cmake"     # project directory structure
    "DocTools.cmake"                # documentation generation
    "ExportTools.cmake"             # exporting targets
    "ExternalData.cmake"            # storing data externally
    "ExternalData_config.cmake.in"  # config of ExternalData.cmake module
    "InterpTools.cmake"             # script interpreter tools
    "InstallationTools.cmake"       # software installation
    "MatlabTools.cmake"             # MATLAB support
    "ProjectTools.cmake"            # main project/module helpers
    "RevisionTools.cmake"           # revision control system support
    "SlicerTools.cmake"             # support build of Slicer Extensions
    "TargetTools.cmake"             # addition of build target
    "TopologicalSort.cmake"         # topological sorting of build dependencies
    "UtilitiesTools.cmake"          # generation/configuration of utilities
)

# ----------------------------------------------------------------------------
# Auxiliary files of CMake BASIS Modules
set (
  BASIS_MODULES_AUXILIARY_FILES
    "Config.cmake.in"               # default <Pkg>Config.cmake template
    "ModuleConfig.cmake.in"         # default <Pkg>Config.cmake template for modules
    "BasisConfigSettings.cmake"     # settings for default <Pkg>Config.cmake
    "ConfigVersion.cmake.in"        # default <Pkg>ConfigVersion.cmake template
    "ConfigUse.cmake.in"            # default <Pkg>Use.cmake template
    "ModuleConfigUse.cmake.in"      # default <Pkg>Use.cmake template for modules
    "GenerateConfig.cmake"          # generation/configuration of <Pkg>Config.cmake et al.
    "buildtimestamp.cmd"            # get build timestamp on Windows
    "CheckPublicHeaders.cmake"      # check if header files were added/removed
    "cmake_uninstall.cmake.in"      # uninstall script
    "configure_script.cmake.in"     # used to configure ("build") script files
    "ConfigureIncludeFiles.cmake"   # used to configure public header files
    "Directories.cmake.in"          # documents project directory structure
    "doxyfilter.bat.in"             # wraps Doxygen filter command on Windows
    "ExecuteProcess.cmake"          # wraps CMake's execute_process() function
    "get_python_lib.py"             # auxiliary script to determine site-packages directory
    "glob.cmake"                    # used to glob for source files
    "ScriptConfig.cmake.in"         # default script configuration file
    "PostprocessChangeLog.cmake"    # used to postprocess generated ChangeLog
    "ProjectSettings.cmake.in"      # default project-dependent settings
)

# ----------------------------------------------------------------------------
# CMake BASIS Find<Pkg>.cmake modules
set (
  BASIS_FIND_MODULES
    "FindAFNI.cmake"
    "FindBASH.cmake"
    "FindBLAS.cmake"
    "FindBoostNumericBindings.cmake"
    "FindITK.cmake"
    "FindJythonInterp.cmake"
    "FindLIBLINEAR.cmake"
    "FindMATLAB.cmake"
    "FindMatlabNiftiTools.cmake"
    "FindMOSEK.cmake"
    "FindNiftiCLib.cmake"
    "FindOpenCV.cmake"
    "FindPerl.cmake"
    "FindPerlLibs.cmake"
    "FindPythonInterp.cmake"
    "FindPythonModules.cmake"
    "FindSparseBayes.cmake"
    "FindSphinx.cmake"
    "FindSVMTorch.cmake"
    "FindWeka.cmake"
)

# ----------------------------------------------------------------------------
# MIRTK CMake modules
set (
  CMAKE_MODULES
    "MIRTKTools.cmake"           # includes MIRTK command definitions
    "MIRTKAddLibrary.cmake"      # defines mirtk_add_library
    "MIRTKAddExecutable.cmake"   # defines mirtk_add_executable
    "MIRTKAddTest.cmake"         # defines mirtk_add_test
    "MIRTKProjectBegin.cmake"    # defines mirtk_project_begin
    "MIRTKProjectEnd.cmake"      # defines mirtk_project_end
    "MIRTKConfigureModule.cmake" # defines mirtk_configure_module
)

# ----------------------------------------------------------------------------
# MIRTK Find<Pkg>.cmake modules
set (
  CMAKE_FIND_MODULES
)

# ----------------------------------------------------------------------------
# install CMake files
set(files)
foreach (fname IN LISTS BASIS_MODULES BASIS_FIND_MODULES BASIS_MODULES_AUXILIARY_FILES)
  list(APPEND files "Basis/${fname}")
endforeach ()
foreach (fname IN LISTS CMAKE_MODULES CMAKE_FIND_MODULES)
  list(APPEND files "Modules/${fname}")
endforeach ()
basis_install(FILES ${files} DESTINATION "${INSTALL_CMAKE_MODULES_DIR}")

# Attention: DO NOT add the uninstall subdirectory here! Instead, the
#            CMakeLists.txt file is installed in a separate directory.
#            This directory is added to the build configuration of each
#            BASIS project (including BASIS) in basis_project_end().
basis_install (
  FILES       "Basis/uninstall/CMakeLists.txt"
  DESTINATION "${INSTALL_CMAKE_MODULES_DIR}/uninstall"
)
