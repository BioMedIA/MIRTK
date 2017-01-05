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
# This CMake configuration file adds installation rules for the MIRTK CMake
# modules and commonly used CMake BASIS Find modules such that MIRTK Packages
# can be built separately outside the MIRTK top-level project.
#
# CMake BASIS Find Modules: https://github.com/cmake-basis/find-modules.git
##############################################################################

basis_install(DESTINATION "${INSTALL_CMAKE_MODULES_DIR}" FILES
  # MIRTK modules
  "mirtkTools.cmake"              # includes MIRTK command definitions
  "mirtkPolicies.cmake"           # common CMake policies
  "mirtkAddLibrary.cmake"         # defines mirtk_add_library
  "mirtkAddExecutable.cmake"      # defines mirtk_add_executable
  "mirtkAddTest.cmake"            # defines mirtk_add_test
  "mirtkTargetDependencies.cmake" # defines mirtk_target_dependencies
  "mirtkGetTargetName.cmake"      # defines mirtk_get_target_name
  "mirtkProjectBegin.cmake"       # defines mirtk_project_begin
  "mirtkProjectEnd.cmake"         # defines mirtk_project_end
  "mirtkConfigureModule.cmake"    # defines mirtk_configure_module
  # Find modules
  "FindARPACK.cmake"
  "FindEigen3.cmake"
  "FindFLANN.cmake"
  "FindMATLAB.cmake"
  "FindPython.cmake"
  "FindPythonInterp.cmake"
  "FindSphinx.cmake"
  "FindSuiteSparse.cmake"
  "FindTBB.cmake"
  "FindUMFPACK.cmake"
)
