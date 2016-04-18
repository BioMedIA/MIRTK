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
# CMake BASIS Modules: https://github.com/cmake-basis/modules.git
##############################################################################

set(INSTALL_BASIS_MODULES_DIR "${INSTALL_CMAKE_MODULES_DIR}")

add_subdirectory(Basis)
add_subdirectory(Modules)
