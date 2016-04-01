# ==============================================================================
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
# ==============================================================================

##############################################################################
# @file  ConfigSettings.cmake
# @brief Sets variables used in CMake package configuration.
#
# It is suggested to use @c _CONFIG as suffix for variable names that are to
# be substituted in the Config.cmake.in template file in order to distinguish
# these variables from the build configuration.
#
# @note The default BasisConfigSettings.cmake file which is part of the BASIS
#       installation is included prior to this file. Hence, the variables are
#       valid even if a custom project-specific configuration is used and
#       default values can further be overwritten in this file.
#
# @ingroup BasisSettings
##############################################################################

# ============================================================================
# build tree configuration settings
# ============================================================================

if (BUILD_CONFIG_SETTINGS)
  set (BASIS_MODULE_PATH_CONFIG "${BASIS_MODULE_PATH}")
  set (MIRTK_MODULE_PATH_CONFIG "${MIRTK_MODULE_PATH}")
  set (MIRTK_TOOLS_DIR_CONFIG  "${MIRTK_TOOLS_PATH}")
  return ()
endif ()

# ============================================================================
# installation configuration settings
# ============================================================================

## @brief Directory of CMake BASIS Modules.
set (BASIS_MODULE_PATH_CONFIG "\${\${NS}INSTALL_PREFIX}/${INSTALL_CMAKE_MODULES_DIR}")
## @brief Directory of MIRTK CMake Modules.
set (MIRTK_MODULE_PATH_CONFIG "\${\${NS}INSTALL_PREFIX}/${INSTALL_CMAKE_MODULES_DIR}")
## @brief Directory of MIRTK command executables.
set (MIRTK_TOOLS_DIR_CONFIG  "\${\${NS}INSTALL_PREFIX}/${INSTALL_TOOLS_DIR}")
