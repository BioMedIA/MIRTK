# ============================================================================
# Medical Image Registration ToolKit (MIRTK)
#
# Copyright 2013-2017 Imperial College London
# Copyright 2013-2018 Andreas Schuh
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

#################################################################################
# @file  BasisProject.cmake
# @brief Sets basic information about a BASIS Project and calls basis_project().
#
# This file defines basic information about a project by calling 
# the basis_project() function. This basic information, also known as metadata, 
# is used by BASIS to setup the project. Moreover, if the project is a module 
# of another BASIS project, the dependencies to other modules have to be specified 
# here such that the top-level project can analyze the inter-module dependencies.
#
# @sa https://cmake-basis.github.io/standard/modules.html
#
# However, not only dependencies to other modules can be specified here,
# but also dependencies on external packages. A more flexible alternative to
# resolve external dependencies is to add the corresponding basis_find_package()
# statements to the Depends.cmake file. This should, however, only be done
# if specifying the dependencies as arguments to the basis_project() function
# cannot be used to resolve the dependencies properly. If you only need to
# make use of additional variables set by the package configuration file
# of the external package or the corresponding Find<Package>.cmake module,
# add the related CMake code to the Settings.cmake file instead.
#
# @ingroup BasisSettings
##############################################################################

# Note: The #<*> dependency patterns are required by the basisproject tool and
#       should be kept on a separate line as last commented argument of the
#       corresponding options of the basis_project() command. The TEMPLATE
#       option and set argument are also required by this tool and should not
#       be changed manually. The argument is updated by basisproject --update.

basis_project (
  # --------------------------------------------------------------------------
  # meta-data
  NAME         "MIRTK"
  VERSION      "0.0.0" # Version of core (+ external) modules
  SOVERSION    "0"     # API yet unstable
  AUTHORS      "Andreas Schuh"
  DESCRIPTION  "Medical Image Registration ToolKit"
  COPYRIGHT    "2013-2018 Imperial College London, Andreas Schuh"
  LICENSE      "Apache License Version 2.0"
  CONTACT      "Andreas Schuh <andreas.schuh.84@gmail.com>"
  TEMPLATE     "with-basis-submodule/1.0"
  LANGUAGES    C CXX-11
 
  # --------------------------------------------------------------------------
  # directories
  CONFIG_DIR  CMake/Config
  DOC_DIR     Documentation
  MODULES_DIR Modules
  TOOLS_DIR   Applications
  OTHER_DIRS  CMake

  MODULE_DIRS
    Packages/Deformable
    Packages/Mapping
    Packages/DrawEM
    Packages/Scripting
    Packages/Viewer

  # --------------------------------------------------------------------------
  # list of modules enabled by default
  DEFAULT_MODULES
    Common
    Numerics
    Image
    IO
    Transformation
    Registration

  # --------------------------------------------------------------------------
  # list of external modules
  EXTERNAL_MODULES
    Deformable
    Mapping
    DrawEM
    Scripting
    Viewer

  # --------------------------------------------------------------------------
  # dependencies
  DEPENDS
    #<dependency>
  OPTIONAL_DEPENDS
    #<optional-dependency>
  TOOLS_DEPENDS
    Python{Interp}
  TEST_DEPENDS
    #<test-dependency>
  OPTIONAL_TEST_DEPENDS
    #<optional-test-dependency>
)
