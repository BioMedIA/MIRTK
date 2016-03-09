# ============================================================================
# Medical Image Registration ToolKit (MIRTK)
#
# Copyright <copyright>
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
# @brief Build configuration of MIRTK <project> module.
##############################################################################

set(HEADERS
  #mirtkClassName.h
)

set(SOURCES
  #mirtkClassName.cc
)

set(DEPENDS
  LibCommon
  LibNumerics
  LibImage
)
if (MIRTK_PointSet_FOUND)
  list(APPEND DEPENDS LibPointSet)
endif ()

mirtk_add_library(AUTO_REGISTER)
