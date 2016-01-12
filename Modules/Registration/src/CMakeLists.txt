# ============================================================================
# Medical Image Registration ToolKit (MIRTK)
#
# Copyright 2013-2015 Imperial College London
# Copyright 2013-2015 Andreas Schuh
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
# @brief Build configuration of MIRTK Registration library.
##############################################################################

set(HEADERS
  ${BINARY_INCLUDE_DIR}/mirtkRegistrationConfig.h
  mirtkCosineOfNormalizedGradientField.h
  mirtkDataFidelity.h
  mirtkGenericRegistrationDebugger.h
  mirtkGenericRegistrationFilter.h
  mirtkGenericRegistrationLogger.h
  mirtkGradientFieldSimilarity.h
  mirtkImageSimilarity.h
  mirtkIntensityCrossCorrelation.h
  mirtkMutualImageInformation.h
  mirtkNormalizedGradientFieldSimilarity.h
  mirtkNormalizedIntensityCrossCorrelation.h
  mirtkNormalizedMutualImageInformation.h
  mirtkProbabilisticImageSimilarity.h
  mirtkRegistration.h
  mirtkRegistrationEnergy.h
  mirtkRegistrationFilter.h
  mirtkSimilarityMeasure.h
  mirtkSumOfSquaredIntensityDifferences.h
)

set(SOURCES
  mirtkCosineOfNormalizedGradientField.cc
  mirtkDataFidelity.cc
  mirtkGenericRegistrationDebugger.cc
  mirtkGenericRegistrationFilter.cc
  mirtkGenericRegistrationLogger.cc
  mirtkGradientFieldSimilarity.cc
  mirtkImageSimilarity.cc
  mirtkIntensityCrossCorrelation.cc
  mirtkMutualImageInformation.cc
  mirtkNormalizedGradientFieldSimilarity.cc
  mirtkNormalizedIntensityCrossCorrelation.cc
  mirtkNormalizedMutualImageInformation.cc
  mirtkProbabilisticImageSimilarity.cc
  mirtkRegistrationConfig.cc
  mirtkRegistrationEnergy.cc
  mirtkRegistrationEnergyParser.h
  mirtkSumOfSquaredIntensityDifferences.cc
)

set(DEPENDS
  LibCommon
  LibNumerics
  LibImage
  LibTransformation
)

if (TARGET LibPointSet AND VTK_FOUND)
  list(APPEND HEADERS
    mirtkCurrentsDistance.h
    mirtkFiducialRegistrationError.h
    mirtkPointCorrespondenceDistance.h
    mirtkPointSetDistance.h
    mirtkPointSetDistanceMeasure.h
    mirtkSurfaceDistance.h
  )
  list(APPEND SOURCES
    mirtkCurrentsDistance.cc
    mirtkFiducialRegistrationError.cc
    mirtkPointCorrespondenceDistance.cc
    mirtkPointSetDistance.cc
    mirtkSurfaceDistance.cc
  )
  list(APPEND DEPENDS LibPointSet ${VTK_LIBRARIES})
  if (TARGET LibDeformable)
    list(APPEND DEPENDS LibDeformable)
  endif ()
endif ()

mirtk_add_library(AUTO_REGISTER)
