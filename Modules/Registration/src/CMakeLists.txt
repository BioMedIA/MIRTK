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
  ${BINARY_INCLUDE_DIR}/mirtk/RegistrationConfig.h
  ${BINARY_INCLUDE_DIR}/mirtk/RegistrationExport.h
  CosineOfNormalizedGradientField.h
  DataFidelity.h
  GenericRegistrationDebugger.h
  GenericRegistrationFilter.h
  GenericRegistrationLogger.h
  GradientFieldSimilarity.h
  HistogramImageSimilarity.h
  ImageCovariance.h
  ImageSimilarity.h
  IntensityCorrelationRatioXY.h
  IntensityCorrelationRatioYX.h
  IntensityCrossCorrelation.h
  JointImageEntropy.h
  LabelConsistency.h
  MeanSquaredDisplacementError.h
  MutualImageInformation.h
  NormalizedGradientFieldSimilarity.h
  NormalizedIntensityCrossCorrelation.h
  NormalizedMutualImageInformation.h
  PeakSignalToNoiseRatio.h
  Registration.h
  RegistrationEnergy.h
  RegistrationFilter.h
  SimilarityMeasure.h
  SumOfSquaredIntensityDifferences.h
)

set(SOURCES
  CosineOfNormalizedGradientField.cc
  DataFidelity.cc
  GenericRegistrationDebugger.cc
  GenericRegistrationFilter.cc
  GenericRegistrationLogger.cc
  GradientFieldSimilarity.cc
  HistogramImageSimilarity.cc
  ImageCovariance.cc
  ImageSimilarity.cc
  IntensityCorrelationRatioXY.cc
  IntensityCorrelationRatioYX.cc
  IntensityCrossCorrelation.cc
  JointImageEntropy.cc
  LabelConsistency.cc
  MeanSquaredDisplacementError.cc
  MutualImageInformation.cc
  NormalizedGradientFieldSimilarity.cc
  NormalizedIntensityCrossCorrelation.cc
  NormalizedMutualImageInformation.cc
  PeakSignalToNoiseRatio.cc
  RegistrationConfig.cc
  RegistrationEnergy.cc
  RegistrationEnergyParser.h
  SumOfSquaredIntensityDifferences.cc
)

set(DEPENDS
  LibCommon
  LibNumerics
  LibImage
  LibTransformation
)

if (TARGET LibPointSet AND VTK_FOUND)
  list(APPEND HEADERS
    CurrentsDistance.h
    FiducialRegistrationError.h
    PointCorrespondenceDistance.h
    PointSetDistance.h
    PointSetDistanceMeasure.h
    SurfaceDistance.h
  )
  list(APPEND SOURCES
    CurrentsDistance.cc
    FiducialRegistrationError.cc
    PointCorrespondenceDistance.cc
    PointSetDistance.cc
    SurfaceDistance.cc
  )
  list(APPEND DEPENDS LibPointSet ${VTK_LIBRARIES})
  if (TARGET LibDeformable)
    list(APPEND DEPENDS LibDeformable)
  endif ()
endif ()

if (TARGET TBB::tbb)
  list(APPEND DEPENDS TBB::tbb)
endif ()

mirtk_add_library(AUTO_REGISTER)
