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
# @brief Build configuration of MIRTK Transformation library.
##############################################################################

set(HEADERS
  ${BINARY_INCLUDE_DIR}/mirtk/TransformationExport.h
  AffineTransformation.h
  BSplineFreeFormTransformation3D.h
  BSplineFreeFormTransformation4D.h
  BSplineFreeFormTransformationStatistical.h
  BSplineFreeFormTransformationSV.h
  BSplineFreeFormTransformationTD.h
  ConstraintMeasure.h
  EnergyTerm.h
  FFDIntegrationMethod.h
  FluidFreeFormTransformation.h
  FreeFormTransformation.h
  FreeFormTransformation3D.h
  FreeFormTransformation4D.h
  HomogeneousTransformation.h
  HomogeneousTransformationIterator.h
  ImageTransformation.h
  InverseAffineTransformation.h
  JacobianConstraint.h
  LinearElasticityConstraint.h
  LinearFreeFormTransformation3D.h
  LinearFreeFormTransformation4D.h
  LinearFreeFormTransformationTD.h
  LogJacobianConstraint.h
  NegJacobianConstraint.h
  MultiLevelFreeFormTransformation.h
  MultiLevelStationaryVelocityTransformation.h
  MultiLevelTransformation.h
  MultipleVoxelTransformation.h
  PartialAffineTransformation.h
  PartialBSplineFreeFormTransformationSV.h
  PartialMultiLevelStationaryVelocityTransformation.h
  RegisteredImage.h
  RigidTransformation.h
  SimilarityTransformation.h
  SmoothnessConstraint.h
  SparsityConstraint.h
  TopologyPreservationConstraint.h
  Transformation.h
  TransformationApproximationError.h
  TransformationConfig.h
  TransformationConstraint.h
  TransformationJacobian.h
  TransformationModel.h
  Transformations.h
  VolumePreservationConstraint.h
)

set(SOURCES
  AffineTransformation.cc
  BSplineFreeFormTransformation3D.cc
  BSplineFreeFormTransformation4D.cc
  BSplineFreeFormTransformationStatistical.cc
  BSplineFreeFormTransformationSV.cc
  BSplineFreeFormTransformationTD.cc
  EnergyTerm.cc
  FluidFreeFormTransformation.cc
  FreeFormTransformation.cc
  FreeFormTransformation3D.cc
  FreeFormTransformation4D.cc
  FreeFormTransformationIntegration.h
  FreeFormTransformationRungeKutta.cc
  FreeFormTransformationRungeKutta.h
  HomogeneousTransformation.cc
  ImageTransformation.cc
  InverseAffineTransformation.cc
  JacobianConstraint.cc
  LinearElasticityConstraint.cc
  LinearFreeFormTransformation3D.cc
  LinearFreeFormTransformation4D.cc
  LinearFreeFormTransformationTD.cc
  LogJacobianConstraint.cc
  NegJacobianConstraint.cc
  MultiLevelFreeFormTransformation.cc
  MultiLevelStationaryVelocityTransformation.cc
  MultiLevelTransformation.cc
  PartialAffineTransformation.cc
  PartialBSplineFreeFormTransformationSV.cc
  PartialMultiLevelStationaryVelocityTransformation.cc
  RegisteredImage.cc
  RigidTransformation.cc
  SimilarityTransformation.cc
  SmoothnessConstraint.cc
  SparsityConstraint.cc
  TopologyPreservationConstraint.cc
  Transformation.cc
  TransformationApproximationError.cc
  TransformationConfig.cc
  TransformationConstraint.cc
  TransformationInverse.cc
  TransformationUtils.h
  VolumePreservationConstraint.cc
)

set(DEPENDS
  LibCommon
  LibNumerics
  LibImage
)

if (VTK_FOUND)
  list(APPEND DEPENDS vtkCommonCore vtkCommonDataModel)
endif ()

if (TARGET TBB::tbb)
  list(APPEND DEPENDS TBB::tbb)
endif ()

mirtk_add_library(AUTO_REGISTER)
