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
  mirtkAffineTransformation.h
  mirtkBSplineFreeFormTransformation3D.h
  mirtkBSplineFreeFormTransformation4D.h
  mirtkBSplineFreeFormTransformationStatistical.h
  mirtkBSplineFreeFormTransformationSV.h
  mirtkBSplineFreeFormTransformationTD.h
  mirtkConstraintMeasure.h
  mirtkEnergyTerm.h
  mirtkFFDIntegrationMethod.h
  mirtkFluidFreeFormTransformation.h
  mirtkFreeFormTransformation.h
  mirtkFreeFormTransformation3D.h
  mirtkFreeFormTransformation4D.h
  mirtkHomogeneousTransformation.h
  mirtkHomogeneousTransformationIterator.h
  mirtkImageTransformation.h
  mirtkInverseAffineTransformation.h
  mirtkJacobianConstraint.h
  mirtkLinearFreeFormTransformation3D.h
  mirtkLinearFreeFormTransformation4D.h
  mirtkLinearFreeFormTransformationTD.h
  mirtkLogJacobianConstraint.h
  mirtkMinJacobianConstraint.h
  mirtkMultiLevelFreeFormTransformation.h
  mirtkMultiLevelStationaryVelocityTransformation.h
  mirtkMultiLevelTransformation.h
  mirtkMultipleVoxelTransformation.h
  mirtkPartialAffineTransformation.h
  mirtkPartialBSplineFreeFormTransformationSV.h
  mirtkPartialMultiLevelStationaryVelocityTransformation.h
  mirtkRegisteredImage.h
  mirtkRigidTransformation.h
  mirtkSimilarityTransformation.h
  mirtkSmoothnessConstraint.h
  mirtkSparsityConstraint.h
  mirtkTopologyPreservationConstraint.h
  mirtkTransformation.h
  mirtkTransformationApproximationError.h
  mirtkTransformationConfig.h
  mirtkTransformationConstraint.h
  mirtkTransformationJacobian.h
  mirtkTransformationModel.h
  mirtkTransformations.h
  mirtkVolumePreservationConstraint.h
)

set(SOURCES
  mirtkAffineTransformation.cc
  mirtkBSplineFreeFormTransformation3D.cc
  mirtkBSplineFreeFormTransformation4D.cc
  mirtkBSplineFreeFormTransformationStatistical.cc
  mirtkBSplineFreeFormTransformationSV.cc
  mirtkBSplineFreeFormTransformationTD.cc
  mirtkEnergyTerm.cc
  mirtkFluidFreeFormTransformation.cc
  mirtkFreeFormTransformation.cc
  mirtkFreeFormTransformation3D.cc
  mirtkFreeFormTransformation4D.cc
  mirtkFreeFormTransformationIntegration.h
  mirtkFreeFormTransformationRungeKutta.cc
  mirtkFreeFormTransformationRungeKutta.h
  mirtkHomogeneousTransformation.cc
  mirtkImageTransformation.cc
  mirtkInverseAffineTransformation.cc
  mirtkJacobianConstraint.cc
  mirtkLinearFreeFormTransformation3D.cc
  mirtkLinearFreeFormTransformation4D.cc
  mirtkLinearFreeFormTransformationTD.cc
  mirtkLogJacobianConstraint.cc
  mirtkMinJacobianConstraint.cc
  mirtkMultiLevelFreeFormTransformation.cc
  mirtkMultiLevelStationaryVelocityTransformation.cc
  mirtkMultiLevelTransformation.cc
  mirtkPartialAffineTransformation.cc
  mirtkPartialBSplineFreeFormTransformationSV.cc
  mirtkPartialMultiLevelStationaryVelocityTransformation.cc
  mirtkRegisteredImage.cc
  mirtkRigidTransformation.cc
  mirtkSimilarityTransformation.cc
  mirtkSmoothnessConstraint.cc
  mirtkSparsityConstraint.cc
  mirtkTopologyPreservationConstraint.cc
  mirtkTransformation.cc
  mirtkTransformationApproximationError.cc
  mirtkTransformationConfig.cc
  mirtkTransformationConstraint.cc
  mirtkTransformationInverse.cc
  mirtkTransformationUtils.h
  mirtkVolumePreservationConstraint.cc
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
