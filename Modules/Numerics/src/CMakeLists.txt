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
# @brief Build configuration of MIRTK Numerics library.
##############################################################################

set(HEADERS
  ${BINARY_INCLUDE_DIR}/mirtkNumericsConfig.h
  mirtkAdaptiveLineSearch.h
  mirtkArith.h
  mirtkBrentLineSearch.h
  mirtkBSpline.h
  mirtkCharbonnierErrorFunction.h
  mirtkConjugateGradientDescent.h
  mirtkDistanceErrorFunction.h
  mirtkEnergyThreshold.h
  mirtkGaussianErrorFunction.h
  mirtkGradientDescent.h
  mirtkInexactLineSearch.h
  mirtkLineSearch.h
  mirtkLocalOptimizer.h
  mirtkMatrix.h
  mirtkMatrix3x3.h
  mirtkMaxStepLineSearch.h
  mirtkNumerics.h
  mirtkObjectiveFunction.h
  mirtkOptimizationMethod.h
  mirtkPeronaMalikErrorFunction.h
  mirtkPoint.h
  mirtkPointSamples.h
  mirtkPointSet.h
  mirtkRadialErrorFunction.h
  mirtkScalarFunction.h
  mirtkScalarGaussian.h
  mirtkSinc.h
  mirtkSparseMatrix.h
  mirtkSquaredErrorFunction.h
  mirtkStoppingCriterion.h
  mirtkVector.h
  mirtkVector3.h
  mirtkVector3D.h
  mirtkVector4D.h
  # private headers -- TODO: Exclude from installation
  mirtkArpack.h
  mirtkUmfpack.h
  mirtkEigen.h
)

set(SOURCES
  mirtkAdaptiveLineSearch.cc
  mirtkArith.cc
  mirtkBrentLineSearch.cc
  mirtkBSpline.cc
  mirtkConjugateGradientDescent.cc
  mirtkEnergyThreshold.cc
  mirtkGradientDescent.cc
  mirtkInexactLineSearch.cc
  mirtkLineSearch.cc
  mirtkLocalOptimizer.cc
  mirtkMatrix.cc
  mirtkMatrix3x3.cc
  mirtkMaxStepLineSearch.cc
  mirtkNumericsConfig.cc
  mirtkPoint.cc
  mirtkPointSamples.cc
  mirtkPointSet.cc
  mirtkRadialErrorFunction.cc
  mirtkScalarGaussian.cc
  mirtkSinc.cc
  mirtkSparseMatrix.cc
  mirtkStoppingCriterion.cc
  mirtkVector.cc
  mirtkVector3.cc
  mirtkVector3D.cc
  mirtkVector4D.cc
)

set(DEPENDS LibCommon)

if (LBFGS_FOUND)
  list(APPEND HEADERS mirtkLimitedMemoryBFGSDescent.h)
  list(APPEND SOURCES mirtkLimitedMemoryBFGSDescent.cc)
  list(APPEND DEPENDS ${LBFGS_LIBRARIES})
endif ()

if (MATLAB_FOUND)
  list(APPEND DEPENDS ${MATLAB_mwmclmcrrt_LIBRARY})
endif ()
if (VTK_FOUND)
  list(APPEND DEPENDS ${VTK_LIBRARIES})
endif ()
if (ARPACK_FOUND AND UMFPACK_FOUND)
  list(APPEND DEPENDS ${ARPACK_LIBRARY} ${UMFPACK_LIBRARIES})
endif ()

mirtk_add_library(AUTO_REGISTER)
