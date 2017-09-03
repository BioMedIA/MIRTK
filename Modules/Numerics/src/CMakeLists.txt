# ============================================================================
# Medical Image Registration ToolKit (MIRTK)
#
# Copyright 2013-2016 Imperial College London
# Copyright 2013-2016 Andreas Schuh
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
  ${BINARY_INCLUDE_DIR}/mirtk/NumericsConfig.h
  ${BINARY_INCLUDE_DIR}/mirtk/NumericsExport.h
  AdaptiveLineSearch.h
  Arith.h
  BrentLineSearch.h
  BSpline.h
  CharbonnierErrorFunction.h
  ConjugateGradientDescent.h
  DistanceErrorFunction.h
  EnergyThreshold.h
  GaussianErrorFunction.h
  GradientDescent.h
  InexactLineSearch.h
  LineSearch.h
  LocalOptimizer.h
  Matrix.h
  Matrix3x3.h
  MaxStepLineSearch.h
  Numerics.h
  ObjectiveFunction.h
  OptimizationMethod.h
  PeronaMalikErrorFunction.h
  Plane.h
  Point.h
  PointSamples.h
  PointSet.h
  Polynomial.h
  PolynomialSolvers.h
  RadialErrorFunction.h
  ScalarFunction.h
  ScalarGaussian.h
  Sinc.h
  SparseMatrix.h
  SquaredErrorFunction.h
  StoppingCriterion.h
  Vector.h
  Vector3.h
  Vector3D.h
  Vector4D.h
  VectorND.h
  # private headers -- TODO: Exclude from installation
  Arpack.h
  Umfpack.h
  Eigen.h
)

set(SOURCES
  AdaptiveLineSearch.cc
  Arith.cc
  BrentLineSearch.cc
  BSpline.cc
  ConjugateGradientDescent.cc
  EnergyThreshold.cc
  GradientDescent.cc
  InexactLineSearch.cc
  LineSearch.cc
  LocalOptimizer.cc
  Matrix.cc
  Matrix3x3.cc
  MaxStepLineSearch.cc
  NumericsConfig.cc
  Plane.cc
  Point.cc
  PointSamples.cc
  PointSet.cc
  Polynomial.cc
  PolynomialSolvers.cc
  RadialErrorFunction.cc
  ScalarGaussian.cc
  Sinc.cc
  SparseMatrix.cc
  StoppingCriterion.cc
  Vector.cc
  Vector3.cc
  Vector3D.cc
  Vector4D.cc
)

set(DEPENDS LibCommon)

if (LBFGS_FOUND)
  list(APPEND HEADERS LimitedMemoryBFGSDescent.h)
  list(APPEND SOURCES LimitedMemoryBFGSDescent.cc)
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
