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
# @brief Build configuration of MIRTK PointSet library.
##############################################################################

set(HEADERS
  ${BINARY_INCLUDE_DIR}/mirtk/PointSetExport.h
  ClosestCell.h
  ClosestPoint.h
  ClosestPointLabel.h
  EdgeConnectivity.h
  EdgeTable.h
  FiducialMatch.h
  FuzzyCorrespondence.h
  ImplicitSurfaceUtils.h
  PointCorrespondence.h
  PointLocator.h
  PointSetUtils.h
  PolyDataCurvature.h
  PolyDataFilter.h
  PolyDataRemeshing.h
  PolyDataSmoothing.h
  Polyhedron.h
  RegisteredPointSet.h
  RegisteredSurface.h
  RobustClosestPoint.h
  RobustPointMatch.h
  SpectralDecomposition.h
  SpectralMatch.h
  SurfaceCollisions.h
  Triangle.h
)

set(SOURCES
  dfsurface.h
  ClosestCell.cc
  ClosestPoint.cc
  ClosestPointLabel.cc
  EdgeConnectivity.cc
  EdgeTable.cc
  FiducialMatch.cc
  FuzzyCorrespondence.cc
  FuzzyCorrespondenceUtils.cc
  FuzzyCorrespondenceUtils.h
  ImplicitSurfaceUtils.cc
  PointCorrespondence.cc
  PointLocator.cc
  PointSetUtils.cc
  PolyDataCurvature.cc
  PolyDataFilter.cc
  PolyDataRemeshing.cc
  PolyDataSmoothing.cc
  Polyhedron.cc
  RegisteredPointSet.cc
  RegisteredSurface.cc
  RobustClosestPoint.cc
  RobustPointMatch.cc
  SpectralDecomposition.cc
  SpectralMatch.cc
  SurfaceCollisions.cc
)

set(DEPENDS
  LibCommon
  LibNumerics
  LibImage
  LibTransformation
  ${VTK_LIBRARIES}
)

if (TARGET TBB::tbb)
  list(APPEND DEPENDS TBB::tbb)
endif ()

if (FLANN_FOUND)
  list(APPEND DEPENDS ${FLANN_CPP_LIBRARY})
endif ()

mirtk_add_library()
