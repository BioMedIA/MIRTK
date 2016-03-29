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
  mirtkClosestCell.h
  mirtkClosestPoint.h
  mirtkClosestPointLabel.h
  mirtkEdgeConnectivity.h
  mirtkEdgeTable.h
  mirtkFiducialMatch.h
  mirtkFuzzyCorrespondence.h
  mirtkImplicitSurfaceUtils.h
  mirtkPointCorrespondence.h
  mirtkPointLocator.h
  mirtkPointSetUtils.h
  mirtkPolyDataCurvature.h
  mirtkPolyDataFilter.h
  mirtkPolyDataRemeshing.h
  mirtkPolyDataSmoothing.h
  mirtkPolyhedron.h
  mirtkRegisteredPointSet.h
  mirtkRegisteredSurface.h
  mirtkRobustClosestPoint.h
  mirtkRobustPointMatch.h
  mirtkSpectralDecomposition.h
  mirtkSpectralMatch.h
  mirtkSurfaceCollisions.h
  mirtkTriangle.h
)

set(SOURCES
  dfsurface.h
  mirtkClosestCell.cc
  mirtkClosestPoint.cc
  mirtkClosestPointLabel.cc
  mirtkEdgeConnectivity.cc
  mirtkEdgeTable.cc
  mirtkFiducialMatch.cc
  mirtkFuzzyCorrespondence.cc
  mirtkFuzzyCorrespondenceUtils.cc
  mirtkFuzzyCorrespondenceUtils.h
  mirtkImplicitSurfaceUtils.cc
  mirtkPointCorrespondence.cc
  mirtkPointLocator.cc
  mirtkPointSetUtils.cc
  mirtkPolyDataCurvature.cc
  mirtkPolyDataFilter.cc
  mirtkPolyDataRemeshing.cc
  mirtkPolyDataSmoothing.cc
  mirtkPolyhedron.cc
  mirtkRegisteredPointSet.cc
  mirtkRegisteredSurface.cc
  mirtkRobustClosestPoint.cc
  mirtkRobustPointMatch.cc
  mirtkSpectralDecomposition.cc
  mirtkSpectralMatch.cc
  mirtkSurfaceCollisions.cc
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
  list(APPEND DEPENDS ${FLANN_LIBRARY})
endif ()

mirtk_add_library()
