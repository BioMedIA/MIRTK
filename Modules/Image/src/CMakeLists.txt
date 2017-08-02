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
# @brief Build configuration of MIRTK Image library.
##############################################################################

# ------------------------------------------------------------------------------
# Source files
set(HEADERS
  ${BINARY_INCLUDE_DIR}/mirtk/ImageConfig.h
  ${BINARY_INCLUDE_DIR}/mirtk/ImageExport.h
  BaseImage.h
  BinaryVoxelFunction.h
  BSplineInterpolateImageFunction.h
  BSplineInterpolateImageFunction.hxx
  BSplineInterpolateImageFunction2D.h
  BSplineInterpolateImageFunction2D.hxx
  BSplineInterpolateImageFunction3D.h
  BSplineInterpolateImageFunction3D.hxx
  BSplineInterpolateImageFunction4D.h
  BSplineInterpolateImageFunction4D.hxx
  CityBlockDistanceTransform.h
  Closing.h
  ConnectedComponents.h
  ConstExtrapolateImageFunction.h
  ConstExtrapolateImageFunctionWithPeriodicTime.h
  ConstGenericImageIterator.h
  ConstImageIterator.h
  ConvolutionFunction.h
  CubicBSplineConvolution.h
  CSplineInterpolateImageFunction.h
  CSplineInterpolateImageFunction.hxx
  CSplineInterpolateImageFunction2D.h
  CSplineInterpolateImageFunction2D.hxx
  CSplineInterpolateImageFunction3D.h
  CSplineInterpolateImageFunction3D.hxx
  CSplineInterpolateImageFunction4D.h
  CSplineInterpolateImageFunction4D.hxx
  CubicBSplineInterpolateImageFunction.h
  CubicBSplineInterpolateImageFunction.hxx
  CubicBSplineInterpolateImageFunction2D.h
  CubicBSplineInterpolateImageFunction2D.hxx
  CubicBSplineInterpolateImageFunction3D.h
  CubicBSplineInterpolateImageFunction3D.hxx
  CubicBSplineInterpolateImageFunction4D.h
  CubicBSplineInterpolateImageFunction4D.hxx
  DataFunctions.h
  DataOp.h
  DataStatistics.h
  DifferenceOfCompositionLieBracketImageFilter3D.h
  Dilation.h
  DisplacementFieldExp.h
  DisplacementToVelocityField.h
  DisplacementToVelocityFieldBCH.h
  Downsampling.h
  Erosion.h
  EuclideanDistanceTransform.h
  ExtrapolateImageFunction.h
  ExtrapolationMode.h
  FastCubicBSplineInterpolateImageFunction.h
  FastCubicBSplineInterpolateImageFunction.hxx
  FastCubicBSplineInterpolateImageFunction2D.h
  FastCubicBSplineInterpolateImageFunction2D.hxx
  FastCubicBSplineInterpolateImageFunction3D.h
  FastCubicBSplineInterpolateImageFunction3D.hxx
  FastCubicBSplineInterpolateImageFunction4D.h
  FastCubicBSplineInterpolateImageFunction4D.hxx
  FastLinearImageGradientFunction.h
  FastLinearImageGradientFunction.hxx
  FastLinearImageGradientFunction2D.h
  FastLinearImageGradientFunction2D.hxx
  FastLinearImageGradientFunction3D.h
  FastLinearImageGradientFunction3D.hxx
  GaussianBlurring.h
  GaussianBlurring2D.h
  GaussianBlurring4D.h
  GaussianBlurringWithPadding.h
  GaussianBlurringWithPadding2D.h
  GaussianInterpolateImageFunction.h
  GaussianInterpolateImageFunction.hxx
  GaussianInterpolateImageFunction2D.h
  GaussianInterpolateImageFunction2D.hxx
  GaussianInterpolateImageFunction3D.h
  GaussianInterpolateImageFunction3D.hxx
  GaussianInterpolateImageFunction4D.h
  GaussianInterpolateImageFunction4D.hxx
  GaussianPyramidFilter.h
  GenericImage.h
  GenericImageIterator.h
  GradientImageFilter.h
  HashImage.h
  HashImage.hxx
  HessianImageFilter.h
  Histogram1D.h
  Histogram2D.h
  HistogramMatching.h
  Image.h
  ImageAttributes.h
  ImageFunction.h
  ImageGradientFunction.h
  ImageGradientFunction.hxx
  ImageIterator.h
  ImageReader.h
  ImageReaderFactory.h
  ImageRegion.h
  ImageSequence.h
  ImageSequence.hh
  ImageToImage.h
  ImageToInterpolationCoefficients.h
  ImageWriter.h
  ImageWriterFactory.h
  InterpolateImageFunction.h
  InterpolateImageFunction.hxx
  InterpolationMode.h
  LieBracketImageFilter.h
  LieBracketImageFilter2D.h
  LieBracketImageFilter3D.h
  LinearImageGradientFunction.h
  LinearImageGradientFunction.hxx
  LinearImageGradientFunction2D.h
  LinearImageGradientFunction2D.hxx
  LinearImageGradientFunction3D.h
  LinearImageGradientFunction3D.hxx
  LinearInterpolateImageFunction.h
  LinearInterpolateImageFunction.hxx
  LinearInterpolateImageFunction2D.h
  LinearInterpolateImageFunction2D.hxx
  LinearInterpolateImageFunction3D.h
  LinearInterpolateImageFunction3D.hxx
  LinearInterpolateImageFunction4D.h
  LinearInterpolateImageFunction4D.hxx
  MirrorExtrapolateImageFunction.h
  NaryVoxelFunction.h
  NearestNeighborExtrapolateImageFunction.h
  NearestNeighborInterpolateImageFunction.h
  NearestNeighborInterpolateImageFunction.hxx
  NeighborhoodOffsets.h
  RepeatExtrapolateImageFunction.h
  Resampling.h
  ResamplingWithPadding.h
  ScalarFunctionToImage.h
  ScalingAndSquaring.h
  SeparableConvolution.h
  ShapeBasedInterpolateImageFunction.h
  SincInterpolateImageFunction.h
  SincInterpolateImageFunction.hxx
  SincInterpolateImageFunction2D.h
  SincInterpolateImageFunction2D.hxx
  SincInterpolateImageFunction3D.h
  SincInterpolateImageFunction3D.hxx
  SincInterpolateImageFunction4D.h
  SincInterpolateImageFunction4D.hxx
  TernaryVoxelFunction.h
  UnaryVoxelFunction.h
  VelocityFieldExp.h
  VelocityToDisplacementField.h
  VelocityToDisplacementFieldEuler.h
  VelocityToDisplacementFieldSS.h
  Voxel.h
  VoxelCast.h
  VoxelDomain.h
  VoxelFunction.h
)

set(SOURCES
  BaseImage.cc
  CityBlockDistanceTransform.cc
  DataOp.cc
  DifferenceOfCompositionLieBracketImageFilter3D.cc
  Dilation.cc
  DisplacementToVelocityFieldBCH.cc
  Closing.cc
  ConnectedComponents.cc
  CubicBSplineConvolution.cc
  Downsampling.cc
  Erosion.cc
  EuclideanDistanceTransform.cc
  ExtrapolateImageFunction.cc
  GaussianBlurring.cc
  GaussianBlurring2D.cc
  GaussianBlurring4D.cc
  GaussianBlurringWithPadding.cc
  GaussianBlurringWithPadding2D.cc
  GaussianPyramidFilter.cc
  GenericImage.cc
  GradientImageFilter.cc
  HessianImageFilter.cc
  Histogram1D.cc
  Histogram2D.cc
  HistogramMatching.cc
  ImageAttributes.cc
  ImageFunction.cc
  ImageGradientFunction.cc
  ImageReader.cc
  ImageReaderFactory.cc
  ImageToImage.cc
  ImageWriter.cc
  ImageWriterFactory.cc
  InterpolateImageFunction.cc
  NeighborhoodOffsets.cc
  Resampling.cc
  ResamplingWithPadding.cc
  ScalarFunctionToImage.cc
  ScalingAndSquaring.cc
  SeparableConvolution.cc
  ShapeBasedInterpolateImageFunction.cc
  VelocityToDisplacementField.cc
  VelocityToDisplacementFieldEuler.cc
  VelocityToDisplacementFieldSS.cc
  Voxel.cc
)

set(DEPENDS
  LibCommon
  LibNumerics
)

if (TARGET TBB::tbb)
  list(APPEND DEPENDS TBB::tbb)
endif ()

if (VTK_FOUND)
  list(APPEND DEPENDS ${VTK_LIBRARIES})
endif ()

# ------------------------------------------------------------------------------
# ForEachVoxel*Function headers
if (PYTHON_EXECUTABLE AND BUILD_FOREACHVOXELFUNCTION_SOURCE)
  set(FOREACHVOXELFUNCTIONS
    ForEachUnaryVoxelFunction
    ForEachBinaryVoxelFunction
    ForEachTernaryVoxelFunction
    ForEachQuaternaryVoxelFunction
    ForEachQuinaryVoxelFunction
    ForEachSenaryVoxelFunction
    ForEachSeptenaryVoxelFunction
    ForEachOctaryVoxelFunction
    ForEachNonaryVoxelFunction
  )
  
  foreach(f IN LISTS FOREACHVOXELFUNCTIONS)
    string(REGEX REPLACE "ForEach|VoxelFunction" "" arity "${f}")
    string(TOLOWER "${arity}" arity)
    add_custom_command(
      OUTPUT  "${PROJECT_INCLUDE_DIR}/mirtk/${f}.h"
      COMMAND "${PYTHON_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/ForEachVoxelFunction.py"
              ${arity} "${PROJECT_INCLUDE_DIR}/mirtk/${f}.h"
      DEPENDS ForEachVoxelFunction.py
      COMMENT "Generating source file mirtk/${f}.h"
    )
    add_custom_target(foreachvoxel_${arity} DEPENDS "${PROJECT_INCLUDE_DIR}/mirtk/${f}.h")
    list(APPEND FOREACHVOXELDEPENDS foreachvoxel_${arity})
    list(APPEND HEADERS ${f}.h)
  endforeach()
  
  add_custom_target(foreachvoxel
    DEPENDS ${FOREACHVOXELDEPENDS}
    SOURCES ForEachVoxelFunction.py
  )
  
  list(APPEND DEPENDS foreachvoxel)
endif ()

# ------------------------------------------------------------------------------
# Add library target
mirtk_add_library()
