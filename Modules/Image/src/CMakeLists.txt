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
  mirtkBaseImage.h
  mirtkBinaryVoxelFunction.h
  mirtkBSplineInterpolateImageFunction.h
  mirtkBSplineInterpolateImageFunction.hxx
  mirtkBSplineInterpolateImageFunction2D.h
  mirtkBSplineInterpolateImageFunction2D.hxx
  mirtkBSplineInterpolateImageFunction3D.h
  mirtkBSplineInterpolateImageFunction3D.hxx
  mirtkBSplineInterpolateImageFunction4D.h
  mirtkBSplineInterpolateImageFunction4D.hxx
  mirtkCityBlockDistanceTransform.h
  mirtkConnectedComponents.h
  mirtkConstExtrapolateImageFunction.h
  mirtkConstExtrapolateImageFunctionWithPeriodicTime.h
  mirtkConstGenericImageIterator.h
  mirtkConstImageIterator.h
  mirtkConvolutionFunction.h
  mirtkCSplineInterpolateImageFunction.h
  mirtkCSplineInterpolateImageFunction.hxx
  mirtkCSplineInterpolateImageFunction2D.h
  mirtkCSplineInterpolateImageFunction2D.hxx
  mirtkCSplineInterpolateImageFunction3D.h
  mirtkCSplineInterpolateImageFunction3D.hxx
  mirtkCSplineInterpolateImageFunction4D.h
  mirtkCSplineInterpolateImageFunction4D.hxx
  mirtkCubicBSplineInterpolateImageFunction.h
  mirtkCubicBSplineInterpolateImageFunction.hxx
  mirtkCubicBSplineInterpolateImageFunction2D.h
  mirtkCubicBSplineInterpolateImageFunction2D.hxx
  mirtkCubicBSplineInterpolateImageFunction3D.h
  mirtkCubicBSplineInterpolateImageFunction3D.hxx
  mirtkCubicBSplineInterpolateImageFunction4D.h
  mirtkCubicBSplineInterpolateImageFunction4D.hxx
  mirtkDataFunctions.h
  mirtkDataOp.h
  mirtkDataStatistics.h
  mirtkDifferenceOfCompositionLieBracketImageFilter3D.h
  mirtkDilation.h
  mirtkDisplacementFieldExp.h
  mirtkDisplacementToVelocityField.h
  mirtkDisplacementToVelocityFieldBCH.h
  mirtkDownsampling.h
  mirtkErosion.h
  mirtkEuclideanDistanceTransform.h
  mirtkExtrapolateImageFunction.h
  mirtkExtrapolationMode.h
  mirtkFastCubicBSplineInterpolateImageFunction.h
  mirtkFastCubicBSplineInterpolateImageFunction.hxx
  mirtkFastCubicBSplineInterpolateImageFunction2D.h
  mirtkFastCubicBSplineInterpolateImageFunction2D.hxx
  mirtkFastCubicBSplineInterpolateImageFunction3D.h
  mirtkFastCubicBSplineInterpolateImageFunction3D.hxx
  mirtkFastCubicBSplineInterpolateImageFunction4D.h
  mirtkFastCubicBSplineInterpolateImageFunction4D.hxx
  mirtkFastLinearImageGradientFunction.h
  mirtkFastLinearImageGradientFunction.hxx
  mirtkFastLinearImageGradientFunction2D.h
  mirtkFastLinearImageGradientFunction2D.hxx
  mirtkFastLinearImageGradientFunction3D.h
  mirtkFastLinearImageGradientFunction3D.hxx
  mirtkGaussianBlurring.h
  mirtkGaussianBlurring2D.h
  mirtkGaussianBlurring4D.h
  mirtkGaussianBlurringWithPadding.h
  mirtkGaussianBlurringWithPadding2D.h
  mirtkGaussianInterpolateImageFunction.h
  mirtkGaussianInterpolateImageFunction.hxx
  mirtkGaussianInterpolateImageFunction2D.h
  mirtkGaussianInterpolateImageFunction2D.hxx
  mirtkGaussianInterpolateImageFunction3D.h
  mirtkGaussianInterpolateImageFunction3D.hxx
  mirtkGaussianInterpolateImageFunction4D.h
  mirtkGaussianInterpolateImageFunction4D.hxx
  mirtkGaussianPyramidFilter.h
  mirtkGenericImage.h
  mirtkGenericImageIterator.h
  mirtkGradientImageFilter.h
  mirtkHashImage.h
  mirtkHessianImageFilter.h
  mirtkHistogram1D.h
  mirtkHistogram2D.h
  mirtkImage.h
  mirtkImageAttributes.h
  mirtkImageFunction.h
  mirtkImageGradientFunction.h
  mirtkImageGradientFunction.hxx
  mirtkImageIterator.h
  mirtkImageReader.h
  mirtkImageReaderFactory.h
  mirtkImageRegion.h
  mirtkImageSequence.h
  mirtkImageSequence.hh
  mirtkImageToImage.h
  mirtkImageToInterpolationCoefficients.h
  mirtkImageWriter.h
  mirtkImageWriterFactory.h
  mirtkInterpolateImageFunction.h
  mirtkInterpolateImageFunction.hxx
  mirtkInterpolationMode.h
  mirtkLieBracketImageFilter.h
  mirtkLieBracketImageFilter2D.h
  mirtkLieBracketImageFilter3D.h
  mirtkLinearImageGradientFunction.h
  mirtkLinearImageGradientFunction.hxx
  mirtkLinearImageGradientFunction2D.h
  mirtkLinearImageGradientFunction2D.hxx
  mirtkLinearImageGradientFunction3D.h
  mirtkLinearImageGradientFunction3D.hxx
  mirtkLinearInterpolateImageFunction.h
  mirtkLinearInterpolateImageFunction.hxx
  mirtkLinearInterpolateImageFunction2D.h
  mirtkLinearInterpolateImageFunction2D.hxx
  mirtkLinearInterpolateImageFunction3D.h
  mirtkLinearInterpolateImageFunction3D.hxx
  mirtkLinearInterpolateImageFunction4D.h
  mirtkLinearInterpolateImageFunction4D.hxx
  mirtkMirrorExtrapolateImageFunction.h
  mirtkNaryVoxelFunction.h
  mirtkNearestNeighborExtrapolateImageFunction.h
  mirtkNearestNeighborInterpolateImageFunction.h
  mirtkNearestNeighborInterpolateImageFunction.hxx
  mirtkNeighborhoodOffsets.h
  mirtkRepeatExtrapolateImageFunction.h
  mirtkResampling.h
  mirtkResamplingWithPadding.h
  mirtkScalarFunctionToImage.h
  mirtkScalingAndSquaring.h
  mirtkShapeBasedInterpolateImageFunction.h
  mirtkSincInterpolateImageFunction.h
  mirtkSincInterpolateImageFunction.hxx
  mirtkSincInterpolateImageFunction2D.h
  mirtkSincInterpolateImageFunction2D.hxx
  mirtkSincInterpolateImageFunction3D.h
  mirtkSincInterpolateImageFunction3D.hxx
  mirtkSincInterpolateImageFunction4D.h
  mirtkSincInterpolateImageFunction4D.hxx
  mirtkTernaryVoxelFunction.h
  mirtkUnaryVoxelFunction.h
  mirtkVelocityFieldExp.h
  mirtkVelocityToDisplacementField.h
  mirtkVelocityToDisplacementFieldEuler.h
  mirtkVelocityToDisplacementFieldSS.h
  mirtkVoxel.h
  mirtkVoxelCast.h
  mirtkVoxelDomain.h
  mirtkVoxelFunction.h
)

set(SOURCES
  mirtkBaseImage.cc
  mirtkCityBlockDistanceTransform.cc
  mirtkDataOp.cc
  mirtkDifferenceOfCompositionLieBracketImageFilter3D.cc
  mirtkDilation.cc
  mirtkDisplacementToVelocityFieldBCH.cc
  mirtkConnectedComponents.cc
  mirtkDownsampling.cc
  mirtkErosion.cc
  mirtkEuclideanDistanceTransform.cc
  mirtkExtrapolateImageFunction.cc
  mirtkGaussianBlurring.cc
  mirtkGaussianBlurring2D.cc
  mirtkGaussianBlurring4D.cc
  mirtkGaussianBlurringWithPadding.cc
  mirtkGaussianBlurringWithPadding2D.cc
  mirtkGaussianPyramidFilter.cc
  mirtkGenericImage.cc
  mirtkGradientImageFilter.cc
  mirtkHashImage.cc
  mirtkHessianImageFilter.cc
  mirtkHistogram1D.cc
  mirtkHistogram2D.cc
  mirtkImageAttributes.cc
  mirtkImageFunction.cc
  mirtkImageGradientFunction.cc
  mirtkImageReader.cc
  mirtkImageReaderFactory.cc
  mirtkImageToImage.cc
  mirtkImageWriter.cc
  mirtkImageWriterFactory.cc
  mirtkInterpolateImageFunction.cc
  mirtkNeighborhoodOffsets.cc
  mirtkResampling.cc
  mirtkResamplingWithPadding.cc
  mirtkScalarFunctionToImage.cc
  mirtkScalingAndSquaring.cc
  mirtkShapeBasedInterpolateImageFunction.cc
  mirtkVelocityToDisplacementField.cc
  mirtkVelocityToDisplacementFieldEuler.cc
  mirtkVelocityToDisplacementFieldSS.cc
  mirtkVoxel.cc
)

set(DEPENDS LibCommon LibNumerics)
if (TBB_FOUND)
  list(APPEND DEPENDS ${TBB_LIBRARIES})
endif ()

# ------------------------------------------------------------------------------
# ForEachVoxel*Function headers
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
  find_package(PythonInterp REQUIRED)
  string(REGEX REPLACE "ForEach|VoxelFunction" "" arity "${f}")
  string(TOLOWER "${arity}" arity)
  add_custom_command(
    OUTPUT  "${PROJECT_INCLUDE_DIR}/mirtk${f}.h"
    COMMAND "${PYTHON_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/mirtkForEachVoxelFunction.py"
            ${arity} "${PROJECT_INCLUDE_DIR}/mirtk${f}.h"
    DEPENDS mirtkForEachVoxelFunction.py
    COMMENT "Generating source file mirtk${f}.h"
  )
  add_custom_target(foreachvoxel_${arity} DEPENDS "${PROJECT_INCLUDE_DIR}/mirtk${f}.h")
  list(APPEND FOREACHVOXELDEPENDS foreachvoxel_${arity})
  list(APPEND HEADERS mirtk${f}.h)
endforeach()

add_custom_target(foreachvoxel
  DEPENDS ${FOREACHVOXELDEPENDS}
  SOURCES mirtkForEachVoxelFunction.py
)

list(APPEND DEPENDS foreachvoxel)

# ------------------------------------------------------------------------------
# Add library target
mirtk_add_library()
