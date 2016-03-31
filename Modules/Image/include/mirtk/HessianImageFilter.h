/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2015 Imperial College London
 * Copyright 2008-2015 Daniel Rueckert, Julia Schnabel
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef MIRTK_HessianImageFilter_H
#define MIRTK_HessianImageFilter_H

#include "mirtk/ImageToImage.h"


namespace mirtk {


/**
 * Class for calculating the second order gradient of an image.
 *
 * The class provides an interface to calculating the second order gradient in the
 * x-x-, x-y-, x-z-, y-y-, y-z- and z-z- directions.
 */
template <class TVoxel>
class HessianImageFilter : public ImageToImage<TVoxel>
{
  mirtkImageFilterMacro(HessianImageFilter, TVoxel);

  // ---------------------------------------------------------------------------
  // Types
public:

  /// Type of gradient vector to compute
  enum OutputType {
    HESSIAN_XX,
    HESSIAN_XY,
    HESSIAN_XZ,
    HESSIAN_YY,
    HESSIAN_YZ,
    HESSIAN_ZZ,
    HESSIAN_VECTOR,
    HESSIAN_MATRIX
  };

  // ---------------------------------------------------------------------------
  // Attributes
protected:

  /// Whether to return gradient in mm or voxel units
  mirtkPublicAttributeMacro(bool, UseVoxelSize);

  /// Whether to return gradient in world orientation
  mirtkPublicAttributeMacro(bool, UseOrientation);

  /// Padding value
  mirtkPublicAttributeMacro(double, PaddingValue);

  /// Which Hessian component(s) to output
  mirtkPublicAttributeMacro(OutputType, Type);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
public:

  /// Constructor
  HessianImageFilter(OutputType type = HESSIAN_MATRIX);

  // ---------------------------------------------------------------------------
  // Execution

  /// Run the Hessian filter
  virtual void Run();

protected:

  /// Initialize the filter
  virtual void Initialize();

};


} // namespace mirtk

#endif // MIRTK_HessianImageFilter_H
