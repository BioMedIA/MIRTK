/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2015 Imperial College London
 * Copyright 2008-2013 Daniel Rueckert, Julia Schnabel
 * Copyright 2013-2015 Andreas Schuh
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

#ifndef MIRTK_GradientImageFilter_H
#define MIRTK_GradientImageFilter_H

#include "mirtk/ImageToImage.h"


namespace mirtk {


/**
 * Class for calculating the gradient of an image.
 *
 * The class provides an interface to calculating the gradient in the
 * x, y and z directions. It also can process vector-valued input images,
 * where each vector component is treated separately. If the output for each
 * component is a 3D gradient vector, the output image contains first the
 * gradient vector for the first input component, then the gradient vector
 * of the second component, and so on. In case of a scalar gradient output
 * such as the gradient vector magnitude, the output image has the same
 * number of components (t dimension) as the input image.
 */
template <class TVoxel>
class GradientImageFilter : public ImageToImage<TVoxel>
{
  mirtkImageFilterMacro(GradientImageFilter, TVoxel);

  // ---------------------------------------------------------------------------
  // Types

public:

  /// Type of image gradient to compute
  enum GradientType {
    GRADIENT_X,
    GRADIENT_Y,
    GRADIENT_Z,
    GRADIENT_MAGNITUDE,
    GRADIENT_VECTOR,
    NORMALISED_GRADIENT_VECTOR,
    GRADIENT_DOT_PRODUCT,
    NUM_GRADIENT_TYPES // Keep as last enumeration entry!
  };

  // ---------------------------------------------------------------------------
  // Attributes

protected:

  /// Type of gradient
  mirtkPublicAttributeMacro(GradientType, Type);

  /// Whether to return gradient in mm or voxel units
  mirtkPublicAttributeMacro(bool, UseVoxelSize);

  /// Whether to return gradient in world orientation
  mirtkPublicAttributeMacro(bool, UseOrientation);

  /// Padding value
  mirtkPublicAttributeMacro(double, PaddingValue);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Constructor
  GradientImageFilter(GradientType type = GRADIENT_MAGNITUDE);

  // ---------------------------------------------------------------------------
  // Execution

  /// Run the convolution filter
  virtual void Run();

protected:

  /// Initialize filter after inputs and parameters are set
  virtual void Initialize();

  /// Finalize filter execution
  virtual void Finalize();

};


} // namespace mirtk

#endif // MIRTK_GradientImageFilter_H
