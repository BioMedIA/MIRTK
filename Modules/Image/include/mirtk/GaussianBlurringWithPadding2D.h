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

#ifndef MIRTK_GaussianBlurringWithPadding2D_H
#define MIRTK_GaussianBlurringWithPadding2D_H

#include "mirtk/GaussianBlurringWithPadding.h"


namespace mirtk {


/**
 * Class for Gaussian blurring of padded images
 *
 * This class defines and implements the Gaussian blurring of padded images.
 * It takes 2D and 3D images but blurres only in the x and y direction.
 * The blurring is implemented by two successive 1D convolutions with a 1D
 * Gaussian kernel. If more than 50% of the voxels used for the convolution
 * have intensities smaller or equal to the padding value, the blurred voxel
 * will be filled with the padding value.
 */
template <class TVoxel>
class GaussianBlurringWithPadding2D : public GaussianBlurringWithPadding<TVoxel>
{
  mirtkInPlaceImageFilterMacro(GaussianBlurringWithPadding2D, TVoxel);

public:

  /// Constructor
  GaussianBlurringWithPadding2D(double = 1.0, VoxelType = -1);

  /// Constructor
  GaussianBlurringWithPadding2D(double, double, VoxelType);

};


} // namespace mirtk

#endif // MIRTK_GaussianBlurringWithPadding2D_H
