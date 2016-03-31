/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2015 Imperial College London
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

#ifndef MIRTK_GaussianPyramidFilter_H
#define MIRTK_GaussianPyramidFilter_H

#include "mirtk/ImageToImage.h"
#include "mirtk/ImageFunction.h"


namespace mirtk {


/**
 * Filter for down-/upsampling of a Gaussian pyramid image from one level to another
 */
template <class TVoxel>
class GaussianPyramidFilter : public ImageToImage<TVoxel>
{
  mirtkImageFilterMacro(GaussianPyramidFilter, TVoxel);

  /// Level of the input image
  mirtkPublicAttributeMacro(int, InputLevel);

  /// Level of the output image
  mirtkPublicAttributeMacro(int, OutputLevel);

protected:

  /// Initialize the filter
  virtual void Initialize();

  /// Downsample image by one level
  virtual void Downsample();

  /// Upsample image by one level
  virtual void Upsample();

public:

  /// Constructor
  GaussianPyramidFilter(int = 1);

  /// Constructor
  GaussianPyramidFilter(int, int);

  /// Run the downsampling filter
  virtual void Run();

};


} // namespace mirtk

#endif // MIRTK_GaussianPyramidFilter_H
