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

#ifndef MIRTK_Downsampling_H
#define MIRTK_Downsampling_H

#include "mirtk/GenericImage.h"
#include "mirtk/ImageToImage.h"
#include "mirtk/ScalarFunction.h"


namespace mirtk {


/**
 * Filter for downsampling of images by a power of two
 *
 * By default, a Gaussian kernel is used to average the voxel intensities within
 * the neighborhood corresponding to each downsampled voxel. This can be disabled
 * by specifying a zero sigma value for the Gaussian kernel.
 */
template <class TVoxel>
class Downsampling : public ImageToImage<TVoxel>
{
  mirtkImageFilterMacro(Downsampling, TVoxel);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Downsampling factor in x dimension
  mirtkPublicAttributeMacro(int, DownsampleFactorX);

  /// Downsampling factor in y dimension
  mirtkPublicAttributeMacro(int, DownsampleFactorY);

  /// Downsampling factor in z dimension
  mirtkPublicAttributeMacro(int, DownsampleFactorZ);

  /// Custom kernel to be applied at each downsampled voxel or NULL
  ScalarFunction *_Kernel;

  /// Extend of discretized kernel
  int _KernelSize;

  /// Whether to normalize the discretized kernel weights
  bool _NormalizeKernel;

protected:

  /// Initialize the filter
  virtual void Initialize();

  /// Discretize continious kernel
  void DiscreteKernel(GenericImage<RealPixel> &, double);

public:

  /// Constructor
  Downsampling(int = 2);

  /// Constructor
  Downsampling(int, int, int = 1);

  /// Set downsampling kernel and desired radius of discretized kernel
  void Kernel(ScalarFunction *, int, bool = true);

  /// Run the downsampling filter
  virtual void Run();

};


} // namespace mirtk

#endif // MIRTK_Downsampling_H
