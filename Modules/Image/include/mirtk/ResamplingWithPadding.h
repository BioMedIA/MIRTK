/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2017 Imperial College London
 * Copyright 2008-2015 Daniel Rueckert, Julia Schnabel
 * Copyright 2017 Andreas Schuh
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

#ifndef MIRTK_ResamplingWithPadding_H
#define MIRTK_ResamplingWithPadding_H

#include "mirtk/Resampling.h"


namespace mirtk {


/**
 * Class for resampling of padded images
 *
 * This class defines and implements the resampling of images with arbitrary
 * voxel dimensions. The new image intensity of the voxels is calculated by
 * interpolation of the old image intensities, using
 * InterpolateImageFunction::EvaluateWithPadding.
 */
template <class TVoxel>
class ResamplingWithPadding : public Resampling<TVoxel>
{
  mirtkImageFilterMacro(ResamplingWithPadding, TVoxel);

  // ---------------------------------------------------------------------------
  // Attributes

protected:

  /// Padding value
  mirtkPublicAttributeMacro(VoxelType, PaddingValue);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Constructor
  ResamplingWithPadding(double, double, double, VoxelType);

  /// Constructor
  ResamplingWithPadding(int, int, int, VoxelType);

  /// Constructor
  ResamplingWithPadding(int, int, int, double, double, double, VoxelType);

  // ---------------------------------------------------------------------------
  // Execution

  /// Run the resampling filter
  virtual void Run();

protected:

  /// Initialize the filter
  virtual void Initialize();

};


} // namespace mirtk

#endif // MIRTK_ResamplingWithPadding_H
