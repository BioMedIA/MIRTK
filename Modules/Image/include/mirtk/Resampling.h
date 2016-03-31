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

#ifndef MIRTK_Resampling_H
#define MIRTK_Resampling_H

#include "mirtk/ImageToImage.h"


namespace mirtk {


// Interpolator base type
class InterpolateImageFunction;


/**
 * Class for resampling of images
 *
 * This class defines and implements the resampling of images with arbitrary
 * voxel dimensions.  The new image intensity of the voxels is calculated by
 * interpolation of the old image intensities. Possible interpolation schemes
 * are nearest neighbor, linear, cubic spline and B-spline interpolation.
 */
template <class TVoxel>
class Resampling : public ImageToImage<TVoxel>
{
  mirtkImageFilterMacro(Resampling, TVoxel);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Size of output after resampling in x dimension
  mirtkPublicAttributeMacro(int, X);

  /// Size of output after resampling in y dimension
  mirtkPublicAttributeMacro(int, Y);

  /// Size of output after resampling in z dimension
  mirtkPublicAttributeMacro(int, Z);

  /// Voxel size of output after resampling in x dimension
  mirtkPublicAttributeMacro(double, XSize);

  /// Voxel size of output after resampling in y dimension
  mirtkPublicAttributeMacro(double, YSize);

  /// Voxel size of output after resampling in z dimension
  mirtkPublicAttributeMacro(double, ZSize);

  /// Image function used to interpolate/extrapolate input
  mirtkPublicAggregateMacro(InterpolateImageFunction, Interpolator);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Constructor
  Resampling(double, double, double);

  /// Constructor
  Resampling(int, int, int);

  /// Constructor
  Resampling(int, int, int, double, double, double);

  // ---------------------------------------------------------------------------
  // Execution

  /// Run the resampling filter
  virtual void Run();

protected:

  /// Initialize the filter
  virtual void Initialize();

  /// Initialize filter output
  virtual void InitializeOutput();

};


} // namespace mirtk

#endif // MIRTK_Resampling_H
