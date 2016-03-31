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

#ifndef MIRTK_VelocityToDisplacementFieldSS_H
#define MIRTK_VelocityToDisplacementFieldSS_H

#include "mirtk/VelocityToDisplacementField.h"
#include "mirtk/GenericImage.h"


namespace mirtk {


// Vector field interpolator
class InterpolateImageFunction;


/**
 * Computes a displacement field from a stationary velocity field.
 *
 * This class implements an image filter which computes the exponential map
 * of a stationary velocity field using the scaling and squaring method.
 * The result is a diffeomorphic displacement field.
 */
template <class TVoxel>
class VelocityToDisplacementFieldSS : public VelocityToDisplacementField<TVoxel>
{
  mirtkInPlaceImageFilterMacro(VelocityToDisplacementFieldSS, TVoxel);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Number of squaring steps
  mirtkPublicAttributeMacro(int, NumberOfSquaringSteps);

  /// Whether to upsample the input vector field before
  mirtkPublicAttributeMacro(bool, Upsample);

  /// Whether to use a Gaussian pyramid downsampler when downsampling the
  /// previously upsample vector field to obey Shannon's sampling theorem.
  /// Otherwise a simple interpolation without smoothing kernel is used.
  mirtkPublicAttributeMacro(bool, SmoothBeforeDownsampling);

  /// Maximum velocity after scaling
  ///
  /// Set to zero in order to scale velocities by exactly pow(2.0, _NumberOfSquaringSteps).
  /// Otherwise, the number of squaring steps is increased in order to ensure that
  /// the maximum velocity in each dimension is less or equal the specified value.
  mirtkPublicAttributeMacro(VoxelType, MaxScaledVelocity);

  /// External memory that can be used for intermedate displacement field
  mirtkPublicAggregateMacro(ImageType, ExternalCache);

  /// Intermediate displacement field
  mirtkAggregateMacro(ImageType, Displacement);

  /// Interpolator of intermediate displacement field
  mirtkAggregateMacro(InterpolateImageFunction, Interpolator);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Constructor
  VelocityToDisplacementFieldSS();

  /// Destructor
  virtual ~VelocityToDisplacementFieldSS();

  // ---------------------------------------------------------------------------
  // Execution

  /// Compute output = log(input)
  virtual void Run();

protected:

  /// Initialize filter
  virtual void Initialize();

  /// Finalize filter
  virtual void Finalize();

};


} // namespace mirtk

#endif // MIRTK_VelocityToDisplacementFieldSS_H
