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

#ifndef MIRTK_VelocityToDisplacementFieldEuler_H
#define MIRTK_VelocityToDisplacementFieldEuler_H

#include "mirtk/VelocityToDisplacementField.h"


namespace mirtk {


// Vector field interpolator
class InterpolateImageFunction;


/**
 * Computes a displacement field from a stationary velocity field.
 *
 * This class implements an image filter which computes the exponential map
 * of a stationary velocity field using the forward Euler integration method.
 * The result is a diffeomorphic displacement field.
 */
template <class TVoxel>
class VelocityToDisplacementFieldEuler : public VelocityToDisplacementField<TVoxel>
{
  mirtkImageFilterMacro(VelocityToDisplacementFieldEuler, TVoxel);

protected:

  /// Image function for interpolation of velocities
  InterpolateImageFunction *_VelocityInterpolator;

  /// Initialize filter
  virtual void Initialize();

public:

  /// Constructor
  VelocityToDisplacementFieldEuler();

  /// Destructor
  virtual ~VelocityToDisplacementFieldEuler();

  /// Compute output = log(input)
  virtual void Run();

};


} // namespace mirtk

#endif // MIRTK_VelocityToDisplacementFieldEuler_H
