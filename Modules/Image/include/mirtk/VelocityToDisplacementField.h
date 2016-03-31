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

#ifndef MIRTK_VelocityToDisplacementField_H
#define MIRTK_VelocityToDisplacementField_H

#include "mirtk/ImageToImage.h"

#include "mirtk/InterpolationMode.h"
#include "mirtk/ExtrapolationMode.h"


namespace mirtk {


/**
 * Computes a displacement field from a stationary velocity field.
 *
 * Base class for image filters which compute the exponential map of a
 * stationary velocity field. The filter output is a diffeomorphic
 * displacement field.
 */
template <class TVoxel>
class VelocityToDisplacementField : public ImageToImage<TVoxel>
{
  mirtkAbstractImageFilterMacro(VelocityToDisplacementField, TVoxel);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Vector field interpolation mode
  mirtkPublicAttributeMacro(InterpolationMode, Interpolation);
  
  /// Vector field extrapolation mode
  mirtkPublicAttributeMacro(ExtrapolationMode, Extrapolation);
  
  /// Whether to compute interpolation coefficients from the given input
  /// or if the input images contain the coefficients already
  mirtkPublicAttributeMacro(bool, ComputeInterpolationCoefficients);

  /// Whether filter should compute the inverse displacement field
  ///
  /// Note that this is equivalent to simply changing the sign of \c _T.
  /// It is in particular used by the DisplacementToVelocityField filter.
  mirtkPublicAttributeMacro(bool, ComputeInverse);

  /// Number of integration steps
  mirtkPublicAttributeMacro(int, NumberOfSteps);

  /// Upper integration limit (negative <=> inverse)
  mirtkPublicAttributeMacro(double, UpperIntegrationLimit);

protected:

  /// Optional input displacement field which is composed with the
  /// exponential of the velocity field, i.e., exp(v) ° phi
  const ImageType *_InputDisplacementField;

  /// Initialize filter
  virtual void Initialize();

  /// Finalize filter
  virtual void Finalize();

  /// Constructor
  VelocityToDisplacementField();

public:  

  /// Destructor
  virtual ~VelocityToDisplacementField();
 
  // Import other overloads
  using Baseclass::Input;

  /// Set n-th input (0: velocity field, 1: input displacement field, optional)
  virtual void Input(int, const ImageType *);

  /// Get n-th input (0: velocity field, 1: input displacement field, optional)
  virtual const ImageType *Input(int);

};


} // namespace mirtk

#endif // MIRTK_VelocityToDisplacementField_H
