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

#ifndef MIRTK_DifferenceOfCompositionLieBracketImageFilter3D_H
#define MIRTK_DifferenceOfCompositionLieBracketImageFilter3D_H

#include "mirtk/LieBracketImageFilter.h"
#include "mirtk/InterpolateImageFunction.h"


namespace mirtk {


/**
 * Image filter for computation of Lie bracket of two 3D vector fields.
 *
 * This filter implements the definition of the Lie bracket as the difference
 * of the composition of the first vector field with the second and vice versa,
 * i.e., [X,Y] = X(Y) - Y(X).
 */
template <class TVoxel>
class DifferenceOfCompositionLieBracketImageFilter3D
: public LieBracketImageFilter<TVoxel>
{
  mirtkImageFilterMacro(DifferenceOfCompositionLieBracketImageFilter3D, TVoxel);

  /// Vector field interpolation mode
  mirtkPublicAttributeMacro(InterpolationMode, Interpolation);

  /// Vector field extrapolation mode
  mirtkPublicAttributeMacro(ExtrapolationMode, Extrapolation);

  /// Whether to compute interpolation coefficients from the given input
  /// or if the input images contain the coefficients already
  mirtkPublicAttributeMacro(bool, ComputeInterpolationCoefficients);

protected:

  typedef GenericInterpolateImageFunction<ImageType> InterpolatorType;

  InterpolatorType *_Interpolator[2]; /// Input vector field interpolators
  double            _Scaling     [2]; /// Scaling of input vector fields

  /// Initialize filter
  virtual void Initialize();

public:

  /// Constructor
  DifferenceOfCompositionLieBracketImageFilter3D();

  /// Destructor
  virtual ~DifferenceOfCompositionLieBracketImageFilter3D();

  /// Set scaling of n-th input vector field
  void Scaling(int, double);

  /// Get scaling of n-th input vector field
  double Scaling(int) const;

  /// Run filter on every voxel
  virtual void Run();

  /// Run filter on single voxel
  virtual void Run(double [3], int, int, int);

  /// Run filter on single voxel and component
  virtual double Run(int, int, int, int);
};


} // namespace mirtk

#endif // MIRTK_DifferenceOfCompositionLieBracketImageFilter3D_H
