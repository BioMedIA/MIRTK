/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2010-2015 Imperial College London
 * Copyright 2010      Wenzhe Shi
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

#ifndef MIRTK_ShapeBasedInterpolateImageFunction_H
#define MIRTK_ShapeBasedInterpolateImageFunction_H

#include "mirtk/InterpolateImageFunction.h"

#include "mirtk/Math.h"
#include "mirtk/GenericImage.h"
#include "mirtk/LinearInterpolateImageFunction.h"
#include "mirtk/NearestNeighborInterpolateImageFunction.h"


namespace mirtk {


/**
 * Shape based interpolation of 3D scalar images
 *
 * This class defines and implements a shape based interpolation of images.
 */
class ShapeBasedInterpolateImageFunction : public InterpolateImageFunction
{
  mirtkInterpolatorMacro(
    ShapeBasedInterpolateImageFunction,
    Interpolation_SBased
  );

  // ---------------------------------------------------------------------------
  // Attributes

protected:

  /// Threshold map for input image
  RealImage _tinput;

  /// Distance map for input image
  RealImage _dmap;

  /// Isotropic resampled distance map for input image with linear interpolation
  RealImage _rdmap;

  /// Isotropic resampled Input image after shape based interpolation
  RealImage _rinput;

  /// Isotropic resampled cache image for refine procedure
  RealImage _rcdmap;

  /// Image function for nearest neighbor interpolation of resampled input
  NearestNeighborInterpolateImageFunction _nn_interpolator;

  /// Image function for linear interpolation of resampled input
  LinearInterpolateImageFunction _linear_interpolator;

  /// Initialize second step, fix the union property
  void Refine();

public:

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Constructor
  ShapeBasedInterpolateImageFunction();

  /// Destructor
  virtual ~ShapeBasedInterpolateImageFunction();

  /// Initialize
  virtual void Initialize(bool = false);

  // ---------------------------------------------------------------------------
  // Domain checks

  /// Returns interval of discrete image indices whose values are needed for
  /// interpolation of the image value at a given continuous coordinate
  virtual void BoundingInterval(double, int &, int &) const;

  // ---------------------------------------------------------------------------
  // Evaluation

  // import non-overridden overloads
  using InterpolateImageFunction::EvaluateInside;
  using InterpolateImageFunction::EvaluateOutside;
  using InterpolateImageFunction::EvaluateWithPaddingInside;
  using InterpolateImageFunction::EvaluateWithPaddingOutside;

  /// Evaluate scalar image without handling boundary conditions
  ///
  /// This version is faster than EvaluateOutside, but is only defined inside
  /// the domain for which all image values required for interpolation are
  /// defined and thus require no extrapolation of the finite image.
  virtual double EvaluateInside(double, double, double, double = 0) const;

  /// Evaluate scalar image at an arbitrary location (in pixels)
  virtual double EvaluateOutside(double, double, double, double = 0) const;

  /// Evaluate scalar image at an arbitrary location (in pixels)
  ///
  /// If the location is partially inside the foreground region of the image,
  /// only the foreground values are interpolated. Otherwise, the _DefaultValue
  /// is returned.
  ///
  /// This version is faster than EvaluateWithPaddingOutside, but is only defined
  /// inside the domain for which all image values required for interpolation are
  /// defined and thus require no extrapolation of the finite image.
  virtual double EvaluateWithPaddingInside(double, double, double, double = 0) const;

  /// Evaluate scalar image at an arbitrary location (in pixels)
  ///
  /// If the location is partially inside the foreground region of the image,
  /// only the foreground values are interpolated. Otherwise, the _DefaultValue
  /// is returned.
  virtual double EvaluateWithPaddingOutside(double, double, double, double = 0) const;

  /// Evaluate vector image without handling boundary conditions
  ///
  /// This version is faster than EvaluateOutside, but is only defined inside
  /// the domain for which all image values required for interpolation are
  /// defined and thus require no extrapolation of the finite image.
  ///
  /// \attention This interpolator is implemented only for scalar images!
  virtual void EvaluateInside(Vector &, double, double, double, double = 0) const;

  /// Evaluate vector image at an arbitrary location (in pixels)
  ///
  /// \attention This interpolator is implemented only for scalar images!
  virtual void EvaluateOutside(Vector &, double, double, double, double = 0) const;

  /// Evaluate vector image at an arbitrary location (in pixels)
  ///
  /// If the location is partially inside the foreground region of the image,
  /// only the foreground values are interpolated. Otherwise, a vector set to
  /// the _DefaultValue is returned.
  ///
  /// This version is faster than EvaluateWithPaddingOutside, but is only defined
  /// inside the domain for which all image values required for interpolation are
  /// defined and thus require no extrapolation of the finite image.
  ///
  /// \attention This interpolator is implemented only for scalar images!
  virtual void EvaluateWithPaddingInside(Vector &, double, double, double, double = 0) const;

  /// Evaluate vector image at an arbitrary location (in pixels)
  ///
  /// If the location is partially inside the foreground region of the image,
  /// only the foreground values are interpolated. Otherwise, a vector set to
  /// the _DefaultValue is returned.
  ///
  /// \attention This interpolator is implemented only for scalar images!
  virtual void EvaluateWithPaddingOutside(Vector &, double, double, double, double = 0) const;

  /// Evaluate the filter at arbitrary image location (in pixels)
  virtual double EvaluateLinear(double, double, double, double = 0) const;

  /// Evaluate the filter at an arbitrary image location (in pixels) without
  /// handling boundary conditions. This version is faster than the method
  /// above, but is only defined inside the image domain.
  virtual double EvaluateInsideLinear(double, double, double, double = 0) const;

  /// Evaluate the filter at a boundary image location (in pixels)
  virtual double EvaluateOutsideLinear(double, double, double, double = 0) const;

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Domain checks
// =============================================================================

// -----------------------------------------------------------------------------
inline void ShapeBasedInterpolateImageFunction
::BoundingInterval(double x, int &i, int &I) const
{
  i = I = iround(x);
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
inline void ShapeBasedInterpolateImageFunction
::EvaluateInside(Vector &v, double x, double y, double z, double t) const
{
  v = this->EvaluateInside(x, y, z, t);
}

// -----------------------------------------------------------------------------
inline void ShapeBasedInterpolateImageFunction
::EvaluateOutside(Vector &v, double x, double y, double z, double t) const
{
  v = this->EvaluateOutside(x, y, z, t);
}

// -----------------------------------------------------------------------------
inline void ShapeBasedInterpolateImageFunction
::EvaluateWithPaddingInside(Vector &v, double x, double y, double z, double t) const
{
  v = this->EvaluateWithPaddingInside(x, y, z, t);
}

// -----------------------------------------------------------------------------
inline void ShapeBasedInterpolateImageFunction
::EvaluateWithPaddingOutside(Vector &v, double x, double y, double z, double t) const
{
  v = this->EvaluateWithPaddingOutside(x, y, z, t);
}


} // namespace mirtk

#endif // MIRTK_ShapeBasedInterpolateImageFunction_H
