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

#ifndef MIRTK_CubicBSplineInterpolateImageFunction4D_H
#define MIRTK_CubicBSplineInterpolateImageFunction4D_H

#include "mirtk/BaseImage.h"
#include "mirtk/CubicBSplineInterpolateImageFunction.h"


namespace mirtk {


/**
 * Cubic B-spline interpolation of generic 4D image
 */
template <class TImage>
class GenericCubicBSplineInterpolateImageFunction4D
: public GenericCubicBSplineInterpolateImageFunction<TImage>
{
  mirtkObjectMacro(GenericCubicBSplineInterpolateImageFunction4D);
  mirtkGenericInterpolatorTypes(GenericCubicBSplineInterpolateImageFunction);

public:

  /// Default constructor
  GenericCubicBSplineInterpolateImageFunction4D();

  /// Get value of given image at arbitrary location (in pixels)
  ///
  /// This function is used to interpolate the image value at arbitrary
  /// locations when no extrapolator was set.
  VoxelType Get(double, double, double, double) const;

  /// Get value of given image at arbitrary location (in pixels)
  ///
  /// This function is used to only interpolate foreground image values.
  /// If fully outside the foreground region, the _DefaultValue is returned.
  VoxelType GetWithPadding(double, double, double, double) const;

  /// Get value of given image at arbitrary location (in pixels)
  ///
  /// If the location is inside the finite domain of the image, an actual image
  /// instance can be passed as first argument directly such as an instance of
  /// GenericImage. Otherwise, an image function which extends the finite
  /// image domain to an infinite lattice is needed, i.e., an instance of a
  /// subclass of ExtrapolateImageFunction.
  template <class TOtherImage> typename TOtherImage::VoxelType
  Get(const TOtherImage *, double, double, double, double) const;

  /// Get value of given image at arbitrary location (in pixels)
  ///
  /// This function is used to only interpolate foreground image values.
  /// If fully outside the foreground region, the _DefaultValue is returned.
  ///
  /// If the location is inside the finite domain of the image, an actual image
  /// instance can be passed as first argument directly such as an instance of
  /// GenericImage. Otherwise, an image function which extends the finite
  /// image domain to an infinite lattice is needed, i.e., an instance of a
  /// subclass of ExtrapolateImageFunction.
  template <class TOtherImage, class TCoefficient> typename TCoefficient::VoxelType
  GetWithPadding(const TOtherImage *, const TCoefficient *,
                 double, double, double, double) const;

  /// Evaluate generic image without handling boundary conditions
  ///
  /// This version is faster than EvaluateOutside, but is only defined inside
  /// the domain for which all image values required for interpolation are
  /// defined and thus require no extrapolation of the finite image.
  virtual VoxelType GetInside(double, double, double, double) const;

  /// Evaluate generic image at an arbitrary location (in pixels)
  virtual VoxelType GetOutside(double, double, double, double) const;

  /// Evaluate generic image without handling boundary conditions
  ///
  /// This function is used to only interpolate foreground image values.
  /// If fully outside the foreground region, the _DefaultValue is returned.
  ///
  /// This version is faster than GetWithPaddingOutside, but is only defined
  /// inside the domain for which all image values required for interpolation
  /// are defined and thus require no extrapolation of the finite image.
  virtual VoxelType GetWithPaddingInside(double, double, double, double) const;

  /// Evaluate generic image at an arbitrary location (in pixels)
  ///
  /// This function is used to only interpolate foreground image values.
  /// If fully outside the foreground region, the _DefaultValue is returned.
  virtual VoxelType GetWithPaddingOutside(double, double, double, double) const;

};

/**
 * Cubic B-spline interpolation of any 4D image
 */
class CubicBSplineInterpolateImageFunction4D
: public GenericCubicBSplineInterpolateImageFunction4D<BaseImage>
{
  mirtkObjectMacro(CubicBSplineInterpolateImageFunction4D);

public:

  /// Constructor
  CubicBSplineInterpolateImageFunction4D() {}

};


} // namespace mirtk

#endif // MIRTK_CubicBSplineInterpolateImageFunction4D_H 
