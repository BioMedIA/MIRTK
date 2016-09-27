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

#ifndef MIRTK_NearestNeighborInterpolateImageFunction_H
#define MIRTK_NearestNeighborInterpolateImageFunction_H

#include "mirtk/InterpolateImageFunction.h"


namespace mirtk {


/**
 * Nearest neighbor interpolation of generic image
 */
template <class TImage>
class GenericNearestNeighborInterpolateImageFunction
: public GenericInterpolateImageFunction<TImage>
{
  mirtkGenericInterpolatorMacro(
    GenericNearestNeighborInterpolateImageFunction,
    Interpolation_NN
  );

public:

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Default constructor
  GenericNearestNeighborInterpolateImageFunction();

  /// Destructor
  virtual ~GenericNearestNeighborInterpolateImageFunction();

  // ---------------------------------------------------------------------------
  // Domain checks

  /// Returns interval of discrete image indices whose values are needed for
  /// interpolation of the image value at a given continuous coordinate
  virtual void BoundingInterval(double, int &, int &) const;

  // ---------------------------------------------------------------------------
  // Evaluation

  /// Get value of given image at arbitrary location (in pixels)
  ///
  /// This function is used to interpolate the image value at arbitrary
  /// locations when no extrapolator was set. If the nearest image location is
  /// outside the finite domain of the image, the _DefaultValue is returned.
  VoxelType Get(double, double, double = 0, double = 0) const;

  /// Get value of given image at arbitrary location (in pixels)
  ///
  /// This function is used to only interpolate foreground image values.
  /// If fully outside the foreground region, the _DefaultValue is returned.
  VoxelType GetWithPadding(double, double, double = 0, double = 0) const;

  /// Get value of given image at arbitrary location (in pixels)
  ///
  /// If the location is inside the finite domain of the image, an actual image
  /// instance can be passed as first argument directly such as an instance of
  /// GenericImage. Otherwise, an image function which extends the finite
  /// image domain to an infinite lattice is needed, i.e., an instance of a
  /// subclass of ExtrapolateImageFunction.
  template <class TOtherImage> typename TOtherImage::VoxelType
  Get(const TOtherImage *, double, double, double = 0, double = 0) const;

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
  template <class TOtherImage> typename TOtherImage::VoxelType
  GetWithPadding(const TOtherImage *, double, double, double = 0, double = 0) const;

  /// Evaluate generic image without handling boundary conditions
  ///
  /// This version is faster than EvaluateOutside, but is only defined inside
  /// the domain for which all image values required for interpolation are
  /// defined and thus require no extrapolation of the finite image.
  virtual VoxelType GetInside(double, double, double = 0, double = 0) const;

  /// Evaluate generic image at an arbitrary location (in pixels)
  virtual VoxelType GetOutside(double, double, double = 0, double = 0) const;

  /// Evaluate generic image without handling boundary conditions
  ///
  /// If the location is partially inside the foreground region of the image,
  /// only the foreground values are interpolated. Otherwise, the _DefaultValue
  /// is returned.
  ///
  /// This version is faster than GetWithPaddingOutside, but is only defined
  /// inside the domain for which all image values required for interpolation
  /// are defined and thus require no extrapolation of the finite image.
  virtual VoxelType GetWithPaddingInside(double, double, double = 0, double = 0) const;

  /// Evaluate generic image at an arbitrary location (in pixels)
  ///
  /// If the location is partially inside the foreground region of the image,
  /// only the foreground values are interpolated. Otherwise, the _DefaultValue
  /// is returned.
  virtual VoxelType GetWithPaddingOutside(double, double, double = 0, double = 0) const;

};


/**
 * Nearest neighbor interpolation of any scalar image
 */
class NearestNeighborInterpolateImageFunction
: public GenericNearestNeighborInterpolateImageFunction<BaseImage>
{
  mirtkObjectMacro(NearestNeighborInterpolateImageFunction);

public:

  /// Constructor
  NearestNeighborInterpolateImageFunction() {}

};


} // namespace mirtk

#endif // MIRTK_NearestNeighorInterpolateImageFunction_H
