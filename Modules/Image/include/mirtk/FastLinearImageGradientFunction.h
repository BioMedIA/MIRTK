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

#ifndef MIRTK_FastLinearImageGradientFunction_H
#define MIRTK_FastLinearImageGradientFunction_H

#include "mirtk/BaseImage.h"
#include "mirtk/ImageGradientFunction.h"


namespace mirtk {


/**
 * Fast linear interpolation of generic image gradient
 *
 * This image gradient evaluation function approximates the image gradient
 * using the derivative of the linear interpolation kernel. The difference
 * to the GenericLinearImageGradientFunction is that the resulting finite
 * differences are not necessarily central. In the extreme case, where the
 * gradient is evaluated at a voxel center, the gradient computation
 * corresponds to a forward difference scheme. When the evaluation point is in
 * between two voxel centers, however, the resulting gradient interpolation is
 * a central difference. Side effects resulting from this are generally only
 * noticeable at coarse image resolutions.
 */
template <class TImage>
class GenericFastLinearImageGradientFunction
: public GenericImageGradientFunction<TImage>
{
  mirtkGenericGradientInterpolatorMacro(
    GenericFastLinearImageGradientFunction,
    Interpolation_FastLinear
  );

public:

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Default constructor
  GenericFastLinearImageGradientFunction();

  /// Destructor
  virtual ~GenericFastLinearImageGradientFunction();

  /// Initialize interpolation function
  virtual void Initialize(bool = false);

  // ---------------------------------------------------------------------------
  // Domain checks

  /// Returns interval of discrete image indices whose values are needed for
  /// interpolation of the image value at a given continuous coordinate
  virtual void BoundingInterval(double, int &, int &) const;

  // ---------------------------------------------------------------------------
  // Linear interpolation weights

  /// Returns truncated integral part of floating point
  ///
  /// \param[in]  x Floating point number.
  /// \param[out] w Linear interpolation weights for closest lattice nodes.
  static int ComputeWeights(double, Real[2]);

  // ---------------------------------------------------------------------------
  // Evaluation

  /// Get gradient of given 2D image at arbitrary location (in pixels)
  ///
  /// This function is used to interpolate the image gradient at arbitrary
  /// locations when no extrapolator was set.
  GradientType Get2D(double, double, double = 0, double = 0) const;

  /// Get gradient of given 2D image at arbitrary location (in pixels)
  ///
  /// This function is used to only interpolate foreground image values.
  /// If fully outside the foreground region, the _DefaultValue is returned.
  GradientType GetWithPadding2D(double, double, double = 0, double = 0) const;

  /// Get gradient of given 2D image at arbitrary location (in pixels)
  ///
  /// If the location is inside the finite domain of the image, an actual image
  /// instance can be passed as first argument directly such as an instance of
  /// GenericImage. Otherwise, an image function which extends the finite
  /// image domain to an infinite lattice is needed, i.e., an instance of a
  /// subclass of ExtrapolateImageFunction.
  template <class TOtherImage>
  GradientType Get2D(const TOtherImage *, double, double, double = 0, double = 0) const;

  /// Get gradient of given 2D image at arbitrary location (in pixels)
  ///
  /// This function is used to only interpolate foreground image values.
  /// If fully outside the foreground region, the _DefaultValue is returned.
  ///
  /// If the location is inside the finite domain of the image, an actual image
  /// instance can be passed as first argument directly such as an instance of
  /// GenericImage. Otherwise, an image function which extends the finite
  /// image domain to an infinite lattice is needed, i.e., an instance of a
  /// subclass of ExtrapolateImageFunction.
  template <class TOtherImage>
  GradientType GetWithPadding2D(const TOtherImage *, double, double, double = 0, double = 0) const;

  /// Get gradient of given 3D image at arbitrary location (in pixels)
  ///
  /// This function is used to interpolate the image value at arbitrary
  /// locations when no extrapolator was set.
  GradientType Get3D(double, double, double = 0, double = 0) const;

  /// Get gradient of given 3D image at arbitrary location (in pixels)
  ///
  /// This function is used to only interpolate foreground image values.
  /// If fully outside the foreground region, the _DefaultValue is returned.
  GradientType GetWithPadding3D(double, double, double = 0, double = 0) const;

  /// Get gradient of given 3D image at arbitrary location (in pixels)
  ///
  /// If the location is inside the finite domain of the image, an actual image
  /// instance can be passed as first argument directly such as an instance of
  /// GenericImage. Otherwise, an image function which extends the finite
  /// image domain to an infinite lattice is needed, i.e., an instance of a
  /// subclass of ExtrapolateImageFunction.
  template <class TOtherImage>
  GradientType Get3D(const TOtherImage *, double, double, double = 0, double = 0) const;

  /// Get gradient of given 3D image at arbitrary location (in pixels)
  ///
  /// This function is used to only interpolate foreground image values.
  /// If fully outside the foreground region, the _DefaultValue is returned.
  ///
  /// If the location is inside the finite domain of the image, an actual image
  /// instance can be passed as first argument directly such as an instance of
  /// GenericImage. Otherwise, an image function which extends the finite
  /// image domain to an infinite lattice is needed, i.e., an instance of a
  /// subclass of ExtrapolateImageFunction.
  template <class TOtherImage>
  GradientType GetWithPadding3D(const TOtherImage *, double, double, double = 0, double = 0) const;

  /// Get gradient of given 3D image at arbitrary location (in pixels)
  ///
  /// This function is used to interpolate the image value at arbitrary
  /// locations when no extrapolator was set.
  GradientType Get4D(double, double, double = 0, double = 0) const;

  /// Get gradient of given 4D image at arbitrary location (in pixels)
  ///
  /// This function is used to only interpolate foreground image values.
  /// If fully outside the foreground region, the _DefaultValue is returned.
  GradientType GetWithPadding4D(double, double, double = 0, double = 0) const;

  /// Get gradient of given 4D image at arbitrary location (in pixels)
  ///
  /// If the location is inside the finite domain of the image, an actual image
  /// instance can be passed as first argument directly such as an instance of
  /// GenericImage. Otherwise, an image function which extends the finite
  /// image domain to an infinite lattice is needed, i.e., an instance of a
  /// subclass of ExtrapolateImageFunction.
  template <class TOtherImage>
  GradientType Get4D(const TOtherImage *, double, double, double = 0, double = 0) const;

  /// Get gradient of given 4D image at arbitrary location (in pixels)
  ///
  /// This function is used to only interpolate foreground image values.
  /// If fully outside the foreground region, the _DefaultValue is returned.
  ///
  /// If the location is inside the finite domain of the image, an actual image
  /// instance can be passed as first argument directly such as an instance of
  /// GenericImage. Otherwise, an image function which extends the finite
  /// image domain to an infinite lattice is needed, i.e., an instance of a
  /// subclass of ExtrapolateImageFunction.
  template <class TOtherImage>
  GradientType GetWithPadding4D(const TOtherImage *, double, double, double = 0, double = 0) const;

  /// Get gradient of given image at arbitrary location (in pixels)
  ///
  /// This function is used to interpolate the image value at arbitrary
  /// locations when no extrapolator was set.
  virtual GradientType Get(double, double, double = 0, double = 0) const;

  /// Get value of given image at arbitrary location (in pixels)
  ///
  /// This function is used to only interpolate foreground image values.
  /// If fully outside the foreground region, the _DefaultValue is returned.
  virtual GradientType GetWithPadding(double, double, double = 0, double = 0) const;

  /// Get gradient of given image at arbitrary location (in pixels)
  ///
  /// If the location is inside the finite domain of the image, an actual image
  /// instance can be passed as first argument directly such as an instance of
  /// GenericImage. Otherwise, an image function which extends the finite
  /// image domain to an infinite lattice is needed, i.e., an instance of a
  /// subclass of ExtrapolateImageFunction.
  template <class TOtherImage>
  GradientType Get(const TOtherImage *, double, double, double = 0, double = 0) const;

  /// Get gradient of given image at arbitrary location (in pixels)
  ///
  /// This function is used to only interpolate foreground image values.
  /// If fully outside the foreground region, the _DefaultValue is returned.
  ///
  /// If the location is inside the finite domain of the image, an actual image
  /// instance can be passed as first argument directly such as an instance of
  /// GenericImage. Otherwise, an image function which extends the finite
  /// image domain to an infinite lattice is needed, i.e., an instance of a
  /// subclass of ExtrapolateImageFunction.
  template <class TOtherImage>
  GradientType GetWithPadding(const TOtherImage *, double, double, double = 0, double = 0) const;

  /// Evaluate generic image gradient without handling boundary conditions
  ///
  /// This version is faster than EvaluateOutside, but is only defined inside
  /// the domain for which all image values required for interpolation are
  /// defined and thus require no extrapolation of the finite image.
  virtual GradientType GetInside(double, double, double = 0, double = 0) const;

  /// Evaluate generic image at an arbitrary location (in pixels)
  virtual GradientType GetOutside(double, double, double = 0, double = 0) const;

  /// Evaluate generic image gradient without handling boundary conditions
  ///
  /// If the location is partially inside the foreground region of the image,
  /// only the foreground values are interpolated. Otherwise, the _DefaultValue
  /// is returned.
  ///
  /// This version is faster than GetWithPaddingOutside, but is only defined
  /// inside the domain for which all image values required for interpolation
  /// are defined and thus require no extrapolation of the finite image.
  virtual GradientType GetWithPaddingInside(double, double, double = 0, double = 0) const;

  /// Evaluate generic image gradient at an arbitrary location (in pixels)
  ///
  /// If the location is partially inside the foreground region of the image,
  /// only the foreground values are interpolated. Otherwise, the _DefaultValue
  /// is returned.
  virtual GradientType GetWithPaddingOutside(double, double, double = 0, double = 0) const;

};

/**
 * Fast linear interpolation of any scalar image gradient
 */
class FastLinearImageGradientFunction
: public GenericFastLinearImageGradientFunction<BaseImage>
{
  mirtkObjectMacro(FastLinearImageGradientFunction);

public:

  /// Constructor
  FastLinearImageGradientFunction() {}

};


} // namespace mirtk

#endif // MIRTK_FastLinearImageGradientFunction_H 
