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

#ifndef MIRTK_LinearInterpolateImageFunction_H
#define MIRTK_LinearInterpolateImageFunction_H

#include "mirtk/InterpolateImageFunction.h"
#include "mirtk/BaseImage.h"


namespace mirtk {


/**
 * Linear interpolation of generic image
 */
template <class TImage>
class GenericLinearInterpolateImageFunction
: public GenericInterpolateImageFunction<TImage>
{
  mirtkGenericInterpolatorMacro(
    GenericLinearInterpolateImageFunction,
    Interpolation_Linear
  );

  /// Offsets for fast pixel access
  int _Offset[16];

public:

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Default constructor
  GenericLinearInterpolateImageFunction();

  /// Destructor
  virtual ~GenericLinearInterpolateImageFunction();

  /// Initialize interpolation function
  virtual void Initialize(bool = false);

  // ---------------------------------------------------------------------------
  // Linear interpolation weights

  /// Returns truncated integral part of floating point
  ///
  /// \param[in]  x Floating point number.
  /// \param[out] w Linear interpolation weights for closest lattice nodes.
  static int ComputeWeights(double, Real[2]);

  // ---------------------------------------------------------------------------
  // Domain checks

  /// Returns interval of discrete image indices whose values are needed for
  /// interpolation of the image value at a given continuous coordinate
  virtual void BoundingInterval(double, int &, int &) const;

  // ---------------------------------------------------------------------------
  // Evaluation

  /// Get value of given 2D image at arbitrary location (in pixels)
  ///
  /// This function is used to interpolate the image value at arbitrary
  /// locations when no extrapolator was set.
  VoxelType Get2D(double, double, double = 0, double = 0) const;

  /// Get value of given 2D image at arbitrary location (in pixels)
  ///
  /// This function is used to only interpolate foreground image values.
  /// If fully outside the foreground region, the _DefaultValue is returned.
  VoxelType GetWithPadding2D(double, double, double = 0, double = 0) const;

  /// Get value of given 2D image at arbitrary location (in pixels)
  ///
  /// If the location is inside the finite domain of the image, an actual image
  /// instance can be passed as first argument directly such as an instance of
  /// GenericImage. Otherwise, an image function which extends the finite
  /// image domain to an infinite lattice is needed, i.e., an instance of a
  /// subclass of ExtrapolateImageFunction.
  template <class TOtherImage> typename TOtherImage::VoxelType
  Get2D(const TOtherImage *, double, double, double = 0, double = 0) const;

  /// Get value of given 2D image at arbitrary location (in pixels)
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
  GetWithPadding2D(const TOtherImage *, double, double, double = 0, double = 0) const;

  /// Get value of given 3D image at arbitrary location (in pixels)
  ///
  /// This function is used to interpolate the image value at arbitrary
  /// locations when no extrapolator was set.
  VoxelType Get3D(double, double, double = 0, double = 0) const;

  /// Get value of given 3D image at arbitrary location (in pixels)
  ///
  /// This function is used to only interpolate foreground image values.
  /// If fully outside the foreground region, the _DefaultValue is returned.
  VoxelType GetWithPadding3D(double, double, double = 0, double = 0) const;

  /// Get value of given 3D image at arbitrary location (in pixels)
  ///
  /// If the location is inside the finite domain of the image, an actual image
  /// instance can be passed as first argument directly such as an instance of
  /// GenericImage. Otherwise, an image function which extends the finite
  /// image domain to an infinite lattice is needed, i.e., an instance of a
  /// subclass of ExtrapolateImageFunction.
  template <class TOtherImage> typename TOtherImage::VoxelType
  Get3D(const TOtherImage *, double, double, double = 0, double = 0) const;

  /// Get value of given 3D image at arbitrary location (in pixels)
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
  GetWithPadding3D(const TOtherImage *, double, double, double = 0, double = 0) const;

  /// Get value of given 3D image at arbitrary location (in pixels)
  ///
  /// This function is used to interpolate the image value at arbitrary
  /// locations when no extrapolator was set.
  VoxelType Get4D(double, double, double = 0, double = 0) const;

  /// Get value of given 4D image at arbitrary location (in pixels)
  ///
  /// This function is used to only interpolate foreground image values.
  /// If fully outside the foreground region, the _DefaultValue is returned.
  VoxelType GetWithPadding4D(double, double, double = 0, double = 0) const;

  /// Get value of given 4D image at arbitrary location (in pixels)
  ///
  /// If the location is inside the finite domain of the image, an actual image
  /// instance can be passed as first argument directly such as an instance of
  /// GenericImage. Otherwise, an image function which extends the finite
  /// image domain to an infinite lattice is needed, i.e., an instance of a
  /// subclass of ExtrapolateImageFunction.
  template <class TOtherImage> typename TOtherImage::VoxelType
  Get4D(const TOtherImage *, double, double, double = 0, double = 0) const;

  /// Get value of given 4D image at arbitrary location (in pixels)
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
  GetWithPadding4D(const TOtherImage *, double, double, double = 0, double = 0) const;

  /// Get value of given image at arbitrary location (in pixels)
  ///
  /// This function is used to interpolate the image value at arbitrary
  /// locations when no extrapolator was set.
  VoxelType Get(double, double, double = 0, double = 0) const;

  /// Get value of given image at arbitrary location (in pixels)
  ///
  /// This function is used to only interpolate foreground image values.
  /// If fully outside the foreground region, the _DefaultValue is returned.
  virtual VoxelType GetWithPadding(double, double, double = 0, double = 0) const;

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

  /// Evaluate generic 2D image without handling boundary conditions
  virtual VoxelType GetInside2D(double, double, double = 0, double = 0) const;

  /// Evaluate generic 3D image without handling boundary conditions
  virtual VoxelType GetInside3D(double, double, double = 0, double = 0) const;

  /// Evaluate generic 4D image without handling boundary conditions
  virtual VoxelType GetInside4D(double, double, double = 0, double = 0) const;

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

  // ---------------------------------------------------------------------------
  // Derivative evaluation

  /// Get 1st order derivatives of given image at arbitrary location (in pixels)
  ///
  /// When the image has scalar data type and stores vector components in the
  /// fourth dimension, the derivatives of all components are evaluated when
  /// the t coordinate is set to NaN. Otherwise, only the derivatives of the
  /// specified t component are evaluated.
  virtual void EvaluateJacobianInside(Matrix &, double, double, double = 0, double = NaN) const;

  /// Get 1st order derivatives of given image at arbitrary location (in pixels)
  ///
  /// When the image has scalar data type and stores vector components in the
  /// fourth dimension, the derivatives of all components are evaluated when
  /// the t coordinate is set to NaN. Otherwise, only the derivatives of the
  /// specified t component are evaluated.
  virtual void EvaluateJacobianOutside(Matrix &, double, double, double = 0, double = NaN) const;

  /// Get 1st order derivatives of given image at arbitrary location (in pixels)
  ///
  /// When the image has scalar data type and stores vector components in the
  /// fourth dimension, the derivatives of all components are evaluated when
  /// the t coordinate is set to NaN. Otherwise, only the derivatives of the
  /// specified t component are evaluated.
  virtual void EvaluateJacobianWithPaddingInside(Matrix &, double, double, double = 0, double = NaN) const;

  /// Get 1st order derivatives of given image at arbitrary location (in pixels)
  ///
  /// When the image has scalar data type and stores vector components in the
  /// fourth dimension, the derivatives of all components are evaluated when
  /// the t coordinate is set to NaN. Otherwise, only the derivatives of the
  /// specified t component are evaluated.
  virtual void EvaluateJacobianWithPaddingOutside(Matrix &, double, double, double = 0, double = NaN) const;

  /// Get 1st order derivatives of given 3D image at arbitrary location (in pixels)
  ///
  /// This function is used when no extrapolator was set.
  ///
  /// When the image has scalar data type and stores vector components in the
  /// fourth dimension, the derivatives of all components are evaluated when
  /// the t coordinate is set to NaN. Otherwise, only the derivatives of the
  /// specified t component are evaluated.
  void Jacobian2D(Matrix &, double, double, double = 0., double = NaN) const;

  /// Get 1st order derivatives of given 3D image at arbitrary location (in pixels)
  ///
  /// This function is used when no extrapolator was set.
  ///
  /// When the image has scalar data type and stores vector components in the
  /// fourth dimension, the derivatives of all components are evaluated when
  /// the t coordinate is set to NaN. Otherwise, only the derivatives of the
  /// specified t component are evaluated.
  void Jacobian3D(Matrix &, double, double, double = 0., double = NaN) const;

  /// Get 1st order derivatives of given 4D image at arbitrary location (in pixels)
  ///
  /// This function is used when no extrapolator was set.
  void Jacobian4D(Matrix &, double, double, double = 0., double = NaN) const;

  /// Get 1st order derivatives of given image at arbitrary location (in pixels)
  void Jacobian(Matrix &, double, double, double = 0., double = NaN) const;


  /// Get 1st order derivatives of given 3D image at arbitrary location (in pixels)
  ///
  /// Set derivatives to zero outside the image foreground.
  ///
  /// When the image has scalar data type and stores vector components in the
  /// fourth dimension, the derivatives of all components are evaluated when
  /// the t coordinate is set to NaN. Otherwise, only the derivatives of the
  /// specified t component are evaluated.
  void JacobianWithPadding2D(Matrix &, double, double, double = 0., double = NaN) const;

  /// Get 1st order derivatives of given 3D image at arbitrary location (in pixels)
  ///
  /// Set derivatives to zero outside the image foreground.
  ///
  /// When the image has scalar data type and stores vector components in the
  /// fourth dimension, the derivatives of all components are evaluated when
  /// the t coordinate is set to NaN. Otherwise, only the derivatives of the
  /// specified t component are evaluated.
  void JacobianWithPadding3D(Matrix &, double, double, double = 0., double = NaN) const;

  /// Get 1st order derivatives of given 4D image at arbitrary location (in pixels)
  ///
  /// Set derivatives to zero outside the image foreground.
  void JacobianWithPadding4D(Matrix &, double, double, double = 0., double = NaN) const;

  /// Get 1st order derivatives of given image at arbitrary location (in pixels)
  void JacobianWithPadding(Matrix &, double, double, double = 0, double = NaN) const;


  /// Get 1st order derivatives of given 3D image at arbitrary location (in pixels)
  ///
  /// If the location is inside the finite domain of the image, an actual image
  /// instance can be passed as first argument directly such as an instance of
  /// GenericImage. Otherwise, an image function which extends the finite
  /// image domain to an infinite lattice is needed, i.e., an instance of a
  /// subclass of ExtrapolateImageFunction.
  ///
  /// When the image has scalar data type and stores vector components in the
  /// fourth dimension, the derivatives of all components are evaluated when
  /// the t coordinate is set to NaN. Otherwise, only the derivatives of the
  /// specified t component are evaluated.
  template <class TOtherImage>
  void Jacobian2D(Matrix &, const TOtherImage *, double, double, double = 0., double = NaN) const;

  /// Get 1st order derivatives of given 3D image at arbitrary location (in pixels)
  ///
  /// If the location is inside the finite domain of the image, an actual image
  /// instance can be passed as first argument directly such as an instance of
  /// GenericImage. Otherwise, an image function which extends the finite
  /// image domain to an infinite lattice is needed, i.e., an instance of a
  /// subclass of ExtrapolateImageFunction.
  ///
  /// When the image has scalar data type and stores vector components in the
  /// fourth dimension, the derivatives of all components are evaluated when
  /// the t coordinate is set to NaN. Otherwise, only the derivatives of the
  /// specified t component are evaluated.
  template <class TOtherImage>
  void Jacobian3D(Matrix &, const TOtherImage *, double, double, double = 0., double = NaN) const;

  /// Get 1st order derivatives of given 4D image at arbitrary location (in pixels)
  ///
  /// If the location is inside the finite domain of the image, an actual image
  /// instance can be passed as first argument directly such as an instance of
  /// GenericImage. Otherwise, an image function which extends the finite
  /// image domain to an infinite lattice is needed, i.e., an instance of a
  /// subclass of ExtrapolateImageFunction.
  template <class TOtherImage>
  void Jacobian4D(Matrix &, const TOtherImage *, double, double, double = 0., double = NaN) const;

  /// Get 1st order derivatives of given image at arbitrary location (in pixels)
  ///
  /// If the location is inside the finite domain of the image, an actual image
  /// instance can be passed as first argument directly such as an instance of
  /// GenericImage. Otherwise, an image function which extends the finite
  /// image domain to an infinite lattice is needed, i.e., an instance of a
  /// subclass of ExtrapolateImageFunction.
  template <class TOtherImage>
  void Jacobian(Matrix &, const TOtherImage *, double, double, double = 0., double = NaN) const;

};


/**
 * Linear interpolation of any scalar image
 */

class LinearInterpolateImageFunction
: public GenericLinearInterpolateImageFunction<BaseImage>
{
  mirtkObjectMacro(LinearInterpolateImageFunction);

public:

  /// Constructor
  LinearInterpolateImageFunction() {}

};


} // namespace mirtk

#endif // MIRTK_LinearInterpolateImageFunction_H
