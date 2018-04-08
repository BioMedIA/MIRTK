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

#ifndef MIRTK_FastCubicBSplineInterpolateImageFunction_H
#define MIRTK_FastCubicBSplineInterpolateImageFunction_H

#include "mirtk/BaseImage.h"
#include "mirtk/InterpolateImageFunction.h"
#include "mirtk/MirrorExtrapolateImageFunction.h"


namespace mirtk {


// Jacobian matrix type
class Matrix;

// Interpolation kernel
template <class T> class BSpline;


/**
 * Fast cubic B-spline interpolation of generic image
 *
 * Note that the coefficient, i.e., input voxel type may also be a vector
 * type such as Vector3D in particular. This is required by the B-spline
 * free-form transformations. If the input image contains the cubic B-spline
 * interpolation coefficients already, make sure that the VoxelType template
 * argument matches the voxel type of the input image to avoid an unnecessary
 * copy of the input image. This image function will use the memory of the
 * input image directly to access the coefficients in this case.
 */
template <class TImage>
class GenericFastCubicBSplineInterpolateImageFunction
: public GenericInterpolateImageFunction<TImage>
{
  mirtkGenericInterpolatorMacro(
    GenericFastCubicBSplineInterpolateImageFunction,
    Interpolation_FastCubicBSpline
  );

  // ---------------------------------------------------------------------------
  // Types

public:

  typedef GenericImage<RealType>                                   CoefficientImage;
  typedef GenericExtrapolateImageFunction<CoefficientImage>        CoefficientExtrapolator;
  typedef GenericMirrorExtrapolateImageFunction<CoefficientImage>  DefaultExtrapolator;
  typedef BSpline<Real>                                            Kernel;

  // ---------------------------------------------------------------------------
  // Attributes

  /// Input image contains spline coefficients
  mirtkAttributeMacro(bool, UseInputCoefficients);

  /// Image of spline coefficients
  mirtkAttributeMacro(CoefficientImage, Coefficient);

  /// Infinite discrete coefficient image obtained by extrapolation
  mirtkAggregateMacro(CoefficientExtrapolator, InfiniteCoefficient);

protected:

  /// Strides for fast iteration over coefficient image
  int _s2, _s3, _s4;

public:

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Constructor
  GenericFastCubicBSplineInterpolateImageFunction();

  /// Destructor
  virtual ~GenericFastCubicBSplineInterpolateImageFunction();

  /// Initialize image function
  virtual void Initialize(bool = false);

  /// Update spline coefficients
  virtual void Update();

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
  template <class TOtherImage, class TCoefficient> typename TCoefficient::VoxelType
  GetWithPadding2D(const TOtherImage *, const TCoefficient *,
                   double, double, double = 0, double = 0) const;

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
  template <class TOtherImage, class TCoefficient> typename TCoefficient::VoxelType
  GetWithPadding3D(const TOtherImage *, const TCoefficient *,
                   double, double, double = 0, double = 0) const;

  /// Get value of given 4D image at arbitrary location (in pixels)
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
  template <class TOtherImage, class TCoefficient> typename TCoefficient::VoxelType
  GetWithPadding4D(const TOtherImage *, const TCoefficient *,
                   double, double, double = 0, double = 0) const;

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
  template <class TOtherImage, class TCoefficient> typename TCoefficient::VoxelType
  GetWithPadding(const TOtherImage *, const TCoefficient *,
                 double, double, double = 0, double = 0) const;

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
  /// This function is used to only interpolate foreground image values.
  /// If fully outside the foreground region, the _DefaultValue is returned.
  ///
  /// This version is faster than GetWithPaddingOutside, but is only defined
  /// inside the domain for which all image values required for interpolation
  /// are defined and thus require no extrapolation of the finite image.
  virtual VoxelType GetWithPaddingInside(double, double, double = 0, double = 0) const;

  /// Evaluate generic image at an arbitrary location (in pixels)
  ///
  /// This function is used to only interpolate foreground image values.
  /// If fully outside the foreground region, the _DefaultValue is returned.
  virtual VoxelType GetWithPaddingOutside(double, double, double = 0, double = 0) const;

  // ---------------------------------------------------------------------------
  // Derivative evaluation

  /// Get 1st order derivatives of given 2D image at arbitrary location (in pixels)
  ///
  /// This function is used when no extrapolator was set.
  void Jacobian2D(Matrix &, double, double, double = 0, double = 0) const;

  /// Get 1st order derivatives of given 2D image at arbitrary location (in pixels)
  ///
  /// If the location is inside the finite domain of the image, an actual image
  /// instance can be passed as first argument directly such as an instance of
  /// GenericImage. Otherwise, an image function which extends the finite
  /// image domain to an infinite lattice is needed, i.e., an instance of a
  /// subclass of ExtrapolateImageFunction.
  template <class TOtherImage>
  void Jacobian2D(Matrix &, const TOtherImage *, double, double, double = 0, double = 0) const;

  /// Get 1st order derivatives of given 3D image at arbitrary location (in pixels)
  ///
  /// This function is used when no extrapolator was set.
  void Jacobian3D(Matrix &, double, double, double = 0, double = 0) const;

  /// Get 1st order derivatives of given 3D image at arbitrary location (in pixels)
  ///
  /// If the location is inside the finite domain of the image, an actual image
  /// instance can be passed as first argument directly such as an instance of
  /// GenericImage. Otherwise, an image function which extends the finite
  /// image domain to an infinite lattice is needed, i.e., an instance of a
  /// subclass of ExtrapolateImageFunction.
  template <class TOtherImage>
  void Jacobian3D(Matrix &, const TOtherImage *, double, double, double = 0, double = 0) const;

  /// Get 1st order derivatives of given 4D image at arbitrary location (in pixels)
  ///
  /// This function is used when no extrapolator was set.
  void Jacobian4D(Matrix &, double, double, double = 0, double = 0) const;

  /// Get 1st order derivatives of given 4D image at arbitrary location (in pixels)
  ///
  /// If the location is inside the finite domain of the image, an actual image
  /// instance can be passed as first argument directly such as an instance of
  /// GenericImage. Otherwise, an image function which extends the finite
  /// image domain to an infinite lattice is needed, i.e., an instance of a
  /// subclass of ExtrapolateImageFunction.
  template <class TOtherImage>
  void Jacobian4D(Matrix &, const TOtherImage *, double, double, double = 0, double = 0) const;

  /// Get 1st order derivatives of given image at arbitrary location (in pixels)
  void Jacobian(Matrix &, double, double, double = 0, double = 0) const;

  /// Get 1st order derivatives of given image at arbitrary location (in pixels)
  ///
  /// If the location is inside the finite domain of the image, an actual image
  /// instance can be passed as first argument directly such as an instance of
  /// GenericImage. Otherwise, an image function which extends the finite
  /// image domain to an infinite lattice is needed, i.e., an instance of a
  /// subclass of ExtrapolateImageFunction.
  template <class TOtherImage>
  void Jacobian(Matrix &, const TOtherImage *, double, double, double = 0, double = 0) const;

  /// Get 1st order derivatives of given image at arbitrary location (in pixels)
  virtual void EvaluateJacobianInside(Matrix &, double, double, double = 0, double = 0) const;

  /// Get 1st order derivatives of given image at arbitrary location (in pixels)
  virtual void EvaluateJacobianOutside(Matrix &, double, double, double = 0, double = 0) const;

};


/**
 * Fast cubic B-spline interpolation of any scalar image
 */
class FastCubicBSplineInterpolateImageFunction
: public GenericFastCubicBSplineInterpolateImageFunction<BaseImage>
{
  mirtkObjectMacro(FastCubicBSplineInterpolateImageFunction);

public:

  /// Constructor
  FastCubicBSplineInterpolateImageFunction() {}

};


} // namespace mirtk

#endif // MIRTK_FastCubicBSplineInterpolateImageFunction_H
