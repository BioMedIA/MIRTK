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

#ifndef MIRTK_InterpolateImageFunction_H
#define MIRTK_InterpolateImageFunction_H

#include "mirtk/InterpolationMode.h"
#include "mirtk/Point.h"
#include "mirtk/Vector.h"
#include "mirtk/Voxel.h"
#include "mirtk/VoxelCast.h"
#include "mirtk/BaseImage.h"
#include "mirtk/GenericImage.h"
#include "mirtk/ImageFunction.h"
#include "mirtk/UnaryVoxelFunction.h"

#include "mirtk/ExtrapolateImageFunction.h"


namespace mirtk {


////////////////////////////////////////////////////////////////////////////////
// Abstract interpolation interface
////////////////////////////////////////////////////////////////////////////////

/**
 * Abstract base class for any general image interpolation function
 */
class InterpolateImageFunction : public ImageFunction
{
  mirtkAbstractMacro(InterpolateImageFunction);

  // ---------------------------------------------------------------------------
  // Data members

  /// Number of dimensions to interpolate
  ///
  /// Determined either from input dimensions by default or set to a fixed
  /// number of dimension by specialized subclasses.
  mirtkPublicAttributeMacro(int, NumberOfDimensions);

protected:

  /// Infinite discrete image obtained by extrapolation of finite input image.
  /// Unused by default, i.e., NULL which corresponds to extrapolation mode
  /// Extrapolation_None. If \c NULL, the interpolator has to deal with
  /// boundary cases itself either by only partially interpolating the
  /// available values or returning the _DefaultValue.
  ExtrapolateImageFunction *_InfiniteInput;

  /// Whether infinite discrete image was instantiated by this image function
  bool _InfiniteInputOwner;

  /// Domain of finite input image for which the interpolation is defined
  /// without requiring any extrapolation: [x1, x2]x[y1, y2]x[z1, z2]x[t1, t2]
  double _x1, _y1, _z1, _t1, _x2, _y2, _z2, _t2;

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Default constructor
  InterpolateImageFunction();

public:

  /// Destructor
  virtual ~InterpolateImageFunction();

  /// Construct interpolator with default infinite extension of input image
  static InterpolateImageFunction *New(enum InterpolationMode = Interpolation_Default,
                                       const BaseImage * = NULL);

  /// Construct extrapolator which is compatible with this interpolator
  virtual ExtrapolateImageFunction *New(enum ExtrapolationMode,
                                        const BaseImage * = NULL);

  /// Construct interpolator with specified infinite extension of input image
  ///
  /// The caller is required to set the input, initialize, and destroy the
  /// interpolator only, the extrapolator is initialized and destroyed by the
  /// interpolator unless the extrapolator has been replaced using the setter.
  static InterpolateImageFunction *New(enum InterpolationMode,
                                       enum ExtrapolationMode, const BaseImage * = NULL);

  // ---------------------------------------------------------------------------
  // Initialization

  /// Set input image
  void Input(const BaseImage *);

  /// Get input image
  const BaseImage *Input() const;

  /// Get interpolation mode corresponding to this interpolator
  virtual enum InterpolationMode InterpolationMode() const = 0;

  /// Get extrapolation mode used by this interpolator
  enum ExtrapolationMode ExtrapolationMode() const;

  /// Set extrapolate image function for evaluation outside of image domain
  virtual void Extrapolator(ExtrapolateImageFunction *, bool = false);

  /// Get extrapolate image function for evaluation outside of image domain
  /// or \c NULL if extrapolation mode is \c Extrapolation_None
  ExtrapolateImageFunction *Extrapolator();

  /// Get extrapolate image function for evaluation outside of image domain
  /// or \c NULL if extrapolation mode is \c Extrapolation_None
  const ExtrapolateImageFunction *Extrapolator() const;

  /// Initialize image function
  ///
  /// \note Override the virtual Initialize(bool) member function in subclasses,
  ///       but not this member function which is required by the abstract
  ///       image function interface (irtkImageFunction::Initialize).
  virtual void Initialize();

  /// Initialize image function
  ///
  /// \param[in] coeff Whether input image contains interpolation coefficients
  ///                  already. Otherwise, the interpolate image function will
  ///                  compute these coefficients from the input intensities.
  virtual void Initialize(bool coeff);

  /// Update internal state when input image content has changed
  ///
  /// When the attributes of the input have changed, call Initialize instead.
  /// This function is used for example by B-spline based interpolation functions
  /// to re-compute the spline coefficients that interpolate the input image.
  virtual void Update();

  // ---------------------------------------------------------------------------
  // Lattice

  /// Lattice attributes
  const ImageAttributes &Attributes() const;

  /// Image size along x axis
  int X() const;

  /// Image size along y axis
  int Y() const;

  /// Image size along z axis
  int Z() const;

  /// Image size along t axis
  int T() const;

  /// Image spacing along x axis
  double XSize() const;

  /// Image spacing along y axis
  double YSize() const;

  /// Image spacing along z axis
  double ZSize() const;

  /// Image spacing along t axis
  double TSize() const;

  /// Convert world coordinates (in mm) to image location (in pixels)
  void WorldToImage(double &, double &) const;

  /// Convert world coordinates (in mm) to image location (in pixels)
  void WorldToImage(double &, double &, double &) const;

  /// Convert world coordinates (in mm) to image location (in pixels)
  void WorldToImage(Point &) const;

  /// Convert world coordinates vector (in mm) to image vector (in pixels)
  void WorldToImage(Vector3 &) const;

  /// Convert image location (in pixels) to world coordinates (in mm)
  void ImageToWorld(double &, double &) const;

  /// Convert image location (in pixels) to world coordinates (in mm)
  void ImageToWorld(double &, double &, double &) const;

  /// Convert image location (in pixels) to world coordinates (in mm)
  void ImageToWorld(Point &) const;

  /// Convert image vector (in pixels) to world coordinates (in mm)
  void ImageToWorld(Vector3 &) const;

  // ---------------------------------------------------------------------------
  // Domain checks

  /// Returns the image domain for which this image interpolation function
  /// can be used without handling any form of boundary conditions
  void Inside(double &, double &,
              double &, double &) const;

  /// Returns the image domain for which this image interpolation function
  /// can be used without handling any form of boundary conditions
  void Inside(double &, double &, double &,
              double &, double &, double &) const;

  /// Returns the image domain for which this image interpolation function
  /// can be used without handling any form of boundary conditions
  void Inside(double &, double &, double &, double &,
              double &, double &, double &, double &) const;

  /// Returns interval of discrete image indices whose values are needed for
  /// interpolation of the image value at a given continuous coordinate
  virtual void BoundingInterval(double, int &, int &) const = 0;

  /// Returns discrete boundaries of local 2D image region needed for interpolation
  virtual void BoundingBox(double, double, int &, int &,
                                           int &, int &) const;

  /// Returns discrete boundaries of local 3D image region needed for interpolation
  virtual void BoundingBox(double, double, double, int &, int &, int &,
                                                   int &, int &, int &) const;

  /// Returns discrete boundaries of local 4D image region needed for interpolation
  virtual void BoundingBox(double, double, double, double,
                           int &, int &, int &, int &,
                           int &, int &, int &, int &) const;

  /// Check if the location (in pixels) is inside the domain for which this image
  /// interpolation can be used without handling any form of boundary condition
  bool IsInside(double, double) const;

  /// Check if the location (in pixels) is inside the domain for which this image
  /// interpolation can be used without handling any form of boundary condition
  bool IsInside(double, double, double) const;

  /// Check if the location (in pixels) is inside the domain for which this image
  /// interpolation can be used without handling any form of boundary condition
  bool IsInside(double, double, double, double) const;

  /// Check if the location (in pixels) is inside the domain for which this image
  /// interpolation can be used without handling any form of boundary condition
  bool IsInside(const Point &) const;

  /// Check if the location (in pixels) is inside the domain for which this image
  /// interpolation can be used without handling any form of boundary condition
  bool IsInside(const Point &, double) const;

  /// Check if the location (in pixels) is outside the domain for which this image
  /// interpolation can be used without handling any form of boundary condition
  bool IsOutside(double, double) const;

  /// Check if the location (in pixels) is outside the domain for which this image
  /// interpolation can be used without handling any form of boundary condition
  bool IsOutside(double, double, double) const;

  /// Check if the location (in pixels) is outside the domain for which this image
  /// interpolation can be used without handling any form of boundary condition
  bool IsOutside(double, double, double, double) const;

  /// Check if the location (in pixels) is outside the domain for which this image
  /// interpolation can be used without handling any form of boundary condition
  bool IsOutside(const Point &) const;

  /// Check if the location (in pixels) is outside the domain for which this image
  /// interpolation can be used without handling any form of boundary condition
  bool IsOutside(const Point &, double) const;

  /// Check if the location is fully inside the foreground of the image, i.e.,
  /// including all discrete image locations required for interpolation
  bool IsForeground(double, double) const;

  /// Check if the location is fully inside the foreground of the image, i.e.,
  /// including all discrete image locations required for interpolation
  bool IsForeground(double, double, double) const;

  /// Check if the location is fully inside the foreground of the image, i.e.,
  /// including all discrete image locations required for interpolation
  bool IsForeground(double, double, double, double) const;

  /// Check if the location is fully inside the foreground of the image, i.e.,
  /// including all discrete image locations required for interpolation
  bool IsForeground(const Point &) const;

  /// Check if the location is fully inside the foreground of the image, i.e.,
  /// including all discrete image locations required for interpolation
  bool IsForeground(const Point &, double) const;

  // ---------------------------------------------------------------------------
  // Evaluation

  /// Evaluate scalar image without handling boundary conditions
  ///
  /// This version is faster than EvaluateOutside, but is only defined inside
  /// the domain for which all image values required for interpolation are
  /// defined and thus require no extrapolation of the finite image.
  virtual double EvaluateInside(double, double, double = 0, double = 0) const = 0;

  /// Evaluate scalar image at an arbitrary location (in pixels)
  virtual double EvaluateOutside(double, double, double = 0, double = 0) const = 0;

  /// Evaluate scalar image at an arbitrary location (in pixels)
  ///
  /// If the location is inside the domain for which the filter can perform
  /// an interpolation without considering a particular boundary condition,
  /// the faster EvaluateInside method is called. Otherwise, the EvaluateOutside
  /// method which makes use of the extrapolation of the discrete image domain
  /// in order to interpolate also at boundary or outside locations is used.
  double Evaluate(double, double, double = 0, double = 0) const;

  /// Evaluate scalar image at an arbitrary location (in pixels)
  ///
  /// If the location is inside the domain for which the filter can perform
  /// an interpolation without considering a particular boundary condition,
  /// the faster EvaluateInside method is called. Otherwise, the EvaluateOutside
  /// method which makes use of the extrapolation of the discrete image domain
  /// in order to interpolate also at boundary or outside locations is used.
  double Evaluate(const Point &, double = 0) const;

  /// Evaluate scalar image at an arbitrary location (in pixels)
  ///
  /// If the location is inside the domain for which the filter can perform
  /// an interpolation without considering a particular boundary condition,
  /// the faster EvaluateInside method is called. Otherwise, the EvaluateOutside
  /// method which makes use of the extrapolation of the discrete image domain
  /// in order to interpolate also at boundary or outside locations is used.
  ///
  /// \note This overloaded function corrects for the const-ness of the
  ///       corresponding virtual base class function ImageFunction::Evaluate.
  double Evaluate(double, double, double = 0, double = 0);

  /// Evaluate scalar image at an arbitrary location (in pixels)
  ///
  /// If the location is partially inside the foreground region of the image,
  /// only the foreground values are interpolated. Otherwise, the _DefaultValue
  /// is returned.
  ///
  /// This version is faster than EvaluateWithPaddingOutside, but is only defined
  /// inside the domain for which all image values required for interpolation are
  /// defined and thus require no extrapolation of the finite image.
  virtual double EvaluateWithPaddingInside(double, double, double = 0, double = 0) const = 0;

  /// Evaluate scalar image at an arbitrary location (in pixels)
  ///
  /// If the location is partially inside the foreground region of the image,
  /// only the foreground values are interpolated. Otherwise, the _DefaultValue
  /// is returned.
  virtual double EvaluateWithPaddingOutside(double, double, double = 0, double = 0) const = 0;

  /// Evaluate scalar image at an arbitrary location (in pixels)
  ///
  /// If the location is partially inside the foreground region of the image,
  /// only the foreground values are interpolated. Otherwise, the _DefaultValue
  /// is returned.
  ///
  /// If the location is inside the domain for which the filter can perform
  /// an interpolation without considering a particular boundary condition,
  /// the faster EvaluateWithPaddingInside method is called. Otherwise, the
  /// EvaluateWithPaddingOutside method which makes use of the extrapolation
  /// of the discrete image domain in order to interpolate also at boundary or
  /// outside locations is used.
  double EvaluateWithPadding(double, double, double = 0, double = 0) const;

  /// Evaluate multi-channel image without handling boundary conditions
  ///
  /// This version is faster than EvaluateOutside, but is only defined inside
  /// the domain for which all image values required for interpolation are
  /// defined and thus require no extrapolation of the finite image.
  virtual void EvaluateInside(double *, double, double, double = 0, int = 1) const;

  /// Evaluate multi-channel image at an arbitrary location (in pixels)
  virtual void EvaluateOutside(double *, double, double, double = 0, int = 1) const;

  /// Evaluate multi-channel image at an arbitrary location (in pixels)
  ///
  /// If the location is inside the domain for which the filter can perform
  /// an interpolation without considering a particular boundary condition,
  /// the faster EvaluateInside method is called. Otherwise, the EvaluateOutside
  /// method which makes use of the extrapolation of the discrete image domain
  /// in order to interpolate also at boundary or outside locations is used.
  void Evaluate(double *, double, double, double = 0, int = 1) const;

  /// Evaluate multi-channel image at an arbitrary location (in pixels)
  ///
  /// If the location is inside the domain for which the filter can perform
  /// an interpolation without considering a particular boundary condition,
  /// the faster EvaluateInside method is called. Otherwise, the EvaluateOutside
  /// method which makes use of the extrapolation of the discrete image domain
  /// in order to interpolate also at boundary or outside locations is used.
  void Evaluate(double *, const Point &, int = 1) const;

  /// Evaluate multi-channel image at an arbitrary location (in pixels)
  ///
  /// If the location is partially inside the foreground region of the image,
  /// only the foreground values are interpolated. Otherwise, the _DefaultValue
  /// is returned.
  ///
  /// This version is faster than EvaluateWithPaddingOutside, but is only defined
  /// inside the domain for which all image values required for interpolation are
  /// defined and thus require no extrapolation of the finite image.
  virtual void EvaluateWithPaddingInside(double *, double, double, double = 0, int = 1) const;

  /// Evaluate multi-channel image at an arbitrary location (in pixels)
  ///
  /// If the location is partially inside the foreground region of the image,
  /// only the foreground values are interpolated. Otherwise, the _DefaultValue
  /// is returned.
  virtual void EvaluateWithPaddingOutside(double *, double, double, double = 0, int = 1) const;

  /// Evaluate multi-channel image at an arbitrary location (in pixels)
  ///
  /// If the location is partially inside the foreground region of the image,
  /// only the foreground values are interpolated. Otherwise, the _DefaultValue
  /// is returned.
  ///
  /// If the location is inside the domain for which the filter can perform
  /// an interpolation without considering a particular boundary condition,
  /// the faster EvaluateWithPaddingInside method is called. Otherwise, the
  /// EvaluateWithPaddingOutside method which makes use of the extrapolation
  /// of the discrete image domain in order to interpolate also at boundary or
  /// outside locations is used.
  void EvaluateWithPadding(double *, double, double, double = 0, int = 1) const;

  /// Evaluate vector image without handling boundary conditions
  ///
  /// This version is faster than EvaluateOutside, but is only defined inside
  /// the domain for which all image values required for interpolation are
  /// defined and thus require no extrapolation of the finite image.
  virtual void EvaluateInside(Vector &, double, double, double = 0, double = 0) const = 0;

  /// Evaluate vector image at an arbitrary location (in pixels)
  virtual void EvaluateOutside(Vector &, double, double, double = 0, double = 0) const = 0;

  /// Evaluate vector image at an arbitrary location (in pixels)
  ///
  /// If the location is inside the domain for which the filter can perform
  /// an interpolation without considering a particular boundary condition,
  /// the faster EvaluateInside method is called. Otherwise, the EvaluateOutside
  /// method which makes use of the extrapolation of the discrete image domain
  /// in order to interpolate also at boundary or outside locations is used.
  void Evaluate(Vector &, double, double, double = 0, double = 0) const;

  /// Evaluate vector image at an arbitrary location (in pixels)
  ///
  /// If the location is partially inside the foreground region of the image,
  /// only the foreground values are interpolated. Otherwise, a vector set to
  /// the _DefaultValue is returned.
  ///
  /// This version is faster than EvaluateWithPaddingOutside, but is only defined
  /// inside the domain for which all image values required for interpolation are
  /// defined and thus require no extrapolation of the finite image.
  virtual void EvaluateWithPaddingInside(Vector &, double, double, double = 0, double = 0) const = 0;

  /// Evaluate vector image at an arbitrary location (in pixels)
  ///
  /// If the location is partially inside the foreground region of the image,
  /// only the foreground values are interpolated. Otherwise, a vector set to
  /// the _DefaultValue is returned.
  virtual void EvaluateWithPaddingOutside(Vector &, double, double, double = 0, double = 0) const = 0;

  /// Evaluate vector image at an arbitrary location (in pixels)
  ///
  /// If the location is partially inside the foreground region of the image,
  /// only the foreground values are interpolated. Otherwise, a vector set to
  /// the _DefaultValue is returned.
  ///
  /// If the location is inside the domain for which the filter can perform
  /// an interpolation without considering a particular boundary condition,
  /// the faster EvaluateWithPaddingInside method is called. Otherwise, the
  /// EvaluateWithPaddingOutside method which makes use of the extrapolation
  /// of the discrete image domain in order to interpolate also at boundary or
  /// outside locations is used.
  void EvaluateWithPadding(Vector &, double, double, double = 0, double = 0) const;

  /// Evaluate image function at all locations of the output image
  template <class TVoxel>
  void Evaluate(GenericImage<TVoxel> &) const;

  // ---------------------------------------------------------------------------
  // Derivatives

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
  void EvaluateJacobian(Matrix &, double, double, double = 0, double = NaN) const;


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

  /// Get 1st order derivatives of given image at arbitrary location (in pixels)
  ///
  /// When the image has scalar data type and stores vector components in the
  /// fourth dimension, the derivatives of all components are evaluated when
  /// the t coordinate is set to NaN. Otherwise, only the derivatives of the
  /// specified t component are evaluated.
  void EvaluateJacobianWithPadding(Matrix &, double, double, double = 0, double = NaN) const;

};

////////////////////////////////////////////////////////////////////////////////
// Generic interpolation interface
////////////////////////////////////////////////////////////////////////////////

/**
 * Abstract base class of generic interpolation functions
 *
 * This interpolation interface is templated over the input image type and
 * thus can access the image data using non-virtual getters which return
 * the image values with the appropriate voxel type. Therefore, it is more
 * efficient and type safe to use this interpolation interface whenever the
 * image type is known. Otherwise, use the abstract InterpolateImageFunction
 * interface instead.
 *
 * \sa InterpolateImageFunction
 */
template <class TImage>
class GenericInterpolateImageFunction : public InterpolateImageFunction
{
  mirtkAbstractMacro(GenericInterpolateImageFunction);

public:

  // ---------------------------------------------------------------------------
  // Types

  typedef TImage                                     ImageType;
  typedef typename ImageType::VoxelType              VoxelType;
  typedef typename ImageType::RealType               RealType;
  typedef typename ImageType::RealScalarType         Real;
  typedef GenericExtrapolateImageFunction<ImageType> ExtrapolatorType;

  // ---------------------------------------------------------------------------
  // Construction/Destruction

protected:

  /// Default constructor
  GenericInterpolateImageFunction();

public:

  /// Destructor
  virtual ~GenericInterpolateImageFunction();

  /// Construct interpolator with default infinite extension of input image
  static GenericInterpolateImageFunction *New(enum InterpolationMode = Interpolation_Default,
                                              const TImage * = NULL);

  /// Construct extrapolator which is compatible with this interpolator
  virtual ExtrapolateImageFunction *New(enum ExtrapolationMode,
                                        const BaseImage * = NULL);

  /// Construct interpolator with specified infinite extension of input image
  ///
  /// The caller is required to set the input, initialize, and destroy the
  /// interpolator only, the extrapolator is initialized and destroyed by the
  /// interpolator unless the extrapolator has been replaced using the setter.
  static GenericInterpolateImageFunction *New(enum InterpolationMode,
                                              enum ExtrapolationMode,
                                              const TImage * = NULL);

  // ---------------------------------------------------------------------------
  // Initialization

  /// Set input image
  virtual void Input(const BaseImage *);

  /// Get input image
  const ImageType *Input() const;

  /// Set extrapolate image function for evaluation outside of image domain
  virtual void Extrapolator(ExtrapolateImageFunction *, bool = false);

  /// Get extrapolate image function for evaluation outside of image domain
  /// or \c NULL if extrapolation mode is \c Extrapolation_None
  ExtrapolatorType *Extrapolator();

  /// Get extrapolate image function for evaluation outside of image domain
  /// or \c NULL if extrapolation mode is \c Extrapolation_None
  const ExtrapolatorType *Extrapolator() const;

  /// Initialize image function
  ///
  /// \param[in] coeff Whether input image contains interpolation coefficients
  ///                  already. Otherwise, the interpolate image function will
  ///                  compute these coefficients from the input intensities.
  virtual void Initialize(bool coeff = false);

  // ---------------------------------------------------------------------------
  // Evaluation (generic type)

  /// Evaluate generic image without handling boundary conditions
  ///
  /// This version is faster than GetOutside, but is only defined inside
  /// the domain for which all image values required for interpolation are
  /// defined and thus require no extrapolation of the finite image.
  virtual VoxelType GetInside(double, double, double = 0, double = 0) const = 0;

  /// Evaluate generic image at an arbitrary location (in pixels)
  virtual VoxelType GetOutside(double, double, double = 0, double = 0) const = 0;

  /// Evaluate generic image without handling boundary conditions
  ///
  /// If the location is partially inside the foreground region of the image,
  /// only the foreground values are interpolated. Otherwise, the _DefaultValue
  /// is returned.
  ///
  /// This version is faster than GetWithPaddingOutside, but is only defined
  /// inside the domain for which all image values required for interpolation
  /// are defined and thus require no extrapolation of the finite image.
  virtual VoxelType GetWithPaddingInside(double, double, double = 0, double = 0) const = 0;

  /// Evaluate generic image at an arbitrary location (in pixels)
  ///
  /// If the location is partially inside the foreground region of the image,
  /// only the foreground values are interpolated. Otherwise, the _DefaultValue
  /// is returned.
  virtual VoxelType GetWithPaddingOutside(double, double, double = 0, double = 0) const = 0;

  /// Evaluate generic image at an arbitrary location (in pixels)
  ///
  /// If the location is partially inside the foreground region of the image,
  /// only the foreground values are interpolated. Otherwise, the _DefaultValue
  /// is returned.
  VoxelType GetWithPadding(double, double, double = 0, double = 0) const;

  /// Evaluate generic image at an arbitrary location (in pixels)
  ///
  /// If the location is inside the domain for which the filter can perform
  /// an interpolation without considering a particular boundary condition,
  /// the faster GetInside method is called. Otherwise, the GetOutside method
  /// which makes use of the extrapolation of the discrete image domain in
  /// order to interpolate also at boundary or outside locations is used.
  VoxelType Get(double, double, double = 0, double = 0) const;

  /// Evaluate generic image at an arbitrary location (in pixels)
  VoxelType operator ()(double, double, double = 0, double = 0) const;

  // ---------------------------------------------------------------------------
  // Evaluation (general type)

  /// Evaluate scalar image without handling boundary conditions
  ///
  /// This version is faster than EvaluateOutside, but is only defined inside
  /// the domain for which all image values required for interpolation are
  /// defined and thus require no extrapolation of the finite image.
  virtual double EvaluateInside(double, double, double = 0, double = 0) const;

  /// Evaluate scalar image at an arbitrary location (in pixels)
  virtual double EvaluateOutside(double, double, double = 0, double = 0) const;

  /// Evaluate scalar image at an arbitrary location (in pixels)
  ///
  /// If the location is partially inside the foreground region of the image,
  /// only the foreground values are interpolated. Otherwise, the _DefaultValue
  /// is returned.
  ///
  /// This version is faster than EvaluateWithPaddingOutside, but is only defined
  /// inside the domain for which all image values required for interpolation are
  /// defined and thus require no extrapolation of the finite image.
  virtual double EvaluateWithPaddingInside(double, double, double = 0, double = 0) const;

  /// Evaluate scalar image at an arbitrary location (in pixels)
  ///
  /// If the location is partially inside the foreground region of the image,
  /// only the foreground values are interpolated. Otherwise, the _DefaultValue
  /// is returned.
  virtual double EvaluateWithPaddingOutside(double, double, double = 0, double = 0) const;

  /// Evaluate multi-channel image without handling boundary conditions
  ///
  /// This version is faster than EvaluateOutside, but is only defined inside
  /// the domain for which all image values required for interpolation are
  /// defined and thus require no extrapolation of the finite image.
  virtual void EvaluateInside(double *, double, double, double = 0, int = 1) const;

  /// Evaluate multi-channel image at an arbitrary location (in pixels)
  virtual void EvaluateOutside(double *, double, double, double = 0, int = 1) const;

  /// Evaluate multi-channel image at an arbitrary location (in pixels)
  ///
  /// If the location is partially inside the foreground region of the image,
  /// only the foreground values are interpolated. Otherwise, the _DefaultValue
  /// is returned.
  ///
  /// This version is faster than EvaluateWithPaddingOutside, but is only defined
  /// inside the domain for which all image values required for interpolation are
  /// defined and thus require no extrapolation of the finite image.
  virtual void EvaluateWithPaddingInside(double *, double, double, double = 0, int = 1) const;

  /// Evaluate multi-channel image at an arbitrary location (in pixels)
  ///
  /// If the location is partially inside the foreground region of the image,
  /// only the foreground values are interpolated. Otherwise, the _DefaultValue
  /// is returned.
  virtual void EvaluateWithPaddingOutside(double *, double, double, double = 0, int = 1) const;

  /// Evaluate vector image without handling boundary conditions
  ///
  /// This version is faster than EvaluateOutside, but is only defined inside
  /// the domain for which all image values required for interpolation are
  /// defined and thus require no extrapolation of the finite image.
  virtual void EvaluateInside(Vector &, double, double, double = 0, double = 0) const;

  /// Evaluate vector image at an arbitrary location (in pixels)
  virtual void EvaluateOutside(Vector &, double, double, double = 0, double = 0) const;

  /// Evaluate vector image at an arbitrary location (in pixels)
  ///
  /// If the location is partially inside the foreground region of the image,
  /// only the foreground values are interpolated. Otherwise, a vector set to
  /// the _DefaultValue is returned.
  ///
  /// This version is faster than EvaluateWithPaddingOutside, but is only defined
  /// inside the domain for which all image values required for interpolation are
  /// defined and thus require no extrapolation of the finite image.
  virtual void EvaluateWithPaddingInside(Vector &, double, double, double = 0, double = 0) const;

  /// Evaluate vector image at an arbitrary location (in pixels)
  ///
  /// If the location is partially inside the foreground region of the image,
  /// only the foreground values are interpolated. Otherwise, a vector set to
  /// the _DefaultValue is returned.
  virtual void EvaluateWithPaddingOutside(Vector &, double, double, double = 0, double = 0) const;

};

////////////////////////////////////////////////////////////////////////////////
// Auxiliary macro for subclass implementation
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
#define mirtkInterpolatorMacro(clsname, mode)                                  \
    mirtkObjectMacro(clsname);                                                 \
  public:                                                                      \
    /** Get interpolation mode implemented by this interpolator */             \
    inline virtual enum InterpolationMode InterpolationMode() const            \
    { return mode; }                                                           \
    /** Get interpolation mode implemented by this class */                    \
    inline static  enum InterpolationMode InterpolationType()                  \
    { return mode; }                                                           \
  private:

// -----------------------------------------------------------------------------
#define mirtkGenericInterpolatorTypes(superclsname)                            \
  public:                                                                      \
    typedef superclsname<TImage>                           Superclass;         \
    typedef typename Superclass::ImageType                 ImageType;          \
    typedef typename Superclass::VoxelType                 VoxelType;          \
    typedef typename Superclass::RealType                  RealType;           \
    typedef typename Superclass::Real                      Real;               \
    typedef typename Superclass::ExtrapolatorType          ExtrapolatorType;   \
  private:

// -----------------------------------------------------------------------------
#define mirtkGenericInterpolatorMacro(clsname, mode)                           \
    mirtkInterpolatorMacro(clsname, mode);                                     \
    mirtkGenericInterpolatorTypes(GenericInterpolateImageFunction);            \
  private:

////////////////////////////////////////////////////////////////////////////////
// Inline definitions -- InterpolateImageFunction
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Initialization
// =============================================================================

// -----------------------------------------------------------------------------
inline void InterpolateImageFunction::Input(const BaseImage *input)
{
  this->_Input = const_cast<BaseImage *>(input);
}

// -----------------------------------------------------------------------------
inline const BaseImage *InterpolateImageFunction::Input() const
{
  return this->_Input;
}

// -----------------------------------------------------------------------------
inline void InterpolateImageFunction::Initialize()
{
  this->Initialize(false);
}

// -----------------------------------------------------------------------------
inline void InterpolateImageFunction::Update()
{
}

// -----------------------------------------------------------------------------
inline void InterpolateImageFunction
::Extrapolator(ExtrapolateImageFunction *input, bool owner)
{
  if (_InfiniteInputOwner) delete _InfiniteInput;
  _InfiniteInput      = input;
  _InfiniteInputOwner = _InfiniteInput && owner;
}

// -----------------------------------------------------------------------------
inline ExtrapolateImageFunction *InterpolateImageFunction::Extrapolator()
{
  return _InfiniteInput;
}

// -----------------------------------------------------------------------------
inline const ExtrapolateImageFunction *InterpolateImageFunction::Extrapolator() const
{
  return _InfiniteInput;
}

// -----------------------------------------------------------------------------
inline enum ExtrapolationMode InterpolateImageFunction::ExtrapolationMode() const
{
  return (_InfiniteInput ? _InfiniteInput->ExtrapolationMode() : Extrapolation_None);
}

// =============================================================================
// Image attributes
// =============================================================================

// -----------------------------------------------------------------------------
inline const ImageAttributes &InterpolateImageFunction::Attributes() const
{
  return this->_Input->Attributes();
}

// -----------------------------------------------------------------------------
inline int InterpolateImageFunction::X() const
{
  return this->_Input->X();
}

// -----------------------------------------------------------------------------
inline int InterpolateImageFunction::Y() const
{
  return this->_Input->Y();
}

// -----------------------------------------------------------------------------
inline int InterpolateImageFunction::Z() const
{
  return this->_Input->Z();
}

// -----------------------------------------------------------------------------
inline int InterpolateImageFunction::T() const
{
  return this->_Input->T();
}

// -----------------------------------------------------------------------------
inline double InterpolateImageFunction::XSize() const
{
  return this->_Input->XSize();
}

// -----------------------------------------------------------------------------
inline double InterpolateImageFunction::YSize() const
{
  return this->_Input->YSize();
}

// -----------------------------------------------------------------------------
inline double InterpolateImageFunction::ZSize() const
{
  return this->_Input->ZSize();
}

// -----------------------------------------------------------------------------
inline double InterpolateImageFunction::TSize() const
{
  return this->_Input->TSize();
}

// ----------------------------------------------------------------------------
inline void InterpolateImageFunction::WorldToImage(double &x, double &y) const
{
  this->_Input->WorldToImage(x, y);
}

// ----------------------------------------------------------------------------
inline void InterpolateImageFunction::WorldToImage(double &x, double &y, double &z) const
{
  this->_Input->WorldToImage(x, y, z);
}

// ----------------------------------------------------------------------------
inline void InterpolateImageFunction::WorldToImage(Point &p) const
{
  this->_Input->WorldToImage(p);
}

// ----------------------------------------------------------------------------
inline void InterpolateImageFunction::WorldToImage(Vector3 &v) const
{
  this->_Input->WorldToImage(v);
}

// ----------------------------------------------------------------------------
inline void InterpolateImageFunction::ImageToWorld(double &x, double &y) const
{
  this->_Input->ImageToWorld(x, y);
}

// ----------------------------------------------------------------------------
inline void InterpolateImageFunction::ImageToWorld(double &x, double &y, double &z) const
{
  this->_Input->ImageToWorld(x, y, z);
}

// ----------------------------------------------------------------------------
inline void InterpolateImageFunction::ImageToWorld(Point &p) const
{
  this->_Input->ImageToWorld(p);
}

// ----------------------------------------------------------------------------
inline void InterpolateImageFunction::ImageToWorld(Vector3 &v) const
{
  this->_Input->ImageToWorld(v);
}

// =============================================================================
// Domain checks
// =============================================================================

// -----------------------------------------------------------------------------
inline void InterpolateImageFunction::Inside(double &x1, double &y1,
                                             double &x2, double &y2) const
{
  x1 = _x1, y1 = _y1;
  x2 = _x2, y2 = _y2;
}

// -----------------------------------------------------------------------------
inline void InterpolateImageFunction::Inside(double &x1, double &y1, double &z1,
                                             double &x2, double &y2, double &z2) const
{
  x1 = _x1, y1 = _y1, z1 = _z1;
  x2 = _x2, y2 = _y2, z2 = _z2;
}

// -----------------------------------------------------------------------------
inline void InterpolateImageFunction::Inside(double &x1, double &y1, double &z1, double &t1,
                                             double &x2, double &y2, double &z2, double &t2) const
{
  x1 = _x1, y1 = _y1, z1 = _z1, t1 = _t1;
  x2 = _x2, y2 = _y2, z2 = _z2, t2 = _t2;
}

// -----------------------------------------------------------------------------
inline void InterpolateImageFunction::BoundingBox(double x, double y,
                                                  int &i1, int &j1,
                                                  int &i2, int &j2) const
{
  this->BoundingInterval(x, i1, i2);
  this->BoundingInterval(y, j1, j2);
}

// -----------------------------------------------------------------------------
inline void InterpolateImageFunction::BoundingBox(double x, double y, double z,
                                                  int &i1, int &j1, int &k1,
                                                  int &i2, int &j2, int &k2) const
{
  if (this->NumberOfDimensions() >= 3) {
    this->BoundingInterval(x, i1, i2);
    this->BoundingInterval(y, j1, j2);
    this->BoundingInterval(z, k1, k2);
  } else {
    this->BoundingInterval(x, i1, i2);
    this->BoundingInterval(y, j1, j2);
    k1 = k2 = iround(z);
  }
}

// -----------------------------------------------------------------------------
inline void InterpolateImageFunction::BoundingBox(double x, double y, double z, double t,
                                                  int &i1, int &j1, int &k1, int &l1,
                                                  int &i2, int &j2, int &k2, int &l2) const
{
  if (this->NumberOfDimensions() >= 4) {
    this->BoundingInterval(x, i1, i2);
    this->BoundingInterval(y, j1, j2);
    this->BoundingInterval(z, k1, k2);
    this->BoundingInterval(t, l1, l2);
  } else if (this->NumberOfDimensions() == 3) {
    this->BoundingInterval(x, i1, i2);
    this->BoundingInterval(y, j1, j2);
    this->BoundingInterval(z, k1, k2);
    l1 = l2 = iround(t);
  } else {
    this->BoundingInterval(x, i1, i2);
    this->BoundingInterval(y, j1, j2);
    k1 = k2 = iround(z);
    l1 = l2 = iround(t);
  }
}

// -----------------------------------------------------------------------------
inline bool InterpolateImageFunction::IsInside(double x, double y) const
{
  return (_x1 <= x && x <= _x2) && (_y1 <= y && y <= _y2);
}

// -----------------------------------------------------------------------------
inline bool InterpolateImageFunction::IsInside(double x, double y, double z) const
{
  return (_x1 <= x && x <= _x2) && (_y1 <= y && y <= _y2) && (_z1 <= z && z <= _z2);
}

// -----------------------------------------------------------------------------
inline bool InterpolateImageFunction::IsInside(double x, double y, double z, double t) const
{
  return (_x1 <= x && x <= _x2) && (_y1 <= y && y <= _y2) && (_z1 <= z && z <= _z2) && (_t1 <= t && t <= _t2);
}

// -----------------------------------------------------------------------------
inline bool InterpolateImageFunction::IsInside(const Point &p) const
{
  return IsInside(p._x, p._y, p._z);
}

// -----------------------------------------------------------------------------
inline bool InterpolateImageFunction::IsInside(const Point &p, double t) const
{
  return IsInside(p._x, p._y, p._z, t);
}

// -----------------------------------------------------------------------------
inline bool InterpolateImageFunction::IsOutside(double x, double y) const
{
  return !IsInside(x, y);
}

// -----------------------------------------------------------------------------
inline bool InterpolateImageFunction::IsOutside(double x, double y, double z) const
{
  return !IsInside(x, y, z);
}

// -----------------------------------------------------------------------------
inline bool InterpolateImageFunction::IsOutside(double x, double y, double z, double t) const
{
  return !IsInside(x, y, z, t);
}

// -----------------------------------------------------------------------------
inline bool InterpolateImageFunction::IsOutside(const Point &p) const
{
  return IsOutside(p._x, p._y, p._z);
}

// -----------------------------------------------------------------------------
inline bool InterpolateImageFunction::IsOutside(const Point &p, double t) const
{
  return IsOutside(p._x, p._y, p._z, t);
}

// -----------------------------------------------------------------------------
inline bool InterpolateImageFunction::IsForeground(double x, double y) const
{
  int i1, j1, i2, j2;
  BoundingBox(x, y,                             i1, j1, i2, j2);
  return Input()->IsBoundingBoxInsideForeground(i1, j1, i2, j2);
}

// -----------------------------------------------------------------------------
inline bool InterpolateImageFunction::IsForeground(double x, double y, double z) const
{
  int i1, j1, k1, i2, j2, k2;
  BoundingBox(x, y, z,                          i1, j1, k1, i2, j2, k2);
  return Input()->IsBoundingBoxInsideForeground(i1, j1, k1, i2, j2, k2);
}

// -----------------------------------------------------------------------------
inline bool InterpolateImageFunction::IsForeground(double x, double y, double z, double t) const
{
  int i1, j1, k1, l1, i2, j2, k2, l2;
  BoundingBox(x, y, z, t,                       i1, j1, k1, l1, i2, j2, k2, l2);
  return Input()->IsBoundingBoxInsideForeground(i1, j1, k1, l1, i2, j2, k2, l2);
}

// -----------------------------------------------------------------------------
inline bool InterpolateImageFunction::IsForeground(const Point &p) const
{
  return IsForeground(p._x, p._y, p._z);
}

// -----------------------------------------------------------------------------
inline bool InterpolateImageFunction::IsForeground(const Point &p, double t) const
{
  return IsForeground(p._x, p._y, p._z, t);
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
inline double InterpolateImageFunction::Evaluate(double x, double y, double z, double t) const
{
  if (IsInside(x, y, z, t)) return this->EvaluateInside (x, y, z, t);
  else                      return this->EvaluateOutside(x, y, z, t);
}

// -----------------------------------------------------------------------------
inline double InterpolateImageFunction::Evaluate(double x, double y, double z, double t)
{
  return const_cast<const InterpolateImageFunction *>(this)->Evaluate(x, y, z, t);
}

// -----------------------------------------------------------------------------
inline double InterpolateImageFunction::Evaluate(const Point &p, double t) const
{
  return Evaluate(p._x, p._y, p._z, t);
}

// -----------------------------------------------------------------------------
inline double InterpolateImageFunction::EvaluateWithPadding(double x, double y, double z, double t) const
{
  if (IsInside(x, y, z, t)) return this->EvaluateWithPaddingInside (x, y, z, t);
  else                      return this->EvaluateWithPaddingOutside(x, y, z, t);
}

// -----------------------------------------------------------------------------
inline void InterpolateImageFunction::EvaluateInside(double *v, double x, double y, double z, int vt) const
{
  for (int t = 0; t < Input()->T(); ++t, v += vt) {
    (*v) = this->EvaluateInside(x, y, z, static_cast<double>(t));
  }
}

// -----------------------------------------------------------------------------
inline void InterpolateImageFunction::EvaluateOutside(double *v, double x, double y, double z, int vt) const
{
  for (int t = 0; t < Input()->T(); ++t, v += vt) {
    (*v) = this->EvaluateOutside(x, y, z, static_cast<double>(t));
  }
}

// -----------------------------------------------------------------------------
inline void InterpolateImageFunction::Evaluate(double *v, double x, double y, double z, int vt) const
{
  if (IsInside(x, y, z)) this->EvaluateInside (v, x, y, z, vt);
  else                   this->EvaluateOutside(v, x, y, z, vt);
}

// -----------------------------------------------------------------------------
inline void InterpolateImageFunction::Evaluate(double *v, const Point &p, int vt) const
{
  Evaluate(v, p._x, p._y, p._z, vt);
}

// -----------------------------------------------------------------------------
inline void InterpolateImageFunction::EvaluateWithPaddingInside(double *v, double x, double y, double z, int vt) const
{
  for (int t = 0; t < Input()->T(); ++t, v += vt) {
    (*v) = this->EvaluateWithPaddingInside(x, y, z, static_cast<double>(t));
  }
}

// -----------------------------------------------------------------------------
inline void InterpolateImageFunction::EvaluateWithPaddingOutside(double *v, double x, double y, double z, int vt) const
{
  for (int t = 0; t < Input()->T(); ++t, v += vt) {
    (*v) = this->EvaluateWithPaddingOutside(x, y, z, static_cast<double>(t));
  }
}

// -----------------------------------------------------------------------------
inline void InterpolateImageFunction::EvaluateWithPadding(double *v, double x, double y, double z, int vt) const
{
  if (IsInside(x, y, z)) this->EvaluateWithPaddingInside (v, x, y, z, vt);
  else                   this->EvaluateWithPaddingOutside(v, x, y, z, vt);
}

// -----------------------------------------------------------------------------
inline void InterpolateImageFunction::Evaluate(Vector &v, double x, double y, double z, double t) const
{
  if (IsInside(x, y, z, t)) this->EvaluateInside (v, x, y, z, t);
  else                      this->EvaluateOutside(v, x, y, z, t);
}

// -----------------------------------------------------------------------------
inline void InterpolateImageFunction::EvaluateWithPadding(Vector &v, double x, double y, double z, double t) const
{
  if (IsInside(x, y, z, t)) this->EvaluateWithPaddingInside (v, x, y, z, t);
  else                      this->EvaluateWithPaddingOutside(v, x, y, z, t);
}

// -----------------------------------------------------------------------------
template <class TVoxel>
inline void InterpolateImageFunction::Evaluate(GenericImage<TVoxel> &output) const
{
  // Interpolate multi-channel image (or 3D+t vector image)
  if (output.T() > 1 && IsZero(output.TSize())) {
    UnaryVoxelFunction::InterpolateMultiChannelImage<InterpolateImageFunction> eval(this, &output);
    ParallelForEachVoxel(output.Attributes(), output, eval);
  // Interpolate scalar image
  } else if (output.N() == 1) {
    UnaryVoxelFunction::InterpolateScalarImage<InterpolateImageFunction> eval(this, &output);
    ParallelForEachVoxel(output.Attributes(), output, eval);
  // Interpolate vector image
  } else {
    UnaryVoxelFunction::InterpolateImage<InterpolateImageFunction> eval(this, &output);
    ParallelForEachVoxel(output.Attributes(), output, eval);
  }
}

// =============================================================================
// Derivatives
// =============================================================================

// -----------------------------------------------------------------------------
inline void InterpolateImageFunction
::EvaluateJacobianInside(Matrix &, double, double, double, double) const
{
  Throw(ERR_NotImplemented, __FUNCTION__, "Not implemented");
}

// -----------------------------------------------------------------------------
inline void InterpolateImageFunction
::EvaluateJacobianOutside(Matrix &, double, double, double, double) const
{
  Throw(ERR_NotImplemented, __FUNCTION__, "Not implemented");
}

// -----------------------------------------------------------------------------
inline void InterpolateImageFunction
::EvaluateJacobian(Matrix &jac, double x, double y, double z, double t) const
{
  if (this->IsInside(x, y, z, t)) this->EvaluateJacobianInside (jac, x, y, z, t);
  else                            this->EvaluateJacobianOutside(jac, x, y, z, t);
}

// -----------------------------------------------------------------------------
inline void InterpolateImageFunction
::EvaluateJacobianWithPaddingInside(Matrix &, double, double, double, double) const
{
  Throw(ERR_NotImplemented, __FUNCTION__, "Not implemented");
}

// -----------------------------------------------------------------------------
inline void InterpolateImageFunction
::EvaluateJacobianWithPaddingOutside(Matrix &, double, double, double, double) const
{
  Throw(ERR_NotImplemented, __FUNCTION__, "Not implemented");
}

// -----------------------------------------------------------------------------
inline void InterpolateImageFunction
::EvaluateJacobianWithPadding(Matrix &jac, double x, double y, double z, double t) const
{
  if (this->IsInside(x, y, z, t)) this->EvaluateJacobianWithPaddingInside (jac, x, y, z, t);
  else                            this->EvaluateJacobianWithPaddingOutside(jac, x, y, z, t);
}

////////////////////////////////////////////////////////////////////////////////
// Inline definitions -- GenericInterpolateImageFunction
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Getters
// =============================================================================

// -----------------------------------------------------------------------------
template <class TImage>
inline const TImage *GenericInterpolateImageFunction<TImage>::Input() const
{
  return reinterpret_cast<const TImage *>(this->_Input);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericInterpolateImageFunction<TImage>::ExtrapolatorType *
GenericInterpolateImageFunction<TImage>::Extrapolator()
{
  return reinterpret_cast<ExtrapolatorType *>(_InfiniteInput);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline const typename GenericInterpolateImageFunction<TImage>::ExtrapolatorType *
GenericInterpolateImageFunction<TImage>::Extrapolator() const
{
  return reinterpret_cast<const ExtrapolatorType *>(_InfiniteInput);
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
template <class TImage>
inline typename TImage::VoxelType
GenericInterpolateImageFunction<TImage>
::Get(double x, double y, double z, double t) const
{
  if (IsInside(x, y, z, t)) return this->GetInside (x, y, z, t);
  else                      return this->GetOutside(x, y, z, t);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename TImage::VoxelType
GenericInterpolateImageFunction<TImage>
::GetWithPadding(double x, double y, double z, double t) const
{
  if (IsInside(x, y, z, t)) return this->GetWithPaddingInside (x, y, z, t);
  else                      return this->GetWithPaddingOutside(x, y, z, t);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename TImage::VoxelType GenericInterpolateImageFunction<TImage>
::operator ()(double x, double y, double z, double t) const
{
  return this->Get(x, y, z, t);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline double GenericInterpolateImageFunction<TImage>
::EvaluateInside(double x, double y, double z, double t) const
{
  return voxel_cast<double>(this->GetInside(x, y, z, t));
}

// -----------------------------------------------------------------------------
template <class TImage>
inline double GenericInterpolateImageFunction<TImage>
::EvaluateOutside(double x, double y, double z, double t) const
{
  return voxel_cast<double>(this->GetOutside(x, y, z, t));
}

// -----------------------------------------------------------------------------
template <class TImage>
inline double GenericInterpolateImageFunction<TImage>
::EvaluateWithPaddingInside(double x, double y, double z, double t) const
{
  return voxel_cast<double>(this->GetWithPaddingInside(x, y, z, t));
}

// -----------------------------------------------------------------------------
template <class TImage>
inline double GenericInterpolateImageFunction<TImage>
::EvaluateWithPaddingOutside(double x, double y, double z, double t) const
{
  return voxel_cast<double>(this->GetWithPaddingOutside(x, y, z, t));
}

// -----------------------------------------------------------------------------
template <class TImage>
inline void GenericInterpolateImageFunction<TImage>
::EvaluateInside(double *v, double x, double y, double z, int vt) const
{
  for (int t = 0; t < Input()->T(); ++t, v += vt) {
    (*v) = voxel_cast<double>(this->GetInside(x, y, z, static_cast<double>(t)));
  }
}

// -----------------------------------------------------------------------------
template <class TImage>
inline void GenericInterpolateImageFunction<TImage>
::EvaluateOutside(double *v, double x, double y, double z, int vt) const
{
  for (int t = 0; t < Input()->T(); ++t, v += vt) {
    (*v) = voxel_cast<double>(this->GetOutside(x, y, z, static_cast<double>(t)));
  }
}

// -----------------------------------------------------------------------------
template <class TImage>
inline void GenericInterpolateImageFunction<TImage>
::EvaluateWithPaddingInside(double *v, double x, double y, double z, int vt) const
{
  for (int t = 0; t < Input()->T(); ++t, v += vt) {
    (*v) = voxel_cast<double>(this->GetWithPaddingInside(x, y, z, static_cast<double>(t)));
  }
}

// -----------------------------------------------------------------------------
template <class TImage>
inline void GenericInterpolateImageFunction<TImage>
::EvaluateWithPaddingOutside(double *v, double x, double y, double z, int vt) const
{
  for (int t = 0; t < Input()->T(); ++t, v += vt) {
    (*v) = voxel_cast<double>(this->GetWithPaddingOutside(x, y, z, static_cast<double>(t)));
  }
}

// -----------------------------------------------------------------------------
template <class TImage>
inline void GenericInterpolateImageFunction<TImage>
::EvaluateInside(Vector &v, double x, double y, double z, double t) const
{
  v = voxel_cast<Vector>(this->GetInside(x, y, z, t));
}

// -----------------------------------------------------------------------------
template <class TImage>
inline void GenericInterpolateImageFunction<TImage>
::EvaluateOutside(Vector &v, double x, double y, double z, double t) const
{
  v = voxel_cast<Vector>(this->GetOutside(x, y, z, t));
}

// -----------------------------------------------------------------------------
template <class TImage>
inline void GenericInterpolateImageFunction<TImage>
::EvaluateWithPaddingInside(Vector &v, double x, double y, double z, double t) const
{
  v = voxel_cast<Vector>(this->GetWithPaddingInside(x, y, z, t));
}

// -----------------------------------------------------------------------------
template <class TImage>
inline void GenericInterpolateImageFunction<TImage>
::EvaluateWithPaddingOutside(Vector &v, double x, double y, double z, double t) const
{
  v = voxel_cast<Vector>(this->GetWithPaddingOutside(x, y, z, t));
}


} // namespace mirtk

#endif // MIRTK_InterpolateImageFunction_H
