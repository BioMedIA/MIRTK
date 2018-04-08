/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2017 Imperial College London
 * Copyright 2013-2017 Andreas Schuh
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

#ifndef MIRTK_RegisteredImage_H
#define MIRTK_RegisteredImage_H

#include "mirtk/GenericImage.h"

#include "mirtk/Parallel.h"
#include "mirtk/Transformation.h"
#include "mirtk/InterpolationMode.h"
#include "mirtk/ExtrapolationMode.h"

#include "mirtk/TestProd.h"


namespace mirtk {


/**
 * Registered image such as fixed target image or transformed source image
 *
 * An instance of this (multi-channel) image class stores the intensities
 * of a registered image after alignment using the current transformation
 * estimate. If no (changing) transformation is set, the image is considered
 * to be fixed and is thus only updated once during the initialization.
 * If a transformation is set, however, upon every update of the input
 * transformation, the new moving image intensities are computed by applying
 * the current transformation and interpolating the input image intensities.
 * For the computation of the similarity measure gradient, the warped image
 * derivatives are stored along with the intensities in the consecutive image
 * channels of the registered image instance (i.e., t dimension). In case of
 * a gradient field similarity measure, these derivatives include not only the
 * 1st order, but also 2nd order derivatives needed for the similarity measure
 * gradient computation. An image similarity specifies which channels are
 * required for its evaluation by calling the Initialize function with a suitable
 * number of channels (i.e., 1, 4, or 10) and request their udpate using the
 * appropriate parameters of the Update function.
 *
 * The assignment of the transformed intensities and derivatives to the fourth
 * dimension (i.e. channels) of the registered image is as follows:
 *
 * - t=0: Transformed intensity
 * - t=1: Transformed 1st order derivative w.r.t x
 * - t=2: Transformed 1st order derivative w.r.t y
 * - t=3: Transformed 1st order derivative w.r.t z
 * - t=4: Transformed 2nd order derivative w.r.t xx
 * - t=5: Transformed 2nd order derivative w.r.t xy
 * - t=6: Transformed 2nd order derivative w.r.t xz
 * - t=7: Transformed 2nd order derivative w.r.t yy
 * - t=8: Transformed 2nd order derivative w.r.t yz
 * - t=9: Transformed 2nd order derivative w.r.t zz
 */
class RegisteredImage : public GenericImage<double>
{
  mirtkObjectMacro(RegisteredImage);

public:

  // Do not override other base class overloads
  using GenericImage<VoxelType>::ImageToWorld;

  // ---------------------------------------------------------------------------
  // Types
public:

  /// Enumeration of registered image channel indices
  enum Channel { I = 0, Dx, Dy, Dz, Dxx, Dxy, Dxz, Dyx, Dyy, Dyz, Dzx, Dzy, Dzz };

  /// Type of untransformed input image
  typedef GenericImage<VoxelType>   InputImageType;

  /// Type of untransformed gradient image
  typedef GenericImage<VoxelType>   InputGradientType;

  /// Type of untransformed Hessian image
  typedef GenericImage<VoxelType>   InputHessianType;

  /// Type of cached displacement fields
  typedef GenericImage<double>   DisplacementImageType;

  // ---------------------------------------------------------------------------
  // Attributes

  /// Untransformed input image
  mirtkPublicAggregateMacro(InputImageType, InputImage);

  /// Untransformed input gradient image
  mirtkPublicComponentMacro(InputGradientType, InputGradient);

  /// Untransformed input Hessian image
  mirtkPublicComponentMacro(InputHessianType, InputHessian);

  /// Current transformation estimate
  mirtkPublicAggregateMacro(const class Transformation, Transformation);

  /// Interpolation mode
  mirtkPublicAttributeMacro(enum InterpolationMode, InterpolationMode);

  /// Extrapolation mode
  mirtkPublicAttributeMacro(enum ExtrapolationMode, ExtrapolationMode);

  /// Pre-computed image to world coordinates
  mirtkPublicAggregateMacro(WorldCoordsImage, WorldCoordinates);

  /// Pre-computed image to world coordinates
  mirtkPublicComponentMacro(WorldCoordsImage, ImageToWorld);

  /// Externally pre-computed displacements to use
  mirtkPublicAggregateMacro(DisplacementImageType, ExternalDisplacement);

  /// Pre-computed fixed displacements
  mirtkComponentMacro(DisplacementImageType, FixedDisplacement);

  /// Pre-computed displacements
  mirtkComponentMacro(DisplacementImageType, Displacement);

  /// Whether to use pre-computed world coordinates
  mirtkPublicAttributeMacro(bool, CacheWorldCoordinates);

  /// Whether to use pre-computed fixed displacements
  mirtkPublicAttributeMacro(bool, CacheFixedDisplacement);

  /// Whether to use pre-computed displacements
  mirtkPublicAttributeMacro(bool, CacheDisplacement);

  /// Whether self-update is enabled
  mirtkPublicAttributeMacro(bool, SelfUpdate);

  /// Minimum foreground intensity of input image
  mirtkReadOnlyAttributeMacro(double, MinInputIntensity);

  /// Maximum foreground intensity of input image
  mirtkReadOnlyAttributeMacro(double, MaxInputIntensity);

  /// Minimum foreground intensity of warped image or NaN
  mirtkPublicAttributeMacro(double, MinIntensity);

  /// Maximum foreground intensity of warped image or NaN
  mirtkPublicAttributeMacro(double, MaxIntensity);

  /// Standard deviation of Gaussian smoothing kernel applied before
  /// 1st order derivative computations in voxel units
  mirtkPublicAttributeMacro(double, GradientSigma);

  /// Standard deviation of Gaussian smoothing kernel applied before
  /// 2nd order derivative computations in voxel units
  mirtkPublicAttributeMacro(double, HessianSigma);

  /// Whether to precompute image derivatives
  /// true:  Compute derivatives of input image and transform these
  /// false: Use derivative of interpolation kernel to evaluate image derivative
  mirtkPublicAttributeMacro(bool, PrecomputeDerivatives);

  /// Maximum gradient magnitude as percentile of all gradient magnitudes
  ///
  /// \note Used only when _PrecomputeDerivatives is enabled.
  mirtkPublicAttributeMacro(int, MaxGradientPercentile);

  /// When positive and not inf, linearly rescale gradient vectors to this maximum magnitude
  mirtkPublicAttributeMacro(double, MaxGradientMagnitude);

protected:

  /// Number of active levels
  int _NumberOfActiveLevels;

  /// Number of passive levels (fixed transformation)
  int _NumberOfPassiveLevels;

  /// Offsets of the different registered image channels
  int _Offset[13];

  /// (Pre-)compute gradient of input image
  /// \param[in] sigma Standard deviation of Gaussian smoothing filter in voxels.
  void ComputeInputGradient(double sigma);

  /// (Pre-)compute Hessian of input image
  /// \param[in] sigma Standard deviation of Gaussian smoothing filter in voxels.
  void ComputeInputHessian(double sigma);

private:

  /// Internal voxel-wise update functor object
  SharedPtr<Object> _Evaluator;

  // ---------------------------------------------------------------------------
  // Construction/Destruction
public:

  /// Constructor
  RegisteredImage();

  /// Copy constructor
  RegisteredImage(const RegisteredImage &);

  /// Assignment operator
  RegisteredImage &operator =(const RegisteredImage &);

  /// Destructor
  ~RegisteredImage();

  /// Get interpolation mode
  ///
  /// When _InterpolationMode is Interpolation_Default, this function returns
  /// the respective default interpolation mode with or without padding.
  enum InterpolationMode GetInterpolationMode() const;

  // ---------------------------------------------------------------------------
  // Channels

  /// Number of registered channels
  int NumberOfChannels() const;

  /// Number of voxels per registered channel
  int NumberOfVoxels() const;

  /// Offset of channel \p c with respect to the start of the image data
  int Offset(int) const;

  // ---------------------------------------------------------------------------
  // Initialization/Update

  /// Initialize image
  ///
  /// This function is called by the image similarity to set the attributes
  /// of the registered image, i.e., the attributes of the discrete image domain
  /// on which similarity is to be evaluated. It initializes the image to zero.
  /// The Update function must be called to initialize the output image
  /// intensities and derivatives before similarity can be evaluated.
  void Initialize(const ImageAttributes &, int = 0);

  /// Update image intensity and transformed derivatives
  ///
  /// This function updates the requested output image channels only within the
  /// specified image region. If the image is not transformed or only warped by
  /// a fixed transformation, this function does nothing. Set \c force to \c true
  /// to initialize this image even if it does not change over the course of the
  /// registration or use the corresponding Recompute function.
  ///
  /// \param[in] region    Image region to update.
  /// \param[in] intensity Request update of scalar intensities.
  /// \param[in] gradient  Request update of 1st order derivatives.
  /// \param[in] hessian   Request update of 2nd order derivatives.
  /// \param[in] force     Force update in any case.
  void Update(const blocked_range3d<int> &region,
              bool intensity = true, bool gradient = false, bool hessian = false,
              bool force     = false);

  /// Update image intensity and transformed derivatives using custom displacements
  ///
  /// \param[in] region    Image region to update.
  /// \param[in] disp      Custom displacement field to use.
  /// \param[in] intensity Request update of scalar intensities.
  /// \param[in] gradient  Request update of 1st order derivatives.
  /// \param[in] hessian   Request update of 2nd order derivatives.
  void Update(const blocked_range3d<int> &region,
              const DisplacementImageType *disp,
              bool intensity = true, bool gradient = false, bool hessian = false);

  /// Update image intensity and transformed derivatives
  ///
  /// This function only updates the specified image channels if the self-update
  /// of this image is enabled and only if a (changing) transformation is set.
  /// If the image is not transformed or only warped by a fixed transformation,
  /// this function does nothing. Use \c force=true to initialize this image
  /// even if does not change over the course of the registration.
  ///
  /// \param[in] intensity Request update of scalar intensities.
  /// \param[in] gradient  Request update of 1st order derivatives.
  /// \param[in] hessian   Request update of 2nd order derivatives.
  /// \param[in] force     Force update in any case.
  void Update(bool intensity = true, bool gradient = false, bool hessian = false,
              bool force     = false);

  /// Force update of all output channels
  ///
  /// This convenience function always recomputes all output channels within
  /// the specified image region. To do so, it calls the Update function with
  /// \c force set to \c true. It is commonly used after Initialize in an
  /// application where the input image and transformation are not being modified
  /// or when the SelfUpdate of the image has been disabled to explicitly
  /// recompute the image content after the input was modified.
  void Recompute(const blocked_range3d<int> &region);

  /// Force update of all output channels
  ///
  /// This convenience function always recomputes all output channels.
  /// To do so, it calls the Update function with \c force set to \c true.
  /// It is commonly used after Initialize in an application where the
  /// input image and transformation are not being modified or when the
  /// SelfUpdate of the image has been disabled to explicitly recompute
  /// the image content after the input was modified.
  void Recompute();

private:

  template <class>                      void Update1(const blocked_range3d<int> &, bool, bool, bool);
  template <class, class, class, class> void Update2(const blocked_range3d<int> &, bool, bool, bool);
  template <class, class>               void Update3(const blocked_range3d<int> &, bool, bool, bool);

  FRIEND_TEST(RegisteredImage, GlobalAndLocalTransformation);
};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
inline int RegisteredImage::NumberOfChannels() const
{
  return _attr._t;
}

// -----------------------------------------------------------------------------
inline int RegisteredImage::NumberOfVoxels() const
{
  return _attr._x * _attr._y * _attr._z;
}

// -----------------------------------------------------------------------------
inline int RegisteredImage::Offset(int c) const
{
  return _Offset[c];
}


} // namespace mirtk

#endif // MIRTK_RegisteredImage_H
