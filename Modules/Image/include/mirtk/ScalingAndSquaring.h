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

#ifndef MIRTK_ScalingAndSquaring_H
#define MIRTK_ScalingAndSquaring_H

#include "mirtk/Object.h"
#include "mirtk/GenericImage.h"
#include "mirtk/InterpolateImageFunction.h"
#include "mirtk/FastCubicBSplineInterpolateImageFunction.h"


namespace mirtk {


/**
 * Computes the exponential map of a SVF and its derivatives
 *
 * This class implements an image filter which computes the exponential map
 * of a stationary velocity field using the scaling and squaring method.
 * The result is a diffeomorphic displacement field. Additionally, this filter
 * computes also the derivatives of the exponential map either with respect to
 * the spatial coordinate or the stationary velocity field itself.
 */
template <class TReal>
class ScalingAndSquaring : public Object
{
  mirtkObjectMacro(ScalingAndSquaring);

  // ---------------------------------------------------------------------------
  // Types

public:

  /// Image type of input, output, and intermediate images
  typedef GenericImage<TReal> ImageType;

  /// Type of vector field interpolator
  typedef GenericInterpolateImageFunction<ImageType> VectorField;

  /// Type of input velocity field interpolator
  typedef GenericFastCubicBSplineInterpolateImageFunction<ImageType>  VelocityField;

private:

  // ---------------------------------------------------------------------------
  // Input

  /// Input velocity field
  mirtkPublicAggregateMacro(const ImageType, InputVelocity);

  /// Input displacement field
  mirtkPublicAggregateMacro(const ImageType, InputDisplacement);

  /// Input deformation field
  mirtkPublicAggregateMacro(const ImageType, InputDeformation);

  // ---------------------------------------------------------------------------
  // Interim

  /// Intermediate displacement field
  mirtkReadOnlyAttributeMacro(UniquePtr<ImageType>, InterimDisplacement);

  /// Intermediate Jacobian w.r.t. x
  mirtkReadOnlyAttributeMacro(UniquePtr<ImageType>, InterimJacobian);

  /// Intermediate determinant of Jacobian w.r.t. x
  mirtkReadOnlyAttributeMacro(UniquePtr<ImageType>, InterimDetJacobian);

  /// Intermediate log of determinant of Jacobian w.r.t. x
  mirtkReadOnlyAttributeMacro(UniquePtr<ImageType>, InterimLogJacobian);

  /// Intermediate Jacobian w.r.t. v
  mirtkReadOnlyAttributeMacro(UniquePtr<ImageType>, InterimJacobianDOFs);

  // ---------------------------------------------------------------------------
  // Output

  /// Output displacement field
  mirtkPublicAggregateMacro(ImageType, OutputDisplacement);

  /// Output deformation field
  mirtkPublicAggregateMacro(ImageType, OutputDeformation);

  /// Jacobian of output deformation field w.r.t. x
  mirtkPublicAggregateMacro(ImageType, OutputJacobian);

  /// Determinant of Jacobian of output deformation field w.r.t. x
  mirtkPublicAggregateMacro(ImageType, OutputDetJacobian);

  /// Log of determinant of Jacobian of output deformation field w.r.t. x
  mirtkPublicAggregateMacro(ImageType, OutputLogJacobian);

  /// Jacobian of output deformation (or displacement) field w.r.t. v
  mirtkPublicAggregateMacro(ImageType, OutputJacobianDOFs);

  // ---------------------------------------------------------------------------
  // Settings

  /// Attributes of intermediate images, defaults to output attributes
  mirtkPublicAttributeMacro(ImageAttributes, InterimAttributes);

  /// Attributes of output images, defaults to input attributes
  mirtkPublicAttributeMacro(ImageAttributes, OutputAttributes);

  /// Interpolation used for each squaring step
  mirtkPublicAttributeMacro(InterpolationMode, Interpolation);

  /// Whether to compute interpolation coefficients from the given input.
  /// Set to \c false if the input is the coefficient image.
  mirtkPublicAttributeMacro(bool, ComputeInterpolationCoefficients);

  /// Whether filter should invert the velocities
  /// This is equivalent to changing the sign of the upper integration limit.
  mirtkPublicAttributeMacro(bool, ComputeInverse);

  /// Upper integration limit (negative <=> inverse)
  mirtkPublicAttributeMacro(double, UpperIntegrationLimit);

  /// Number of integration steps
  mirtkPublicAttributeMacro(int, NumberOfSteps);

  /// Number of squaring steps
  mirtkPublicAttributeMacro(int, NumberOfSquaringSteps);

  /// Maximum velocity after scaling
  ///
  /// Set to zero in order to scale velocities by exactly
  /// pow(2.0, _NumberOfSquaringSteps). Otherwise, the number of squaring steps
  /// is increased in order to ensure that the maximum velocity in each
  /// dimension is less or equal the specified value.
  mirtkPublicAttributeMacro(double, MaxScaledVelocity);

  /// Whether to upsample the input velocity field
  mirtkPublicAttributeMacro(bool, Upsample);

  /// Whether to use a Gaussian pyramid downsampler when downsampling the
  /// previously upsampled output fields to obey Shannon's sampling theorem.
  /// Otherwise a simple interpolation without smoothing kernel is used.
  mirtkPublicAttributeMacro(bool, SmoothBeforeDownsampling);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Free allocated memory
  void Clear();

public:

  /// Constructor
  ScalingAndSquaring();

  /// Destructor
  virtual ~ScalingAndSquaring();

  // ---------------------------------------------------------------------------
  // Evaluation

protected:

  /// Initialize filter
  virtual void Initialize();

  /// Finalize filter
  virtual void Finalize();

  /// Resample intermediate filter output
  void Resample(const ImageType *, ImageType *);

public:

  /// Compute output deformation/displacement field and its derivatives
  virtual void Run();

};


} // namespace mirtk

#endif // MIRTK_ScalingAndSquaring_H
