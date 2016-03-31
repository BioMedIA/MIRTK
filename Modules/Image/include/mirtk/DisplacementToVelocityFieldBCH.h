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

#ifndef MIRTK_DisplacementToVelocityFieldBCH_H
#define MIRTK_DisplacementToVelocityFieldBCH_H

#include "mirtk/GenericImage.h"
#include "mirtk/DisplacementToVelocityField.h"
#include "mirtk/VelocityToDisplacementField.h"


namespace mirtk {


/**
 * Computes a stationary velocity field from a given displacement field.
 *
 * This class implements an image filter which computes the group logarithm
 * of a (diffeomorphic) displacement field using the Baker-Campbell-Hausdorff (BCH)
 * formula. The result is a stationary velocity field. Integrating this velocity
 * field from 0 to _T, e.g., using N forward Euler integration steps or the
 * scaling and squaring (SS) method, yields the original displacement field
 * (with an approximation error).
 *
 * M. Bossa and S. Olmos, "A new algorithm for the computation of the group
 * logarithm of diffeomorphisms", MFCA'08
 */
template <class VoxelType>
class DisplacementToVelocityFieldBCH : public DisplacementToVelocityField<VoxelType>
{
  mirtkObjectMacro(DisplacementToVelocityFieldBCH);

  // ---------------------------------------------------------------------------
  // Types
public:

  /// Type of internal images
  typedef GenericImage<VoxelType> ImageType;

  /// Type of exponential filter
  typedef VelocityToDisplacementField<VoxelType> ExponentialFilterType;

  // ---------------------------------------------------------------------------
  // Attributes
private:

  /// Image filter for computation of exponential map of inverse velocities
  mirtkReadOnlyAggregateMacro(ExponentialFilterType, ExponentialFilter);

  /// Whether exponential image filter was set by user
  bool _CustomExponentialFilter;

  /// Number of update steps
  mirtkPublicAttributeMacro(int, NumberOfIterations);

  /// Number of terms to use
  mirtkPublicAttributeMacro(int, NumberOfTerms);

  /// Whether to use the Jacobian for the Lie bracket computation
  mirtkPublicAttributeMacro(bool, UseJacobian);

  /// Whether to smooth the velocity fields to stabilize the computation
  mirtkPublicAttributeMacro(bool, SmoothVelocities);

  /// Result of composition: exp(-v) ° phi
  ImageType *_dv;

  /// Result of Lie bracket: [v, dv]
  ImageType *_l1;

  /// Result of Lie bracket: [v, [v, dv]]
  ImageType *_l2;

  /// Result of Lie bracket: [dv, [v, dv]]
  ImageType *_l3;

  /// Result of Lie bracket: [dv, [v, [v, dv]]]
  ImageType *_l4;

  // ---------------------------------------------------------------------------
  // Construction/Destruction
public:

  /// Constructor
  DisplacementToVelocityFieldBCH();

  /// Destructor
  virtual ~DisplacementToVelocityFieldBCH();

  // ---------------------------------------------------------------------------
  // Setter/Getters

  /// Set image filter for computation of exponential map of inverse velocities
  void ExponentialFilter(ExponentialFilterType *);

  /// Set upper integration limit of exponential filter
  void UpperIntegrationLimit(double);

  /// Get upper integration limit of exponential filter
  double UpperIntegrationLimit();

  /// Set number of integration steps
  void NumberOfSteps(int);

  /// Get number of integration steps
  int NumberOfSteps();

  // ---------------------------------------------------------------------------
  // Execution

  /// Compute output = log(input)
  virtual void Run();

protected:

  /// Initialize filter
  virtual void Initialize();

  /// Finalize filter
  virtual void Finalize();

};


} // namespace mirtk

#endif // MIRTK_DisplacementToVelocityFieldBCH_H
