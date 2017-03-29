/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2017 Imperial College London
 * Copyright 2017 Andreas Schuh
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

#ifndef MIRTK_MeanSquaredDisplacementError_H
#define MIRTK_MeanSquaredDisplacementError_H

#include "mirtk/DataFidelity.h"
#include "mirtk/RegisteredImage.h"


namespace mirtk {


/**
 * Constrains displacement of each target voxel to be close to another transformation
 *
 * This constraint can be used to enforce some consistency of the optimized
 * displacements with a given mean or pre-composed target displacement.
 * The constraint may be relaxed for voxels with a stronger misalignment.
 */
class MeanSquaredDisplacementError : public DataFidelity
{
  mirtkEnergyTermMacro(DisplacementConstraint, EM_MeanSquaredDisplacementError);

  // ---------------------------------------------------------------------------
  // Types
public:

  /// Type of non-parametric gradient components
  typedef double   GradientType;

  /// Type of non-parametric gradient image
  typedef GenericImage<GradientType>   GradientImageType;

  /// Type of cached displacement fields
  typedef RegisteredImage::DisplacementImageType   DisplacementImageType;

  // ---------------------------------------------------------------------------
  // Attributes
protected:

  /// Finite regular grid over which to integrate error
  mirtkPublicAttributeMacro(ImageAttributes, Domain);

  /// Target transformation
  mirtkPublicAggregateMacro(const class Transformation, TargetTransformation);

  /// Target displacements
  mirtkAttributeMacro(DisplacementImageType, TargetDisplacement);

  /// Current displacements
  ///
  /// This displacement field is used only when the transformation requires
  /// the caching of the entire dense displacements for efficienty reasons
  /// such as the scaling and squaring of a SV FFD. When an externally
  /// computed/updated displacement field is given, it is used instead.
  mirtkAttributeMacro(DisplacementImageType, CurrentDisplacement);

  /// Pre-computed displacements which are updated by an external process
  ///
  /// When this displacement field is given, it is used instead of the
  /// current transformation. The external process managing the optimization
  /// has to ensure that the displacement field is updated whenever the
  /// transformation changes.
  mirtkPublicAggregateMacro(DisplacementImageType, ExternalDisplacement);

  /// Non-parametric gradient
  mirtkAttributeMacro(GradientImageType, NonParametricGradient);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
public:

  /// Constructor
  MeanSquaredDisplacementError(const char * = "", double = 1.);

  // ---------------------------------------------------------------------------
  // Evaluation

  /// Initialize energy term once input and parameters have been set
  virtual void Initialize();

  /// Update internal state after change of DoFs
  ///
  /// \param[in] gradient Whether to also update internal state for evaluation
  ///                     of energy gradient. If \c false, only the internal state
  ///                     required for the energy evaluation need to be updated.
  virtual void Update(bool gradient = true);

protected:

  /// Compute penalty for current transformation estimate
  virtual double Evaluate();

  /// Compute gradient of penalty term w.r.t transformation parameters
  virtual void EvaluateGradient(double *, double, double);

};


} // namespace mirtk

#endif // MIRTK_MeanSquaredDisplacementError_H
