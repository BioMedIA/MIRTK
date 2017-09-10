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

#ifndef MIRTK_LinearElasticityConstraint_H
#define MIRTK_LinearElasticityConstraint_H

#include "mirtk/TransformationConstraint.h"


namespace mirtk {


/**
 * Linear elasticity constraint
 *
 * \todo Generalize implementation and extend FreeFormTransformation API such that this
 *       constraint can be used with any deformable transformation model. The current
 *       implementation is specific to the 3D cubic B-spline FFD.
 */
class LinearElasticityConstraint : public TransformationConstraint
{
  mirtkEnergyTermMacro(SmoothnessConstraint, EM_LinearElasticity);

  /// Whether to integrate linear elasticity energy over control point lattice only
  mirtkPublicAttributeMacro(bool, Approximate);

  /// Whether to evaluate derivatives of smoothness term w.r.t. world coordinates.
  ///
  /// When \c false, the smoothness penalty is evaluated w.r.t the local lattice coordinates
  mirtkPublicAttributeMacro(bool, WithRespectToWorld);

  /// Whether to use control point spacing when derivatives are computed w.r.t. world coordinates
  mirtkPublicAttributeMacro(bool, UseLatticeSpacing);

  /// Whether to include rotation component (default: true)
  mirtkPublicAttributeMacro(bool, ConstrainRotation);

  /// Lame's first parameter (default: 0)
  ///
  /// This parameter is multiplied with the constraint weight.
  mirtkPublicAttributeMacro(double, Lambda);

  /// Lame's second parameter (default: 1)
  ///
  /// This parameter is multiplied with the constraint weight.
  mirtkPublicAttributeMacro(double, Mu);

  /// Jacobian matrices of transformation
  mirtkAttributeMacro(Array<Matrix>, Jacobian);

public:

  /// Constructor
  LinearElasticityConstraint(const char * = "", double = 1.);

  // ---------------------------------------------------------------------------
  // Parameters

protected:

  /// Set parameter value from string
  virtual bool SetWithoutPrefix(const char *, const char *);

public:

  // Import other overloads
  using TransformationConstraint::Parameter;

  /// Get parameter key/value as string map
  virtual ParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Evaluation

  /// Initialize energy term once input and parameters have been set
  virtual void Initialize();

  /// Update internal state upon change of input
  virtual void Update(bool = true);

protected:

  /// Compute penalty for current transformation estimate
  virtual double Evaluate();

  /// Compute gradient of penalty term w.r.t transformation parameters
  virtual void EvaluateGradient(double *, double, double);

public:

  // ---------------------------------------------------------------------------
  // Debugging

  /// Write gradient of penalty term
  virtual void WriteGradient(const char *, const char *) const;

};


} // namespace mirtk

#endif // MIRTK_LinearElasticityConstraint_H
