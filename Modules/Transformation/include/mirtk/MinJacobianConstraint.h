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

#ifndef MIRTK_MinJacobianConstraint_H
#define MIRTK_MinJacobianConstraint_H

#include <mirtkJacobianConstraint.h>


namespace mirtk {


/**
 * Penalizes non-diffeomorphic transformations by enforcing a minimum Jacobian determinant
 *
 * Rueckert et al., Diffeomorphic Registration Using B-Splines, MICCAI 2006.
 */
class MinJacobianConstraint : public JacobianConstraint
{
  mirtkEnergyTermMacro(MinJacobianConstraint, EM_MinDetJac);

  /// Lower Jacobian determinant threshold, only determinant values less or
  /// equal this threshold are penalized
  mirtkPublicAttributeMacro(double, Gamma);

public:

  /// Constructor
  MinJacobianConstraint(const char * = "");

  /// Destructor
  virtual ~MinJacobianConstraint();

  // ---------------------------------------------------------------------------
  // Parameters

protected:

  /// Set parameter value from string
  virtual bool SetWithoutPrefix(const char *, const char *);

public:

  // Import other overloads
  using JacobianConstraint::Parameter;

  /// Get parameter key/value as string map
  virtual ParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Evaluation

protected:

  /// Compute penalty for current transformation estimate
  virtual double Evaluate();

  /// Compute gradient of penalty term w.r.t transformation parameters
  virtual void EvaluateGradient(double *, double, double);

};


} // namespace mirtk

#endif // MIRTK_MinJacobianConstraint_H
