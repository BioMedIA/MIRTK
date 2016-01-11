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

#ifndef MIRTK_SmoothnessConstraint_H
#define MIRTK_SmoothnessConstraint_H

#include <mirtkTransformationConstraint.h>


namespace mirtk {


/**
 * Smoothness term based on bending energy of transformation
 *
 * Daniel Rueckert et al., Nonrigid registration using free-form deformations:
 * Application to breast MR images, IEEE Transactions on Medical Imaging, 18(8),
 * 712–21 (1999; doi:10.1109/42.796284)
 */
class SmoothnessConstraint : public TransformationConstraint
{
  mirtkEnergyTermMacro(SmoothnessConstraint, EM_BendingEnergy);

  /// Whether to evaluate smoothness w.r.t world coordinates.
  /// Otherwise, the smoothness penalty is evaluated w.r.t the local lattice
  /// coordinates of the free-form deformation.
  ///
  /// \attention The use of this option is experimental!
  mirtkPublicAttributeMacro(bool, WithRespectToWorld);

  /// Smoothness penalty annealing rate.
  /// Used in conjuction with RobustPointMatch.
  ///
  /// \attention The use of this option is experimental and it may be removed.
  mirtkPublicAttributeMacro(double, AnnealingRate);

  /// Current additional weight factor due to annealing process
  double _AnnealingWeight;

public:

  /// Constructor
  SmoothnessConstraint(const char * = "", double = 1.0);

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

  /// Update energy term after convergence
  virtual bool Upgrade();

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

protected:

  /// Write gradient of penalty term w.r.t. control point parameters
  void WriteFFDGradient(const char *, const FreeFormTransformation *, const double *) const;

};


} // namespace mirtk

#endif // MIRTK_SmoothnessConstraint_H 
