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

#ifndef MIRTK_LogJacobianConstraint_H
#define MIRTK_LogJacobianConstraint_H

#include "mirtk/JacobianConstraint.h"


namespace mirtk {


class FreeFormTransformation;


/**
 * Constrains log of Jacobian determinant of FFD parameterization
 *
 * This constraint is based on a penalty imposed by the squared logarithm of
 * the Jacobian determinant. It strongly penalizes small and negative
 * Jacobian determinant values, and slowly increases for large Jacobian
 * determinant values. It therefore locally preserves volume because the
 * penalty is zero only for a Jacobian determinant value of one.
 * Volume increase is however less penalized than volume reduction.
 *
 * By default, the Jacobian determinant of the spline function is penalised
 * by this class. For a classic FFD model parameterized by control point
 * displacements, this is identical to the VolumePreservationConstraint.
 * For a velocity based model, the Jacobian of the velocity field is constraint.
 *
 * Torsten Rohlﬁng and Calvin R. Maurer, Jr., Intensity-Based Non-rigid
 * Registration Using Adaptive Multilevel Free-Form Deformation with an
 * Incompressibility Constraint, MICCAI 2001.
 *
 * @sa VolumePreservationConstraint
 */
class LogJacobianConstraint : public JacobianConstraint
{
  mirtkEnergyTermMacro(LogJacobianConstraint, EM_SqLogDetJac);

  /// Small value below which Jacobian determinant is penalised linearly
  /// with increasing smaller (i.e., negative) value
  ///
  /// The penalty implemented in IRTK can be used with _Epsilon set to a very
  /// negative value, i.e., negative with great magnitude, -inf, or NaN.
  mirtkPublicAttributeMacro(double, Epsilon);

public:

  /// Constructor
  LogJacobianConstraint(const char * = "", bool = true);

  /// Destructor
  virtual ~LogJacobianConstraint();

protected:

  /// Set parameter value from string
  virtual bool SetWithPrefix(const char *, const char *);

  /// Set parameter value from string
  virtual bool SetWithoutPrefix(const char *, const char *);

public:

  // Import other overloads
  using JacobianConstraint::Parameter;

  /// Get parameter key/value as string map
  virtual ParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Penalty

public:

  /// Evaluate penalty at control point location given Jacobian determinant value
  virtual double Penalty(double det) const
  {
    if (det < _Epsilon) {
      double l = log(_Epsilon);
      double m = 2. * log(_Epsilon) / _Epsilon;
      double t = l * l - m * _Epsilon;
      return m * det + t;
    } else {
      double l = log(det);
      return l * l;
    }
  }

  /// Evaluate penalty derivative at control point location w.r.t. Jacobian determinant value
  virtual double DerivativeWrtJacobianDet(double det) const
  {
    if (det < _Epsilon) det = _Epsilon;
    return 2. * log(det) / det;
  }

};


} // namespace mirtk

#endif // MIRTK_LogJacobianConstraint_H
