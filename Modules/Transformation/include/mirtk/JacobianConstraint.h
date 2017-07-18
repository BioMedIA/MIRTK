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

#ifndef MIRTK_JacobianConstraint_H
#define MIRTK_JacobianConstraint_H

#include "mirtk/TransformationConstraint.h"

#include "mirtk/Matrix.h"


namespace mirtk {


class FreeFormTransformation;


/**
 * Base class of soft transformation constraints penalizing the Jacobian determinant
 *
 * \note The Jacobian, Penalty, and Derivative functions are public such that
 *       these can be called by a TBB parallel_for/parall_reduce body in a subclass
 *       implementation. These should not be considered as public API and may be
 *       subject to change.
 */
class JacobianConstraint : public TransformationConstraint
{
  mirtkAbstractMacro(JacobianConstraint);

  // ---------------------------------------------------------------------------
  // Attributes

protected:

  double *_DetJacobian; ///< Determinant of Jacobian at each control point
  Matrix *_AdjJacobian; ///< Adjugate of Jacobian at each control point
  int     _NumberOfCPs; ///< Number of control points

  // ---------------------------------------------------------------------------
  // Construction/destruction

  /// Constructor
  JacobianConstraint(const char * = "");

public:

  /// Destructor
  virtual ~JacobianConstraint();

  // ---------------------------------------------------------------------------
  // Evaluation

  /// Update internal state upon change of input
  virtual void Update(bool = true);

  /// Compute determinant and adjugate of Jacobian of transformation
  ///
  /// This function is called by Update for each (control) point at which the
  /// Jacobian determinant is evaluated. Can be overridden in subclasses to
  /// either compute a different Jacobian or to threshold the determinant value.
  virtual double Jacobian(const FreeFormTransformation *ffd,
                          double x, double y, double z, double t,
                          Matrix &adj) const
  {
    double det;
    ffd->Jacobian(adj, x, y, z, t, t);
    adj.Adjugate(det);
    return det;
  }

  /// Evaluate penalty at control point location given Jacobian determinant value
  ///
  /// \param[in] det Jacobian determinant value.
  virtual double Penalty(double det) const = 0;

  /// Evaluate derivative of penalty w.r.t. Jacobian determinant value
  ///
  /// \param[in] det Jacobian determinant value.
  virtual double DerivativeWrtJacobianDet(double det) const = 0;

protected:

  /// Compute penalty for current transformation estimate
  virtual double Evaluate();

  /// Compute gradient of penalty term w.r.t transformation parameters
  virtual void EvaluateGradient(double *, double, double);

  // ---------------------------------------------------------------------------
  // Debugging

public:

  /// Write Jacobian determinant
  virtual void WriteDataSets(const char *, const char *, bool = true) const;

protected:

  /// Write Jacobian determinant values to scalar image
  void WriteJacobian(const char *, const FreeFormTransformation *, const double *) const;

};


} // namespace mirtk

#endif // MIRTK_JacobianConstraint_H
