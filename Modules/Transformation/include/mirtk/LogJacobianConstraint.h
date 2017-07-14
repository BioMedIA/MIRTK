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
 * This constraint prevents folding of the transformation parameterization,
 * i.e., either of the control point displacements or velocities. It preserves
 * volume in case of a classical FFD model and is in this case equivalent to
 * the VolumePreservationConstraint.
 *
 * Torsten Rohlﬁng and Calvin R. Maurer, Jr., Intensity-Based Non-rigid
 * Registration Using Adaptive Multilevel Free-Form Deformation with an
 * Incompressibility Constraint, MICCAI 2001.
 *
 * Modat et al., Log-Euclidean free-form deformation, SPIE 2011.
 *
 * @sa VolumePreservationConstraint
 */
class LogJacobianConstraint : public JacobianConstraint
{
  mirtkEnergyTermMacro(LogJacobianConstraint, EM_SqLogDetJac);

public:

  /// Constructor
  LogJacobianConstraint(const char * = "");

  /// Destructor
  virtual ~LogJacobianConstraint();

  /// Compute determinant and adjugate of Jacobian of transformation
  virtual double Jacobian(const FreeFormTransformation *ffd,
                          double x, double y, double z, double t,
                          Matrix &adj) const
  {
    double det;
    ffd->FFDJacobianWorld(adj, x, y, z, t, t);
    adj.Adjugate(det);
    if (det < 1e-7) det = 1e-7;
    return det;
  }

  /// Evaluate penalty at control point location given Jacobian determinant value
  virtual double Penalty(double det) const
  {
    double logdet = log(det);
    return logdet * logdet;
  }

  /// Evaluate penalty derivative at control point location w.r.t. Jacobian determinant value
  virtual double DerivativeWrtJacobianDet(double det) const
  {
    return 2.0 * log(det) / det;
  }

};


} // namespace mirtk

#endif // MIRTK_LogJacobianConstraint_H
