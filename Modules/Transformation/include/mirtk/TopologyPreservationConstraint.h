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

#ifndef MIRTK_TopologyPreservationConstraint_H
#define MIRTK_TopologyPreservationConstraint_H

#include "mirtk/JacobianConstraint.h"


namespace mirtk {


/**
 * Topology preservation constraint for deformable image registration
 *
 * Rueckert et al. (2006). Diffeomorphic registration using B-splines. MICCAI, 702–709.
 */
class TopologyPreservationConstraint : public JacobianConstraint
{
  mirtkEnergyTermMacro(TopologyPreservationConstraint, EM_TopologyPreservation);

  /// Jacobian determinant threshold
  mirtkPublicAttributeMacro(double, Gamma);

public:

  /// Constructor
  TopologyPreservationConstraint(const char * = "");

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
  // Penalty

public:

  /// Evaluate penalty at control point location given Jacobian determinant value
  virtual double Penalty(double det) const
  {
    if (det >= _Gamma) return 0.;
    det *= 10. * det;
    return det + (1. / det) - 2.;
  }

  /// Evaluate derivative of penalty w.r.t. Jacobian determinant value
  virtual double DerivativeWrtJacobianDet(double det) const
  {
    if (det >= _Gamma) return 0.;
    return 20. * det - 2. / (10. * det * det * det);
  }

};


} // namespace mirtk

#endif // MIRTK_TopologyPreservationConstraint_H
