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

#ifndef MIRTK_ConstraintMeasure_H
#define MIRTK_ConstraintMeasure_H

#include "mirtk/EnergyMeasure.h"


namespace mirtk {


// -----------------------------------------------------------------------------
/// Enumeration of available transformation regularization terms
///
/// \note This enumeration constains only a subset of all EnergyMeasure
///       enumeration values, whereby the integer value of corresponding
///       enumeration values is equal.
///
/// \sa EnergyMeasure
enum ConstraintMeasure
{
  CM_Unknown              = CM_Begin,                ///< Unknown/invalid regularizer
  CM_VolumePreservation   = EM_VolumePreservation,   ///< Volume preservation constraint
  CM_TopologyPreservation = EM_TopologyPreservation, ///< Topology preservation constraint
  CM_Sparsity             = EM_Sparsity,             ///< Default sparsity constraint
  CM_BendingEnergy        = EM_BendingEnergy,        ///< Thin-plate spline bending energy
  CM_LinearElasticity     = EM_LinearElasticity,     ///< Linear elastic energy
  CM_L0Norm               = EM_L0Norm,               ///< Sparsity constraint based on l0-norm
  CM_L1Norm               = EM_L1Norm,               ///< Sparsity constraint based on l1-norm
  CM_L2Norm               = EM_L2Norm,               ///< Sparsity constraint based on l2-norm
  CM_SqLogDetJac          = EM_SqLogDetJac,          ///< Squared logarithm of the Jacobian determinant
  CM_NegDetJac            = EM_NegDetJac             ///< Penalize negative Jacobian determinant
};

// -----------------------------------------------------------------------------
template <>
inline string ToString(const ConstraintMeasure &cm, int w, char c, bool left)
{
  EnergyMeasure em = static_cast<EnergyMeasure>(cm);
  if (em <= CM_Begin || em >= CM_End) return ToString("Unknown", w, c, left);
  return ToString(em, w, c, left);
}

// -----------------------------------------------------------------------------
template <>
inline bool FromString(const char *str, ConstraintMeasure &cm)
{
  EnergyMeasure em = EM_Unknown;
  if (FromString(str, em) && CM_Begin < em && em < CM_End) {
    cm = static_cast<ConstraintMeasure>(em);
    return true;
  } else {
    cm = CM_Unknown;
    return false;
  }
}


} // namespace mirtk

#endif // MIRTK_ConstraintMeasure_H
