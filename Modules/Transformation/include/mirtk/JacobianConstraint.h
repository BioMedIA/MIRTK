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
#include "mirtk/GenericImage.h"


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
  // Types

public:

  /// Enumeration of target image sub-domains on which to evaluate penalty
  enum SubDomainEnum
  {
    SD_Image,   ///< Evaluate penalty at image voxels
    SD_Lattice, ///< Evaluate penalty at control points
    SD_SubDiv,  ///< Evaluate penalty at subdivided control point lattice
    SD_Shifted  ///< Evaluate penalty in-between control points
  };

  // ---------------------------------------------------------------------------
  // Attributes

  /// Sub-domain on which to evaluate penalty
  mirtkPublicAttributeMacro(SubDomainEnum, SubDomain);

  /// Whether to always apply constraint to Jacobian of FFD spline function
  ///
  /// This option has no effect when the FFD parameters are displacements.
  /// In case of a velocity based FFD model, it applies the constraint to
  /// the velocity field instead of the corresponding displacement field.
  mirtkReadOnlyAttributeMacro(bool, ConstrainParameterization);

  /// Whether to evaluate derivatives of smoothness term w.r.t. world coordinates.
  ///
  /// When \c false, the smoothness penalty is evaluated w.r.t the local lattice coordinates
  mirtkPublicAttributeMacro(bool, WithRespectToWorld);

  /// Whether to use control point spacing when derivatives are computed w.r.t. world coordinates
  mirtkPublicAttributeMacro(bool, UseLatticeSpacing);

  /// Whether to apply symmetric penalty in case of SVFFD model (default: true)
  mirtkPublicAttributeMacro(bool, Symmetric);

protected:

  int     _NumJacobian; ///< Number of allocated Jacobian matrices
  double *_DetJacobian; ///< Determinant of Jacobian at each control point
  Matrix *_AdjJacobian; ///< Adjugate of Jacobian at each control point
  Array<Matrix> _MatW2L; ///< World to sub-domain lattice coordinates
  Array<Matrix> _MatL2W; ///< Sub-domain lattice coordinates to world
  Array<ImageAttributes> _SubDomains; ///< Discrete sub-domain over which to integrate penalty

  // ---------------------------------------------------------------------------
  // Construction/destruction

  /// Constructor
  JacobianConstraint(const char *, bool = false);

public:

  /// Destructor
  virtual ~JacobianConstraint();

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

  /// Write gradient of penalty term
  virtual void WriteGradient(const char *, const char *) const;

};

///////////////////////////////////////////////////////////////////////////////
// Inline definitions
///////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
template <>
inline string ToString(const JacobianConstraint::SubDomainEnum &value, int w, char c, bool left)
{
  const char *str;
  switch (value) {
    case JacobianConstraint::SD_Image:   str = "Image";   break;
    case JacobianConstraint::SD_Lattice: str = "Lattice"; break;
    case JacobianConstraint::SD_Shifted: str = "Shifted"; break;
    case JacobianConstraint::SD_SubDiv:  str = "SubDiv";  break;
    default:                             str = "Unknown"; break;
  }
  return ToString(str, w, c, left);
}

// ----------------------------------------------------------------------------
template <>
inline bool FromString(const char *str, JacobianConstraint::SubDomainEnum &value)
{
  string lstr = ToLower(str);
  if (lstr == "image" || lstr == "target" || lstr == "imagelattice") {
    value = JacobianConstraint::SD_Image;
  } else if (lstr == "lattice" || lstr == "ffd" || lstr == "ffdlattice" || lstr == "cp" || lstr == "cps" || lstr == "controlpoints") {
    value = JacobianConstraint::SD_Lattice;
  } else if (lstr == "shifted" || lstr == "betweencontrolpoints") {
    value = JacobianConstraint::SD_Shifted;
  } else if (lstr == "subdiv" || lstr == "subdivided" || lstr == "subdividedlattice") {
    value = JacobianConstraint::SD_SubDiv;
  } else {
    return false;
  }
  return true;
}


} // namespace mirtk

#endif // MIRTK_JacobianConstraint_H
