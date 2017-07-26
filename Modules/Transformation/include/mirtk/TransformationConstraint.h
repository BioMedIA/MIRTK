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

#ifndef MIRTK_TransformationConstraint_H
#define MIRTK_TransformationConstraint_H

#include "mirtk/EnergyTerm.h"
#include "mirtk/ConstraintMeasure.h"
#include "mirtk/ImageAttributes.h"
#include "mirtk/FreeFormTransformation.h"
#include "mirtk/MultiLevelTransformation.h"


namespace mirtk {


/**
 * Base class for a penalty term imposed on a transformation
 *
 * The penalty is minimized by the registration using an instance of
 * RegistrationEnergy with set data similarity and regularization terms.
 * Higher penalties lead to stronger enforcement of the constraint.
 */
class TransformationConstraint : public EnergyTerm
{
  mirtkAbstractMacro(TransformationConstraint);

  /// Image domain on which penalty is applied
  mirtkPublicAttributeMacro(ImageAttributes, Domain);

  /// Whether to apply constraint also at passive DoFs (control points)
  mirtkPublicAttributeMacro(bool, ConstrainPassiveDoFs);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
protected:

  /// Constructor
  TransformationConstraint(const char * = "", double = 1.0);

  /// Copy constructor
  TransformationConstraint(const TransformationConstraint &);

  /// Assignment operator
  TransformationConstraint &operator =(const TransformationConstraint &);

public:

  /// Instantiate new transformation constraint of given kind
  static TransformationConstraint *New(ConstraintMeasure,
                                       const char * = "", double = 1.0);

  /// Destructor
  virtual ~TransformationConstraint();

  // ---------------------------------------------------------------------------
  // Configuration

protected:

  /// Set parameter value from string
  virtual bool SetWithPrefix(const char *, const char *);

  /// Set parameter value from string
  virtual bool SetWithoutPrefix(const char *, const char *);

public:

  // Import other overloads
  using EnergyTerm::Parameter;

  /// Get parameter name/value pairs
  virtual ParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Subclass helper
protected:

  /// Get transformation as specific type or NULL if dynamic cast fails
  template <class TransformationType>
  const TransformationType *TransformationAs() const;

  /// Get transformation as free-form deformation or NULL if it is none
  const FreeFormTransformation *FFD() const;

  /// Get transformation as multi-level transformation or NULL if it is none
  const MultiLevelTransformation *MFFD() const;

  /// Write gradient of penalty term w.r.t. control point parameters
  void WriteFFDGradient(const char *, const FreeFormTransformation *, const double *) const;

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Subclass helper
// =============================================================================

// -----------------------------------------------------------------------------
template <class TransformationType>
inline const TransformationType *TransformationConstraint::TransformationAs() const
{
  return dynamic_cast<const TransformationType *>(Transformation());
}

// -----------------------------------------------------------------------------
inline const FreeFormTransformation *TransformationConstraint::FFD() const
{
  return TransformationAs<FreeFormTransformation>();
}

// -----------------------------------------------------------------------------
inline const MultiLevelTransformation *TransformationConstraint::MFFD() const
{
  return TransformationAs<MultiLevelTransformation>();
}


} // namespace mirtk

#endif // MIRTK_TransformationConstraint_H
