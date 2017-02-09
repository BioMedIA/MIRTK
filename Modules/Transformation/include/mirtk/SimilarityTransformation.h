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

#ifndef MIRTK_SimilarityTransformation_H
#define MIRTK_SimilarityTransformation_H

#include "mirtk/RigidTransformation.h"


namespace mirtk {


/**
 * Class for similarity transformations.
 *
 * This class defines and implements similarity transformations. In addition to
 * the rigid body transformation parameters, similarity transformations are
 * parameterized by a global scaling parameter. The scaling parameter defines
 * the scaling along all axis of the coordinate transformations.
 */
class SimilarityTransformation : public RigidTransformation
{
  mirtkTransformationMacro(SimilarityTransformation);

protected:

  /// Update transformation matrix after change of parameter
  virtual void UpdateMatrix();

  /// Update transformation parameters after change of matrix
  virtual void UpdateDOFs();

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Default constructor with given number of parameters
  SimilarityTransformation(int);

  /// Copy constructor with given number of parameters
  SimilarityTransformation(const RigidTransformation &, int);

  /// Copy constructor with given number of parameters
  SimilarityTransformation(const SimilarityTransformation &, int);

public:

  /// Default constructor
  SimilarityTransformation();

  /// Copy constructor
  SimilarityTransformation(const RigidTransformation &);

  /// Copy constructor
  SimilarityTransformation(const SimilarityTransformation &);

  /// Destructor
  virtual ~SimilarityTransformation();

  // ---------------------------------------------------------------------------
  // Approximation

  /// Approximate displacements: This function takes a set of points and a set
  /// of displacements and finds !new! parameters such that the resulting
  /// transformation approximates the displacements as good as possible.
  virtual void ApproximateDOFs(const double *, const double *, const double *, const double *,
                               const double *, const double *, const double *, int);

  // ---------------------------------------------------------------------------
  // Transformation parameters

  /// Copy active transformation parameters (DoFs) from given
  /// transformation if possible and return \c false, otherwise
  virtual bool CopyFrom(const Transformation *);

  /// Puts scaling factor
  virtual void PutScale(double);

  /// Gets scaling factor
  virtual double GetScale() const;

  /// Construct a matrix based on parameters passed in the array
  static Matrix DOFs2Matrix(const double *);

  // ---------------------------------------------------------------------------
  // Derivatives

  // Do not overwrite other base class overloads
  using Transformation::JacobianDOFs;

  /// Calculates the Jacobian of the transformation w.r.t the parameters
  virtual void JacobianDOFs(double [3], int, double, double, double, double = 0, double = -1) const;

  /// Calculates the derivative of the Jacobian of the transformation (w.r.t. world coordinates) w.r.t. a transformation parameter
  virtual void DeriveJacobianWrtDOF(Matrix &, int, double, double, double, double = 0, double = -1) const;

  // ---------------------------------------------------------------------------
  // I/O

  // Do not hide methods of base class
  using RigidTransformation::Print;
  using RigidTransformation::Write;

  /// Prints the parameters of the transformation
  virtual void Print(ostream &, Indent = 0) const;

  /// Whether this transformation can read a file of specified type (i.e. format)
  virtual bool CanRead(TransformationType) const;

  /// Writes transformation to a file stream
  virtual Cofstream &Write(Cofstream &) const;

protected:

  /// Reads transformation from a file stream
  virtual Cifstream &ReadDOFs(Cifstream &, TransformationType);

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Transformation parameters
// =============================================================================

// -----------------------------------------------------------------------------
inline void SimilarityTransformation::PutScale(double s)
{
  Put(SG, s);
}

// -----------------------------------------------------------------------------
inline double SimilarityTransformation::GetScale() const
{
  return Get(SG);
}


} // namespace mirtk

#endif // MIRTK_SimilarityTransformation_H
