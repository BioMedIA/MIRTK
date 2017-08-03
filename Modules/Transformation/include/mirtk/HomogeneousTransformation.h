/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2015 Imperial College London
 * Copyright 2008-2013 Daniel Rueckert, Julia Schnabel
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

#ifndef MIRTK_HomogeneousTransformation_H
#define MIRTK_HomogeneousTransformation_H

#include "mirtk/Transformation.h"

#include "mirtk/Matrix.h"

#include <cmath>


namespace mirtk {


class MultiThreadedImageHomogeneousTransformation;

enum HomogeneousTransformationParameterIndex
{
  TX, TY, TZ, RX, RY, RZ, SX, SY, SZ, SXY, SYZ, SXZ, SYX, SZY, SZX,
  SG = SX // global/isotropic scaling for similarity transformation
};


/**
 * Base class for homogeneous transformations.
 *
 * This class defines and implements homogeneous transformations which can
 * be represented by a 4 x 4 transformation matrix. The transformation of
 * a point is implemented by post multiplying the transformation matrix with
 * the point in homogeneous coordinates. The transformation is parameterized
 * by twelve degrees of freedom.
 *
 */

class HomogeneousTransformation : public Transformation
{
  mirtkTransformationMacro(HomogeneousTransformation);

protected:

  // ---------------------------------------------------------------------------
  // Types

  enum AttributeSelector { MATRIX, DOFS, DOFS_MATRIX };

  // ---------------------------------------------------------------------------
  // Data members

  /// 4x4 transformation matrix for homogeneous coordinates
  Matrix _matrix;

  /// Inverse 4x4 transformation matrix for homogeneous coordinates
  Matrix _inverse;

  /// Updates transformation matrix after change of parameter
  virtual void UpdateMatrix();

  /// Updates transformation parameters after change of matrix
  virtual void UpdateDOFs();

  /// Update transformation parameters and/or matrices
  void Update(AttributeSelector);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Default constructor with given number of parameters
  HomogeneousTransformation(int);

  /// Copy constructor with given number of parameters
  HomogeneousTransformation(const HomogeneousTransformation &, int);

public:

  /// Default constructor
  HomogeneousTransformation();

  /// Construct from 4x4 transformation matrix (without checks)
  HomogeneousTransformation(const Matrix &);

  /// Copy Constructor
  HomogeneousTransformation(const HomogeneousTransformation &);

  /// Destructor
  virtual ~HomogeneousTransformation();

  // ---------------------------------------------------------------------------
  // Approximation

  // Import other overloads
  using Transformation::Approximate;
  using Transformation::ApproximateAsNew;

  /// Approximate displacements: This function takes a set of points and a set
  /// of displacements and finds a transformation which approximates these
  /// displacements. After approximation, the displacements are replaced by the
  /// residual displacement errors at the points.
  virtual double Approximate(const ImageAttributes &, double *, double *, double *,
                             int = 1, double = .0);

  /// Approximate displacements: This function takes a set of points and a set
  /// of displacements and finds a transformation which approximates these
  /// displacements. After approximation, the displacements are replaced by the
  /// residual displacement errors at the points.
  virtual double Approximate(const double *, const double *, const double *,
                             double *,       double *,       double *, int,
                             int = 1, double = .0);

  /// Approximate displacements: This function takes a set of points and a set
  /// of displacements and finds a transformation which approximates these
  /// displacements. After approximation, the displacements are replaced by the
  /// residual displacement errors at the points.
  virtual double Approximate(const double *, const double *, const double *, const double *,
                             double *,       double *,       double *,       int,
                             int = 1, double = .0);

  /// Approximate given homogeneous coordinate transformation matrix
  virtual double Approximate(const Matrix &);

  /// Approximate given homogeneous coordinate transformation matrix
  virtual double ApproximateAsNew(const Matrix &);

  /// Approximate displacements: This function takes a set of points and a set
  /// of displacements and finds !new! parameters such that the resulting
  /// transformation approximates the displacements as good as possible.
  virtual void ApproximateDOFs(const double *, const double *, const double *, const double *,
                               const double *, const double *, const double *, int);

  // ---------------------------------------------------------------------------
  // Transformation parameters (DOFs)

  /// Copy active transformation parameters (DoFs) from given
  /// transformation if possible and return \c false, otherwise
  virtual bool CopyFrom(const Transformation *);

  /// Puts a transformation parameter
  virtual void Put(int, DOFValue);

  /// Puts transformation parameters
  virtual void Put(const DOFValue *);

  /// Add change to transformation parameters
  /// \note Also passive parameters are modified by this function.
  virtual void Add(const DOFValue *);

  /// Update transformation parameters given parametric gradient
  /// \note Only active parameters are modified by this function.
  virtual double Update(const DOFValue *);

  /// Set transformation parameters from a homogeneous transformation matrix
  ///
  /// The transformation sets its parameters to match the given homogeneous
  /// transformation matrix as good as possible and then updates its
  /// internal transformation matrix from these parameters. The resulting
  /// homogeneous transformations may therefore differ from the given
  /// transformation matrix if the transformation is not a full 12 DoF
  /// affine transformation, but has fewer free parameters.
  ///
  /// \note Use one of the Approximate or ApproximateAsNew functions to
  ///       minimize the mean squared error of the approximation at given
  ///       evaluation points instead for a more accurate approximation.
  void PutMatrix(const Matrix &);

  /// Gets the transformation matrix
  const Matrix &GetMatrix() const;

  /// Gets the inverse transformation matrix
  const Matrix &GetInverseMatrix() const;

  /// Reset transformation
  virtual void Reset();

  /// Inverts the transformation
  void Invert();

  // ---------------------------------------------------------------------------
  // Point transformation

  // Import other overloads from base class
  using Transformation::Transform;
  using Transformation::Inverse;

  /// Transforms a single point using the global transformation component only
  virtual void GlobalTransform(double &, double &, double &, double = 0, double = -1) const;

  /// Transforms a single point using the local transformation component only
  virtual void LocalTransform(double &, double &, double &, double = 0, double = -1) const;

  /// Transforms a single point
  virtual void Transform(double &, double &, double &, double = 0, double = -1) const;

  /// Transforms a single point using the inverse of the global transformation only
  virtual void GlobalInverse(double &, double &, double &, double = 0, double = -1) const;

  /// Transforms a single point using the inverse of the local transformation only
  virtual bool LocalInverse(double &, double &, double &, double = 0, double = -1) const;

  /// Transforms a single point using the inverse of the transformation
  virtual bool Inverse(double &, double &, double &, double = 0, double = -1) const;

  // ---------------------------------------------------------------------------
  // Derivatives

  /// Calculates the Jacobian of the global transformation w.r.t world coordinates
  virtual void GlobalJacobian(Matrix &, double, double, double, double = 0, double = -1) const;

  /// Calculates the Jacobian of the local transformation w.r.t world coordinates
  virtual void LocalJacobian(Matrix &, double, double, double, double = 0, double = -1) const;

  /// Calculates the Jacobian of the transformation w.r.t world coordinates
  virtual void Jacobian(Matrix &, double, double, double, double = 0, double = -1) const;

  // ---------------------------------------------------------------------------
  // Properties

  /// Checks whether transformation is an identity mapping
  virtual bool IsIdentity() const;

  // ---------------------------------------------------------------------------
  // I/O

  // Do not hide methods of base class
  using Transformation::Print;

  /// Prints the parameters of the transformation
  virtual void Print(ostream &, Indent = 0) const;

protected:
  
  /// Reads transformation parameters from a file stream
  virtual Cifstream &ReadDOFs(Cifstream &, TransformationType);

  // ---------------------------------------------------------------------------
  // Backwards compatibility

public:

  /// \deprecated Replaced by UpdateDOFs, which, however, does not have to be
  ///             called explicitly. The parameters are automatically updated
  ///             when the transformation matrix changed. Some older code may
  ///             yet call this deprecated method after PutMatrix, which no
  ///             longer does anything and the call should be optimized away.
  void UpdateParameter() {}

  friend class MultiThreadedImageHomogeneousTransformation;
};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Transformation parameters
// =============================================================================

// -----------------------------------------------------------------------------
inline void HomogeneousTransformation::Update(AttributeSelector what)
{
  // Update attributes
  switch (what) {
    case MATRIX:
      this->UpdateMatrix();
      break;
    case DOFS:
      this->UpdateDOFs();
      break;
    case DOFS_MATRIX:
      this->UpdateDOFs();
      this->UpdateMatrix();
      break;
  }
  // Mark transformation as changed
  this->Changed(true);
  // Notify observers, keeping changed state intact
  Broadcast(ModifiedEvent, &_matrix);
}

// -----------------------------------------------------------------------------
inline void HomogeneousTransformation::Put(int i, DOFValue x)
{
  Transformation::Put(i, x);
  this->Update(MATRIX);
}

// -----------------------------------------------------------------------------
inline void HomogeneousTransformation::Put(const DOFValue *x)
{
  Transformation::Put(x);
  this->Update(MATRIX);
}

// -----------------------------------------------------------------------------
inline void HomogeneousTransformation::Add(const DOFValue *dx)
{
  double delta, max_delta = .0;
  for (int idx = 0; idx < _NumberOfDOFs; ++idx) {
    // Add change also to passive parameters
    _Param[idx] += dx[idx];
    delta = abs(dx[idx]);
    if (delta > max_delta) max_delta = delta;
  }
  if (max_delta > .0) this->Update(MATRIX);
}

// -----------------------------------------------------------------------------
inline double HomogeneousTransformation::Update(const DOFValue *dx)
{
  double delta, max_delta = .0;
  for (int idx = 0; idx < _NumberOfDOFs; ++idx) {
    if (_Status[idx] == Active) {
      _Param[idx] += dx[idx];
      delta = abs(dx[idx]);
      if (delta > max_delta) max_delta = delta;
    }
  }
  if (max_delta > .0) this->Update(MATRIX);
  return max_delta;
}

// -----------------------------------------------------------------------------
inline void HomogeneousTransformation::PutMatrix(const Matrix &matrix)
{
  if (matrix.Rows() == 3 && matrix.Cols() == 3) {
    for (int c = 0; c < 3; ++c)
    for (int r = 0; r < 3; ++r) {
      _matrix(r, c) = matrix(r, c);
    }
    for (int r = 0; r < 3; ++r) {
      _matrix(r, 3) = 0.;
    }
  } else if (matrix.Rows() == 4 && matrix.Cols() == 4) {
    _matrix = matrix;
  } else {
    Throw(ERR_InvalidArgument, __FUNCTION__, "Transformation matrix must be 3x3 or 4x4");
  }
  _matrix(3, 0) = _matrix(3, 1) = _matrix(3, 2) = 0., _matrix(3, 3) = 1.;
  this->Update(DOFS_MATRIX); // also updates _inverse
}

// -----------------------------------------------------------------------------
inline const Matrix &HomogeneousTransformation::GetMatrix() const
{
  return _matrix;
}

// -----------------------------------------------------------------------------
inline const Matrix &HomogeneousTransformation::GetInverseMatrix() const
{
  return _inverse;
}

// =============================================================================
// Virtual Transformation "Impl" methods
// =============================================================================

// -----------------------------------------------------------------------------
inline void HomogeneousTransformation::GlobalTransform(double &x, double &y, double &z, double, double) const
{
  double a = _matrix(0, 0) * x + _matrix(0, 1) * y + _matrix(0, 2) * z + _matrix(0, 3);
  double b = _matrix(1, 0) * x + _matrix(1, 1) * y + _matrix(1, 2) * z + _matrix(1, 3);
  double c = _matrix(2, 0) * x + _matrix(2, 1) * y + _matrix(2, 2) * z + _matrix(2, 3);

  x = a;
  y = b;
  z = c;
}

// -----------------------------------------------------------------------------
inline void HomogeneousTransformation::LocalTransform(double &, double &, double &, double, double) const
{
  // No local component
}

// -----------------------------------------------------------------------------
inline void HomogeneousTransformation::Transform(double &x, double &y, double &z, double t, double t0) const
{
  this->GlobalTransform(x, y, z, t, t0);
}

// -----------------------------------------------------------------------------
inline void HomogeneousTransformation::GlobalInverse(double &x, double &y, double &z, double, double) const
{
  double a = _inverse(0, 0) * x + _inverse(0, 1) * y + _inverse(0, 2) * z + _inverse(0, 3);
  double b = _inverse(1, 0) * x + _inverse(1, 1) * y + _inverse(1, 2) * z + _inverse(1, 3);
  double c = _inverse(2, 0) * x + _inverse(2, 1) * y + _inverse(2, 2) * z + _inverse(2, 3);

  x = a;
  y = b;
  z = c;
}

// -----------------------------------------------------------------------------
inline bool HomogeneousTransformation::LocalInverse(double &, double &, double &, double, double) const
{
  // No local component
  return true;
}

// -----------------------------------------------------------------------------
inline bool HomogeneousTransformation::Inverse(double &x, double &y, double &z, double t, double t0) const
{
  this->GlobalInverse(x, y, z, t, t0);
  return true;
}

// -----------------------------------------------------------------------------
inline void HomogeneousTransformation::GlobalJacobian(Matrix &jac, double, double, double, double, double) const
{
  jac.Initialize(3, 3);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      jac(i, j) = _matrix(i, j);
    }
  }
}

// -----------------------------------------------------------------------------
inline void HomogeneousTransformation::LocalJacobian(Matrix &jac, double, double, double, double, double) const
{
  jac.Initialize(3, 3);
  jac(0, 0) = 1;
  jac(1, 1) = 1;
  jac(2, 2) = 1;
}

// -----------------------------------------------------------------------------
inline void HomogeneousTransformation::Jacobian(Matrix &jac, double x, double y, double z, double t, double t0) const
{
  this->GlobalJacobian(jac, x, y, z, t, t0);
}


} // namespace mirtk

#endif // MIRTK_HomogeneousTransformation_H
