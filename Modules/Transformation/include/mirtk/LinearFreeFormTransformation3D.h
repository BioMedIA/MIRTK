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

#ifndef MIRTK_LinearFreeFormTransformation3D_H
#define MIRTK_LinearFreeFormTransformation3D_H

#include "mirtk/FreeFormTransformation3D.h"
#include "mirtk/LinearInterpolateImageFunction.h"

#include <cmath>


namespace mirtk {


// Forward declaration
class BSplineFreeFormTransformation3D;


/**
 * Free-form transformation with linear 3D interpolation kernel
 */
class LinearFreeFormTransformation3D : public FreeFormTransformation3D
{
  mirtkTransformationMacro(LinearFreeFormTransformation3D);

  // ---------------------------------------------------------------------------
  // Types

  typedef GenericLinearInterpolateImageFunction<CPImage> Interpolator;

  // ---------------------------------------------------------------------------
  // Data members

protected:

  /// Interpolates control point values at arbitrary lattice locations
  Interpolator _FFD;

public:

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Default constructor
  LinearFreeFormTransformation3D();

  /// Construct free-form transformation for given image domain and lattice spacing
  LinearFreeFormTransformation3D(double, double, double,
                                 double, double, double,
                                 double, double, double,
                                 double *, double *, double *);

  /// Construct free-form transformation for given image domain and lattice spacing
  explicit LinearFreeFormTransformation3D(const ImageAttributes &,
                                          double = -1, double = -1, double = -1);

  /// Construct free-form transformation for given target image and lattice spacing
  explicit LinearFreeFormTransformation3D(const BaseImage &,
                                          double, double, double);

  /// Construct free-form transformation from existing 3D+t displacement field
  explicit LinearFreeFormTransformation3D(const GenericImage<double> &, bool = true);

  /// Construct free-form transformation from B-spline FFD
  explicit LinearFreeFormTransformation3D(const BSplineFreeFormTransformation3D &);

  /// Copy Constructor
  LinearFreeFormTransformation3D(const LinearFreeFormTransformation3D &);

  /// Destructor
  virtual ~LinearFreeFormTransformation3D();

  // ---------------------------------------------------------------------------
  // Approximation/Interpolation

  /// Approximate displacements: This function takes a set of points and a set
  /// of displacements and finds !new! parameters such that the resulting
  /// transformation approximates the displacements as good as possible.
  virtual void ApproximateDOFs(const double *, const double *, const double *, const double *,
                               const double *, const double *, const double *, int);

  /// Finds gradient of approximation error: This function takes a set of points
  /// and a set of errors. It finds a gradient w.r.t. the transformation parameters
  /// which minimizes the L2 norm of the approximation error and adds it to the
  /// input gradient with the given weight.
  virtual void ApproximateDOFsGradient(const double *, const double *, const double *, const double *,
                                       const double *, const double *, const double *, int,
                                       double *, double = 1.0) const;

  /// Interpolates displacements: This function takes a set of displacements defined
  /// at the control points and finds a FFD which interpolates these displacements.
  virtual void Interpolate(const double *, const double *, const double *);

  // ---------------------------------------------------------------------------
  // Lattice
  using FreeFormTransformation::BoundingBox;

  /// Returns the bounding box for a control point (in mm). The last
  /// parameter specifies what fraction of the bounding box to return.
  /// The default is 1 which equals 100% of the bounding box.
  virtual void BoundingBox(int, double &, double &, double &,
                                double &, double &, double &, double = 1) const;

  // ---------------------------------------------------------------------------
  // Evaluation

  /// Evaluates the free-form transformation at a lattice point
  void Evaluate(double &, double &, double &, int, int, int) const;

  /// Evaluates the free-form transformation at a point in lattice coordinates
  void Evaluate(double &, double &, double &) const;

  /// Calculates the Jacobian of the local transformation w.r.t world coordinates
  void EvaluateJacobian(Matrix &, double, double, double) const;

  // ---------------------------------------------------------------------------
  // Point transformation

  /// Transforms a single point using the local transformation only
  virtual void LocalTransform(double &, double &, double &, double = 0, double = NaN) const;

  /// Transforms a single point using the inverse of the local transformation only
  virtual bool LocalInverse(double &, double &, double &, double = 0, double = NaN) const;

  // ---------------------------------------------------------------------------
  // Derivatives
  using FreeFormTransformation3D::LocalJacobian;
  using FreeFormTransformation3D::JacobianDOFs;

  /// Calculates the Jacobian of the local transformation w.r.t world coordinates
  virtual void LocalJacobian(Matrix &, double, double, double, double = 0, double = NaN) const;

  /// Calculates the Jacobian of the transformation w.r.t the transformation parameters
  virtual void JacobianDOFs(double [3], int, int, int, double, double, double) const;

  // ---------------------------------------------------------------------------
  // Properties
  using FreeFormTransformation3D::BendingEnergy;

  /// Size of support region of the used kernel
  virtual int KernelSize() const;

protected:

  /// Calculate the bending energy of the transformation at control points
  double BendingEnergy(int i, int j, int k) const;

public:

  /// Approximates the bending energy on the control point lattice
  virtual double BendingEnergy(bool = false, bool = true) const;

  /// Approximates the gradient of the bending energy on the control point
  /// lattice w.r.t the transformation parameters and adds it with the given weight
  virtual void BendingEnergyGradient(double *, double = 1, bool = false, bool = true, bool = true) const;

  // ---------------------------------------------------------------------------
  // I/O

  // Do not hide methods of base class
  using FreeFormTransformation3D::Print;

  /// Prints the parameters of the transformation
  virtual void Print(ostream &, Indent = 0) const;

  /// Whether this transformation can read a file of specified type (i.e. format)
  virtual bool CanRead(TransformationType) const;

  // ---------------------------------------------------------------------------
  // Others

  /// Compose this transformation (T1) with second transformation (T2).
  /// The result is defined as T = T1 o T2.
  virtual void Compose(const Transformation *);

};

// Backwards compatibility
typedef LinearFreeFormTransformation3D LinearFreeFormTransformation;

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
inline void LinearFreeFormTransformation3D
::Evaluate(double &x, double &y, double &z, int i, int j, int k) const
{
  Vector d = _CPValue->Get(i, j, k);
  x = d._x, y = d._y, z = d._z;
}

// -----------------------------------------------------------------------------
inline void LinearFreeFormTransformation3D
::Evaluate(double &x, double &y, double &z) const
{
  Vector d = _FFD(x, y, z);
  x = d._x, y = d._y, z = d._z;
}

// =============================================================================
// Point transformation
// =============================================================================

// -----------------------------------------------------------------------------
inline void LinearFreeFormTransformation3D
::LocalTransform(double &x, double &y, double &z, double, double) const
{
  // Convert to lattice coordinates
  double dx = x, dy = y, dz = z;
  this->WorldToLattice(dx, dy, dz);
  // Evaluate displacement
  this->Evaluate(dx, dy, dz);
  // Transform point
  x += dx, y += dy, z += dz;
}

// -----------------------------------------------------------------------------
inline bool LinearFreeFormTransformation3D
::LocalInverse(double &x, double &y, double &z, double, double) const
{
  // Convert to lattice coordinates
  double dx = x, dy = y, dz = z;
  this->WorldToLattice(dx, dy, dz);
  // Evaluate displacement
  this->Evaluate(dx, dy, dz);
  // Transform point
  x -= dx, y -= dy, z -= dz;
  return true;
}

// =============================================================================
// Derivatives
// =============================================================================

// -----------------------------------------------------------------------------
inline void LinearFreeFormTransformation3D
::LocalJacobian(Matrix &jac, double x, double y, double z, double, double) const
{
  // Convert to lattice coordinates
  this->WorldToLattice(x, y, z);
  // Compute 1st order derivatives
  this->EvaluateJacobian(jac, x, y, z);
  // Convert derivatives to world coordinates
  JacobianToWorld(jac);
  // Add derivatives of "x" term in T(x) = x + this->Evaluate(x)
  jac(0, 0) += 1.0;
  jac(1, 1) += 1.0;
  jac(2, 2) += 1.0;
}

// -----------------------------------------------------------------------------
inline void LinearFreeFormTransformation3D
::JacobianDOFs(double jac[3], int ci, int cj, int ck, double x, double y, double z) const
{
  // Convert point to lattice coordinates
  this->WorldToLattice(x, y, z);
  // Calculate tensor product of 1D linear interpolation kernels
  const double dx = abs(ci - x);
  const double dy = abs(cj - y);
  const double dz = abs(ck - z);
  if (dx < 1.0 && dy < 1.0 && dz < 1.0) {
    jac[0] = jac[1] = jac[2] = (1.0 - dx) * (1.0 - dy) * (1.0 - dz);
  } else {
    jac[0] = jac[1] = jac[3] = .0;
  }
}


} // namespace mirtk

#endif // MIRTK_LinearFreeFormTransformation3D_H
