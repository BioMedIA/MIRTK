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

#ifndef MIRTK_BSplineFreeFormTransformation4D_H
#define MIRTK_BSplineFreeFormTransformation4D_H

#include "mirtk/FreeFormTransformation4D.h"

#include "mirtk/FastCubicBSplineInterpolateImageFunction4D.h"
#include "mirtk/TransformationJacobian.h"


namespace mirtk {


/**
 * Class for free-form transformations based on tensor product B-splines.
 *
 * This class implements 4D free-form transformation using B-splines.
 *
 * For more details about the implementation see Lee, Wolberg and Shin, IEEE
 * Transactions on Visualization and Computer Graphics, Vol. 3, No. 3, 1997.
 */
class BSplineFreeFormTransformation4D : public FreeFormTransformation4D
{
  mirtkTransformationMacro(BSplineFreeFormTransformation4D);

  // ---------------------------------------------------------------------------
  // Types

public:

  typedef GenericFastCubicBSplineInterpolateImageFunction4D<CPImage> Interpolator;
  typedef Interpolator::Kernel                                       Kernel;

  // -------------------------------------------------------------------------
  // Data members

protected:

  /// Interpolates control point values at arbitrary lattice locations
  Interpolator _FFD;

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Default constructor
  BSplineFreeFormTransformation4D();

  /// Construct free-form transformation for given image domain and lattice spacing
  BSplineFreeFormTransformation4D(double, double, double, double,
                                  double, double, double, double,
                                  double, double, double, double,
                                  double *, double *, double *);

  /// Construct free-form transformation for given image domain and lattice spacing
  explicit BSplineFreeFormTransformation4D(const ImageAttributes &,
                                           double = -1, double = -1, double = -1, double = -1);

  /// Construct free-form transformation for given target image and lattice spacing
  explicit BSplineFreeFormTransformation4D(const BaseImage &,
                                           double, double, double, double);

  /// Copy Constructor
  BSplineFreeFormTransformation4D(const BSplineFreeFormTransformation4D &);

  /// Destructor
  virtual ~BSplineFreeFormTransformation4D();

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

  using FreeFormTransformation4D::BoundingBox;

  /// Number of control points in x after subdivision
  virtual int GetXAfterSubdivision() const;

  /// Number of control points in y after subdivision
  virtual int GetYAfterSubdivision() const;

  /// Number of control points in z after subdivision
  virtual int GetZAfterSubdivision() const;

  /// Number of control points in t after subdivision
  virtual int GetTAfterSubdivision() const;

  /// Returns the control point spacing in x after the subdivision
  virtual double GetXSpacingAfterSubdivision() const;

  /// Returns the control point spacing in y after the subdivision
  virtual double GetYSpacingAfterSubdivision() const;

  /// Returns the control point spacing in z after the subdivision
  virtual double GetZSpacingAfterSubdivision() const;

  /// Returns the control point spacing in t after the subdivision
  virtual double GetTSpacingAfterSubdivision() const;

  /// Subdivide FFD lattice in the specified dimensions
  virtual void Subdivide(bool = true, bool = true, bool = true, bool = true);

  /// Gets the temporal bounding box for a control point. The last parameter
  /// specifies what fraction of the bounding box to return. The default
  /// is 1 which equals 100% of the bounding box.
  virtual void BoundingBox(int, double &, double &, double = 1) const;

  /// Gets the spatial bounding box for a control point. The last parameter
  /// specifies what fraction of the bounding box to return. The default
  /// is 1 which equals 100% of the bounding box.
  virtual void BoundingBox(int, double &, double &, double &,
                                double &, double &, double &, double = 1) const;

  /// Gets the spatio-temporal bounding box for a control point. The last parameter
  /// specifies what fraction of the bounding box to return. The default
  /// is 1 which equals 100% of the bounding box.
  virtual void BoundingBox(int, double &, double &, double &, double &,
                                double &, double &, double &, double &, double = 1) const;

  // ---------------------------------------------------------------------------
  // Evaluation

  /// Calculates the FFD (for a point in FFD coordinates)
  void Evaluate(double &, double &, double &, double) const;

  /// Evaluates the FFD at a point in lattice coordinates inside the FFD domain
  void EvaluateInside(double &, double &, double &, double) const;

  /// Calculates the spatial Jacobian of the FFD at a point in lattice coordinates
  void EvaluateJacobian(Matrix &, int, int, int, int) const;

  /// Calculates the spatial Jacobian of the FFD at a point in lattice coordinates
  void EvaluateJacobian(Matrix &, double, double, double, double) const;

  /// Calculates the spatial Jacobian of the FFD at a point in lattice coordinates
  /// and converts the resulting Jacobian to derivatives w.r.t world coordinates
  void EvaluateJacobianWorld(Matrix &, double, double, double, double) const;

  /// Calculates the Jacobian of the FFD at a point in lattice coordinates
  /// w.r.t the control point with lattice coordinates (i, j, k, l)
  void EvaluateJacobianDOFs(double [3], int,    int,    int,    int,
                                        double, double, double, double) const;

  /// Calculates the Jacobian of the FFD at a point in lattice coordinates
  /// w.r.t the parameters. The resulting Jacobian matrix contains one column
  /// for each parameter, but only non-zero columns, i.e., of control points
  /// for whom the point is in the local support region, are actually stored.
  void EvaluateJacobianDOFs(TransformationJacobian &jac, double, double, double, double) const;

  /// Calculates the Hessian of the FFD at a point in lattice coordinates
  void EvaluateHessian(Matrix [3], int, int, int, int) const;

  /// Calculates the Hessian of the FFD at a point in lattice coordinates
  void EvaluateHessian(Matrix [3], double, double, double, double) const;

  /// Calculates the Hessian of the FFD at a point in lattice coordinates
  /// and converts the resulting Hessian to derivatives w.r.t world coordinates
  void EvaluateHessianWorld(Matrix [3], double, double, double, double) const;

  /// Calculates the Laplacian of the FFD at a point in lattice coordinates
  void EvaluateLaplacian(double [3], int, int, int, int) const;

  /// Calculates the Laplacian of the FFD at a point in lattice coordinates
  void EvaluateLaplacian(double [3], int, int, int, double) const;

  /// Calculates the Laplacian of the FFD at a point in lattice coordinates
  void EvaluateLaplacian(double [3], double, double, double, double) const;

  /// Calculates the Laplacian of the FFD at a point in lattice coordinates
  void EvaluateLaplacian(double &, double &, double &, double) const;

  // ---------------------------------------------------------------------------
  // Point transformation

  /// Transforms a single point using the local transformation component only
  virtual void LocalTransform(double &, double &, double &, double, double = NaN) const;

  // ---------------------------------------------------------------------------
  // Derivatives

  using FreeFormTransformation4D::LocalJacobian;
  using FreeFormTransformation4D::JacobianDOFs;
  using FreeFormTransformation4D::ParametricGradient;

  /// Calculates the Jacobian of the local transformation w.r.t world coordinates
  virtual void LocalJacobian(Matrix &, double, double, double, double, double = NaN) const;

  /// Calculates the Hessian for each component of the local transformation w.r.t world coordinates
  virtual void LocalHessian(Matrix [3], double, double, double, double, double = NaN) const;

  /// Calculates the Jacobian of the transformation w.r.t the transformation parameters
  virtual void JacobianDOFs(Matrix &, int, int, int, int, double, double, double, double, double = NaN) const;

  /// Calculates the Jacobian of the transformation w.r.t the transformation parameters
  virtual void JacobianDOFs(double [3], int, int, int, int, double, double, double, double, double = NaN) const;

  /// Calculates the Jacobian of the transformation w.r.t the transformation parameters
  virtual void JacobianDOFs(TransformationJacobian &, double, double, double, double, double = NaN) const;

  /// Applies the chain rule to convert spatial non-parametric gradient
  /// to a gradient w.r.t the parameters of this transformation
  virtual void ParametricGradient(const GenericImage<double> *, double *,
                                  const WorldCoordsImage *,
                                  const WorldCoordsImage *,
                                  double = NaN, double = 1.0) const;

  // ---------------------------------------------------------------------------
  // Properties
  using FreeFormTransformation4D::BendingEnergy;

  /// Size of support region of the used kernel
  virtual int KernelSize() const;

  /// Calculates the bending of the transformation
  virtual double BendingEnergy(double, double, double, double = 0, double = NaN, bool = true) const;

  /// Approximates the bending energy on the control point lattice
  virtual double BendingEnergy(bool = false, bool = true) const;

  /// Approximates the bending energy on the specified discrete domain
  virtual double BendingEnergy(const ImageAttributes &, double = NaN, bool = true) const;

  /// Approximates the gradient of the bending energy on the control point
  /// lattice w.r.t the transformation parameters and adds it with the given weight
  virtual void BendingEnergyGradient(double *, double = 1.0, bool = false, bool = true, bool = true) const;

  // ---------------------------------------------------------------------------
  // I/O

  // Do not hide methods of base class
  using FreeFormTransformation4D::Print;

  /// Prints the parameters of the transformation
  virtual void Print(ostream &, Indent = 0) const;

  /// Whether this transformation can read a file of specified type (i.e. format)
  virtual bool CanRead(TransformationType) const;

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
inline void BSplineFreeFormTransformation4D
::Evaluate(double &x, double &y, double &z, double t) const
{
  Vector d = _FFD(x, y, z, t);
  x = d._x, y = d._y, z = d._z;
}

// -----------------------------------------------------------------------------
inline void BSplineFreeFormTransformation4D
::EvaluateInside(double &x, double &y, double &z, double t) const
{
  Vector d = _FFD.GetInside(x, y, z, t);
  x = d._x, y = d._y, z = d._z;
}

// -----------------------------------------------------------------------------
inline void BSplineFreeFormTransformation4D
::EvaluateJacobianWorld(Matrix &jac, double x, double y, double z, double t) const
{
  this->EvaluateJacobian(jac, x, y, z, t);
  JacobianToWorld(jac);
}

// -----------------------------------------------------------------------------
inline void BSplineFreeFormTransformation4D
::EvaluateJacobianDOFs(double jac[3], int    i, int    j, int    k, int    l,
                                      double x, double y, double z, double t) const
{
  jac[0] = jac[1] = jac[2] = Kernel::Weight(x - i) *
                             Kernel::Weight(y - j) *
                             Kernel::Weight(z - k) *
                             Kernel::Weight(t - l);
}

// -----------------------------------------------------------------------------
inline void BSplineFreeFormTransformation4D
::EvaluateJacobianDOFs(TransformationJacobian &jac, double x, double y, double z, double t) const
{
  int i = static_cast<int>(floor(x));
  int j = static_cast<int>(floor(y));
  int k = static_cast<int>(floor(z));
  int l = static_cast<int>(floor(t));

  const int A = Kernel::VariableToIndex(x - i);
  const int B = Kernel::VariableToIndex(y - j);
  const int C = Kernel::VariableToIndex(z - k);
  const int D = Kernel::VariableToIndex(t - l);

  --i, --j, --k, --l;

  double basis[4];
  int    ci, cj, ck, cl, xdof, ydof, zdof;

  for (int d = 0; d <= 3; ++d) {
    cl = l + d;
    if (0 <= cl && cl < _t) {
      basis[3] = Kernel::LookupTable[D][d];
      for (int c = 0; c <= 3; ++c) {
        ck = k + c;
        if (0 <= ck && ck < _z) {
          basis[2] = Kernel::LookupTable[C][c] * basis[3];
          for (int b = 0; b <= 3; ++b) {
            cj = j + b;
            if (0 <= cj && cj < _y) {
              basis[1] = Kernel::LookupTable[B][b] * basis[2];
              for (int a = 0; a <= 3; ++a) {
                ci = i + a;
                if (0 <= ci && ci < _x) {
                  basis[0] = Kernel::LookupTable[A][a] * basis[1];
                  IndexToDOFs(LatticeToIndex(ci, cj, ck, cl), xdof, ydof, zdof);
                  jac(xdof)._x = basis[0];
                  jac(ydof)._y = basis[0];
                  jac(zdof)._z = basis[0];
                }
              }
            }
          }
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------
inline void BSplineFreeFormTransformation4D
::EvaluateHessianWorld(Matrix hessian[3], double x, double y, double z, double t) const
{
  this->EvaluateHessian(hessian, x, y, z, t);
  HessianToWorld(hessian);
}

// =============================================================================
// Point transformation
// =============================================================================

// -----------------------------------------------------------------------------
inline void BSplineFreeFormTransformation4D
::LocalTransform(double &x, double &y, double &z, double t, double) const
{
  // Convert to lattice coordinates
  double dx = x, dy = y, dz = z;
  this->WorldToLattice(dx, dy, dz);
  t = this->TimeToLattice(t);
  // Evaluate displacement
  this->Evaluate(dx, dy, dz, t);
  // Transform point
  x += dx, y += dy, z += dz;
}

// =============================================================================
// Derivatives
// =============================================================================

// -----------------------------------------------------------------------------
inline void BSplineFreeFormTransformation4D
::LocalJacobian(Matrix &jac, double x, double y, double z, double t, double) const
{
  // Convert to lattice coordinates
  this->WorldToLattice(x, y, z);
  t = this->TimeToLattice(t);
  // Compute 1st order derivatives
  this->EvaluateJacobianWorld(jac, x, y, z, t);
  // Add derivatives of "x" term in T(x) = x + FFD(x)
  jac(0, 0) += 1.0;
  jac(1, 1) += 1.0;
  jac(2, 2) += 1.0;
}

// -----------------------------------------------------------------------------
inline void BSplineFreeFormTransformation4D
::LocalHessian(Matrix hessian[3], double x, double y, double z, double t, double) const
{
  // Convert to lattice coordinates
  this->WorldToLattice(x, y, z);
  t = this->TimeToLattice(t);
  // Compute 2nd order derivatives
  this->EvaluateHessianWorld(hessian, x, y, z, t);
}

// -----------------------------------------------------------------------------
inline void BSplineFreeFormTransformation4D
::JacobianDOFs(double jac[3], int    i, int    j, int    k, int    l,
                              double x, double y, double z, double t, double) const
{
  // Convert to lattice coordinates
  this->WorldToLattice(x, y, z);
  t = this->TimeToLattice(t);
  // Calculate derivatives w.r.t. control point coordinates
  this->EvaluateJacobianDOFs(jac, i, j, k, l, x, y, z, t);
}

// -----------------------------------------------------------------------------
inline void BSplineFreeFormTransformation4D
::JacobianDOFs(Matrix &jac, int    i, int    j, int    k, int    l,
                                double x, double y, double z, double t, double t0) const
{
  double tmp[3];
  this->JacobianDOFs(tmp, i, j, k, l, x, y, z, t, t0);
  jac.Initialize(3, 3);
  jac(0, 0) = tmp[0];
  jac(1, 1) = tmp[1];
  jac(2, 2) = tmp[2];
}

// -----------------------------------------------------------------------------
inline void BSplineFreeFormTransformation4D
::JacobianDOFs(TransformationJacobian &jac, double x, double y, double z, double t, double) const
{
  // Convert to lattice coordinates
  this->WorldToLattice(x, y, z);
  t = this->TimeToLattice(t);
  // Calculate derivatives w.r.t. control point coordinates
  this->EvaluateJacobianDOFs(jac, x, y, z, t);
}


} // namespace mirtk

#endif // MIRTK_BSplineFreeFormTransformation4D_H
