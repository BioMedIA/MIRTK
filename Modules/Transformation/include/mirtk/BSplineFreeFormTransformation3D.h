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

#ifndef MIRTK_BSplineFreeFormTransformation3D_H
#define MIRTK_BSplineFreeFormTransformation3D_H

#include "mirtk/FreeFormTransformation3D.h"

#include "mirtk/FastCubicBSplineInterpolateImageFunction2D.h"
#include "mirtk/FastCubicBSplineInterpolateImageFunction3D.h"


namespace mirtk {


/**
 * Class for free-form transformations based on tensor product B-splines.
 *
 * This class implements 3D free-form transformation using B-splines.
 *
 * For more details about the implementation see Lee, Wolberg and Shin, IEEE
 * Transactions on Visualization and Computer Graphics, Vol. 3, No. 3, 1997.
 */
class BSplineFreeFormTransformation3D : public FreeFormTransformation3D
{
  mirtkTransformationMacro(BSplineFreeFormTransformation3D);

  // ---------------------------------------------------------------------------
  // Types

  /// \todo Remove once BSplineFreeFormTransformation2D is implemented.
  typedef GenericFastCubicBSplineInterpolateImageFunction2D<CPImage> Interpolator2D;

public:

  typedef GenericFastCubicBSplineInterpolateImageFunction3D<CPImage> Interpolator;
  typedef Interpolator::Kernel                                       Kernel;

  /// Options for parametric gradient calculation
  enum ParametricGradientType
  {
    PG_Default,       ///< Default gradient computation.
    PG_Analytic,      ///< Analytic derivation w.r.t. DoFs.
    PG_Convolution,   ///< Convolution with cubic B-spline filter.
    PG_Approximation  ///< Approximate voxel-wise non-parametric gradient
                      ///< field with cubic B-spline function, also known
                      ///< as "directly manipulated FFD (DMFFD)".
  };

  // ---------------------------------------------------------------------------
  // Attributes

protected:

  /// Whether to compute parametric gradient using cubic B-spline convolution
  mirtkPublicAttributeMacro(ParametricGradientType, ParametricGradientCalculation);

  /// Interpolates control point values at arbitrary lattice locations
  Interpolator _FFD;

  /// Interpolates control point values at arbitrary 2D lattice locations
  /// \todo Remove once BSplineFreeFormTransformation2D is implemented.
  Interpolator2D _FFD2D;

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Default constructor
  BSplineFreeFormTransformation3D();

  /// Construct free-form transformation for given image domain and lattice spacing
  BSplineFreeFormTransformation3D(double, double, double,
                                  double, double, double,
                                  double, double, double,
                                  double *, double *, double *);

  /// Construct free-form transformation for given image domain and lattice spacing
  explicit BSplineFreeFormTransformation3D(const ImageAttributes &,
                                           double = -1, double = -1, double = -1);

  /// Construct free-form transformation for given target image and lattice spacing
  explicit BSplineFreeFormTransformation3D(const BaseImage &,
                                           double, double, double);

  /// Construct free-form transformation from existing 3D+t deformation field
  explicit BSplineFreeFormTransformation3D(const GenericImage<double> &, bool = false);

  /// Copy Constructor
  BSplineFreeFormTransformation3D(const BSplineFreeFormTransformation3D &);

  /// Destructor
  virtual ~BSplineFreeFormTransformation3D();

  // ---------------------------------------------------------------------------
  // Parameters (non-DoFs)

  // Import other Parameter overloads
  using FreeFormTransformation3D::Parameter;

  /// Set named (non-DoF) parameter from value as string
  virtual bool Set(const char *, const char *);

  /// Get (non-DoF) parameters as key/value as string map
  virtual ParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Approximation/Interpolation

  /// Approximate displacements: This function takes a set of points and a set
  /// of displacements and finds !new! parameters such that the resulting
  /// transformation approximates the displacements as good as possible.
  /// These parameters are added to the given coefficients using specified weight.
  void AddApproximateSplineCoefficients(const double *, const double *, const double *,
                                        const double *, const double *, const double *,
                                        int, double *, double = 1., bool = false) const;

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
  virtual void Interpolate(const double *, const double *, const double * = NULL);

  // ---------------------------------------------------------------------------
  // Lattice

  using FreeFormTransformation3D::BoundingBox;

  /// Number of control points in x after subdivision
  virtual int GetXAfterSubdivision() const;

  /// Number of control points in y after subdivision
  virtual int GetYAfterSubdivision() const;

  /// Number of control points in z after subdivision
  virtual int GetZAfterSubdivision() const;

  /// Number of control points in t after subdivision
  virtual int GetTAfterSubdivision() const;

  /// Control point spacing in x after subdivision
  virtual double GetXSpacingAfterSubdivision() const;

  /// Control point spacing in y after subdivision
  virtual double GetYSpacingAfterSubdivision() const;

  /// Control point spacing in z after subdivision
  virtual double GetZSpacingAfterSubdivision() const;

  /// Control point spacing in t after subdivision
  virtual double GetTSpacingAfterSubdivision() const;

  /// Subdivide FFD lattice
  virtual void Subdivide(bool = true, bool = true, bool = true, bool = true);

  /// Returns the bounding box for a control point (in mm). The last
  /// parameter specifies what fraction of the bounding box to return.
  /// The default is 1 which equals 100% of the bounding box.
  virtual void BoundingBox(int, double &, double &, double &,
                                double &, double &, double &, double = 1) const;

  // ---------------------------------------------------------------------------
  // Evaluation

  /// Evaluates the FFD at a lattice point
  void Evaluate(double &, double &, double &, int, int) const;

  /// Evaluates the FFD at a lattice point
  void Evaluate(double &, double &, double &, int, int, int) const;

  /// Evaluates the FFD at a point in lattice coordinates
  void Evaluate(double &, double &, double &) const;

  /// Evaluates the FFD at a point in lattice coordinates inside the FFD domain
  void EvaluateInside(double &, double &, double &) const;

  /// Calculates the Jacobian of the FFD at a point in lattice coordinates
  void EvaluateJacobian(Matrix &, int, int) const;

  /// Calculates the Jacobian of the FFD at a point in lattice coordinates
  void EvaluateJacobian(Matrix &, int, int, int) const;

  /// Calculates the Jacobian of the FFD at a point in lattice coordinates
  void EvaluateJacobian(Matrix &, double, double) const;

  /// Calculates the Jacobian of the FFD at a point in lattice coordinates
  void EvaluateJacobian(Matrix &, double, double, double) const;

  /// Calculates the Jacobian of the FFD at a point in lattice coordinates
  /// and converts the resulting Jacobian to derivatives w.r.t world coordinates
  void EvaluateJacobianWorld(Matrix &, double, double) const;

  /// Calculates the Jacobian of the FFD at a point in lattice coordinates
  /// and converts the resulting Jacobian to derivatives w.r.t world coordinates
  void EvaluateJacobianWorld(Matrix &, double, double, double) const;

  /// Calculates the Jacobian of the FFD at a point in lattice coordinates
  /// w.r.t the control point with lattice coordinates (i, j)
  void EvaluateJacobianDOFs(double [3], int, int, double, double) const;

  /// Calculates the Jacobian of the FFD at a point in lattice coordinates
  /// w.r.t the control point with lattice coordinates (i, j, k)
  void EvaluateJacobianDOFs(double [3], int,    int,    int,
                                        double, double, double) const;

  /// Calculates the Jacobian of the FFD at a point in lattice coordinates
  /// w.r.t the control point with lattice coordinates (i, j, k)
  ///
  /// \note The temporal coordinates are required by the templated integration
  ///       methods. These arguments are ignored by this function.
  void EvaluateJacobianDOFs(double [3], int,    int,    int,    int,
                                        double, double, double, double) const;

  /// Calculates the derivative of the Jacobian of the FFD at a point in lattice coordinates
  /// w.r.t a transformation parameter
  void EvaluateDerivativeOfJacobianWrtDOF(Matrix &, int, double, double) const;

  /// Calculates the derivative of the Jacobian of the FFD at a point in lattice coordinates
  /// w.r.t a transformation parameter
  void EvaluateDerivativeOfJacobianWrtDOF(Matrix &, int, double, double, double) const;

  /// Calculates the Hessian of the 2D FFD at a point in lattice coordinates
  void EvaluateHessian(Matrix [3], int, int) const;

  /// Calculates the Hessian of the 2D FFD at a point in lattice coordinates
  void EvaluateHessian(Matrix [3], double, double) const;

  /// Calculates the Hessian of the FFD at a point in lattice coordinates
  void EvaluateHessian(Matrix [3], int, int, int) const;

  /// Calculates the Hessian of the FFD at a point in lattice coordinates
  void EvaluateHessian(Matrix [3], double, double, double) const;

  /// Calculates the Laplacian of the FFD at a point in lattice coordinates
  void EvaluateLaplacian(double [3], int, int, int) const;

  /// Calculates the Laplacian of the FFD at a point in lattice coordinates
  void EvaluateLaplacian(double [3], double, double, double) const;

  /// Calculates the Laplacian of the FFD at a point in lattice coordinates
  void EvaluateLaplacian(double &, double &, double &) const;

  /// Calculate derivatives of Jacobian determinant w.r.t. DoFs of control point
  ///
  /// \param[out] dJ           Partial derivatives of Jacobian determinant w.r.t. DoFs of control point.
  /// \param[in]  adj          Adjugate of Jacobian matrix evaluated at (x, y, z).
  /// \param[in]  a            Distance from control point along x axis of lattice in lattice units.
  /// \param[in]  b            Distance from control point along y axis of lattice in lattice units.
  /// \param[in]  c            Distance from control point along z axis of lattice in lattice units.
  /// \param[in]  wrt_world    Whether derivatives are computed w.r.t. world coordinate system.
  /// \param[in]  use_spacing  Whether to use grid spacing when \p wrt_world is \c true.
  void EvaluateJacobianDetDerivative(double dJ[3], const Matrix &adj, double a, double b, double c,
                                     bool wrt_world = true, bool use_spacing = true) const;

  /// Calculate derivatives of Jacobian determinant w.r.t. DoFs of control point
  ///
  /// \param[out] dJ           Partial derivatives of Jacobian determinant w.r.t. DoFs of control point.
  /// \param[in]  adj          Adjugate of Jacobian matrix evaluated at (x, y, z).
  /// \param[in]  a            Distance from control point along x axis of lattice in lattice units.
  /// \param[in]  b            Distance from control point along y axis of lattice in lattice units.
  /// \param[in]  c            Distance from control point along z axis of lattice in lattice units.
  /// \param[in]  wrt_world    Whether derivatives are computed w.r.t. world coordinate system.
  /// \param[in]  use_spacing  Whether to use grid spacing when \p wrt_world is \c true.
  void EvaluateJacobianDetDerivative(double dJ[3], const Matrix &adj, int a, int b, int c,
                                     bool wrt_world = true, bool use_spacing = true) const;

  /// Calculate derivatives of Jacobian determinant w.r.t. DoFs of control point
  ///
  /// \param[out] dJ           Partial derivatives of Jacobian determinant w.r.t. DoFs of control point.
  /// \param[in]  adj          Adjugate of Jacobian matrix evaluated at (x, y, z).
  /// \param[in]  i            Index of control point along x axis of lattice.
  /// \param[in]  j            Index of control point along y axis of lattice.
  /// \param[in]  k            Index of control point along z axis of lattice.
  /// \param[in]  x            Point coordinate along x axis of lattice in lattice units.
  /// \param[in]  y            Point coordinate along y axis of lattice in lattice units.
  /// \param[in]  z            Point coordinate along z axis of lattice in lattice units.
  /// \param[in]  wrt_world    Whether derivatives are computed w.r.t. world coordinate system.
  /// \param[in]  use_spacing  Whether to use grid spacing when \p wrt_world is \c true.
  void EvaluateJacobianDetDerivative(double dJ[3], const Matrix &adj,
                                     int i, int j, int k,
                                     double x, double y, double z,
                                     bool wrt_world = true, bool use_spacing = true) const;

  /// Calculate derivatives of Jacobian determinant w.r.t. DoFs of control point
  ///
  /// \param[out] dJ           Partial derivatives of Jacobian determinant w.r.t. DoFs of control point.
  /// \param[in]  adj          Adjugate of Jacobian matrix evaluated at (x, y, z).
  /// \param[in]  cp           Linear index of control point.
  /// \param[in]  x            Point coordinate along x axis of lattice in lattice units.
  /// \param[in]  y            Point coordinate along y axis of lattice in lattice units.
  /// \param[in]  z            Point coordinate along z axis of lattice in lattice units.
  /// \param[in]  wrt_world    Whether derivatives are computed w.r.t. world coordinate system.
  /// \param[in]  use_spacing  Whether to use grid spacing when \p wrt_world is \c true.
  void EvaluateJacobianDetDerivative(double dJ[3], const Matrix &adj, int cp,
                                     double x, double y, double z,
                                     bool wrt_world = true, bool use_spacing = true) const;

  // ---------------------------------------------------------------------------
  // Point transformation

  /// Transforms a single point using the local transformation component only
  virtual void LocalTransform(double &, double &, double &, double = 0, double = NaN) const;

  /// Whether this transformation implements a more efficient update of a given
  /// displacement field given the desired change of a transformation parameter
  virtual bool CanModifyDisplacement(int = -1) const;

  /// Updates the displacement vectors for a whole image domain
  ///
  /// \param[in]     dof Transformation parameter.
  /// \param[in]     dv  Change of transformation parameter value.
  /// \param[in,out] dx  Displacement field to be updated.
  /// \param[in]     t   Time point of start point.
  /// \param[in]     t0  Time point of end point.
  /// \param[in]     i2w Pre-computed world coordinates.
  virtual void DisplacementAfterDOFChange(int dof, double dv,
                                          GenericImage<double> &dx,
                                          double t, double t0 = -1,
                                          const WorldCoordsImage *i2w = NULL) const;

  // ---------------------------------------------------------------------------
  // Derivatives

  using FreeFormTransformation3D::LocalJacobian;
  using FreeFormTransformation3D::LocalHessian;
  using FreeFormTransformation3D::JacobianDOFs;
  using FreeFormTransformation3D::ParametricGradient;

  /// Calculates the Jacobian of the transformation w.r.t either control point displacements or velocities
  virtual void FFDJacobianWorld(Matrix &, double, double, double, double = 0, double = NaN) const;

  /// Calculates the Jacobian of the local transformation w.r.t world coordinates
  virtual void LocalJacobian(Matrix &, double, double, double, double = 0, double = NaN) const;

  /// Calculates the Hessian for each component of the local transformation w.r.t world coordinates
  virtual void LocalHessian(Matrix [3], double, double, double, double = 0, double = NaN) const;

  /// Calculates the Jacobian of the transformation w.r.t. the parameters of a control point
  virtual void JacobianDOFs(double [3], int, int, int, double, double, double) const;

  /// Calculates derivatives of the Jacobian determinant of spline function w.r.t. DoFs of a control point
  ///
  /// This function is identical to JacobianDetDerivative when the DoFs of the control points are displacements.
  /// When the DoFs are velocities, however, this function computes the derivatives of the Jacobian determinant
  /// of the velocity field instead.
  ///
  /// \param[out] dJ  Partial derivatives of Jacobian determinant at (x, y, z) w.r.t. DoFs of control point.
  /// \param[in]  cp  Index of control point w.r.t. whose DoFs the derivatives are computed.
  /// \param[in]  x   World coordinate along x axis at which to evaluate derivatives.
  /// \param[in]  y   World coordinate along y axis at which to evaluate derivatives.
  /// \param[in]  z   World coordinate along z axis at which to evaluate derivatives.
  /// \param[in]  adj Adjugate of Jacobian matrix evaluated at (x, y, z).
  /// \param[in]  wrt_world    Whether derivatives are computed w.r.t. world coordinate system.
  /// \param[in]  use_spacing  Whether to use grid spacing when \p wrt_world is \c true.
  virtual void FFDJacobianDetDerivative(double dJ[3], const Matrix &adj,
                                        int cp, double x, double y, double z, double = 0, double = NaN,
                                        bool wrt_world = true, bool use_spacing = true) const;

  /// Calculates the derivative of the Jacobian of the transformation (w.r.t. world coordinates) w.r.t. a transformation parameter
  virtual void DeriveJacobianWrtDOF(Matrix &, int, double, double, double, double = 0, double = NaN) const;

    /// Applies the chain rule to convert spatial non-parametric gradient
  /// to a gradient w.r.t the parameters of this transformation
  virtual void ParametricGradient(const GenericImage<double> *, double *,
                                  const WorldCoordsImage *,
                                  const WorldCoordsImage *,
                                  double = NaN, double = 1) const;

  // ---------------------------------------------------------------------------
  // Properties

  /// Size of support region of the used kernel
  virtual int KernelSize() const;

  /// Calculates the bending energy of the transformation
  virtual double BendingEnergy(double, double, double, double = 0, double = NaN, bool = true) const;

  /// Approximates the bending energy on the control point lattice
  virtual double BendingEnergy(bool = false, bool = true) const;

  /// Approximates the bending energy on the specified discrete domain
  virtual double BendingEnergy(const ImageAttributes &, double = NaN, bool = true) const;

  /// Approximates and adds the gradient of the bending energy on the control point
  /// lattice w.r.t the transformation parameters using the given weight
  virtual void BendingEnergyGradient(double *, double = 1, bool = false, bool = true, bool = true) const;

  // ---------------------------------------------------------------------------
  // I/O

  // Do not hide methods of base class
  using FreeFormTransformation3D::Print;

  /// Prints the parameters of the transformation
  virtual void Print(ostream &, Indent = 0) const;

  /// Whether this transformation can read a file of specified type (i.e. format)
  virtual bool CanRead(TransformationType) const;

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// String conversion of enum values
// =============================================================================

// ----------------------------------------------------------------------------
template <>
inline string ToString(const BSplineFreeFormTransformation3D::ParametricGradientType &value, int w, char c, bool left)
{
  const char *str;
  switch (value) {
    case BSplineFreeFormTransformation3D::PG_Default:       str = "Default";       break;
    case BSplineFreeFormTransformation3D::PG_Analytic:      str = "Analytic";      break;
    case BSplineFreeFormTransformation3D::PG_Convolution:   str = "Convolution";   break;
    case BSplineFreeFormTransformation3D::PG_Approximation: str = "Approximation"; break;
    default:                                                str = "Unknown";       break;
  }
  return ToString(str, w, c, left);
}

// ----------------------------------------------------------------------------
template <>
inline bool FromString(const char *str, BSplineFreeFormTransformation3D::ParametricGradientType &value)
{
  string lstr = ToLower(str);
  if (lstr == "default") {
    value = BSplineFreeFormTransformation3D::PG_Default;
  } else if (lstr == "analytic" || lstr == "exact") {
    value = BSplineFreeFormTransformation3D::PG_Analytic;
  } else if (lstr == "convolution") {
    value = BSplineFreeFormTransformation3D::PG_Convolution;
  } else if (lstr == "dmffd" || lstr == "directlymanipulated" || lstr == "directly manipulated" || lstr == "approximation") {
    value = BSplineFreeFormTransformation3D::PG_Approximation;
  } else {
    return false;
  }
  return true;
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
inline void BSplineFreeFormTransformation3D
::Evaluate(double &x, double &y, double &z) const
{
  Vector d = ((_z == 1) ? _FFD2D(x, y) : _FFD(x, y, z));
  x = d._x, y = d._y, z = d._z;
}

// -----------------------------------------------------------------------------
inline void BSplineFreeFormTransformation3D
::EvaluateInside(double &x, double &y, double &z) const
{
  Vector d = ((_z == 1) ? _FFD2D.GetInside(x, y) : _FFD.GetInside(x, y, z));
  x = d._x, y = d._y, z = d._z;
}

// -----------------------------------------------------------------------------
inline void BSplineFreeFormTransformation3D
::EvaluateJacobianWorld(Matrix &jac, double x, double y) const
{
  // Compute 1st order derivatives
  EvaluateJacobian(jac, x, y);
  // Convert derivatives to world coordinates
  JacobianToWorld(jac);
}

// -----------------------------------------------------------------------------
inline void BSplineFreeFormTransformation3D
::EvaluateJacobianWorld(Matrix &jac, double x, double y, double z) const
{
  // Compute 1st order derivatives
  EvaluateJacobian(jac, x, y, z);
  // Convert derivatives to world coordinates
  JacobianToWorld(jac);
}

// -----------------------------------------------------------------------------
inline void BSplineFreeFormTransformation3D
::EvaluateJacobianDOFs(double jac[3], int i, int j, double x, double y) const
{
  jac[0] = jac[1] = jac[2] = Kernel::Weight(x - i) * Kernel::Weight(y - j);
}

// -----------------------------------------------------------------------------
inline void BSplineFreeFormTransformation3D
::EvaluateJacobianDOFs(double jac[3], int    i, int    j, int    k,
                                      double x, double y, double z) const
{
  jac[0] = jac[1] = jac[2] = Kernel::Weight(x - i) *
                             Kernel::Weight(y - j) *
                             Kernel::Weight(z - k);
}

// -----------------------------------------------------------------------------
inline void BSplineFreeFormTransformation3D
::EvaluateJacobianDOFs(double jac[3], int    i, int    j, int    k, int,
                                      double x, double y, double z, double) const
{
  EvaluateJacobianDOFs(jac, i, j, k, x, y, z);
}

// -----------------------------------------------------------------------------
inline void BSplineFreeFormTransformation3D
::EvaluateDerivativeOfJacobianWrtDOF(Matrix &dJdp, int dof, double x, double y) const
{
  const int dim = this->DOFToDimension(dof);

  int i, j;
  this->IndexToLattice(this->DOFToIndex(dof), i, j);

  double val;
  if (dim == 0) val = Kernel::B_I(x - i) * Kernel::B  (y - j);
  else          val = Kernel::B  (x - i) * Kernel::B_I(y - j);

  dJdp.Initialize(3, 3);
  dJdp(0, dim) = dJdp(1, dim) = dJdp(2, dim) = val;
}

// -----------------------------------------------------------------------------
inline void BSplineFreeFormTransformation3D
::EvaluateDerivativeOfJacobianWrtDOF(Matrix &dJdp, int dof, double x, double y, double z) const
{
  const int dim = this->DOFToDimension(dof);

  int i, j, k;
  this->IndexToLattice(this->DOFToIndex(dof), i, j, k);

  double val;
  if      (dim == 0) val = Kernel::B_I(x - i) * Kernel::B  (y - j) * Kernel::B  (z - k);
  else if (dim == 1) val = Kernel::B  (x - i) * Kernel::B_I(y - j) * Kernel::B  (z - k);
  else               val = Kernel::B  (x - i) * Kernel::B  (y - j) * Kernel::B_I(z - k);

  dJdp.Initialize(3, 3);
  dJdp(0, dim) = dJdp(1, dim) = dJdp(2, dim) = val;
}

// =============================================================================
// Point transformation
// =============================================================================

// -----------------------------------------------------------------------------
inline void BSplineFreeFormTransformation3D
::LocalTransform(double &x, double &y, double &z, double, double) const
{
  // Convert to lattice coordinates
  double dx = x, dy = y, dz = z;
  this->WorldToLattice(dx, dy, dz);
  // Evaluate displacement
  Evaluate(dx, dy,  dz);
  // Transform point
  x += dx, y += dy, z += dz;
}

// =============================================================================
// Derivatives
// =============================================================================

// -----------------------------------------------------------------------------
inline void BSplineFreeFormTransformation3D::FFDJacobianWorld(Matrix &jac, double x, double y, double z, double, double) const
{
  // Convert to lattice coordinates
  this->WorldToLattice(x, y, z);
  // Compute 1st order derivatives
  if (_z == 1) EvaluateJacobianWorld(jac, x, y);
  else         EvaluateJacobianWorld(jac, x, y, z);
  // Add derivatives of "x" term in T(x) = x + FFD(x)
  jac(0, 0) += 1.0;
  jac(1, 1) += 1.0;
  jac(2, 2) += 1.0;
}

// -----------------------------------------------------------------------------
inline void BSplineFreeFormTransformation3D::LocalJacobian(Matrix &jac, double x, double y, double z, double t, double t0) const
{
  BSplineFreeFormTransformation3D::FFDJacobianWorld(jac, x, y, z, t, t0);
}

// -----------------------------------------------------------------------------
inline void BSplineFreeFormTransformation3D::LocalHessian(Matrix hessian[3], double x, double y, double z, double, double) const
{
  // Convert to lattice coordinates
  this->WorldToLattice(x, y, z);
  // Compute 2nd order derivatives
  if (_z == 1) EvaluateHessian(hessian, x, y);
  else         EvaluateHessian(hessian, x, y, z);
  // Convert derivatives to world coordinates
  HessianToWorld(hessian);
}

// -----------------------------------------------------------------------------
inline void BSplineFreeFormTransformation3D
::JacobianDOFs(double jac[3], int ci, int cj, int ck, double x, double y, double z) const
{
  // Convert point to lattice coordinates
  this->WorldToLattice(x, y, z);
  // Evaluate derivatives w.r.t. transformation parameters
  if (_z == 1) EvaluateJacobianDOFs(jac, ci, cj,     x, y);
  else         EvaluateJacobianDOFs(jac, ci, cj, ck, x, y, z);
}

// -----------------------------------------------------------------------------
inline void BSplineFreeFormTransformation3D
::DeriveJacobianWrtDOF(Matrix &dJdp, int dof, double x, double y, double z, double, double) const
{
  // Convert point to lattice coordinates
  this->WorldToLattice(x, y, z);
  // Evaluate derivatives w.r.t. transformation parameters
  if (_z == 1) EvaluateDerivativeOfJacobianWrtDOF(dJdp, dof, x, y);
  else         EvaluateDerivativeOfJacobianWrtDOF(dJdp, dof, x, y, z);
  // Convert derivatives to world coordinates
  JacobianToWorld(dJdp);
}

// -----------------------------------------------------------------------------
inline void BSplineFreeFormTransformation3D
::EvaluateJacobianDetDerivative(double dJ[3], const Matrix &adj,
                                int i, int j, int k, double x, double y, double z,
                                bool wrt_world, bool use_spacing) const
{
  EvaluateJacobianDetDerivative(dJ, adj, x - i, y - j, z - k, wrt_world, use_spacing);
}

// -----------------------------------------------------------------------------
inline void BSplineFreeFormTransformation3D
::EvaluateJacobianDetDerivative(double dJ[3], const Matrix &adj,
                                int cp, double x, double y, double z,
                                bool wrt_world, bool use_spacing) const
{
  int i, j, k;
  this->IndexToLattice(cp, i, j, k);
  EvaluateJacobianDetDerivative(dJ, adj, x - i, y - j, z - k, wrt_world, use_spacing);
}

// -----------------------------------------------------------------------------
inline void BSplineFreeFormTransformation3D
::FFDJacobianDetDerivative(double dJ[3], const Matrix &adj,
                           int cp, double x, double y, double z, double, double,
                           bool wrt_world, bool use_spacing) const
{
  int i, j, k;
  this->IndexToLattice(cp, i, j, k);
  this->WorldToLattice(x, y, z);
  EvaluateJacobianDetDerivative(dJ, adj, x - i, y - j, z - k, wrt_world, use_spacing);
}


} // namespace mirtk

#endif // MIRTK_BSplineFreeFormTransformation3D_H
