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

#ifndef MIRTK_BSplineFreeFormTransformationTD_H
#define MIRTK_BSplineFreeFormTransformationTD_H

#include "mirtk/BSplineFreeFormTransformation4D.h"
#include "mirtk/FFDIntegrationMethod.h"


namespace mirtk {


/**
 * Temporal diffeomorphic free-form transformation.
 *
 * This class implements a free-form transformation which is represented
 * by a time-varying velocity field (3D+t). The 3D displacement field at a
 * specific time is obtained by integrating the velocity field starting at the
 * reference time point. The integration steps are adjusted if necessary in
 * order to ensure that the resulting spatial transformation is diffeomorphic.
 *
 * For more details about the implementation see De Craene et al. (2012).
 * Temporal diffeomorphic free-form deformation: application to motion and
 * strain estimation from 3D echocardiography.
 * Medical image analysis, 16(2), 427, 2012. doi:10.1016/j.media.2011.10.006
 */

class BSplineFreeFormTransformationTD : public BSplineFreeFormTransformation4D
{
  mirtkTransformationMacro(BSplineFreeFormTransformationTD);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Numerical integration method
  mirtkPublicAttributeMacro(FFDIntegrationMethod, IntegrationMethod);

  /// Minimum length of temporal steps (in ms)
  mirtkPublicAttributeMacro(double, MinTimeStep);

  /// Maximum length of temporal steps (in ms)
  mirtkPublicAttributeMacro(double, MaxTimeStep);

  /// Local integration error tolerance (in mm)
  mirtkPublicAttributeMacro(double, Tolerance);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Default constructor
  BSplineFreeFormTransformationTD();

  /// Construct free-form transformation for given image domain and lattice spacing
  explicit BSplineFreeFormTransformationTD(const ImageAttributes &,
                                           double = -1, double = -1, double = -1, double = -1);

  /// Construct free-form transformation for given target image and lattice spacing
  explicit BSplineFreeFormTransformationTD(const BaseImage &,
                                           double, double, double, double);

  /// Copy constructor
  BSplineFreeFormTransformationTD(const BSplineFreeFormTransformationTD &);

  /// Destructor
  virtual ~BSplineFreeFormTransformationTD();

  // ---------------------------------------------------------------------------
  // Approximation/Interpolation

  using BSplineFreeFormTransformation4D::ApproximateAsNew;

  /// Approximate displacements: This function takes a set of 3D displacement fields
  /// and corresponding time point and temporal interval. Given these inputs,
  /// it finds a time-varying velocity field which approximates these displacements.
  virtual void ApproximateDOFs(const GenericImage<double> * const *,
                               const double *, const double *, int,
                               bool = false, int = 3, int = 8);

  /// Approximates displacements: This function takes a set of points and a set
  /// of displacements and finds a !new! FFD which approximates these displacements.
  virtual void ApproximateDOFs(const double *, const double *, const double *, const double *,
                               const double *, const double *, const double *, int);

  /// Finds gradient of approximation error: This function takes a set of points
  /// and a set of errors. It finds a gradient w.r.t. the transformation parameters
  /// which minimizes the L2 norm of the approximation error and adds it to the
  /// input gradient with the given weight.
  virtual void ApproximateDOFsGradient(const double *, const double *, const double *, const double *,
                                       const double *, const double *, const double *, int,
                                       double *, double = 1.0) const;

  /// Approximate displacements: This function takes a set of 3D displacement fields
  /// and corresponding time point and temporal interval. Given these inputs,
  /// it finds a time-varying velocity field which approximates these displacements.
  /// The displacements are replaced by the residual displacements of the newly
  /// approximated transformation. Returns the approximation error of the resulting FFD.
  virtual double ApproximateAsNew(GenericImage<double> **,
                                  const double *, const double *, int,
                                  bool = false, int = 3, int = 8);

  /// Interpolates displacements: This function takes a set of displacements defined at
  /// the control points and finds a time-varying velocity field such that the
  /// resulting transformation interpolates these displacements.
  virtual void Interpolate(const double *, const double *, const double *);

  // ---------------------------------------------------------------------------
  // Parameters (non-DoFs)

  /// Set named (non-DoF) parameter from value as string
  virtual bool Set(const char *, const char *);

  /// Get (non-DoF) parameters as key/value as string map
  virtual ParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Point transformation
  using BSplineFreeFormTransformation4D::Displacement;

  /// Transforms a single point
  virtual void LocalTransform(double &, double &, double &, double, double) const;

  /// Transforms a single point using the inverse transformation
  virtual bool LocalInverse(double &, double &, double &, double, double) const;

  // ---------------------------------------------------------------------------
  // Derivatives
  using BSplineFreeFormTransformation4D::JacobianDOFs;
  using BSplineFreeFormTransformation4D::ParametricGradient;

  /// Calculates the Jacobian of the (local) transformation w.r.t world coordinates
  /// and transforms the given point at the same time
  virtual void TransformAndJacobian(Matrix &, double &, double &, double &, double, double) const;

  /// Calculates the Jacobian of the transformation w.r.t the transformation parameters
  /// of the specified control point and transforms the given point at the same time
  virtual void TransformAndJacobianDOFs(Matrix &, int, int, int, int, double &, double &, double &, double, double) const;

  /// Calculates the Jacobian of the transformation w.r.t the transformation parameters
  /// of the specified control point and transforms the given point at the same time
  virtual void TransformAndJacobianDOFs(Matrix &, int, double &, double &, double &, double, double) const;

  /// Calculates the Jacobian of the transformation w.r.t all transformation parameters
  /// and transforms the given point at the same time
  virtual void TransformAndJacobianDOFs(TransformationJacobian &, double &, double &, double &, double, double) const;

  /// Calculates the Jacobian of the local transformation w.r.t world coordinates
  virtual void LocalJacobian(Matrix &, double, double, double, double, double) const;

  /// Calculates the Hessian for each component of the local transformation w.r.t world coordinates
  virtual void LocalHessian(Matrix [3], double, double, double, double, double) const;

  /// Calculates the Jacobian of the transformation w.r.t the transformation parameters
  virtual void JacobianDOFs(Matrix &, int, int, int, int, double, double, double, double, double) const;

  /// Calculates the Jacobian of the transformation w.r.t the transformation parameters
  virtual void JacobianDOFs(double [3], int, int, int, int, double, double, double, double, double) const;

  /// Calculates the Jacobian of the transformation w.r.t the transformation parameters
  virtual void JacobianDOFs(TransformationJacobian &, double, double, double, double, double) const;

  /// Applies the chain rule to convert spatial non-parametric gradient
  /// to a gradient w.r.t the parameters of this transformation
  virtual void ParametricGradient(const GenericImage<double> *, double *,
                                  const WorldCoordsImage *,
                                  const WorldCoordsImage *,
                                  double = 1, double = 1) const;

  /// Applies the chain rule to convert point-wise non-parametric gradient
  /// to a gradient w.r.t the parameters of this transformation.
  virtual void ParametricGradient(const PointSet &, const Vector3D<double> *,
                                  double *, double = 0, double = -1, double = 1) const;

  // ---------------------------------------------------------------------------
  // I/O

  // Do not hide methods of base class
  using BSplineFreeFormTransformation4D::Print;

  /// Prints the parameters of the transformation
  virtual void Print(ostream &, Indent = 0) const;

  /// Whether this transformation can read a file of specified type (i.e. format)
  virtual bool CanRead(TransformationType) const;

protected:

  /// Reads transformation parameters from a file stream
  virtual Cifstream &ReadDOFs(Cifstream &, TransformationType);

  /// Writes transformation parameters to a file stream
  virtual Cofstream &WriteDOFs(Cofstream &) const;

public:

  // ---------------------------------------------------------------------------
  // Others

  /// Invert the transformation
  virtual void Invert();

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////
  
// =============================================================================
// Point transformation
// =============================================================================

// -----------------------------------------------------------------------------
inline bool BSplineFreeFormTransformationTD
::LocalInverse(double &x, double &y, double &z, double t, double t0) const
{
  this->LocalTransform(x, y, z, t, t0);
  return true;
}

// =============================================================================
// Derivatives
// =============================================================================

// -----------------------------------------------------------------------------
inline void BSplineFreeFormTransformationTD
::LocalJacobian(Matrix &jac, double x, double y, double z, double t, double t0) const
{
  this->TransformAndJacobian(jac, x, y, z, t, t0);
}

// -----------------------------------------------------------------------------
inline void BSplineFreeFormTransformationTD
::LocalHessian(Matrix [3], double, double, double, double, double) const
{
  cerr << this->NameOfClass() << "::LocalHessian: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
inline void BSplineFreeFormTransformationTD
::JacobianDOFs(Matrix &jac, int    i, int    j, int    k, int    l,
                                double x, double y, double z, double t, double t0) const
{
  this->TransformAndJacobianDOFs(jac, i, j, k, l, x, y, z, t, t0);
}

// -----------------------------------------------------------------------------
inline void BSplineFreeFormTransformationTD
::JacobianDOFs(TransformationJacobian &jac, double x, double y, double z, double t, double t0) const
{
  jac.Clear();
  this->TransformAndJacobianDOFs(jac, x, y, z, t, t0);
}

// -----------------------------------------------------------------------------
inline void BSplineFreeFormTransformationTD
::JacobianDOFs(double [3], int, int, int, int, double, double, double, double, double) const
{
  cerr << this->NameOfClass() << "::JacobianDOFs: Jacobian is full symmetric 3x3 matrix, not only a diagonal matrix" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
inline void BSplineFreeFormTransformationTD
::ParametricGradient(const GenericImage<double> *in, double *out,
                     const GenericImage<double> *i2w, const WorldCoordsImage *wc,
                     double t0, double w) const
{
  // Use general implementation provided by FFD base class, not the one
  // of FreeFormTransformation4D which is only valid for FFDs that are
  // parameterized by displacements instead of velocities.
  FreeFormTransformation::ParametricGradient(in, out, i2w, wc, t0, w);
}

// -----------------------------------------------------------------------------
/*
void BSplineFreeFormTransformationTD
::ParametricGradient(const GenericImage<double> **in, int n, double *out,
                     const GenericImage<double> *i2w, const WorldCoordsImage *wc,
                     double *t0, double w) const
{
  // TODO: Can be more efficient by re-using already computed trajectories
  //       especially when different source images are frames of a temporal
  //       sequence. See ImageTDFFDRegistration.
}
*/


} // namespace mirtk

#endif // MIRTK_BSplineFreeFormTransformationTD_H
