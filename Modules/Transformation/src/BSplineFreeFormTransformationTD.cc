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

#include "mirtk/BSplineFreeFormTransformationTD.h"

#include "mirtk/Math.h"
#include "mirtk/Memory.h"
#include "mirtk/DisplacementToVelocityFieldBCH.h"
#include "mirtk/ImageToInterpolationCoefficients.h"

#include "FreeFormTransformationIntegration.h"


namespace mirtk {


// =============================================================================
// Integration methods
// =============================================================================

MIRTK_FFDIM2(RKE1,   BSplineFreeFormTransformationTD);
MIRTK_FFDIM2(RKEH12, BSplineFreeFormTransformationTD);
MIRTK_FFDIM2(RKE2,   BSplineFreeFormTransformationTD);
MIRTK_FFDIM2(RKH2,   BSplineFreeFormTransformationTD);
MIRTK_FFDIM2(RKBS23, BSplineFreeFormTransformationTD);
MIRTK_FFDIM2(RK4,    BSplineFreeFormTransformationTD);
MIRTK_FFDIM2(RKF45,  BSplineFreeFormTransformationTD);
MIRTK_FFDIM2(RKCK45, BSplineFreeFormTransformationTD);
MIRTK_FFDIM2(RKDP45, BSplineFreeFormTransformationTD);

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
BSplineFreeFormTransformationTD::BSplineFreeFormTransformationTD()
:
  _IntegrationMethod(FFDIM_RKE2),
  _MinTimeStep      (0.01),
  _MaxTimeStep      (0.1),
  _Tolerance        (1.e-3)
{
  _ExtrapolationMode = Extrapolation_NN;
}

// -----------------------------------------------------------------------------
BSplineFreeFormTransformationTD
::BSplineFreeFormTransformationTD(const ImageAttributes &attr,
                                  double dx, double dy, double dz, double dt)
:
  _IntegrationMethod(FFDIM_RKE2),
  _MinTimeStep      (0.01),
  _MaxTimeStep      (0.1),
  _Tolerance        (1.e-3)
{
  _ExtrapolationMode = Extrapolation_NN;
  Initialize(attr, dx, dy, dz, dt);
}

// -----------------------------------------------------------------------------
BSplineFreeFormTransformationTD
::BSplineFreeFormTransformationTD(const BaseImage &target,
                                  double dx, double dy, double dz, double dt)
:
  _IntegrationMethod(FFDIM_RKE2),
  _MinTimeStep      (0.01),
  _MaxTimeStep      (0.1),
  _Tolerance        (1.e-3)
{
  _ExtrapolationMode = Extrapolation_NN;
  Initialize(target.Attributes(), dx, dy, dz, dt);
}

// -----------------------------------------------------------------------------
BSplineFreeFormTransformationTD
::BSplineFreeFormTransformationTD(const BSplineFreeFormTransformationTD &ffd)
:
  BSplineFreeFormTransformation4D(ffd),
  _IntegrationMethod(ffd._IntegrationMethod),
  _MinTimeStep      (ffd._MinTimeStep),
  _MaxTimeStep      (ffd._MaxTimeStep),
  _Tolerance        (ffd._Tolerance)
{
}

// -----------------------------------------------------------------------------
BSplineFreeFormTransformationTD
::~BSplineFreeFormTransformationTD()
{
}

// =============================================================================
// Approximation/Interpolation
// =============================================================================

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformationTD
::ApproximateDOFs(const GenericImage<double> * const *disp,
                  const double *t1, const double *t2, int no,
                  bool smooth, int nterms, int niter)
{
  // TODO: Refactor this method taking advantage of recent changes related
  //       to velocity field based transformation parameterization.
  //       - Use Vector as voxel type for 4D vector field image.
  //       - Use fast inter-/extrapolators which also handle vector images.
  //       - Possibly use cubic B-spline interpolator for BCH filter,
  //         similar to the BSplineFreeFormTransformationSV::Add method.
  InterpolateImageFunction *f;

  DisplacementToVelocityFieldBCH<double> dtov;
  dtov.SmoothVelocities  (smooth);
  dtov.NumberOfTerms     (nterms);
  dtov.NumberOfIterations(niter);

  GenericImage<double> vx(_attr);
  GenericImage<double> vy(_attr);
  GenericImage<double> vz(_attr);

  double *vw = new double[_t];
  memset(vw, 0, _t * sizeof(double));

  GenericImage<double> d(_attr, 3);
  GenericImage<double> v(_attr, 3);
 
  for (int n = 0; n < no; n++) {
    // Get control point frames covered by n-th displacement field
    int l1 = ifloor(this->TimeToLattice(t1[n]));
    int l2 = iceil (this->TimeToLattice(t2[n]));
    if (l1 > l2) swap(l1, l2);

    // Check if within temporal control point domain
    if ((l1 < 0 || l1 >= _t) && (l2 < 0 || l2 >= _t)) continue;

    // Sample displacement field at control points using linear interpolation
    if (disp[n]->GetImageAttributes() == _attr) {
      for (int k = 0; k < _z; ++k)
      for (int j = 0; j < _y; ++j)
      for (int i = 0; i < _x; ++i) {
        d.Put(i, j, k, 0, disp[n]->Get(i, j, k, 0));
        d.Put(i, j, k, 1, disp[n]->Get(i, j, k, 1));
        d.Put(i, j, k, 2, disp[n]->Get(i, j, k, 2));
      }
    } else {
      double x, y, z, vec[3] = {.0, .0, .0};
      f = InterpolateImageFunction::New(Interpolation_Linear, Extrapolation_NN, disp[n]);
      f->Input(disp[n]);
      f->Initialize();
      for (int k = 0; k < _z; ++k)
      for (int j = 0; j < _y; ++j)
      for (int i = 0; i < _x; ++i) {
        x = i, y = j, z = k;
        d.ImageToWorld(x, y, z);
        disp[n]->WorldToImage(x, y, z);
        f->Evaluate(vec, x, y, z);
        d.Put(i, j, k, 0, vec[0]);
        d.Put(i, j, k, 1, vec[1]);
        d.Put(i, j, k, 2, vec[2]);
      }
    }

    // Compute stationary velocity field
    const double T = t2[n] - t1[n];

    dtov.NumberOfSteps(iround(T / _MinTimeStep));
    dtov.Input (&d);
    dtov.Output(&v);
    dtov.Run();

    // Add stationary velocities to frames in [l1, l2]
    for (int l = l1; l <= l2; l++) {
      if (0 <= l && l < _t) {
        for (int k = 0; k < _z; ++k)
        for (int j = 0; j < _y; ++j)
        for (int i = 0; i < _x; ++i) {
          vx(i, j, k, l) += v(i, j, k, 0);
          vy(i, j, k, l) += v(i, j, k, 1);
          vz(i, j, k, l) += v(i, j, k, 2);
        }
        vw[l] += T;
      }
    }
  }
 
  // Compute average velocities at control points
  for (int l = 0; l < _t; l++) {
    if (vw[l] == .0) vw[l] = 1.0;
    for (int k = 0; k < _z; ++k)
    for (int j = 0; j < _y; ++j)
    for (int i = 0; i < _x; ++i) {
      _CPImage(i, j, k, l) = vx(i, j, k, l) / vw[l];
      _CPImage(i, j, k, l) = vy(i, j, k, l) / vw[l];
      _CPImage(i, j, k, l) = vz(i, j, k, l) / vw[l];
    }
  }

  // Convert to B-spline coefficients
  ConvertToSplineCoefficients(3, _CPImage);

  // Clean up
  delete[] vw;
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformationTD
::ApproximateDOFs(const double *, const double *, const double *, const double *,
                  const double *, const double *, const double *, int)
{
  cerr << this->NameOfClass() << "::ApproximateDOFs: Not implemented" << endl;
  cerr << "  --> Try the other overloaded approximation function" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformationTD
::ApproximateDOFsGradient(const double *, const double *, const double *, const double *,
                          const double *, const double *, const double *, int,
                          double *, double) const
{
  cerr << this->NameOfClass() << "::ApproximateDOFsGradient: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
double BSplineFreeFormTransformationTD
::ApproximateAsNew(GenericImage<double> **disp,
                   const double *t1, const double *t2, int no,
                   bool smooth, int nterms, int niter)
{
  // Approximate 3D+t velocity field
  this->ApproximateDOFs(disp, t1, t2, no, smooth, nterms, niter);

  // Evaluate RMS of approximation error
  double error  = .0;
  int    ntotal = 0;
  double dx, dy, dz;

  for (int n = 0; n < no; ++n) {
    GenericImage<double> &d = *disp[n];
    for (int k = 0; k < d.Z(); ++k)
    for (int j = 0; j < d.Y(); ++j)
    for (int i = 0; i < d.X(); ++i) {
      dx = i, dy = j, dz = k;
      d.ImageToWorld(dx, dy, dz);
      this->Displacement(dx, dy, dz, t1[n], t2[n]);
      d(i, j, k, 0) -= dx;
      d(i, j, k, 1) -= dy;
      d(i, j, k, 2) -= dz;
      error += sqrt(d(i, j, k, 0) * d(i, j, k, 0) +
                    d(i, j, k, 1) * d(i, j, k, 1) +
                    d(i, j, k, 2) * d(i, j, k, 2));
      ++ntotal;
    }
  }
  if (ntotal > 0) error /= ntotal;
  return error;
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformationTD::Interpolate(const double *, const double *, const double *)
{
  cerr << this->NameOfClass() << "::Interpolate: Not implemented" << endl;
  exit(1);
}

// =============================================================================
// Parameters (non-DoFs)
// =============================================================================

// -----------------------------------------------------------------------------
bool BSplineFreeFormTransformationTD::Set(const char *name, const char *value)
{
  if (strcmp(name, "Integration method")          == 0 ||
      strcmp(name, "Velocity integration method") == 0) {
    return FromString(value, _IntegrationMethod);
  }
  if (strcmp(name, "Length of integration steps") == 0) {
    double dt;
    if (!FromString(value, dt) || dt <= .0) return false;
    _MinTimeStep = _MaxTimeStep = dt;
    return true;
  }
  if (strcmp(name, "Minimum length of integration steps") == 0) {
    return FromString(value, _MinTimeStep) && _MinTimeStep > .0;
  }
  if (strcmp(name, "Maximum length of integration steps") == 0) {
    return FromString(value, _MaxTimeStep) && _MaxTimeStep > .0;
  }
  if (strcmp(name, "Integration tolerance")          == 0 ||
             strcmp(name, "Velocity integration tolerance") == 0) {
    return FromString(value, _Tolerance) && _Tolerance >= .0;
  }
  return BSplineFreeFormTransformation4D::Set(name, value);
}

// -----------------------------------------------------------------------------
ParameterList BSplineFreeFormTransformationTD::Parameter() const
{
  ParameterList params = BSplineFreeFormTransformation4D::Parameter();
  Insert(params, "Integration method",                  ToString(_IntegrationMethod));
  Insert(params, "Minimum length of integration steps", ToString(_MinTimeStep));
  Insert(params, "Maximum length of integration steps", ToString(_MaxTimeStep));
  Insert(params, "Integration tolerance",               ToString(_Tolerance));
  return params;
}

// =============================================================================
// Point transformation
// =============================================================================

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformationTD::LocalTransform(double &x, double &y, double &z, double t, double t0) const
{
  if      (_IntegrationMethod == FFDIM_RKE1)   RKE1  ::Transform(this, x, y, z, t0, t, _MinTimeStep);
  else if (_IntegrationMethod == FFDIM_RKE2)   RKE2  ::Transform(this, x, y, z, t0, t, _MinTimeStep);
  else if (_IntegrationMethod == FFDIM_RKH2)   RKH2  ::Transform(this, x, y, z, t0, t, _MinTimeStep);
  else if (_IntegrationMethod == FFDIM_RK4)    RK4   ::Transform(this, x, y, z, t0, t, _MinTimeStep);
  else if (_IntegrationMethod == FFDIM_RKEH12) RKEH12::Transform(this, x, y, z, t0, t, _MinTimeStep, _MaxTimeStep, _Tolerance);
  else if (_IntegrationMethod == FFDIM_RKBS23) RKBS23::Transform(this, x, y, z, t0, t, _MinTimeStep, _MaxTimeStep, _Tolerance);
  else if (_IntegrationMethod == FFDIM_RKF45)  RKF45 ::Transform(this, x, y, z, t0, t, _MinTimeStep, _MaxTimeStep, _Tolerance);
  else if (_IntegrationMethod == FFDIM_RKCK45) RKCK45::Transform(this, x, y, z, t0, t, _MinTimeStep, _MaxTimeStep, _Tolerance);
  else if (_IntegrationMethod == FFDIM_RKDP45) RKDP45::Transform(this, x, y, z, t0, t, _MinTimeStep, _MaxTimeStep, _Tolerance);
  else {
    cerr << "BSplineFreeFormTransformationTD::LocalTransform: Unknown integration method: " << _IntegrationMethod << endl;
    exit(1);
  }
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformationTD::TransformAndJacobian(Matrix &jac, double &x, double &y, double &z, double t, double t0) const
{
  if      (_IntegrationMethod == FFDIM_RKE1)   RKE1  ::Jacobian(this, jac, x, y, z, t0, t, _MinTimeStep);
  else if (_IntegrationMethod == FFDIM_RKE2)   RKE2  ::Jacobian(this, jac, x, y, z, t0, t, _MinTimeStep);
  else if (_IntegrationMethod == FFDIM_RKH2)   RKH2  ::Jacobian(this, jac, x, y, z, t0, t, _MinTimeStep);
  else if (_IntegrationMethod == FFDIM_RK4)    RK4   ::Jacobian(this, jac, x, y, z, t0, t, _MinTimeStep);
  else if (_IntegrationMethod == FFDIM_RKEH12) RKEH12::Jacobian(this, jac, x, y, z, t0, t, _MinTimeStep, _MaxTimeStep, _Tolerance);
  else if (_IntegrationMethod == FFDIM_RKBS23) RKBS23::Jacobian(this, jac, x, y, z, t0, t, _MinTimeStep, _MaxTimeStep, _Tolerance);
  else if (_IntegrationMethod == FFDIM_RKF45)  RKF45 ::Jacobian(this, jac, x, y, z, t0, t, _MinTimeStep, _MaxTimeStep, _Tolerance);
  else if (_IntegrationMethod == FFDIM_RKCK45) RKCK45::Jacobian(this, jac, x, y, z, t0, t, _MinTimeStep, _MaxTimeStep, _Tolerance);
  else if (_IntegrationMethod == FFDIM_RKDP45) RKDP45::Jacobian(this, jac, x, y, z, t0, t, _MinTimeStep, _MaxTimeStep, _Tolerance);
  else {
    cerr << "BSplineFreeFormTransformationTD::TransformAndJacobian: Unknown integration method: " << _IntegrationMethod << endl;
    exit(1);
  }
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformationTD::TransformAndJacobianDOFs(Matrix &jac, int ci, int cj, int ck, int cl, double &x, double &y, double &z, double t, double t0) const
{
  if      (_IntegrationMethod == FFDIM_RKE1)   RKE1  ::JacobianDOFs(this, jac, ci, cj, ck, cl, x, y, z, t0, t, _MinTimeStep);
  else if (_IntegrationMethod == FFDIM_RKE2)   RKE2  ::JacobianDOFs(this, jac, ci, cj, ck, cl, x, y, z, t0, t, _MinTimeStep);
  else if (_IntegrationMethod == FFDIM_RKH2)   RKH2  ::JacobianDOFs(this, jac, ci, cj, ck, cl, x, y, z, t0, t, _MinTimeStep);
  else if (_IntegrationMethod == FFDIM_RK4)    RK4   ::JacobianDOFs(this, jac, ci, cj, ck, cl, x, y, z, t0, t, _MinTimeStep);
  else if (_IntegrationMethod == FFDIM_RKEH12) RKEH12::JacobianDOFs(this, jac, ci, cj, ck, cl, x, y, z, t0, t, _MinTimeStep, _MaxTimeStep, _Tolerance);
  else if (_IntegrationMethod == FFDIM_RKBS23) RKBS23::JacobianDOFs(this, jac, ci, cj, ck, cl, x, y, z, t0, t, _MinTimeStep, _MaxTimeStep, _Tolerance);
  else if (_IntegrationMethod == FFDIM_RKF45)  RKF45 ::JacobianDOFs(this, jac, ci, cj, ck, cl, x, y, z, t0, t, _MinTimeStep, _MaxTimeStep, _Tolerance);
  else if (_IntegrationMethod == FFDIM_RKCK45) RKCK45::JacobianDOFs(this, jac, ci, cj, ck, cl, x, y, z, t0, t, _MinTimeStep, _MaxTimeStep, _Tolerance);
  else if (_IntegrationMethod == FFDIM_RKDP45) RKDP45::JacobianDOFs(this, jac, ci, cj, ck, cl, x, y, z, t0, t, _MinTimeStep, _MaxTimeStep, _Tolerance);
  else {
    cerr << "BSplineFreeFormTransformationTD::TransformAndJacobianDOFs: Unknown integration method: " << _IntegrationMethod << endl;
    exit(1);
  }
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformationTD::TransformAndJacobianDOFs(Matrix &jac, int cp, double &x, double &y, double &z, double t, double t0) const
{
  int ci, cj, ck, cl;
  this->IndexToLattice(cp, ci, cj, ck, cl);
  this->TransformAndJacobianDOFs(jac, ci, cj, ck, cl, x, y, z, t, t0);
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformationTD::TransformAndJacobianDOFs(TransformationJacobian &jac, double &x, double &y, double &z, double t, double t0) const
{
  if      (_IntegrationMethod == FFDIM_RKE1)   RKE1  ::JacobianDOFs(this, jac, x, y, z, t0, t, _MinTimeStep);
  else if (_IntegrationMethod == FFDIM_RKE2)   RKE2  ::JacobianDOFs(this, jac, x, y, z, t0, t, _MinTimeStep);
  else if (_IntegrationMethod == FFDIM_RKH2)   RKH2  ::JacobianDOFs(this, jac, x, y, z, t0, t, _MinTimeStep);
  else if (_IntegrationMethod == FFDIM_RK4)    RK4   ::JacobianDOFs(this, jac, x, y, z, t0, t, _MinTimeStep);
  else if (_IntegrationMethod == FFDIM_RKEH12) RKEH12::JacobianDOFs(this, jac, x, y, z, t0, t, _MinTimeStep, _MaxTimeStep, _Tolerance);
  else if (_IntegrationMethod == FFDIM_RKBS23) RKBS23::JacobianDOFs(this, jac, x, y, z, t0, t, _MinTimeStep, _MaxTimeStep, _Tolerance);
  else if (_IntegrationMethod == FFDIM_RKF45)  RKF45 ::JacobianDOFs(this, jac, x, y, z, t0, t, _MinTimeStep, _MaxTimeStep, _Tolerance);
  else if (_IntegrationMethod == FFDIM_RKCK45) RKCK45::JacobianDOFs(this, jac, x, y, z, t0, t, _MinTimeStep, _MaxTimeStep, _Tolerance);
  else if (_IntegrationMethod == FFDIM_RKDP45) RKDP45::JacobianDOFs(this, jac, x, y, z, t0, t, _MinTimeStep, _MaxTimeStep, _Tolerance);
  else {
    cerr << "BSplineFreeFormTransformationTD::TransformAndJacobianDOFs: Unknown integration method: " << _IntegrationMethod << endl;
    exit(1);
  }
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformationTD
::ParametricGradient(const PointSet &pos, const Vector3D<double> *in,
                     double *out, double t, double t0, double w) const
{
  // Cannot use irtk(BSpline)FreeFormTransformation3D::ParametricGradient
  // computation which assumes a local support region of the interpolation kernel
  // because point trajectories may enter support regions of more than one kernel
  FreeFormTransformation::ParametricGradient(pos, in, out, t, t0, w);
}

// =============================================================================
// I/O
// =============================================================================

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformationTD::Print(ostream &os, Indent indent) const
{
  os << indent << "B-spline TD FFD:" << endl;
  ++indent;
  FreeFormTransformation4D::Print(os, indent);
  os << indent << "Numerical integration:" << endl;
  ++indent;
  os << indent << "Method:      " << ToString(_IntegrationMethod) << endl;
  if (_MinTimeStep < _MaxTimeStep && _IntegrationMethod >= FFDIM_RKEH12) {
    os << indent << "Step size:   [" << _MinTimeStep << ", " << _MaxTimeStep << "]" << endl;
    os << indent << "Local error: " << _Tolerance << endl;
  } else {
    os << indent << "Step size:   " << _MinTimeStep << endl;
  }
}

// -----------------------------------------------------------------------------
bool BSplineFreeFormTransformationTD::CanRead(TransformationType format) const
{
  switch (format) {
    case TRANSFORMATION_BSPLINE_FFD_TD_v1:
    case TRANSFORMATION_BSPLINE_FFD_TD_v2:
    case TRANSFORMATION_BSPLINE_FFD_TD_v3:
    case TRANSFORMATION_BSPLINE_FFD_TD_v4:
      return true;
    default:
      return false;
  }
}

// -----------------------------------------------------------------------------
Cifstream &BSplineFreeFormTransformationTD::ReadDOFs(Cifstream &from, TransformationType format)
{
  // Read FFD data
  switch (format) {
    case TRANSFORMATION_BSPLINE_FFD_TD_v1:
    case TRANSFORMATION_BSPLINE_FFD_TD_v2:
      BSplineFreeFormTransformation4D::ReadDOFs(from, TRANSFORMATION_BSPLINE_FFD_4D_v1);
      break;
    case TRANSFORMATION_BSPLINE_FFD_TD_v3:
      BSplineFreeFormTransformation4D::ReadDOFs(from, TRANSFORMATION_BSPLINE_FFD_4D_v2);
      break;
    default:
      BSplineFreeFormTransformation4D::ReadDOFs(from, TRANSFORMATION_BSPLINE_FFD_4D);
  }

  // Read integration method
  if (format == TRANSFORMATION_BSPLINE_FFD_4D_v1) {
    _IntegrationMethod = FFDIM_RKE1;
  } else {
    unsigned int tmp = FFDIM_Unknown;
    if (!from.ReadAsUInt(&tmp, 1)) {
      cerr << "BSplineFreeFormTransformationTD::ReadDOFs: Failed to read integration method" << endl;
      exit(1);
    }
    _IntegrationMethod = static_cast<FFDIntegrationMethod>(tmp);
  }

  // Read minimum/maximum time step length
  if (!from.ReadAsDouble(&_MinTimeStep, 1) || !from.ReadAsDouble(&_MaxTimeStep, 1)) {
    cerr << "BSplineFreeFormTransformationTD::ReadDOFs: Failed to read min/max time step" << endl;
    exit(1);
  }

  // Read local error tolerance
  if (format == TRANSFORMATION_BSPLINE_FFD_TD_v1) {
    _Tolerance = 1.e-3;
  } else {
    if (!from.ReadAsDouble(&_Tolerance, 1)) {
      cerr << "BSplineFreeFormTransformationTD::ReadDOFs: Failed to read integration tolerance" << endl;
      exit(1);
    }
  }

  return from;
}

// -----------------------------------------------------------------------------
Cofstream &BSplineFreeFormTransformationTD::WriteDOFs(Cofstream &to) const
{
  // Write FFD data
  BSplineFreeFormTransformation4D::WriteDOFs(to);

  // Write integration method
  unsigned int integration_method = _IntegrationMethod;
  to.WriteAsUInt(&integration_method, 1);

  // Write minimum/maximum time step length
  to.WriteAsDouble(&_MinTimeStep, 1);
  to.WriteAsDouble(&_MaxTimeStep, 1);

  // Write local error tolerance
  to.WriteAsDouble(&_Tolerance, 1);

  return to;
}

// =============================================================================
// Others
// =============================================================================

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformationTD::Invert()
{
  Vector tmp;
  for (int l = 0; l < _t / 2; ++l) {
    int L = _t - l - 1;
    for (int k = 0; k < _z; ++k)
    for (int j = 0; j < _y; ++j)
    for (int i = 0; i < _x; ++i) {
      tmp = _CPImage(i, j, k, l);
      _CPImage(i, j, k, l) = - _CPImage(i, j, k, L);
      _CPImage(i, j, k, L) = - tmp;
    }
  }
}


} // namespace mirtk
