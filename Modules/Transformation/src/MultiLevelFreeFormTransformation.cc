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

#include "mirtk/MultiLevelFreeFormTransformation.h"

#include "mirtk/Memory.h"

#include "TransformationUtils.h"


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
MultiLevelFreeFormTransformation::MultiLevelFreeFormTransformation()
{
}

// -----------------------------------------------------------------------------
MultiLevelFreeFormTransformation::MultiLevelFreeFormTransformation(const RigidTransformation &t)
:
  MultiLevelTransformation(t)
{
}

// -----------------------------------------------------------------------------
MultiLevelFreeFormTransformation::MultiLevelFreeFormTransformation(const AffineTransformation &t)
:
  MultiLevelTransformation(t)
{
}

// -----------------------------------------------------------------------------
MultiLevelFreeFormTransformation::MultiLevelFreeFormTransformation(const MultiLevelFreeFormTransformation &t)
:
  MultiLevelTransformation(t)
{
}

// -----------------------------------------------------------------------------
MultiLevelFreeFormTransformation::~MultiLevelFreeFormTransformation()
{
}

// =============================================================================
// Levels
// =============================================================================

// -----------------------------------------------------------------------------
void MultiLevelFreeFormTransformation::CombineLocalTransformation()
{
  FreeFormTransformation *first = NULL, *second = NULL;
  while (_NumberOfLevels > 1) {
    first  = this->PopLocalTransformation();
    second = this->PopLocalTransformation();
    if (first->NumberOfDOFs() == second->NumberOfDOFs()) {
      for (int i = 0; i < first->NumberOfDOFs(); ++i) {
        second->Put(i, first->Get(i) + second->Get(i));
      }
    } else {
      cerr << this->NameOfClass() << "::CombineLocalTransformation: Only implemented for transformations with equal number of DOFs" << endl;
      exit(1);
    }
    this->PushLocalTransformation(second);
    delete first;
  }
}

// -----------------------------------------------------------------------------
void MultiLevelFreeFormTransformation::MergeGlobalIntoLocalDisplacement()
{
  // Do nothing if global transformation is the identity
  if (_GlobalTransformation.IsIdentity()) return;

  FreeFormTransformation *ffd = this->GetLocalTransformation(0);

  // Get a copy for making a FFD interpolation of the global affine component
  Transformation *t = Transformation::New(ffd);
  FreeFormTransformation *ffdCopy = dynamic_cast<FreeFormTransformation *>(t);
  if (ffdCopy == NULL) {
    cerr << this->NameOfClass() << "::MergeGlobalIntoLocalDisplacement: Failed to copy local transformation at level 0" << endl;
    exit(1);
  }

  // Interpolate global transformation by FFD
  InterpolateGlobalDisplacement(ffdCopy);

  // Add the calculated coefficients to the coefficients of the first FFD
  for (int i = 0; i < ffd->NumberOfDOFs(); ++i) {
    ffd->Put(i, ffd->Get(i) + ffdCopy->Get(i));
  }

  // Reset matrix previously used for global transformation to identity
  _GlobalTransformation.Reset();

  // Clean up
  delete ffdCopy;
}

// =============================================================================
// Approximation
// =============================================================================

// -----------------------------------------------------------------------------
double MultiLevelFreeFormTransformation
::Approximate(const ImageAttributes &domain, double *dx, double *dy, double *dz,
              int niter, double max_error)
{
  // Check input arguments
  const int no = domain.NumberOfPoints();
  if (no <= 0) return .0;
  double error = numeric_limits<double>::infinity();
  if (niter < 1) return error;

  // Compute world coordinates of lattice points
  double *x = Allocate<double>(no);
  double *y = Allocate<double>(no);
  double *z = Allocate<double>(no);
  double *t = Allocate<double>(no);
  domain.LatticeToWorld(x, y, z, t);

  // Transform world coordinates by global transformation and passive levels
  double *wx = Allocate<double>(no);
  double *wy = Allocate<double>(no);
  double *wz = Allocate<double>(no);
  memcpy(wx, x, no * sizeof(double));
  memcpy(wy, y, no * sizeof(double));
  memcpy(wz, z, no * sizeof(double));

  GetGlobalTransformation()->Transform(no, wx, wy, wz, t);
  for (int lvl = 0; lvl < NumberOfLevels(); ++lvl) {
    if (!LocalTransformationIsActive(lvl)) {
      GetLocalTransformation(lvl)->Transform(no, wx, wy, wz, t);
    }
  }

  // Compute residual displacements
  for (int idx = 0; idx < no; ++idx) {
    dx[idx] -= (wx[idx] - x[idx]);
    dy[idx] -= (wy[idx] - y[idx]);
    dz[idx] -= (wz[idx] - z[idx]);
  }

  // Free memory
  Deallocate(wx);
  Deallocate(wy);
  Deallocate(wz);
  Deallocate(x);
  Deallocate(y);
  Deallocate(z);
  Deallocate(t);

  // Approximate residual displacements by active levels
  FreeFormTransformation *ffd;
  for (int lvl = 0; lvl < NumberOfLevels(); ++lvl) {
    ffd = GetLocalTransformation(lvl);
    if (LocalTransformationIsActive(lvl)) {
      if (error > max_error) {
        error = ffd->Approximate(domain, dx, dy, dz, niter, max_error);
      } else {
        ffd->Reset();
      }
    }
    ++lvl;
  }

  return error;
}

// -----------------------------------------------------------------------------
double MultiLevelFreeFormTransformation
::Approximate(const double *x,  const double *y,  const double *z,
              double       *dx, double       *dy, double       *dz, int no,
              int niter, double max_error)
{
  // Check input arguments
  if (no <= 0) return .0;
  double error = numeric_limits<double>::infinity();
  if (niter < 1) return error;

  // Transform world coordinates by global transformation and passive levels
  double *wx = Allocate<double>(no);
  double *wy = Allocate<double>(no);
  double *wz = Allocate<double>(no);
  memcpy(wx, x, no * sizeof(double));
  memcpy(wy, y, no * sizeof(double));
  memcpy(wz, z, no * sizeof(double));

  GetGlobalTransformation()->Transform(no, wx, wy, wz);
  for (int lvl = 0; lvl < NumberOfLevels(); ++lvl) {
    if (!LocalTransformationIsActive(lvl)) {
      GetLocalTransformation(lvl)->Transform(no, wx, wy, wz);
    }
  }

  // Compute residual displacements
  for (int idx = 0; idx < no; ++idx) {
    dx[idx] -= (wx[idx] - x[idx]);
    dy[idx] -= (wy[idx] - y[idx]);
    dz[idx] -= (wz[idx] - z[idx]);
  }

  // Free memory
  Deallocate(wx);
  Deallocate(wy);
  Deallocate(wz);

  // Approximate residual displacements by active levels
  FreeFormTransformation *ffd;
  for (int lvl = 0; lvl < NumberOfLevels(); ++lvl) {
    ffd = GetLocalTransformation(lvl);
    if (LocalTransformationIsActive(lvl)) {
      if (error > max_error) {
        error = ffd->Approximate(x, y, z, dx, dy, dz, no, niter, max_error);
      } else {
        ffd->Reset();
      }
    }
    ++lvl;
  }

  return error;
}

// -----------------------------------------------------------------------------
double MultiLevelFreeFormTransformation
::Approximate(const double *x,  const double *y,  const double *z,  const double *t,
              double       *dx, double       *dy, double       *dz, int no,
              int niter, double max_error)
{
  // Check input arguments
  if (no <= 0) return .0;
  double error = numeric_limits<double>::infinity();
  if (niter < 1) return error;

  // Transform world coordinates by global transformation and passive levels
  double *wx = Allocate<double>(no);
  double *wy = Allocate<double>(no);
  double *wz = Allocate<double>(no);
  memcpy(wx, x, no * sizeof(double));
  memcpy(wy, y, no * sizeof(double));
  memcpy(wz, z, no * sizeof(double));

  GetGlobalTransformation()->Transform(no, wx, wy, wz, t);
  for (int lvl = 0; lvl < NumberOfLevels(); ++lvl) {
    if (!LocalTransformationIsActive(lvl)) {
      GetLocalTransformation(lvl)->Transform(no, wx, wy, wz, t);
    }
  }

  // Compute residual displacements
  for (int idx = 0; idx < no; ++idx) {
    dx[idx] -= (wx[idx] - x[idx]);
    dy[idx] -= (wy[idx] - y[idx]);
    dz[idx] -= (wz[idx] - z[idx]);
  }

  // Free memory
  Deallocate(wx);
  Deallocate(wy);
  Deallocate(wz);

  // Approximate residual displacements by active levels
  FreeFormTransformation *ffd;
  for (int lvl = 0; lvl < NumberOfLevels(); ++lvl) {
    ffd = GetLocalTransformation(lvl);
    if (LocalTransformationIsActive(lvl)) {
      if (error > max_error) {
        error = ffd->Approximate(x, y, z, t, dx, dy, dz, no, niter, max_error);
      } else {
        ffd->Reset();
      }
    }
    ++lvl;
  }

  return error;
}

// -----------------------------------------------------------------------------
double MultiLevelFreeFormTransformation
::ApproximateAsNew(const ImageAttributes &domain,
                   double *dx, double *dy, double *dz,
                   int niter, double max_error)
{
  // Reset transformation
  this->Reset();

  // Check input arguments
  const int no = domain.NumberOfPoints();
  if (no <= 0) return .0;
  double error = numeric_limits<double>::infinity();
  if (niter < 1) return error;

  // Approximate global transformation
  AffineTransformation *aff = GetGlobalTransformation();
  error = aff->ApproximateAsNew(domain, dx, dy, dz, niter, max_error);

  // Approximate residual displacements by consecutive active levels
  FreeFormTransformation *ffd;
  for (int lvl = 0; lvl < NumberOfLevels() && error > max_error; ++lvl) {
    if (LocalTransformationIsActive(lvl)) {
      ffd   = GetLocalTransformation(lvl);
      error = ffd->ApproximateAsNew(domain, dx, dy, dz, niter, max_error);
    }
  }

  return error;
}

// -----------------------------------------------------------------------------
double MultiLevelFreeFormTransformation
::ApproximateAsNew(const double *x,  const double *y,  const double *z,
                   double       *dx, double       *dy, double       *dz, int no,
                   int niter, double max_error)
{
  // Reset transformation
  this->Reset();

  // Check input arguments
  if (no <= 0) return .0;
  double error = numeric_limits<double>::infinity();
  if (niter < 1) return error;

  // Approximate global transformation
  AffineTransformation *aff = GetGlobalTransformation();
  error = aff->ApproximateAsNew(x, y, z, dx, dy, dz, no, niter, max_error);

  // Approximate residual displacements by consecutive active levels
  FreeFormTransformation *ffd;
  for (int lvl = 0; lvl < NumberOfLevels() && error > max_error; ++lvl) {
    if (LocalTransformationIsActive(lvl)) {
      ffd   = GetLocalTransformation(lvl);
      error = ffd->ApproximateAsNew(x, y, z, dx, dy, dz, no, niter, max_error);
    }
  }

  return error;
}

// -----------------------------------------------------------------------------
double MultiLevelFreeFormTransformation
::ApproximateAsNew(const double *x,  const double *y,  const double *z,  const double *t,
                   double       *dx, double       *dy, double       *dz, int no,
                   int niter, double max_error)
{
  // Reset transformation
  this->Reset();

  // Check input arguments
  if (no <= 0) return .0;
  double error = numeric_limits<double>::infinity();
  if (niter < 1) return error;

  // Approximate global transformation
  AffineTransformation *aff = GetGlobalTransformation();
  error = aff->ApproximateAsNew(x, y, z, t, dx, dy, dz, no, niter, max_error);

  // Approximate residual displacements by consecutive active levels
  FreeFormTransformation *ffd;
  for (int lvl = 0; lvl < NumberOfLevels() && error > max_error; ++lvl) {
    if (LocalTransformationIsActive(lvl)) {
      ffd   = GetLocalTransformation(lvl);
      error = ffd->ApproximateAsNew(x, y, z, t, dx, dy, dz, no, niter, max_error);
    }
  }

  return error;
}

// -----------------------------------------------------------------------------
void MultiLevelFreeFormTransformation
::ApproximateDOFs(const double *x,  const double *y,  const double *z,  const double *t,
                  const double *dx, const double *dy, const double *dz, int no)
{
  // Reset active levels
  for (int lvl = 0; lvl < _NumberOfLevels; ++lvl) {
    if (LocalTransformationIsActive(lvl)) {
      GetLocalTransformation(lvl)->Reset();
    }
  }

  // Check input arguments
  if (no <= 0) return;

  // Transform world coordinates by global transformation and passive levels
  double *wx = Allocate<double>(no);
  double *wy = Allocate<double>(no);
  double *wz = Allocate<double>(no);
  memcpy(wx, x, no * sizeof(double));
  memcpy(wy, y, no * sizeof(double));
  memcpy(wz, z, no * sizeof(double));

  GetGlobalTransformation()->Transform(no, wx, wy, wz, t);
  for (int lvl = 0; lvl < NumberOfLevels(); ++lvl) {
    if (!LocalTransformationIsActive(lvl)) {
      GetLocalTransformation(lvl)->Transform(no, wx, wy, wz, t);
    }
  }

  // Compute residual displacements
  double *rx = Allocate<double>(no);
  double *ry = Allocate<double>(no);
  double *rz = Allocate<double>(no);
  for (int idx = 0; idx < no; ++idx) {
    rx[idx] = dx[idx] - (wx[idx] - x[idx]);
    ry[idx] = dy[idx] - (wy[idx] - y[idx]);
    rz[idx] = dz[idx] - (wz[idx] - z[idx]);
  }

  // Free memory
  Deallocate(wx);
  Deallocate(wy);
  Deallocate(wz);

  // Approximate residual displacements by active levels
  FreeFormTransformation *ffd;
  for (int lvl = 0; lvl < NumberOfLevels(); ++lvl) {
    ffd = GetLocalTransformation(lvl);
    if (LocalTransformationIsActive(lvl)) {
      ffd->ApproximateDOFs(x, y, z, t, rx, ry, rz, no);
      TransformationUtils::SubDisplacements sub(ffd, x, y, z, t, rx, ry, rz);
      parallel_for(blocked_range<int>(0, no), sub);
    }
  }

  // Free memory
  Deallocate(rx);
  Deallocate(ry);
  Deallocate(rz);
}

// -----------------------------------------------------------------------------
void MultiLevelFreeFormTransformation
::ApproximateDOFsGradient(const double *, const double *, const double *, const double *,
                          const double *, const double *, const double *, int,
                          double *, double) const
{
  cerr << this->NameOfClass() << "::ApproximateDOFsGradient: Not implemented" << endl;
  exit(1);
}

// =============================================================================
// Point transformation
// =============================================================================

// -----------------------------------------------------------------------------
void MultiLevelFreeFormTransformation
::LocalTransform(int m, int n, double &x, double &y, double &z, double t, double t0) const
{
  if (n < 0 || n > _NumberOfLevels) n = _NumberOfLevels;
  const int l1 = (m < 0 ? 0 : m);

  // Sum displacements of local transformations
  double u, v, w, dx = 0, dy = 0, dz = 0;

  for (int l = l1; l < n; ++l) {
    u = x, v = y, w = z;
    _LocalTransformation[l]->Transform(u, v, w, t, t0);
    dx += (u - x);
    dy += (v - y);
    dz += (w - z);
  }

  // Apply local displacements
  x += dx, y += dy, z += dz;
}

// -----------------------------------------------------------------------------
void MultiLevelFreeFormTransformation
::Transform(int m, int n, double &x, double &y, double &z, double t, double t0) const
{
  if (n < 0 || n > _NumberOfLevels) n = _NumberOfLevels;
  const int l1 = (m < 0 ? 0 : m);

  // Sum displacements of local transformations
  double u, v, w, dx = 0, dy = 0, dz = 0;

  for (int l = l1; l < n; ++l) {
    u = x, v = y, w = z;
    _LocalTransformation[l]->Transform(u, v, w, t, t0);
    dx += (u - x);
    dy += (v - y);
    dz += (w - z);
  }

  // Global transformation
  if (m < 0) _GlobalTransformation.Transform(x, y, z, t, t0);

  // Apply local displacements
  x += dx, y += dy, z += dz;
}

// -----------------------------------------------------------------------------
void MultiLevelFreeFormTransformation
::Displacement(int m, int n, GenericImage<double> &disp, double t, double t0, const WorldCoordsImage *wc) const
{
  if (!this->RequiresCachingOfDisplacements()) {
    MultiLevelTransformation::Displacement(m, n, disp, t, t0, wc);
    return;
  }

  if (n < 0 || n > _NumberOfLevels) n = _NumberOfLevels;
  const int l1 = (m < 0 ? 0 : m);

  if (disp.T() < 2 || disp.T() > 3) {
    cerr << "MultiLevelFreeFormTransformation::Displacement: Input/output image must have either 2 or 3 vector components (_t)" << endl;
    exit(1);
  }

  // Sum displacements of local transformations
  GenericImage<double> tmp, local(disp.GetImageAttributes());
  for (int l = l1; l < n; ++l) {
    tmp    = disp;
    _LocalTransformation[l]->Displacement(tmp, t, t0, wc);
    local += tmp;
  }

  // Displacement of global transformation
  if (m < 0) _GlobalTransformation.Displacement(disp, t, t0, wc);

  // Add local displacements
  disp += local;
}

// -----------------------------------------------------------------------------
void MultiLevelFreeFormTransformation
::Displacement(int m, int n, GenericImage<float> &disp, double t, double t0, const WorldCoordsImage *wc) const
{
  if (!this->RequiresCachingOfDisplacements()) {
    MultiLevelTransformation::Displacement(m, n, disp, t, t0, wc);
    return;
  }

  if (n < 0 || n > _NumberOfLevels) n = _NumberOfLevels;
  const int l1 = (m < 0 ? 0 : m);

  if (disp.T() < 2 || disp.T() > 3) {
    cerr << "MultiLevelFreeFormTransformation::Displacement: Input/output image must have either 2 or 3 vector components (_t)" << endl;
    exit(1);
  }

  // Sum displacements of local transformations
  GenericImage<float> tmp, local(disp.GetImageAttributes());
  for (int l = l1; l < n; ++l) {
    tmp    = disp;
    _LocalTransformation[l]->Displacement(tmp, t, t0, wc);
    local += tmp;
  }

  // Displacement of global transformation
  if (m < 0) _GlobalTransformation.Displacement(disp, t, t0, wc);

  // Add local displacements
  disp += local;
}

// -----------------------------------------------------------------------------
bool MultiLevelFreeFormTransformation::CanModifyDisplacement(int dof) const
{
  if (dof < 0) {
    for (int l = 0; l < _NumberOfLevels; ++l) {
      if (!_LocalTransformation[l]->CanModifyDisplacement(-1)) return false;
    }
    return true;
  } else {
    const FreeFormTransformation *ffd;
    DOFIndexToLocalTransformation(this, dof, ffd, dof);
    return ffd->CanModifyDisplacement();
  }
}

// -----------------------------------------------------------------------------
void MultiLevelFreeFormTransformation
::DisplacementAfterDOFChange(int dof, double dv, GenericImage<double> &dx,
                             double t, double t0, const WorldCoordsImage *i2w) const
{
  const FreeFormTransformation *ffd;
  DOFIndexToLocalTransformation(this, dof, ffd, dof);
  ffd->DisplacementAfterDOFChange(dof, dv, dx, t, t0, i2w);
}

// -----------------------------------------------------------------------------
int MultiLevelFreeFormTransformation
::InverseDisplacement(int m, int n, GenericImage<double> &disp, double t, double t0, const WorldCoordsImage *wc) const
{
  if ((m < 0 && !_GlobalTransformation.IsIdentity()) || !this->RequiresCachingOfDisplacements()) {
    return MultiLevelTransformation::InverseDisplacement(m, n, disp, t, t0, wc);
  }

  if (n < 0 || n > _NumberOfLevels) n = _NumberOfLevels;
  const int l1 = (m < 0 ? 0 : m);

  if (disp.T() < 2 || disp.T() > 3) {
    cerr << "MultiLevelFreeFormTransformation::InverseDisplacement: Input/output image must have either 2 or 3 vector components (_t)" << endl;
    exit(1);
  }

  int ninv = 0;

  GenericImage<double> local;
  for (int l = l1; l < n; ++l) {
    local = disp;
    ninv += _LocalTransformation[l]->InverseDisplacement(local, t, t0, wc);
    disp += local;
  }

  return ninv;
}

// -----------------------------------------------------------------------------
int MultiLevelFreeFormTransformation
::InverseDisplacement(int m, int n, GenericImage<float> &disp, double t, double t0, const WorldCoordsImage *wc) const
{
  if ((m < 0 && !_GlobalTransformation.IsIdentity()) || !this->RequiresCachingOfDisplacements()) {
    return MultiLevelTransformation::InverseDisplacement(m, n, disp, t, t0, wc);
  }

  if (n < 0 || n > _NumberOfLevels) n = _NumberOfLevels;
  const int l1 = (m < 0 ? 0 : m);

  if (disp.T() < 2 || disp.T() > 3) {
    cerr << "MultiLevelFreeFormTransformation::InverseDisplacement: Input/output image must have either 2 or 3 vector components (_t)" << endl;
    exit(1);
  }

  int ninv = 0;

  GenericImage<float> local;
  for (int l = l1; l < n; ++l) {
    local = disp;
    ninv += _LocalTransformation[l]->InverseDisplacement(local, t, t0, wc);
    disp += local;
  }

  return ninv;
}

// =============================================================================
// Derivatives
// =============================================================================

// -----------------------------------------------------------------------------
void MultiLevelFreeFormTransformation
::Jacobian(int m, int n, Matrix &jac, double x, double y, double z, double t, double t0) const
{
  Matrix tmp(3, 3);

  if (n < 0 || n > _NumberOfLevels) n = _NumberOfLevels;
  const int l1 = (m < 0 ? 0 : m);

  // Compute global jacobian
  if (m < 0) {
    _GlobalTransformation.Jacobian(jac, x, y, z, t, t0);
  } else {
    jac.Initialize(3, 3);
    jac.Ident();
  }

  // Compute local jacobian
  for (int l = l1; l < n; ++l) {

    // Calculate jacobian
    _LocalTransformation[l]->Jacobian(tmp, x, y, z, t, t0);

    // Subtract identity matrix
    tmp(0, 0) -= 1;
    tmp(1, 1) -= 1;
    tmp(2, 2) -= 1;

    // Add jacobian
    jac += tmp;
  }
}

// -----------------------------------------------------------------------------
void MultiLevelFreeFormTransformation
::Hessian(int m, int n, Matrix hessian[3], double x, double y, double z, double t, double t0) const
{
  Matrix tmp[3];

  if (n < 0 || n > _NumberOfLevels) n = _NumberOfLevels;
  const int l1 = (m < 0 ? 0 : m);

  // Compute 2nd order derivatives of global transformation
  if (m < 0) {
    _GlobalTransformation.Hessian(hessian, x, y, z, t, t0);
  } else {
    hessian[0].Initialize(3, 3);
    hessian[1].Initialize(3, 3);
    hessian[2].Initialize(3, 3);
  }

  // Compute 2nd order derivatives of local transformations
  for (int l = l1; l < n; ++l) {

    // Calculate 2nd order derivatives
    _LocalTransformation[l]->Hessian(tmp, x, y, z, t, t0);

    // Add derivatives
    hessian[0] += tmp[0];
    hessian[1] += tmp[1];
    hessian[2] += tmp[2];
  }
}

// -----------------------------------------------------------------------------
void MultiLevelFreeFormTransformation
::DeriveJacobianWrtDOF(Matrix &jac, int dof, double x, double y, double z, double t, double t0) const
{
  const FreeFormTransformation *ffd;
  int                           ffd_dof;
  DOFIndexToLocalTransformation(this, dof, ffd, ffd_dof);
  return ffd->DeriveJacobianWrtDOF(jac, ffd_dof, x, y, z, t, t0);
}

// -----------------------------------------------------------------------------
void MultiLevelFreeFormTransformation
::ParametricGradient(const GenericImage<double> *in, double *out,
                     const WorldCoordsImage *i2w, const WorldCoordsImage *wc,
                     double t0, double w) const
{
  for (int l = 0; l < this->NumberOfLevels(); ++l) {
    if (this->LocalTransformationIsActive(l)) {
      const FreeFormTransformation *ffd = this->GetLocalTransformation(l);
      ffd->ParametricGradient(in, out, i2w, wc, t0, w);
      out += ffd->NumberOfDOFs();
    }
  }
}

// -----------------------------------------------------------------------------
void MultiLevelFreeFormTransformation
::ParametricGradient(const PointSet &pos, const Vector3D<double> *in,
                     double *out, double t, double t0, double w) const
{
  for (int l = 0; l < this->NumberOfLevels(); ++l) {
    if (this->LocalTransformationIsActive(l)) {
      const FreeFormTransformation *ffd = this->GetLocalTransformation(l);
      ffd->ParametricGradient(pos, in, out, t, t0, w);
      out += ffd->NumberOfDOFs();
    }
  }
}

// =============================================================================
// I/O
// =============================================================================

// -----------------------------------------------------------------------------
void MultiLevelFreeFormTransformation::Print(ostream &os, Indent indent) const
{
  os << indent << "Multi-level FFD:" << endl;
  MultiLevelTransformation::Print(os, indent + 1);
}


} // namespace mirtk
