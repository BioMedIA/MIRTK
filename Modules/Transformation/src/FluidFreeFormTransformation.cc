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

#include "mirtk/FluidFreeFormTransformation.h"

#include "mirtk/Memory.h"


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
FluidFreeFormTransformation::FluidFreeFormTransformation()
{
}

// -----------------------------------------------------------------------------
FluidFreeFormTransformation::FluidFreeFormTransformation(const RigidTransformation &t)
:
  MultiLevelTransformation(t)
{
}

// -----------------------------------------------------------------------------
FluidFreeFormTransformation::FluidFreeFormTransformation(const AffineTransformation &t)
:
  MultiLevelTransformation(t)
{
}

// -----------------------------------------------------------------------------
FluidFreeFormTransformation::FluidFreeFormTransformation(const FluidFreeFormTransformation &t)
:
  MultiLevelTransformation(t)
{
}

// -----------------------------------------------------------------------------
FluidFreeFormTransformation::~FluidFreeFormTransformation()
{
}

// =============================================================================
// Levels
// =============================================================================

// -----------------------------------------------------------------------------
void FluidFreeFormTransformation::CheckTransformation(FreeFormTransformation *transformation) const
{
  MultiLevelTransformation::CheckTransformation(transformation);
}

// -----------------------------------------------------------------------------
void FluidFreeFormTransformation::CombineLocalTransformation()
{
  cerr << "FluidFreeFormTransformation::CombineLocalTransformation: Does not make sense" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
void FluidFreeFormTransformation::MergeGlobalIntoLocalDisplacement()
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

  // Shift all existing FFDs in the stack down by one place
  for (int l = this->NumberOfLevels(); l > 0; --l) {
  	this->PutLocalTransformation(this->GetLocalTransformation(l-1), l);
  }

  // Insert the new FFD for the global transformation into the first position
  this->PutLocalTransformation(ffdCopy, 0);

  // Reset matrix previously used for global transform to identity
  _GlobalTransformation.Reset();
}

// =============================================================================
// Approximation
// =============================================================================

// -----------------------------------------------------------------------------
double FluidFreeFormTransformation
::Approximate(const ImageAttributes &domain,
              double *dx, double *dy, double *dz,
              int niter, double max_error)
{
  // Check input arguments
  const int no = domain.NumberOfPoints();
  if (no <= 0) return .0;
  double error = numeric_limits<double>::infinity();
  if (niter < 1) return error;

  int lvl = 0;

  // Compute world coordinates of lattice points
  double *x = Allocate<double>(no);
  double *y = Allocate<double>(no);
  double *z = Allocate<double>(no);
  double *t = Allocate<double>(no);
  domain.LatticeToWorld(x, y, z, t);

  // Initialize transformed world coordinates
  double *wx = Allocate<double>(no);
  double *wy = Allocate<double>(no);
  double *wz = Allocate<double>(no);
  memcpy(wx, x, no * sizeof(double));
  memcpy(wy, y, no * sizeof(double));
  memcpy(wz, z, no * sizeof(double));

  // Skip global transformation and initial passive levels
  GetGlobalTransformation()->Transform(no, wx, wy, wz, t);
  while (lvl < NumberOfLevels()) {
    if (LocalTransformationIsActive(lvl)) break;
    GetLocalTransformation(lvl)->Transform(no, wx, wy, wz, t);
    ++lvl;
  }

  // Compute residual displacements
  for (int idx = 0; idx < no; ++idx) {
    dx[idx] -= (wx[idx] - x[idx]);
    dy[idx] -= (wy[idx] - y[idx]);
    dz[idx] -= (wz[idx] - z[idx]);
  }

  // Free memory
  Deallocate(x);
  Deallocate(y);
  Deallocate(z);

  // Approximate residual displacements by active levels,
  // resetting any passive transformations in between active levels
  FreeFormTransformation *ffd;
  while (lvl < NumberOfLevels()) {
    ffd = GetLocalTransformation(lvl);
    if (LocalTransformationIsActive(lvl)) {
      error = ffd->Approximate(wx, wy, wz, t, dx, dy, dz, no, niter, max_error);
      if (error <= max_error) break;
      ffd->Transform(no, wx, wy, wz, t);
    } else {
      ffd->Reset();
    }
    ++lvl;
  }

  // Free memory
  Deallocate(wx);
  Deallocate(wy);
  Deallocate(wz);
  Deallocate(t);

  return error;
}

// -----------------------------------------------------------------------------
double FluidFreeFormTransformation
::Approximate(const double *x,  const double *y,  const double *z,
              double       *dx, double       *dy, double       *dz, int no,
              int niter, double max_error)
{
  // Check input arguments
  if (no <= 0) return .0;
  double error = numeric_limits<double>::infinity();
  if (niter < 1) return error;

  int lvl = 0;

  // Initialize transformed world coordinates
  double *wx = Allocate<double>(no);
  double *wy = Allocate<double>(no);
  double *wz = Allocate<double>(no);
  memcpy(wx, x, no * sizeof(double));
  memcpy(wy, y, no * sizeof(double));
  memcpy(wz, z, no * sizeof(double));

  // Skip global transformation and initial passive levels
  GetGlobalTransformation()->Transform(no, wx, wy, wz);
  while (lvl < NumberOfLevels()) {
    if (LocalTransformationIsActive(lvl)) break;
    GetLocalTransformation(lvl)->Transform(no, wx, wy, wz);
    ++lvl;
  }

  // Compute residual displacements
  for (int idx = 0; idx < no; ++idx) {
    dx[idx] -= (wx[idx] - x[idx]);
    dy[idx] -= (wy[idx] - y[idx]);
    dz[idx] -= (wz[idx] - z[idx]);
  }

  // Approximate residual displacements by active levels,
  // resetting any passive transformations in between active levels
  FreeFormTransformation *ffd;
  while (lvl < NumberOfLevels()) {
    ffd = GetLocalTransformation(lvl);
    if (LocalTransformationIsActive(lvl)) {
      error = ffd->Approximate(wx, wy, wz, dx, dy, dz, no, niter, max_error);
      if (error <= max_error) break;
      ffd->Transform(no, wx, wy, wz);
    } else {
      ffd->Reset();
    }
    ++lvl;
  }

  // Free memory
  Deallocate(wx);
  Deallocate(wy);
  Deallocate(wz);

  return error;
}

// -----------------------------------------------------------------------------
double FluidFreeFormTransformation
::Approximate(const double *x,  const double *y,  const double *z,  const double *t,
              double       *dx, double       *dy, double       *dz, int no,
              int niter, double max_error)
{
  // Check input arguments
  if (no <= 0) return .0;
  double error = numeric_limits<double>::infinity();
  if (niter < 1) return error;

  int lvl = 0;

  // Initialize transformed world coordinates
  double *wx = Allocate<double>(no);
  double *wy = Allocate<double>(no);
  double *wz = Allocate<double>(no);
  memcpy(wx, x, no * sizeof(double));
  memcpy(wy, y, no * sizeof(double));
  memcpy(wz, z, no * sizeof(double));

  // Skip global transformation and initial passive levels
  GetGlobalTransformation()->Transform(no, wx, wy, wz, t);
  while (lvl < NumberOfLevels()) {
    if (LocalTransformationIsActive(lvl)) break;
    GetLocalTransformation(lvl)->Transform(no, wx, wy, wz, t);
    ++lvl;
  }

  // Compute residual displacements
  for (int idx = 0; idx < no; ++idx) {
    dx[idx] -= (wx[idx] - x[idx]);
    dy[idx] -= (wy[idx] - y[idx]);
    dz[idx] -= (wz[idx] - z[idx]);
  }

  // Approximate residual displacements by active levels,
  // resetting any passive transformations in between active levels
  FreeFormTransformation *ffd;
  while (lvl < NumberOfLevels()) {
    ffd = GetLocalTransformation(lvl);
    if (LocalTransformationIsActive(lvl)) {
      error = ffd->Approximate(wx, wy, wz, t, dx, dy, dz, no, niter, max_error);
      if (error <= max_error) break;
      ffd->Transform(no, wx, wy, wz, t);
    } else {
      ffd->Reset();
    }
    ++lvl;
  }

  // Free memory
  Deallocate(wx);
  Deallocate(wy);
  Deallocate(wz);

  return error;
}

// -----------------------------------------------------------------------------
double FluidFreeFormTransformation
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

  // Compute world coordinates of lattice points
  double *wx = Allocate<double>(no);
  double *wy = Allocate<double>(no);
  double *wz = Allocate<double>(no);
  double *t  = Allocate<double>(no);
  domain.LatticeToWorld(wx, wy, wz, t);

  // Approximate global transformation
  AffineTransformation *aff = GetGlobalTransformation();
  error = aff->ApproximateAsNew(wx, wy, wz, t, dx, dy, dz, no, niter, max_error);

  // Approximate residual displacements by consecutive active levels
  if (error > max_error) {
    // Transform world coordinates
    aff->Transform(no, wx, wy, wz, t);
    FreeFormTransformation *ffd;
    for (int lvl = 0; lvl < NumberOfLevels(); ++lvl) {
      if (LocalTransformationIsActive(lvl)) {
        ffd = GetLocalTransformation(lvl);
        // Transform residual displacements
        error = ffd->ApproximateAsNew(wx, wy, wz, t, dx, dy, dz, no, niter, max_error);
        if (error <= max_error) break;
        // Transform world coordinates
        ffd->Transform(no, wx, wy, wz, t);
      }
    }
  }

  // Free memory
  Deallocate(wx);
  Deallocate(wy);
  Deallocate(wz);
  Deallocate(t);

  return error;
}

// -----------------------------------------------------------------------------
double FluidFreeFormTransformation
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
  if (error <= max_error) return error;

  // Transform world coordinates
  double *wx = Allocate<double>(no);
  double *wy = Allocate<double>(no);
  double *wz = Allocate<double>(no);
  memcpy(wx, x, no * sizeof(double));
  memcpy(wy, y, no * sizeof(double));
  memcpy(wz, z, no * sizeof(double));
  aff->Transform(no, wx, wy, wz);

  // Approximate residual displacements by consecutive active levels
  FreeFormTransformation *ffd;
  for (int lvl = 0; lvl < NumberOfLevels(); ++lvl) {
    if (LocalTransformationIsActive(lvl)) {
      ffd = GetLocalTransformation(lvl);
      // Transform residual displacements
      error = ffd->ApproximateAsNew(wx, wy, wz, dx, dy, dz, no, niter, max_error);
      if (error <= max_error) break;
      // Transform world coordinates
      ffd->Transform(no, wx, wy, wz);
    }
  }

  // Free memory
  Deallocate(wx);
  Deallocate(wy);
  Deallocate(wz);

  return error;
}

// -----------------------------------------------------------------------------
double FluidFreeFormTransformation
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
  if (error <= max_error) return error;

  // Transform world coordinates
  double *wx = Allocate<double>(no);
  double *wy = Allocate<double>(no);
  double *wz = Allocate<double>(no);
  memcpy(wx, x, no * sizeof(double));
  memcpy(wy, y, no * sizeof(double));
  memcpy(wz, z, no * sizeof(double));
  aff->Transform(no, wx, wy, wz, t);

  // Approximate residual displacements by consecutive active levels
  FreeFormTransformation *ffd;
  for (int lvl = 0; lvl < NumberOfLevels(); ++lvl) {
    if (LocalTransformationIsActive(lvl)) {
      ffd = GetLocalTransformation(lvl);
      // Transform residual displacements
      error = ffd->ApproximateAsNew(wx, wy, wz, t, dx, dy, dz, no, niter, max_error);
      if (error <= max_error) break;
      // Transform world coordinates
      ffd->Transform(no, wx, wy, wz, t);
    }
  }

  // Free memory
  Deallocate(wx);
  Deallocate(wy);
  Deallocate(wz);

  return error;
}

// -----------------------------------------------------------------------------
void FluidFreeFormTransformation
::ApproximateDOFs(const double *x,  const double *y,  const double *z,  const double *t,
                  const double *dx, const double *dy, const double *dz, int no)
{
  // Initialize transformed world coordinates
  double *wx = Allocate<double>(no);
  double *wy = Allocate<double>(no);
  double *wz = Allocate<double>(no);
  memcpy(wx, x, no * sizeof(double));
  memcpy(wy, y, no * sizeof(double));
  memcpy(wz, z, no * sizeof(double));

  // Allocate memory for residual displacements
  double *rx = Allocate<double>(no);
  double *ry = Allocate<double>(no);
  double *rz = Allocate<double>(no);

  // Transform points by global transformation and initial local passive
  // transformations whose parameters are fixed
  int lvl = 0;
  GetGlobalTransformation()->Transform(no, wx, wy, wz, t);
  while (lvl < NumberOfLevels()) {
    if (!LocalTransformationIsActive(lvl)) {
      GetLocalTransformation(lvl)->Transform(no, wx, wy, wz, t);
    }
    ++lvl;
  }

  // Compute residual displacements
  for (int idx = 0; idx < no; ++idx) {
    rx[idx] = dx[idx] - (wx[idx] - x[idx]);
    ry[idx] = dy[idx] - (wy[idx] - y[idx]);
    rz[idx] = dz[idx] - (wz[idx] - z[idx]);
  }

  // Approximate displacements by consecutive levels,
  // resetting any passive transformations in between active levels
  FreeFormTransformation *ffd;
  while (lvl < NumberOfLevels()) {
    ffd = GetLocalTransformation(lvl);
    if (LocalTransformationIsActive(lvl)) {
      ffd->ApproximateDOFs(wx, wy, wz, t, rx, ry, rz, no);
      ffd->Transform(no, wx, wy, wz, t);
      for (int idx = 0; idx < no; ++idx) {
        rx[idx] = dx[idx] - (wx[idx] - x[idx]);
        ry[idx] = dy[idx] - (wy[idx] - y[idx]);
        rz[idx] = dz[idx] - (wz[idx] - z[idx]);
      }
    } else {
      ffd->Reset();
    }
    ++lvl;
  }

  // Free memory
  Deallocate(wx);
  Deallocate(wy);
  Deallocate(wz);
  Deallocate(rx);
  Deallocate(ry);
  Deallocate(rz);
}

// -----------------------------------------------------------------------------
void FluidFreeFormTransformation
::ApproximateDOFsGradient(const double *, const double *, const double *, const double *,
                          const double *, const double *, const double *, int,
                          double *, double) const
{
  cerr << this->NameOfClass() << "::ApproximateDOFsGradient: Not implemented" << endl;
  exit(1);
}

// =============================================================================
// Transformation parameters (DoFs)
// =============================================================================

// -----------------------------------------------------------------------------
bool FluidFreeFormTransformation::CopyFrom(const Transformation *other)
{
  bool ok = MultiLevelTransformation::CopyFrom(other);
  const FluidFreeFormTransformation *fluid;
  if (ok && (fluid = dynamic_cast<const FluidFreeFormTransformation *>(other))) {
    ok = _AffineTransformation.CopyFrom(fluid->GetAffineTransformation());
  }
  return ok;
}

// -----------------------------------------------------------------------------
bool FluidFreeFormTransformation::IsIdentity() const
{
  return _AffineTransformation.IsIdentity() && MultiLevelTransformation::IsIdentity();
}

// -----------------------------------------------------------------------------
void FluidFreeFormTransformation::Reset()
{
  MultiLevelTransformation::Reset();
  _AffineTransformation.Reset();
}

// -----------------------------------------------------------------------------
void FluidFreeFormTransformation::Clear()
{
  MultiLevelTransformation::Clear();
  _AffineTransformation.Reset();
}

// =============================================================================
// Point transformation
// =============================================================================

// -----------------------------------------------------------------------------
void FluidFreeFormTransformation::Transform(int m, int n, double &x, double &y, double &z, double, double) const
{
  // Global transformation
  if (m < 0) _GlobalTransformation.Transform(x, y, z);

  // Local transformations
  const int l1 = (m < 0 ? 0 : m);
  const int l2 = (n < 0 || n > _NumberOfLevels ? _NumberOfLevels : n);
  for (int l = l1; l < l2; ++l) _LocalTransformation[l]->Transform(x, y, z);

  // Affine post-transformation
  if (n < 0 || n > _NumberOfLevels) _AffineTransformation.Transform(x, y, z);
}

// -----------------------------------------------------------------------------
bool FluidFreeFormTransformation::Inverse(int m, int n, double &x, double &y, double &z, double, double) const
{
  bool ok = true;

  // Inverse affine post-transformation
  if (n < 0 || n > _NumberOfLevels) {
    _AffineTransformation.Inverse(x, y, z);
  }

  // Inverse local transformations
  const int l1 = (m < 0 ? 0 : m);
  const int l2 = (n < 0 || n > _NumberOfLevels ? _NumberOfLevels : n);

  for (int l = l2 - 1; l >= l1; --l) {
    ok = _LocalTransformation[l]->Inverse(x, y, z) && ok;
  }

  // Inverse global transformation
  if (m < 0) _GlobalTransformation.Inverse(x, y, z);
  return ok;
}

// -----------------------------------------------------------------------------
void FluidFreeFormTransformation::Displacement(int m, int n, GenericImage<double> &disp, double t, double t0, const WorldCoordsImage *wc) const
{
  // Call base class function if caching not required
  if (!this->RequiresCachingOfDisplacements()) {
    MultiLevelTransformation::Displacement(m, n, disp, t, t0, wc);
    return;
  }

  // Check arguments
  if (disp.GetT() < 2 || disp.GetT() > 3) {
    cerr << "FluidFreeFormTransformation::Displacement: Input/output image must have either 2 or 3 vector components (_t)" << endl;
    exit(1);
  }

  // Displacement of global transformation
  if (m < 0) _GlobalTransformation.Displacement(disp, t, t0, wc);

  // Compose with displacements of local transformations
  const int l1 = (m < 0 ? 0 : m);
  const int l2 = (n < 0 || n > _NumberOfLevels ? _NumberOfLevels : n);

  for (int l = l1; l < l2; ++l) {
    _LocalTransformation[l]->Displacement(disp, t, t0, wc);
  }

  // Compose with affine post-transformation
  if (n < 0 || n > _NumberOfLevels) {
    _AffineTransformation.Displacement(disp, t, t0, wc);
  }
}

// -----------------------------------------------------------------------------
void FluidFreeFormTransformation::Displacement(int m, int n, GenericImage<float> &disp, double t, double t0, const WorldCoordsImage *wc) const
{
  // Call base class function if caching not required
  if (!this->RequiresCachingOfDisplacements()) {
    MultiLevelTransformation::Displacement(m, n, disp, t, t0, wc);
    return;
  }

  // Check arguments
  if (disp.GetT() < 2 || disp.GetT() > 3) {
    cerr << "FluidFreeFormTransformation::Displacement: Input/output image must have either 2 or 3 vector components (_t)" << endl;
    exit(1);
  }

  // Displacement of global transformation
  if (m < 0) _GlobalTransformation.Displacement(disp, t, t0, wc);

  // Compose with displacements of local transformations
  const int l1 = (m < 0 ? 0 : m);
  const int l2 = (n < 0 || n > _NumberOfLevels ? _NumberOfLevels : n);

  for (int l = l1; l < l2; ++l) {
    _LocalTransformation[l]->Displacement(disp, t, t0, wc);
  }

  // Compose with affine post-transformation
  if (n < 0 || n > _NumberOfLevels) {
    _AffineTransformation.Displacement(disp, t, t0, wc);
  }
}

// -----------------------------------------------------------------------------
int FluidFreeFormTransformation::InverseDisplacement(int m, int n, GenericImage<double> &disp, double t, double t0, const WorldCoordsImage *wc) const
{
  // Call base class function if caching not required
  if (!this->RequiresCachingOfDisplacements()) {
    return MultiLevelTransformation::InverseDisplacement(m, n, disp, t, t0, wc);
  }

  // Check arguments
  if (disp.GetT() < 2 || disp.GetT() > 3) {
    cerr << "FluidFreeFormTransformation::InverseDisplacement: Input/output image must have either 2 or 3 vector components (_t)" << endl;
    exit(1);
  }

  // Evaluate displacements of inverse affine post-transformation
  if (n < 0 || n > _NumberOfLevels) {
    _AffineTransformation.InverseDisplacement(disp, t, t0, wc);
  }

  // Compose inverse displacements of local transformations
  const int l1 = (m < 0 ? 0 : m);
  const int l2 = (n < 0 || n > _NumberOfLevels ? _NumberOfLevels : n);

  int ninv = 0;
  for (int l = l2-1; l >= l1; --l) {
    ninv += _LocalTransformation[l]->InverseDisplacement(disp, t, t0, wc);
  }

  // Compose with inverse displacement of global transformation
  if (m < 0) _GlobalTransformation.InverseDisplacement(disp, t, t0, wc);

  return ninv;
}

// -----------------------------------------------------------------------------
int FluidFreeFormTransformation::InverseDisplacement(int m, int n, GenericImage<float> &disp, double t, double t0, const WorldCoordsImage *wc) const
{
  // Call base class function if caching not required
  if (!this->RequiresCachingOfDisplacements()) {
    return MultiLevelTransformation::InverseDisplacement(m, n, disp, t, t0, wc);
  }

  // Check arguments
  if (disp.GetT() < 2 || disp.GetT() > 3) {
    cerr << "FluidFreeFormTransformation::InverseDisplacement: Input/output image must have either 2 or 3 vector components (_t)" << endl;
    exit(1);
  }

  // Evaluate displacements of inverse affine post-transformation
  if (n < 0 || n > _NumberOfLevels) {
    _AffineTransformation.InverseDisplacement(disp, t, t0, wc);
  }

  // Compose inverse displacements of local transformations
  const int l1 = (m < 0 ? 0 : m);
  const int l2 = (n < 0 || n > _NumberOfLevels ? _NumberOfLevels : n);

  int ninv = 0;
  for (int l = l2-1; l >= l1; --l) {
    ninv += _LocalTransformation[l]->InverseDisplacement(disp, t, t0, wc);
  }

  // Compose with inverse displacement of global transformation
  if (m < 0) _GlobalTransformation.InverseDisplacement(disp, t, t0, wc);

  return ninv;
}

// =============================================================================
// Derivatives
// =============================================================================

// -----------------------------------------------------------------------------
void FluidFreeFormTransformation::Jacobian(int m, int n, Matrix &jac, double x, double y, double z, double, double) const
{
  // Initialize to identity matrix
  jac.Initialize(3, 3);
  jac.Ident();
  if (0 <= n && n <= m) return;

  // Compute 1st order derivatives of global transformation
  if (m < 0) {
    _GlobalTransformation.Jacobian(jac, x, y, z);
    _GlobalTransformation.Transform(x, y, z);
  }
  if (n == 0) return;

  // Recursively compute Jacobian of composite transformation (n levels)
  Matrix tmp(3, 3);

  const int l1 = (m < 0 ? 0 : m);
  const int l2 = (n < 0 || n > _NumberOfLevels ? _NumberOfLevels : n);

  for (int l = l1; l < l2; ++l) {
    _LocalTransformation[l]->Jacobian(tmp, x, y, z);
    _LocalTransformation[l]->Transform(x, y, z);
    jac = tmp * jac;
  }

  // Multiply by 1st order derivatives of affine post-transformation
  if (n < 0 || n > _NumberOfLevels) {
    _AffineTransformation.Jacobian(tmp, x, y, z);
    jac = tmp * jac;
  }
}

// -----------------------------------------------------------------------------
void FluidFreeFormTransformation::Hessian(int m, int n, Matrix hessian[3], double x, double y, double z, double, double) const
{
  // Initialize matrices
  hessian[0].Initialize(3, 3);
  hessian[1].Initialize(3, 3);
  hessian[2].Initialize(3, 3);
  if (0 <= n && n <= m) return;

  Matrix jac(3, 3);
  jac.Ident();

  // Compute 1st and 2nd order derivatives of global transformation
  if (m < 0) {
    _GlobalTransformation.Hessian(hessian, x, y, z);
    _GlobalTransformation.Jacobian(jac, x, y, z);
    _GlobalTransformation.Transform(x, y, z);
  }
  if (n == 0) return;

  // Recursively compute Hessian of composite transformation
  Matrix localjac, localhessian[3];
  double dxx[3], dxy[3], dxz[3], dyy[3], dyz[3], dzz[3];

  const int l1 = (m < 0 ? 0 : m);
  const int l2 = (n < 0 || n > _NumberOfLevels ? _NumberOfLevels : n);

  for (int l = l1; l < l2; ++l) {

    // Compute 1st and 2nd order derivatives of i-th level
    _LocalTransformation[l]->Jacobian(localjac,     x, y, z);
    _LocalTransformation[l]->Hessian (localhessian, x, y, z);

    // Compute 2nd order derivatives of composite transformation
    for (int dim = 0; dim < 3; ++dim) {
      // Compute 1st term of the product rule...
      localhessian[dim] = localhessian[dim] * jac; localhessian[dim].Transpose();
      localhessian[dim] = localhessian[dim] * jac;
      // ...and add 2nd term of the product rule
      dxx[dim] = localhessian[dim](0, 0) + localjac(dim, 0) * hessian[0](0, 0)
                                         + localjac(dim, 1) * hessian[1](0, 0)
                                         + localjac(dim, 2) * hessian[2](0, 0);
      dxy[dim] = localhessian[dim](0, 1) + localjac(dim, 0) * hessian[0](0, 1)
                                         + localjac(dim, 1) * hessian[1](0, 1)
                                         + localjac(dim, 2) * hessian[2](0, 1);
      dxz[dim] = localhessian[dim](0, 2) + localjac(dim, 0) * hessian[0](0, 2)
                                         + localjac(dim, 1) * hessian[1](0, 2)
                                         + localjac(dim, 2) * hessian[2](0, 2);
      dyy[dim] = localhessian[dim](1, 1) + localjac(dim, 0) * hessian[0](1, 1)
                                         + localjac(dim, 1) * hessian[1](1, 1)
                                         + localjac(dim, 2) * hessian[2](1, 1);
      dyz[dim] = localhessian[dim](1, 2) + localjac(dim, 0) * hessian[0](1, 2)
                                         + localjac(dim, 1) * hessian[1](1, 2)
                                         + localjac(dim, 2) * hessian[2](1, 2);
      dzz[dim] = localhessian[dim](2, 2) + localjac(dim, 0) * hessian[0](2, 2)
                                         + localjac(dim, 1) * hessian[1](2, 2)
                                         + localjac(dim, 2) * hessian[2](2, 2);
    }

    // Update derivatives of composite transformation
    jac = localjac * jac;
    for (int dim = 0; dim < 3; ++dim) {
      hessian[dim](0, 0) = dxx[dim];
      hessian[dim](0, 1) = hessian[dim](1, 0) = dxy[dim];
      hessian[dim](0, 2) = hessian[dim](2, 0) = dxz[dim];
      hessian[dim](1, 1) = dyy[dim];
      hessian[dim](1, 2) = hessian[dim](2, 1) = dyz[dim];
      hessian[dim](2, 2) = dzz[dim];
    }

    // Transform point
    _LocalTransformation[l]->Transform(x, y, z);
  }

  // Compute Hessian of composition with affine post-transformation
  if (n < 0 || n > _NumberOfLevels) {

    // Compute 1st and 2nd order derivatives of affine post-transformation
    _AffineTransformation.Jacobian(localjac,     x, y, z);
    _AffineTransformation.Hessian (localhessian, x, y, z);

    // Compute 2nd order derivatives of composite transformation
    for (int dim = 0; dim < 3; ++dim) {
      // Compute 1st term of the product rule...
      localhessian[dim] = localhessian[dim] * jac; localhessian[dim].Transpose();
      localhessian[dim] = localhessian[dim] * jac;
      // ...and add 2nd term of the product rule
      dxx[dim] = localhessian[dim](0, 0) + localjac(dim, 0) * hessian[0](0, 0)
                                         + localjac(dim, 1) * hessian[1](0, 0)
                                         + localjac(dim, 2) * hessian[2](0, 0);
      dxy[dim] = localhessian[dim](0, 1) + localjac(dim, 0) * hessian[0](0, 1)
                                         + localjac(dim, 1) * hessian[1](0, 1)
                                         + localjac(dim, 2) * hessian[2](0, 1);
      dxz[dim] = localhessian[dim](0, 2) + localjac(dim, 0) * hessian[0](0, 2)
                                         + localjac(dim, 1) * hessian[1](0, 2)
                                         + localjac(dim, 2) * hessian[2](0, 2);
      dyy[dim] = localhessian[dim](1, 1) + localjac(dim, 0) * hessian[0](1, 1)
                                         + localjac(dim, 1) * hessian[1](1, 1)
                                         + localjac(dim, 2) * hessian[2](1, 1);
      dyz[dim] = localhessian[dim](1, 2) + localjac(dim, 0) * hessian[0](1, 2)
                                         + localjac(dim, 1) * hessian[1](1, 2)
                                         + localjac(dim, 2) * hessian[2](1, 2);
      dzz[dim] = localhessian[dim](2, 2) + localjac(dim, 0) * hessian[0](2, 2)
                                         + localjac(dim, 1) * hessian[1](2, 2)
                                         + localjac(dim, 2) * hessian[2](2, 2);
    }

    // Update derivatives of composite transformation
    //jac = localjac * jac;
    for (int dim = 0; dim < 3; ++dim) {
      hessian[dim](0, 0) = dxx[dim];
      hessian[dim](0, 1) = hessian[dim](1, 0) = dxy[dim];
      hessian[dim](0, 2) = hessian[dim](2, 0) = dxz[dim];
      hessian[dim](1, 1) = dyy[dim];
      hessian[dim](1, 2) = hessian[dim](2, 1) = dyz[dim];
      hessian[dim](2, 2) = dzz[dim];
    }
  }
}

// -----------------------------------------------------------------------------
void FluidFreeFormTransformation
::ParametricGradient(const GenericImage<double> *in, double *out,
                     const WorldCoordsImage *i2w, const WorldCoordsImage *wc,
                     double t0, double w) const
{
  // Compute gradient only if any level is active
  if (this->NumberOfDOFs() == 0) return;

  // Initialize copy of world coordinates
  WorldCoordsImage pos;
  if      (wc ) pos = *wc;
  else if (i2w) pos = *i2w;
  else in->ImageToWorld(pos);

  // Transform world coordinates by global transformation
  this->GetGlobalTransformation()->Transform(pos, t0);

  int previous = -1;
  for (;;) {
    // Determine next active level
    int active = previous + 1;
    while (active < this->NumberOfLevels() && !this->LocalTransformationIsActive(active)) {
      ++active;
    }
    if (active == this->NumberOfLevels()) break;

    // Transform world coordinates by skipped passive local transformation
    for (int l = previous; l < active; ++l) {
      if (l == -1) continue;
      this->GetLocalTransformation(l)->Transform(pos, t0);
    }

    // Add gradient w.r.t. DoFs of this active level
    const FreeFormTransformation *ffd = this->GetLocalTransformation(active);
    ffd->ParametricGradient(in, out, i2w, &pos, t0, w);
    out += ffd->NumberOfDOFs();

    // Continue with next active level
    previous = active;
  }
}

// -----------------------------------------------------------------------------
void FluidFreeFormTransformation::
ParametricGradient(const PointSet &pos, const Vector3D<double> *in,
                   double *out, double t, double t0, double w) const
{
  // Compute gradient only if any level is active
  if (this->NumberOfDOFs() == 0) return;

  PointSet w_pos;
  w_pos.Reserve(pos.Size());

  int previous = -1;
  for (;;) {
    // Determine next active level
    int active = previous + 1;
    while (active < this->NumberOfLevels() && !this->LocalTransformationIsActive(active)) {
      ++active;
    }
    if (active == this->NumberOfLevels()) break;

    // Transform non-parametric gradient
    w_pos = pos;
    this->Transform(active, w_pos);

    // Add gradient w.r.t. DoFs of this active level
    const FreeFormTransformation *ffd = this->GetLocalTransformation(active);
    ffd->ParametricGradient(w_pos, in, out, t, t0, w);
    out += ffd->NumberOfDOFs();

    // Continue with next active level
    previous = active;
  }
}

// =============================================================================
// I/O
// =============================================================================

// -----------------------------------------------------------------------------
void FluidFreeFormTransformation::Print(Indent indent) const
{
  cout << indent << "Fluid FFD:" << endl;
  MultiLevelTransformation::Print(indent + 1);
  cout << (indent + 1) << "Affine post-transformation:" << endl;
  _AffineTransformation.Print(indent + 2);
}

// -----------------------------------------------------------------------------
bool FluidFreeFormTransformation::CanRead(TransformationType format) const
{
  switch (format) {
    case TRANSFORMATION_FLUID_v1:
    case TRANSFORMATION_FLUID_v2:
      return true;
    default:
      return false;
  }
}

// -----------------------------------------------------------------------------
Cifstream &FluidFreeFormTransformation::ReadDOFs(Cifstream &from, TransformationType format)
{
  // Read global and local transformations
  MultiLevelTransformation::ReadDOFs(from, format);

  // Read additional affine transformation
  return _AffineTransformation.Read(from);
}

// -----------------------------------------------------------------------------
Cofstream &FluidFreeFormTransformation::WriteDOFs(Cofstream &to) const
{
  // Write global and local transformations
  MultiLevelTransformation::WriteDOFs(to);

  // Write additional affine transformation
  return _AffineTransformation.Write(to);
}


} // namespace mirtk
