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

#include "mirtk/LinearFreeFormTransformation4D.h"

#include "mirtk/Math.h"
#include "mirtk/Memory.h"
#include "mirtk/BSplineFreeFormTransformation4D.h"


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
LinearFreeFormTransformation4D::LinearFreeFormTransformation4D()
:
  FreeFormTransformation4D(_FFD)
{
}

// -----------------------------------------------------------------------------
LinearFreeFormTransformation4D
::LinearFreeFormTransformation4D(double x1, double y1, double z1, double t1,
                                 double x2, double y2, double z2, double t2,
                                 double dx, double dy, double dz, double dt,
                                 double *xaxis, double *yaxis, double *zaxis)
:
  FreeFormTransformation4D(_FFD)
{
  Initialize(DefaultAttributes(x1, y1, z1, t1,
                               x2, y2, z2, t2,
                               dx, dy, dz, dt,
                               xaxis, yaxis, zaxis));
}

// -----------------------------------------------------------------------------
LinearFreeFormTransformation4D
::LinearFreeFormTransformation4D(const ImageAttributes &attr,
                                 double dx, double dy, double dz, double dt)
:
  FreeFormTransformation4D(_FFD)
{
  Initialize(attr, dx, dy, dz, dt);
}

// -----------------------------------------------------------------------------
LinearFreeFormTransformation4D
::LinearFreeFormTransformation4D(const BaseImage &target,
                                 double dx, double dy, double dz, double dt)
:
  FreeFormTransformation4D(_FFD)
{
  Initialize(target.Attributes(), dx, dy, dz, dt);
}

// -----------------------------------------------------------------------------
LinearFreeFormTransformation4D
::LinearFreeFormTransformation4D(const LinearFreeFormTransformation4D &ffd)
:
  FreeFormTransformation4D(ffd, _FFD)
{
  if (_NumberOfDOFs > 0) InitializeInterpolator();
}

// -----------------------------------------------------------------------------
LinearFreeFormTransformation4D
::LinearFreeFormTransformation4D(const BSplineFreeFormTransformation4D &ffd)
:
  FreeFormTransformation4D(ffd, _FFD)
{
  if (_NumberOfDOFs > 0) InitializeInterpolator();
  // Convert B-spline coefficients to linear interpolation coefficients
  double x, y, z;
  Vector *data = _CPImage.Data();
  for (int l = 0; l < _t; ++l)
  for (int k = 0; k < _z; ++k)
  for (int j = 0; j < _y; ++j)
  for (int i = 0; i < _x; ++i, ++data) {
    x = i, y = j, z = k;
    ffd.Evaluate(x, y, z, l);
    data->_x = x;
    data->_y = y;
    data->_z = z;
  }
}

// -----------------------------------------------------------------------------
LinearFreeFormTransformation4D::~LinearFreeFormTransformation4D()
{
}

// =============================================================================
// Approximation/Interpolation
// =============================================================================

// -----------------------------------------------------------------------------
void LinearFreeFormTransformation4D
::ApproximateDOFs(const double *wx, const double *wy, const double *wz, const double *wt,
                  const double *dx, const double *dy, const double *dz, int no)
{
  int    i, j, k, l, ci, cj, ck, cl;
  double x, y, z, t, w[4], B[4][2];

  // Allocate memory
  Vector ****data = CAllocate<Vector>(_x, _y, _z, _t);
  double ****norm = CAllocate<double>(_x, _y, _z, _t);

  // Initial loop: Calculate change of control points using NN extrapolation
  for (int idx = 0; idx < no; ++idx) {
    x = wx[idx], y = wy[idx], z = wz[idx];
    this->WorldToLattice(x, y, z);
    t = this->TimeToLattice(wt[idx]);

    i = ifloor(x);
    j = ifloor(y);
    k = ifloor(z);
    l = ifloor(t);
    B[0][1] = x - i, B[0][0] = 1.0 - B[0][1];
    B[1][1] = y - j, B[1][0] = 1.0 - B[1][1];
    B[2][1] = z - k, B[2][0] = 1.0 - B[2][1];
    B[3][1] = t - l, B[3][0] = 1.0 - B[3][1];

    for (int d = 0; d <= 1; ++d) {
      cl = l + d;
      if (0 <= cl && cl < _t) {
        w[3] = B[3][d];
        for (int c = 0; c <= 1; ++c) {
          ck = k + c;
          if (0 <= ck && ck < _z) {
            w[2] = B[2][c] * w[3];
            for (int b = 0; b <= 1; ++b) {
              cj = j + b;
              if (0 <= cj && cj < _y) {
                w[1] = B[1][b] * w[2];
                for (int a = 0; a <= 1; ++a) {
                  ci = i + a;
                  if (0 <= ci && ci < _x) {
                    w[0] = B[0][a] * w[1];
                    data[cl][ck][cj][ci]._x += w[0] * dx[idx];
                    data[cl][ck][cj][ci]._y += w[0] * dy[idx];
                    data[cl][ck][cj][ci]._z += w[0] * dz[idx];
                    norm[cl][ck][cj][ci]    += w[0];
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  // Final loop: Calculate new control points
  const Vector zero(.0);

  Vector *in  = data[0][0][0];
  double *div = norm[0][0][0];
  Vector *out = _CPImage.Data();

  for (int cp = 0; cp < NumberOfCPs(); ++cp, ++in, ++div, ++out) {
    (*out) = ((*div) ? (*in) / (*div) : zero);
  }

  // Deallocate memory
  Deallocate(data);
  Deallocate(norm);
}

// -----------------------------------------------------------------------------
void LinearFreeFormTransformation4D
::ApproximateDOFsGradient(const double *wx, const double *wy, const double *wz, const double *wt,
                          const double *dx, const double *dy, const double *dz,
                          int no, double *gradient, double weight) const
{
  int    i, j, k, l, ci, cj, ck, cl;
  double x, y, z, t, w[4], B[4][2];

  // Allocate memory
  Vector ****data = CAllocate<Vector>(_x, _y, _z, _t);

  // Initial loop: Calculate change of control points using NN extrapolation
  for (int idx = 0; idx < no; ++idx) {
    x = wx[idx], y = wy[idx], z = wz[idx];
    this->WorldToLattice(x, y, z);
    t = this->TimeToLattice(wt[idx]);

    i = ifloor(x);
    j = ifloor(y);
    k = ifloor(z);
    l = ifloor(t);
    B[0][1] = x - i, B[0][0] = 1.0 - B[0][1];
    B[1][1] = y - j, B[1][0] = 1.0 - B[1][1];
    B[2][1] = z - k, B[2][0] = 1.0 - B[2][1];
    B[3][1] = t - l, B[3][0] = 1.0 - B[3][1];

    for (int d = 0; d <= 1; ++d) {
      cl = l + d;
      if (0 <= cl && cl < _t) {
        w[3] = B[3][d];
        for (int c = 0; c <= 1; ++c) {
          ck = k + c;
          if (0 <= ck && ck < _z) {
            w[2] = B[2][c] * w[3];
            for (int b = 0; b <= 1; ++b) {
              cj = j + b;
              if (0 <= cj && cj < _y) {
                w[1] = B[1][b] * w[2];
                for (int a = 0; a <= 1; ++a) {
                  ci = i + a;
                  if (0 <= ci && ci < _x) {
                    w[0] = B[0][a] * w[1];
                    data[cl][ck][cj][ci]._x += w[0] * dx[idx];
                    data[cl][ck][cj][ci]._y += w[0] * dy[idx];
                    data[cl][ck][cj][ci]._z += w[0] * dz[idx];
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  // Final loop
  int           xdof, ydof, zdof;
  const Vector *grad = data[0][0][0];

  for (int cp = 0; cp < NumberOfCPs(); ++cp, ++grad) {
    this->IndexToDOFs(cp, xdof, ydof, zdof);
    gradient[xdof] += weight * grad->_x;
    gradient[ydof] += weight * grad->_y;
    gradient[zdof] += weight * grad->_z;
  }

  // Deallocate memory
  Deallocate(data);
}

// -----------------------------------------------------------------------------
void LinearFreeFormTransformation4D
::Interpolate(const double *dx, const double *dy, const double *dz)
{
  Vector *param = _CPImage.Data();
  for (int cp = 0; cp < NumberOfCPs(); ++cp, ++param) {
    param->_x = dx[cp];
    param->_y = dy[cp];
    param->_z = dz[cp];
  }
}

// =============================================================================
// Lattice
// =============================================================================

// -----------------------------------------------------------------------------
void LinearFreeFormTransformation4D
::BoundingBox(int cp, double &t1, double &t2, double fraction) const
{
  int n = _x * _y * _z * _t;
  if (cp >= n) {
    cp -= n;
    if (cp >= n) cp -= n;
  }
  n = cp / (_x * _y * _z);
  t1 = this->LatticeToTime(n - fraction);
  t2 = this->LatticeToTime(n + fraction);
  if (t1 > t2) swap(t1, t2);
}

// -----------------------------------------------------------------------------
void LinearFreeFormTransformation4D
::BoundingBox(int cp, double &x1, double &y1, double &z1,
                      double &x2, double &y2, double &z2, double fraction) const
{
  int i, j, k, l;
  IndexToLattice(cp, i, j, k, l);

  x1 = i - fraction;
  y1 = j - fraction;
  z1 = k - fraction;
  x2 = i + fraction;
  y2 = j + fraction;
  z2 = k + fraction;

  this->LatticeToWorld(x1, y1, z1);
  this->LatticeToWorld(x2, y2, z2);

  if (x1 > x2) swap(x1, x2);
  if (y1 > y2) swap(y1, y2);
  if (z1 > z2) swap(z1, z2);
}

// -----------------------------------------------------------------------------
void LinearFreeFormTransformation4D
::BoundingBox(int cp, double &x1, double &y1, double &z1, double &t1,
                      double &x2, double &y2, double &z2, double &t2, double fraction) const
{
  BoundingBox(cp, x1, y1, z1, x2, y2, z2, fraction);
  BoundingBox(cp, t1, t2, fraction);
  int i, j, k, l;
  IndexToLattice(cp, i, j, k, l);

  x1 = i - fraction;
  y1 = j - fraction;
  z1 = k - fraction;
  t1 = l - fraction;
  x2 = i + fraction;
  y2 = j + fraction;
  z2 = k + fraction;
  t2 = l + fraction;

  this->LatticeToWorld(x1, y1, z1);
  this->LatticeToWorld(x2, y2, z2);
  t1 = this->LatticeToTime(t1);
  t2 = this->LatticeToTime(t2);

  if (x1 > x2) swap(x1, x2);
  if (y1 > y2) swap(y1, y2);
  if (z1 > z2) swap(z1, z2);
  if (t1 > t2) swap(t1, t2);
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
void LinearFreeFormTransformation4D
::EvaluateJacobian(Matrix &jac, double x, double y, double z, double t) const
{
  double x1, y1, z1, x2, y2, z2;

  // Jacobian matrix is 3 x 3
  jac.Initialize(3, 3);

  // Compute derivative in x
  x1 = x - 0.5;
  y1 = y;
  z1 = z;
  x2 = x + 0.5;
  y2 = y;
  z2 = z;
  this->Evaluate(x1, y1, z1, t);
  this->Evaluate(x2, y2, z2, t);
  jac(0, 0) = x2 - x1;
  jac(0, 1) = y2 - y1;
  jac(0, 2) = z2 - z1;

  // Compute derivative in y
  x1 = x;
  y1 = y - 0.5;
  z1 = z;
  x2 = x;
  y2 = y + 0.5;
  z2 = z;
  this->Evaluate(x1, y1, z1, t);
  this->Evaluate(x2, y2, z2, t);
  jac(1, 0) = x2 - x1;
  jac(1, 1) = y2 - y1;
  jac(1, 2) = z2 - z1;

  // Compute derivative in z
  x1 = x;
  y1 = y;
  z1 = z - 0.5;
  x2 = x;
  y2 = y;
  z2 = z + 0.5;
  this->Evaluate(x1, y1, z1, t);
  this->Evaluate(x2, y2, z2, t);
  jac(2, 0) = x2 - x1;
  jac(2, 1) = y2 - y1;
  jac(2, 2) = z2 - z1;
}

// =============================================================================
// Properties
// =============================================================================

// -----------------------------------------------------------------------------
int LinearFreeFormTransformation4D::KernelSize() const
{
  return 2;
}

// -----------------------------------------------------------------------------
// classic pairwise bending energy
double LinearFreeFormTransformation4D::BendingEnergy(int i, int j, int k, int l) const
{
  const Vector &disp    = _CPImage(i,  j, k, l);
  double        bending = .0;

  int I = i + 1;
  int J = j + 1;
  int K = k + 1;
  int L = l + 1;

  if (I >= _x) I = _x - 1;
  if (J >= _y) J = _y - 1;
  if (K >= _z) K = _z - 1;
  if (L >= _t) L = _t - 1;

  if (i != I) {
    Vector d = (_CPImage(I, j, k, l) - disp) / double(I - i) / _dx;
    bending += d._x * d._x + d._y * d._y + d._z * d._z;
  }
  if (j != J) {
    Vector d = (_CPImage(i, J, k, l) - disp) / double(J - j) / _dy;
    bending += d._x * d._x + d._y * d._y + d._z * d._z;
  }
  if (k != K) {
    Vector d = (_CPImage(i, j, K, l) - disp) / double(K - k) / _dy;
    bending += d._x * d._x + d._y * d._y + d._z * d._z;
  }
  if (l != L) {
    Vector d = (_CPImage(i, j, k, L) - disp) / double(L - l) / _dt;
    bending += d._x * d._x + d._y * d._y + d._z * d._z;
  }

  return bending;
}

// -----------------------------------------------------------------------------
double LinearFreeFormTransformation4D::BendingEnergy(bool incl_passive, bool) const
{
  int    nactive = 0;
  double bending = .0;

  for (int l = 0; l < _t; ++l)
  for (int k = 0; k < _z; ++k)
  for (int j = 0; j < _y; ++j)
  for (int i = 0; i < _x; ++i) {
    if (incl_passive || this->IsActive(i, j, k, l)) {
      bending += this->BendingEnergy(i, j, k, l);
      ++nactive;
    }
  }

  if (nactive > 0) bending /= nactive;
  return bending;
}

// -----------------------------------------------------------------------------
// derivative of classic pairwise bending energy
void LinearFreeFormTransformation4D
::BendingEnergyGradient(double *gradient, double weight, bool incl_passive, bool, bool) const
{
  int    i1, j1, k1, l1, i2, j2, k2, l2;
  Vector tmp;

  weight *= 2.0; // pre-multiply weight by constant factor

  int xdof, ydof, zdof;
  for (int l = 0; l < _t; ++l)
  for (int k = 0; k < _z; ++k)
  for (int j = 0; j < _y; ++j)
  for (int i = 0; i < _x; ++i) {
    if (incl_passive || IsActive(i, j, k, l)) {
      i1 = i - 1;
      j1 = j - 1;
      k1 = k - 1;
      l1 = l - 1;
      if (i1 < 0) i1 = 0;
      if (j1 < 0) j1 = 0;
      if (k1 < 0) k1 = 0;
      if (l1 < 0) l1 = 0;

      i2 = i + 1;
      j2 = j + 1;
      k2 = k + 1;
      l2 = l + 1;
      if (i2 >= _x) i2 = _x - 1;
      if (j2 >= _y) j2 = _y - 1;
      if (k2 >= _z) k2 = _z - 1;
      if (l2 >= _t) l2 = _t - 1;

      const Vector ijkl = _CPImage(i,  j,  k,  l) * 2.0;
      const Vector &dx1 = _CPImage(i1, j,  k,  l);
      const Vector &dx2 = _CPImage(i2, j,  k,  l);
      const Vector &dy1 = _CPImage(i,  j1, k,  l);
      const Vector &dy2 = _CPImage(i,  j2, k,  l);
      const Vector &dz1 = _CPImage(i,  j,  k1, l);
      const Vector &dz2 = _CPImage(i,  j,  k2, l);
      const Vector &dt1 = _CPImage(i,  j,  k,  l1);
      const Vector &dt2 = _CPImage(i,  j,  k,  l2);

      tmp = (ijkl - dt1 - dt2) / _dt / _dt
          + (ijkl - dz1 - dz2) / _dz / _dz
          + (ijkl - dy1 - dy2) / _dy / _dy
          + (ijkl - dx1 - dx2) / _dx / _dx;

      IndexToDOFs(LatticeToIndex(i, j, k, l), xdof, ydof, zdof);
      gradient[xdof] = weight * tmp._x;
      gradient[ydof] = weight * tmp._y;
      gradient[zdof] = weight * tmp._z;
    }
  }
}

// =============================================================================
// I/O
// =============================================================================

// -----------------------------------------------------------------------------
void LinearFreeFormTransformation4D::Print(ostream &os, Indent indent) const
{
  os << indent << "4D Linear FFD:" << endl;
  FreeFormTransformation4D::Print(os, indent + 1);
}

// -----------------------------------------------------------------------------
bool LinearFreeFormTransformation4D::CanRead(TransformationType format) const
{
  switch (format) {
    case TRANSFORMATION_LINEAR_FFD_4D_v1:
    case TRANSFORMATION_LINEAR_FFD_4D_v2:
    case TRANSFORMATION_LINEAR_FFD_4D_v3:
      return true;
    default:
      return false;
  }
}


} // namespace mirtk
