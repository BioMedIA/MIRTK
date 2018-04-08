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

#include "mirtk/LinearFreeFormTransformation3D.h"

#include "mirtk/Math.h"
#include "mirtk/Memory.h"
#include "mirtk/BSplineFreeFormTransformation3D.h"


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
LinearFreeFormTransformation3D::LinearFreeFormTransformation3D()
:
  FreeFormTransformation3D(_FFD)
{
}

// -----------------------------------------------------------------------------
LinearFreeFormTransformation3D
::LinearFreeFormTransformation3D(double x1, double y1, double z1,
                                 double x2, double y2, double z2,
                                 double dx, double dy, double dz,
                                 double *xaxis, double *yaxis, double *zaxis)
:
  FreeFormTransformation3D(_FFD)
{
  Initialize(DefaultAttributes(x1, y1, z1, x2, y2, z2, dx, dy, dz, xaxis, yaxis, zaxis));
}

// -----------------------------------------------------------------------------
LinearFreeFormTransformation3D
::LinearFreeFormTransformation3D(const ImageAttributes &attr,
                                 double dx, double dy, double dz)
:
  FreeFormTransformation3D(_FFD)
{
  Initialize(attr, dx, dy, dz);
}

// -----------------------------------------------------------------------------
LinearFreeFormTransformation3D
::LinearFreeFormTransformation3D(const BaseImage &target,
                                 double dx, double dy, double dz)
:
  FreeFormTransformation3D(_FFD)
{
  Initialize(target.Attributes(), dx, dy, dz);
}

// -----------------------------------------------------------------------------
LinearFreeFormTransformation3D
::LinearFreeFormTransformation3D(const GenericImage<double> &image, bool disp)
:
  FreeFormTransformation3D(_FFD)
{
  Initialize(image, disp);
}

// -----------------------------------------------------------------------------
LinearFreeFormTransformation3D
::LinearFreeFormTransformation3D(const LinearFreeFormTransformation3D &ffd)
:
  FreeFormTransformation3D(ffd, _FFD)
{
  if (_NumberOfDOFs > 0) InitializeInterpolator();
}

// -----------------------------------------------------------------------------
LinearFreeFormTransformation3D
::LinearFreeFormTransformation3D(const BSplineFreeFormTransformation3D &ffd)
:
  FreeFormTransformation3D(ffd, _FFD)
{
  if (_NumberOfDOFs > 0) InitializeInterpolator();
  // Convert B-spline coefficients to linear interpolation coefficients
  double dx, dy, dz;
  Vector *data = _CPImage.Data();
  for (int k = 0; k < _z; ++k)
  for (int j = 0; j < _y; ++j)
  for (int i = 0; i < _x; ++i, ++data) {
    dx = i, dy = j, dz = k;
    ffd.Evaluate(dx, dy, dz);
    data->_x = dx;
    data->_y = dy;
    data->_z = dz;
  }
}

// -----------------------------------------------------------------------------
LinearFreeFormTransformation3D::~LinearFreeFormTransformation3D()
{
}

// =============================================================================
// Approximation/Interpolation
// =============================================================================

// -----------------------------------------------------------------------------
void LinearFreeFormTransformation3D
::ApproximateDOFs(const double *wx, const double *wy, const double *wz, const double *,
                  const double *dx, const double *dy, const double *dz, int no)
{
  int    i, j, k, ci, cj, ck;
  double x, y, z, w[3], B[3][2];

  // Allocate memory
  Vector ***data = CAllocate<Vector>(_x, _y, _z);
  double ***norm = CAllocate<double>(_x, _y, _z);

  // Initial loop: Calculate change of control points
  for (int idx = 0; idx < no; ++idx) {
    x = wx[idx], y = wy[idx], z = wz[idx];
    this->WorldToLattice(x, y, z);

    i = ifloor(x);
    j = ifloor(y);
    k = ifloor(z);
    B[0][1] = x - i, B[0][0] = 1.0 - B[0][1];
    B[1][1] = y - j, B[1][0] = 1.0 - B[1][1];
    B[2][1] = z - k, B[2][0] = 1.0 - B[2][1];

    for (int c = 0; c <= 1; ++c) {
      ck = k + c;
      if (0 <= ck && ck < _z) {
        w[2] = B[2][c];
        for (int b = 0; b <= 1; ++b) {
          cj = j + b;
          if (0 <= cj && cj < _y) {
            w[1] = B[1][b] * w[2];
            for (int a = 0; a <= 1; ++a) {
              ci = i + a;
              if (0 <= ci && ci < _x) {
                w[0] = B[0][a] * w[1];
                data[ck][cj][ci]._x += w[0] * dx[idx];
                data[ck][cj][ci]._y += w[0] * dy[idx];
                data[ck][cj][ci]._z += w[0] * dz[idx];
                norm[ck][cj][ci]    += w[0];
              }
            }
          }
        }
      }
    }
  }

  // Final loop: Calculate new control points
  const Vector zero(.0);

  Vector *in  = data[0][0];
  double *div = norm[0][0];
  Vector *out = _CPImage.Data();

  for (int cp = 0; cp < NumberOfCPs(); ++cp, ++in, ++div, ++out) {
    (*out) = (((*div) > .0) ? ((*in) / (*div)) : zero);
  }

  // Deallocate memory
  Deallocate(data);
  Deallocate(norm);
}

// -----------------------------------------------------------------------------
void LinearFreeFormTransformation3D
::ApproximateDOFsGradient(const double *wx, const double *wy, const double *wz, const double *,
                          const double *dx, const double *dy, const double *dz,
                          int no, double *gradient, double weight) const
{
  int    i, j, k, ci, cj, ck;
  double x, y, z, w[3], B[3][2];

  // Allocate memory
  Vector ***data = CAllocate<Vector>(_x, _y, _z);

  // Initial loop: Calculate change of control points using NN extrapolation
  for (int idx = 0; idx < no; ++idx) {
    x = wx[idx], y = wy[idx], z = wz[idx];
    this->WorldToLattice(x, y, z);

    i = ifloor(x);
    j = ifloor(y);
    k = ifloor(z);
    B[0][1] = x - i, B[0][0] = 1.0 - B[0][1];
    B[1][1] = y - j, B[1][0] = 1.0 - B[1][1];
    B[2][1] = z - k, B[2][0] = 1.0 - B[2][1];

    for (int c = 0; c <= 1; ++c) {
      ck = k + c;
      if (0 <= ck && ck < _z) {
        w[2] = B[2][c];
        for (int b = 0; b <= 1; ++b) {
          cj = j + b;
          if (0 <= cj && cj < _y) {
            w[1] = B[1][b] * w[2];
            for (int a = 0; a <= 1; ++a) {
              ci = i + a;
              if (0 <= ci && ci < _x) {
                w[0] = B[0][a] * w[1];
                data[ck][cj][ci]._x += w[0] * dx[idx];
                data[ck][cj][ci]._y += w[0] * dy[idx];
                data[ck][cj][ci]._z += w[0] * dz[idx];
              }
            }
          }
        }
      }
    }
  }

  // Final loop
  int           xdof, ydof, zdof;
  const Vector *grad = data[0][0];

  for (int cp = 0; cp < NumberOfCPs(); ++cp, ++grad) {
    this->IndexToDOFs(cp, xdof, ydof, zdof);
    gradient[xdof] = weight * grad->_x;
    gradient[ydof] = weight * grad->_y;
    gradient[zdof] = weight * grad->_z;
  }

  // Deallocate memory
  Deallocate(data);
}

// -----------------------------------------------------------------------------
void LinearFreeFormTransformation3D
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
void LinearFreeFormTransformation3D
::BoundingBox(int cp, double &x1, double &y1, double &z1,
                      double &x2, double &y2, double &z2, double fraction) const
{
  int i, j, k;
  IndexToLattice(cp, i, j, k);

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

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
void LinearFreeFormTransformation3D
::EvaluateJacobian(Matrix &jac, double x, double y, double z) const
{
  double x1, y1, z1, x2, y2, z2;

  // Jacobian matrix is 3x3
  jac.Initialize(3, 3);

  // Compute derivative in x
  x1 = x - 0.5;
  y1 = y;
  z1 = z;
  x2 = x + 0.5;
  y2 = y;
  z2 = z;
  this->Evaluate(x1, y1, z1);
  this->Evaluate(x2, y2, z2);
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
  this->Evaluate(x1, y1, z1);
  this->Evaluate(x2, y2, z2);
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
  this->Evaluate(x1, y1, z1);
  this->Evaluate(x2, y2, z2);
  jac(2, 0) = x2 - x1;
  jac(2, 1) = y2 - y1;
  jac(2, 2) = z2 - z1;
}

// =============================================================================
// Properties
// =============================================================================

// -----------------------------------------------------------------------------
int LinearFreeFormTransformation3D::KernelSize() const
{
  return 2;
}

// -----------------------------------------------------------------------------
// classic pairwise bending energy
double LinearFreeFormTransformation3D::BendingEnergy(int i, int j, int k) const
{
  const Vector &disp    = _CPImage(i,  j, k);
  double        bending = .0;

  int I = i + 1;
  int J = j + 1;
  int K = k + 1;

  if (I >= _x) I = _x - 1;
  if (J >= _y) J = _y - 1;
  if (K >= _z) K = _z - 1;

  if (i != I) {
    Vector d = (_CPImage(I, j, k) - disp) / double(I - i) / _dx;
    bending += d._x * d._x + d._y * d._y + d._z * d._z;
  }
  if (j != J) {
    Vector d = (_CPImage(i, J, k) - disp) / double(J - j) / _dy;
    bending += d._x * d._x + d._y * d._y + d._z * d._z;
  }
  if (k != K) {
    Vector d = (_CPImage(i, j, K) - disp) / double(K - k) / _dy;
    bending += d._x * d._x + d._y * d._y + d._z * d._z;
  }

  return bending;
}

// -----------------------------------------------------------------------------
double LinearFreeFormTransformation3D::BendingEnergy(bool incl_passive, bool) const
{
  int    nactive = 0;
  double bending = .0;

  for (int k = 0; k < _z; ++k)
  for (int j = 0; j < _y; ++j)
  for (int i = 0; i < _x; ++i) {
    if (incl_passive || this->IsActive(i, j, k)) {
      bending += this->BendingEnergy(i, j, k);
      ++nactive;
    }
  }

  if (nactive > 0) bending /= nactive;
  return bending;
}

// -----------------------------------------------------------------------------
// derivative of classic pairwise bending energy
void LinearFreeFormTransformation3D
::BendingEnergyGradient(double *gradient, double weight, bool incl_passive, bool, bool) const
{
  int    i1, j1, k1, i2, j2, k2;
  Vector tmp;

  weight *= 2.0; // pre-multiply weight with constant factor

  int xdof, ydof, zdof;
  for (int k = 0; k < _z; ++k)
  for (int j = 0; j < _y; ++j)
  for (int i = 0; i < _x; ++i) {
    if (incl_passive || IsActive(i, j, k)) {
      i1 = i - 1;
      j1 = j - 1;
      k1 = k - 1;
      if (i1 < 0) i1 = 0;
      if (j1 < 0) j1 = 0;
      if (k1 < 0) k1 = 0;

      i2 = i + 1;
      j2 = j + 1;
      k2 = k + 1;
      if (i2 >= _x) i2 = _x - 1;
      if (j2 >= _y) j2 = _y - 1;
      if (k2 >= _z) k2 = _z - 1;

      const Vector  ijk = _CPImage(i,  j,  k) * 2.0;
      const Vector &dx1 = _CPImage(i1, j,  k);
      const Vector &dx2 = _CPImage(i2, j,  k);
      const Vector &dy1 = _CPImage(i,  j1, k);
      const Vector &dy2 = _CPImage(i,  j2, k);
      const Vector &dz1 = _CPImage(i,  j,  k1);
      const Vector &dz2 = _CPImage(i,  j,  k2);

      tmp = (ijk - dz1 - dz2) / _dz / _dz
          + (ijk - dy1 - dy2) / _dy / _dy
          + (ijk - dx1 - dx2) / _dx / _dx;

      IndexToDOFs(LatticeToIndex(i, j, k), xdof, ydof, zdof);
      gradient[xdof] += weight * tmp._x;
      gradient[ydof] += weight * tmp._y;
      gradient[zdof] += weight * tmp._z;
    }
  }
}

// =============================================================================
// I/O
// =============================================================================

// -----------------------------------------------------------------------------
void LinearFreeFormTransformation3D::Print(ostream &os, Indent indent) const
{
  os << indent << "3D Linear FFD:" << endl;
  FreeFormTransformation3D::Print(os, indent + 1);
}

// -----------------------------------------------------------------------------
bool LinearFreeFormTransformation3D::CanRead(TransformationType format) const
{
  switch (format) {
    case TRANSFORMATION_LINEAR_FFD_2D_v1:
    case TRANSFORMATION_LINEAR_FFD_3D_v1:
    case TRANSFORMATION_LINEAR_FFD_3D_v2:
    case TRANSFORMATION_LINEAR_FFD_3D_v3:
    case TRANSFORMATION_LINEAR_FFD_3D_v4:
      return true;
    default:
      return false;
  }
}

// =============================================================================
// Others
// =============================================================================

// -----------------------------------------------------------------------------
void LinearFreeFormTransformation3D::Compose(const Transformation *t1)
{
  // Copy transformation
  LinearFreeFormTransformation3D t2(*this);

  // Compose t2 o t1
  double x1, y1, z1, x2, y2, z2;
  for (int k = 0; k < _z; ++k) {
    for (int j = 0; j < _y; ++j) {
      for (int i = 0; i < _x; ++i) {
        x1 = i, y1 = j, z1 = k;
        this->LatticeToWorld(x1, y1, z1);
        x2 = x1, y2 = y1, z2 = z1;
        // Apply first transformation
        t1->Transform(x2, y2, z2);
        // Apply second transformation
        t2.Transform(x2, y2, z2);
        // Update displacement
        _CPImage(i, j, k) = Vector(x2 - x1, y2 - y1, z2 - z2);
      }
    }
  }
}


} // namespace mirtk
