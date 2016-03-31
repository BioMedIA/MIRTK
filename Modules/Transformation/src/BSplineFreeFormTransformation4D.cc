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

#include <mirtkBSplineFreeFormTransformation4D.h>

#include <mirtkMath.h>
#include <mirtkMemory.h>
#include <mirtkParallel.h>
#include <mirtkProfiling.h>

#include <mirtkImageToInterpolationCoefficients.h>


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
BSplineFreeFormTransformation4D::BSplineFreeFormTransformation4D()
:
  FreeFormTransformation4D(_FFD)
{
}

// -----------------------------------------------------------------------------
BSplineFreeFormTransformation4D
::BSplineFreeFormTransformation4D(double x1, double y1, double z1, double t1,
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
BSplineFreeFormTransformation4D
::BSplineFreeFormTransformation4D(const ImageAttributes &attr,
                                  double dx, double dy, double dz, double dt)
:
  FreeFormTransformation4D(_FFD)
{
  Initialize(attr, dx, dy, dz, dt);
}

// -----------------------------------------------------------------------------
BSplineFreeFormTransformation4D
::BSplineFreeFormTransformation4D(const BaseImage &target,
                                  double dx, double dy, double dz, double dt)
:
  FreeFormTransformation4D(_FFD)
{
  Initialize(target.Attributes(), dx, dy, dz, dt);
}

// -----------------------------------------------------------------------------
BSplineFreeFormTransformation4D
::BSplineFreeFormTransformation4D(const BSplineFreeFormTransformation4D &ffd)
:
  FreeFormTransformation4D(ffd, _FFD)
{
  if (_NumberOfDOFs > 0) InitializeInterpolator();
}

// -----------------------------------------------------------------------------
BSplineFreeFormTransformation4D::~BSplineFreeFormTransformation4D()
{
}

// =============================================================================
// Approximation/Interpolation
// =============================================================================

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformation4D
::ApproximateDOFs(const double *wx, const double *wy, const double *wz, const double *wt,
                  const double *dx, const double *dy, const double *dz, int no)
{
  int    i, j, k, l, ci, cj, ck, cl, A, B, C, D;
  double x, y, z, t, w[4], basis, sum;

  // Allocate memory
  Vector ****data = CAllocate<Vector>(_x, _y, _z, _t);
  double ****norm = CAllocate<double>(_x, _y, _z, _t);

  // Initial loop: Calculate change of control points
  for (int idx = 0; idx < no; idx++) {
    x = wx[idx], y = wy[idx], z = wz[idx];
    this->WorldToLattice(x, y, z);
    t = this->TimeToLattice(wt[idx]);

    i = ifloor(x);
    j = ifloor(y);
    k = ifloor(z);
    l = ifloor(t);
    A = Kernel::VariableToIndex(x - i);
    B = Kernel::VariableToIndex(y - j);
    C = Kernel::VariableToIndex(z - k);
    D = Kernel::VariableToIndex(t - l);
    --i, --j, --k, --l;

    sum = .0;
    for (int d = 0; d <= 3; ++d) {
      w[3] = Kernel::LookupTable[D][d];
      for (int c = 0; c <= 3; ++c) {
        w[2] = Kernel::LookupTable[C][c] * w[3];
        for (int b = 0; b <= 3; ++b) {
          w[1] = Kernel::LookupTable[B][b] * w[2];
          for (int a = 0; a <= 3; ++a) {
            w[0] = Kernel::LookupTable[A][a] * w[1];
            sum += w[0] * w[0];
          }
        }
      }
    }

    for (int d = 0; d <= 3; ++d) {
      cl = l + d;
      if (0 <= cl && cl < _t) {
        w[3] = Kernel::LookupTable[D][d];
        for (int c = 0; c <= 3; ++c) {
          ck = k + c;
          if (0 <= ck && ck < _z) {
            w[2] = Kernel::LookupTable[C][c] * w[3];
            for (int b = 0; b <= 3; ++b) {
              cj = j + b;
              if (0 <= cj && cj < _y) {
                w[1] = Kernel::LookupTable[B][b] * w[2];
                for (int a = 0; a <= 3; ++a) {
                  ci = i + a;
                  if (0 <= ci && ci < _x) {
                    w[0]  = Kernel::LookupTable[A][a] * w[1];
                    basis = w[0] * w[0];
                    norm[cl][ck][cj][ci] += basis;
                    basis *= w[0] / sum;
                    data[cl][ck][cj][ci]._x += basis * dx[idx];
                    data[cl][ck][cj][ci]._y += basis * dy[idx];
                    data[cl][ck][cj][ci]._z += basis * dz[idx];
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

  for (int cp = 0; cp < NumberOfCPs(); ++cp, ++out, ++in, ++div) {
    (*out) = ((*div) ? ((*in) / (*div)) : zero);
  }

  // Deallocate memory
  Deallocate(data);
  Deallocate(norm);
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformation4D
::ApproximateDOFsGradient(const double *wx, const double *wy, const double *wz, const double *wt,
                          const double *dx, const double *dy, const double *dz,
                          int no, double *gradient, double weight) const
{
  int    i, j, k, l, ci, cj, ck, cl, A, B, C, D;
  double x, y, z, t, w[4];

  // Allocate memory
  Vector ****data = CAllocate<Vector>(_x, _y, _z, _t);

  // Initial loop: Calculate change of control points
  for (int idx = 0; idx < no; ++idx) {
    x = wx[idx], y = wy[idx], z = wz[idx];
    this->WorldToLattice(x, y, z);
    t = this->TimeToLattice(wt[idx]);

    i = ifloor(x);
    j = ifloor(y);
    k = ifloor(z);
    l = ifloor(t);
    A = Kernel::VariableToIndex(x - i);
    B = Kernel::VariableToIndex(y - j);
    C = Kernel::VariableToIndex(z - k);
    D = Kernel::VariableToIndex(t - l);
    --i, --j, --k, --l;

    for (int d = 0; d <= 3; ++d) {
      cl = l + d;
      if (0 <= cl && cl < _t) {
        w[3] = Kernel::LookupTable[D][d];
        for (int c = 0; c <= 3; ++c) {
          ck = k + c;
          if (0 <= ck && ck < _z) {
            w[2] = Kernel::LookupTable[C][c] * w[3];
            for (int b = 0; b <= 3; ++b) {
              cj = j + b;
              if (0 <= cj && cj < _y) {
                w[1] = Kernel::LookupTable[B][b] * w[2];
                for (int a = 0; a <= 3; ++a) {
                  ci = i + a;
                  if (0 <= ci && ci < _x) {
                    w[0] = Kernel::LookupTable[A][a] * w[1];
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
void BSplineFreeFormTransformation4D::Interpolate(const double *dx, const double *dy, const double *dz)
{
  for (int idx = 0; idx < NumberOfCPs(); ++idx) {
    _CPImage(idx) = Vector(dx[idx], dy[idx], dz[idx]);
  }
  ConvertToSplineCoefficients(3, _CPImage);
}

// =============================================================================
// Lattice
// =============================================================================

// -----------------------------------------------------------------------------
int BSplineFreeFormTransformation4D::GetXAfterSubdivision() const
{
  return (_x == 1 ? _x : 2 * _x - 1);
}

// -----------------------------------------------------------------------------
int BSplineFreeFormTransformation4D::GetYAfterSubdivision() const
{
  return (_y == 1 ? _y : 2 * _y - 1);
}

// -----------------------------------------------------------------------------
int BSplineFreeFormTransformation4D::GetZAfterSubdivision() const
{
  return (_z == 1 ? _z : 2 * _z - 1);
}

// -----------------------------------------------------------------------------
int BSplineFreeFormTransformation4D::GetTAfterSubdivision() const
{
  return (_t == 1 ? _t : 2 * _t - 1);
}

// -----------------------------------------------------------------------------
double BSplineFreeFormTransformation4D::GetXSpacingAfterSubdivision() const
{
  return (_x == 1 ? _dx : 0.5 * _dx);
}

// -----------------------------------------------------------------------------
double BSplineFreeFormTransformation4D::GetYSpacingAfterSubdivision() const
{
  return (_y == 1 ? _dy : 0.5 * _dy);
}

// -----------------------------------------------------------------------------
double BSplineFreeFormTransformation4D::GetZSpacingAfterSubdivision() const
{
  return (_z == 1 ? _dz : 0.5 * _dz);
}

// -----------------------------------------------------------------------------
double BSplineFreeFormTransformation4D::GetTSpacingAfterSubdivision() const
{
  return (_t == 1 ? _dt : 0.5 * _dt);
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformation4D::Subdivide(bool subdivide_x, bool subdivide_y, bool subdivide_z, bool subdivide_t)
{
  if (!subdivide_x && !subdivide_y && !subdivide_z && !subdivide_t) return;

  if (_x == 1) subdivide_x = false;
  if (_y == 1) subdivide_y = false;
  if (_z == 1) subdivide_z = false;
  if (_t == 1) subdivide_t = false;

  // Weights for subdivision
  const double w1[2][3] = {{0.125, 0.75, 0.125}, {0.0, 0.5, 0.5}};
  const double w2[2][3] = {{0.0,   1.0,  0.0},   {0.0, 0.0, 0.0}};

  const double (*wx)[3] = (subdivide_x ? w1 : w2);
  const double (*wy)[3] = (subdivide_y ? w1 : w2);
  const double (*wz)[3] = (subdivide_z ? w1 : w2);
  const double (*wt)[3] = (subdivide_t ? w1 : w2);

  // Size of new control point grid
  const int X = (subdivide_x ? (2 * _x - 1) : _x);
  const int Y = (subdivide_y ? (2 * _y - 1) : _y);
  const int Z = (subdivide_z ? (2 * _z - 1) : _z);
  const int T = (subdivide_t ? (2 * _t - 1) : _t);

  // Allocate memory for new control points
  Vector ****data = Allocate<Vector>(X, Y, Z, T);

  // Limits for inner loops
  const int I2 = (subdivide_x ? 2 : 1); // s.t. i+i2-1 == i if no subdivision
  const int J2 = (subdivide_y ? 2 : 1);
  const int K2 = (subdivide_z ? 2 : 1);
  const int L2 = (subdivide_t ? 2 : 1);

  // Compute control point values on subdivided grid
  int si, sj, sk, sl, I1, J1, K1, L1;
  for (int l = 0; l < _t; ++l) {
    if (subdivide_t) sl = 2 * l, L1 = int(l < _t-1);
    else             sl =     l, L1 = 0;
    for (int k = 0; k < _z; ++k) {
      if (subdivide_z) sk = 2 * k, K1 = int(k < _z-1);
      else             sk =     k, K1 = 0;
      for (int j = 0; j < _y; ++j) {
        if (subdivide_y) sj = 2 * j, J1 = int(j < _y-1);
        else             sj =     j, J1 = 0;
        for (int i = 0; i < _x; ++i) {
          if (subdivide_x) si = 2 * i, I1 = int(i < _x-1);
          else             si =     i, I1 = 0;
          for (int l1 = 0; l1 <= L1; ++l1)
          for (int k1 = 0; k1 <= K1; ++k1)
          for (int j1 = 0; j1 <= J1; ++j1)
          for (int i1 = 0; i1 <= I1; ++i1) {
            for (int l2 = (subdivide_t ? 0 : L2); l2 <= L2; ++l2)
            for (int k2 = (subdivide_z ? 0 : K2); k2 <= K2; ++k2)
            for (int j2 = (subdivide_y ? 0 : J2); j2 <= J2; ++j2)
            for (int i2 = (subdivide_x ? 0 : I2); i2 <= I2; ++i2) {
              // If a dimension is not subdivided,
              // 1) si + i1    == i
              // 2) i + i2 - 1 == i
              // 3) w[i1][i2]  == 1.0
              data[sl+l1][sk+k1][sj+j1][si+i1] += wx[i1][i2] * wy[j1][j2] * wz[k1][k2] * wt[l1][l2]
                                                  * _FFD.Extrapolator()->Get(i+i2-1, j+j2-1, k+k2-1, l+l2-1);
            }
          }
        }
      }
    }
  }

  // Initialize subdivided control points
  ImageAttributes attr = this->Attributes();
  attr._x = X;
  attr._y = Y;
  attr._z = Z;
  attr._t = T;
  if (subdivide_x) attr._dx *= 0.5;
  if (subdivide_y) attr._dy *= 0.5;
  if (subdivide_z) attr._dz *= 0.5;
  if (subdivide_t) attr._dt *= 0.5;
  this->Initialize(attr);

  // Copy subdivided control point data
  memcpy(_CPImage.Data(), data[0][0][0], X * Y * Z * T * sizeof(Vector));

  // Deallocate temporary memory
  Deallocate(data);
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformation4D
::BoundingBox(int cp, double &t1, double &t2, double fraction) const
{
  int n = _x * _y * _z * _t;
  if (cp >= n) {
    cp -= n;
    if (cp >= n) cp -= n;
  }
  n = cp / (_x * _y * _z);
  fraction *= 2;
  t1 = this->LatticeToTime(n - fraction);
  t2 = this->LatticeToTime(n + fraction);

  if (t1 > t2) swap(t1, t2);
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformation4D
::BoundingBox(int cp, double &x1, double &y1, double &z1,
                      double &x2, double &y2, double &z2, double fraction) const
{
  int i, j, k, l;
  IndexToLattice(cp, i, j, k, l);

  fraction *= 2;
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
void BSplineFreeFormTransformation4D
::BoundingBox(int cp, double &x1, double &y1, double &z1, double &t1,
                      double &x2, double &y2, double &z2, double &t2, double fraction) const
{
  int i, j, k, l;
  IndexToLattice(cp, i, j, k, l);

  fraction *= 2;
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
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
// Note: We are only returning the first three columns of the Jacobian
//       (the full Jacobian is a 3x4 matrix and we return a 3x3 one)
template <class CPImage>
void EvaluateJacobian(const CPImage *coeff, Matrix &jac, int i, int j, int k, int l)
{
  typedef BSplineFreeFormTransformation4D::Kernel Kernel;

  const double *w[2] = {
    Kernel::LatticeWeights,
    Kernel::LatticeWeights_I
  };

  typename CPImage::VoxelType dx, dy, dz;
  int                         ia, jb, kc, ld;

  --i, --j, --k, --l;
  for (int d = 0; d < 3; ++d) {
    ld = l + d;
    for (int c = 0; c < 3; ++c) {
      kc = k + c;
      for (int b = 0; b < 3; ++b) {
        jb = j + b;
        for (int a = 0; a < 3; ++a) {
          ia = i + a;
          dx += (w[1][a] * w[0][b] * w[0][c] * w[0][d]) * coeff->Get(ia, jb, kc, ld);
          dy += (w[0][a] * w[1][b] * w[0][c] * w[0][d]) * coeff->Get(ia, jb, kc, ld);
          dz += (w[0][a] * w[0][b] * w[1][c] * w[0][d]) * coeff->Get(ia, jb, kc, ld);
        }
      }
    }
  }

  jac.Initialize(3, 3);
  jac(0, 0) = dx._x; jac(0, 1) = dy._x; jac(0, 2) = dz._x;
  jac(1, 0) = dx._y; jac(1, 1) = dy._y; jac(1, 2) = dz._y;
  jac(2, 0) = dx._z; jac(2, 1) = dy._z; jac(2, 2) = dz._z;
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformation4D
::EvaluateJacobian(Matrix &jac, int i, int j, int k, int l) const
{
  if (_FFD.IsInside(i, j, k, l)) mirtk::EvaluateJacobian(&_CPImage, jac, i, j, k, l);
  else                           mirtk::EvaluateJacobian( _CPValue, jac, i, j, k, l);
}

// -----------------------------------------------------------------------------
// Note: We are only returning the first three columns of the Jacobian
//       (the full Jacobian is a 3x4 matrix and we return a 3x3 one)
template <class CPImage>
void EvaluateJacobian(const CPImage *coeff, Matrix &jac, double x, double y, double z, double t)
{
  typedef BSplineFreeFormTransformation4D::Kernel Kernel;

  int i = ifloor(x);
  int j = ifloor(y);
  int k = ifloor(z);
  int l = ifloor(t);

  const int A = Kernel::VariableToIndex(x - i);
  const int B = Kernel::VariableToIndex(y - j);
  const int C = Kernel::VariableToIndex(z - k);
  const int D = Kernel::VariableToIndex(t - l);

  typename CPImage::VoxelType dx, dy, dz;
  double                      wx[2], wy[2], wz[2], wt;
  int                         ia, jb, kc, ld;

  --i, --j, --k, --l;
  for (int d = 0; d < 4; ++d) {
    ld = l + d;
    wt = Kernel::LookupTable[D][d];
    for (int c = 0; c < 4; ++c) {
      kc = k + c;
      wz[0] = Kernel::LookupTable  [C][c];
      wz[1] = Kernel::LookupTable_I[C][c];
      for (int b = 0; b < 4; ++b) {
        jb = j + b;
        wy[0] = Kernel::LookupTable  [B][b];
        wy[1] = Kernel::LookupTable_I[B][b];
        for (int a = 0; a < 4; ++a) {
          ia = i + a;
          wx[0] = Kernel::LookupTable  [A][a];
          wx[1] = Kernel::LookupTable_I[A][a];
          dx += (wx[1] * wy[0] * wz[0] * wt) * coeff->Get(ia, jb, kc, ld);
          dy += (wx[0] * wy[1] * wz[0] * wt) * coeff->Get(ia, jb, kc, ld);
          dz += (wx[0] * wy[0] * wz[1] * wt) * coeff->Get(ia, jb, kc, ld);
        }
      }
    }
  }

  jac.Initialize(3, 3);
  jac(0, 0) = dx._x; jac(0, 1) = dy._x; jac(0, 2) = dz._x;
  jac(1, 0) = dx._y; jac(1, 1) = dy._y; jac(1, 2) = dz._y;
  jac(2, 0) = dx._z; jac(2, 1) = dy._z; jac(2, 2) = dz._z;
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformation4D
::EvaluateJacobian(Matrix &jac, double x, double y, double z, double t) const
{
  if (_FFD.IsInside(x, y, z, t)) mirtk::EvaluateJacobian(&_CPImage, jac, x, y, z, t);
  else                           mirtk::EvaluateJacobian( _CPValue, jac, x, y, z, t);
}

// -----------------------------------------------------------------------------
template <class CPImage>
void EvaluateHessian(const CPImage *coeff, Matrix hessian[3], int i, int j, int k, int l)
{
  typedef BSplineFreeFormTransformation4D::Kernel Kernel;

  const double *w[3] = {
    Kernel::LatticeWeights,
    Kernel::LatticeWeights_I,
    Kernel::LatticeWeights_II
  };

  typename CPImage::VoxelType dxx, dxy, dxz, dyy, dyz, dzz;
  int                         ia, jb, kc, ld;

  --i, --j, --k, --l;
  for (int d = 0; d < 3; ++d) {
    ld = l + d;
    for (int c = 0; c < 3; ++c) {
      kc = k + c;
      for (int b = 0; b < 3; ++b) {
        jb = j + b;
        for (int a = 0; a < 3; ++a) {
          ia = i + a;
          dxx += (w[2][a] * w[0][b] * w[0][c] * w[0][d]) * coeff->Get(ia, jb, kc, ld);
          dxy += (w[1][a] * w[1][b] * w[0][c] * w[0][d]) * coeff->Get(ia, jb, kc, ld);
          dxz += (w[1][a] * w[0][b] * w[1][c] * w[0][d]) * coeff->Get(ia, jb, kc, ld);
          dyy += (w[0][a] * w[2][b] * w[0][c] * w[0][d]) * coeff->Get(ia, jb, kc, ld);
          dyz += (w[0][a] * w[1][b] * w[1][c] * w[0][d]) * coeff->Get(ia, jb, kc, ld);
          dzz += (w[0][a] * w[0][b] * w[2][c] * w[0][d]) * coeff->Get(ia, jb, kc, ld);
        }
      }
    }
  }

  Matrix &hx = hessian[0];
  hx.Initialize(3, 3);
  hx(0, 0) = dxx._x; hx(0, 1) = dxy._x; hx(0, 2) = dxz._x;
  hx(1, 0) = dxy._x; hx(1, 1) = dyy._x; hx(1, 2) = dyz._x;
  hx(2, 0) = dxz._x; hx(2, 1) = dyz._x; hx(2, 2) = dzz._x;

  Matrix &hy = hessian[1];
  hy.Initialize(3, 3);
  hy(0, 0) = dxx._y; hy(0, 1) = dxy._y; hy(0, 2) = dxz._y;
  hy(1, 0) = dxy._y; hy(1, 1) = dyy._y; hy(1, 2) = dyz._y;
  hy(2, 0) = dxz._y; hy(2, 1) = dyz._y; hy(2, 2) = dzz._y;

  Matrix &hz = hessian[2];
  hz.Initialize(3, 3);
  hz(0, 0) = dxx._z; hz(0, 1) = dxy._z; hz(0, 2) = dxz._z;
  hz(1, 0) = dxy._z; hz(1, 1) = dyy._z; hz(1, 2) = dyz._z;
  hz(2, 0) = dxz._z; hz(2, 1) = dyz._z; hz(2, 2) = dzz._z;
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformation4D
::EvaluateHessian(Matrix hessian[3], int i, int j, int k, int l) const
{
  if (_FFD.IsInside(i, j, k, l)) mirtk::EvaluateHessian(&_CPImage, hessian, i, j, k, l);
  else                           mirtk::EvaluateHessian( _CPValue, hessian, i, j, k, l);
}

// -----------------------------------------------------------------------------
template <class CPImage>
void EvaluateHessian(const CPImage *coeff, Matrix hessian[3], double x, double y, double z, double t)
{
  typedef BSplineFreeFormTransformation4D::Kernel Kernel;

  int i = ifloor(x);
  int j = ifloor(y);
  int k = ifloor(z);
  int l = ifloor(t);

  const int A = Kernel::VariableToIndex(x - i);
  const int B = Kernel::VariableToIndex(y - j);
  const int C = Kernel::VariableToIndex(z - k);
  const int D = Kernel::VariableToIndex(t - l);

  typename CPImage::VoxelType dxx, dxy, dxz, dyy, dyz, dzz;
  double                      wx[3], wy[3], wz[3], wt;
  int                         ia, jb, kc, ld;

  --i, --j, --k, --l;
  for (int d = 0; d < 4; ++d) {
    ld = l + d;
    wt = Kernel::LookupTable[D][d];
    for (int c = 0; c < 4; ++c) {
      kc = k + c;
      wz[0] = Kernel::LookupTable   [C][c];
      wz[1] = Kernel::LookupTable_I [C][c];
      wz[2] = Kernel::LookupTable_II[C][c];
      for (int b = 0; b < 4; ++b) {
        jb = j + b;
        wy[0] = Kernel::LookupTable   [B][b];
        wy[1] = Kernel::LookupTable_I [B][b];
        wy[2] = Kernel::LookupTable_II[B][b];
        for (int a = 0; a < 4; ++a) {
          ia = i + a;
          wx[0] = Kernel::LookupTable   [A][a];
          wx[1] = Kernel::LookupTable_I [A][a];
          wx[2] = Kernel::LookupTable_II[A][a];
          dxx += (wx[2] * wy[0] * wz[0] * wt) * coeff->Get(ia, jb, kc, ld);
          dxy += (wx[1] * wy[1] * wz[0] * wt) * coeff->Get(ia, jb, kc, ld);
          dxz += (wx[1] * wy[0] * wz[1] * wt) * coeff->Get(ia, jb, kc, ld);
          dyy += (wx[0] * wy[2] * wz[0] * wt) * coeff->Get(ia, jb, kc, ld);
          dyz += (wx[0] * wy[1] * wz[1] * wt) * coeff->Get(ia, jb, kc, ld);
          dzz += (wx[0] * wy[0] * wz[2] * wt) * coeff->Get(ia, jb, kc, ld);
        }
      }
    }
  }

  Matrix &hx = hessian[0];
  hx.Initialize(3, 3);
  hx(0, 0) = dxx._x; hx(0, 1) = dxy._x; hx(0, 2) = dxz._x;
  hx(1, 0) = dxy._x; hx(1, 1) = dyy._x; hx(1, 2) = dyz._x;
  hx(2, 0) = dxz._x; hx(2, 1) = dyz._x; hx(2, 2) = dzz._x;

  Matrix &hy = hessian[1];
  hy.Initialize(3, 3);
  hy(0, 0) = dxx._y; hy(0, 1) = dxy._y; hy(0, 2) = dxz._y;
  hy(1, 0) = dxy._y; hy(1, 1) = dyy._y; hy(1, 2) = dyz._y;
  hy(2, 0) = dxz._y; hy(2, 1) = dyz._y; hy(2, 2) = dzz._y;

  Matrix &hz = hessian[2];
  hz.Initialize(3, 3);
  hz(0, 0) = dxx._z; hz(0, 1) = dxy._z; hz(0, 2) = dxz._z;
  hz(1, 0) = dxy._z; hz(1, 1) = dyy._z; hz(1, 2) = dyz._z;
  hz(2, 0) = dxz._z; hz(2, 1) = dyz._z; hz(2, 2) = dzz._z;
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformation4D
::EvaluateHessian(Matrix hessian[3], double x, double y, double z, double t) const
{
  if (_FFD.IsInside(x, y, z, t)) mirtk::EvaluateHessian(&_CPImage, hessian, x, y, z, t);
  else                           mirtk::EvaluateHessian( _CPValue, hessian, x, y, z, t);
}

// -----------------------------------------------------------------------------
template <class CPImage>
void EvaluateLaplacian(const CPImage *coeff, double laplacian[3], int i, int j, int k, int l)
{
  typedef BSplineFreeFormTransformation4D::Kernel Kernel;

  const double *w[3] = {
    Kernel::LatticeWeights,
    Kernel::LatticeWeights_I,
    Kernel::LatticeWeights_II
  };

  typename CPImage::VoxelType v = .0;
  int                         ia, jb, kc, ld;

  --i, --j, --k, --l;
  for (int d = 0; d < 3; ++d) {
    ld = l + d;
    for (int c = 0; c < 3; ++c) {
      kc = k + c;
      for (int b = 0; b < 3; ++b) {
        jb = j + b;
        for (int a = 0; a < 3; ++a) {
          ia = i + a;
          v += (w[2][a] * w[0][b] * w[0][c] * w[0][d]) * coeff->Get(ia, jb, kc, ld);
          v += (w[0][a] * w[2][b] * w[0][c] * w[0][d]) * coeff->Get(ia, jb, kc, ld);
          v += (w[0][a] * w[0][b] * w[2][c] * w[0][d]) * coeff->Get(ia, jb, kc, ld);
        }
      }
    }
  }

  laplacian[0] = v._x;
  laplacian[1] = v._y;
  laplacian[2] = v._z;
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformation4D
::EvaluateLaplacian(double laplacian[3], int i, int j, int k, int l) const
{
  if (_FFD.IsInside(i, j, k, l)) mirtk::EvaluateLaplacian(&_CPImage, laplacian, i, j, k, l);
  else                           mirtk::EvaluateLaplacian( _CPValue, laplacian, i, j, k, l);
}

// -----------------------------------------------------------------------------
template <class CPImage>
void EvaluateLaplacian(const CPImage *coeff, double laplacian[3], int i, int j, int k, double t)
{
  typedef BSplineFreeFormTransformation4D::Kernel Kernel;

  const double *w[3] = {
    Kernel::LatticeWeights,
    Kernel::LatticeWeights_I,
    Kernel::LatticeWeights_II
  };

  int       l = ifloor(t);
  const int D = Kernel::VariableToIndex(t - l);

  typename CPImage::VoxelType v = .0;
  double                      wt;
  int                         ia, jb, kc, ld;

  --i, --j, --k, --l;
  for (int d = 0; d < 4; ++d) {
    ld = l + d;
    wt = Kernel::LookupTable[D][d];
    for (int c = 0; c < 3; ++c) {
      kc = k + c;
      for (int b = 0; b < 3; ++b) {
        jb = j + b;
        for (int a = 0; a < 3; ++a) {
          ia = i + a;
          v += (w[2][a] * w[0][b] * w[0][c] * wt) * coeff->Get(ia, jb, kc, ld);
          v += (w[0][a] * w[2][b] * w[0][c] * wt) * coeff->Get(ia, jb, kc, ld);
          v += (w[0][a] * w[0][b] * w[2][c] * wt) * coeff->Get(ia, jb, kc, ld);
        }
      }
    }
  }

  laplacian[0] = v._x;
  laplacian[1] = v._y;
  laplacian[2] = v._z;
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformation4D
::EvaluateLaplacian(double laplacian[3], int i, int j, int k, double t) const
{
  if (_FFD.IsInside(i, j, k, t)) mirtk::EvaluateLaplacian(&_CPImage, laplacian, i, j, k, t);
  else                           mirtk::EvaluateLaplacian( _CPValue, laplacian, i, j, k, t);
}

// -----------------------------------------------------------------------------
template <class CPImage>
void EvaluateLaplacian(const CPImage *coeff, double &x, double &y, double &z, double t)
{
  typedef BSplineFreeFormTransformation4D::Kernel Kernel;

  int i = ifloor(x);
  int j = ifloor(y);
  int k = ifloor(z);
  int l = ifloor(t);

  const int A = Kernel::VariableToIndex(x - i);
  const int B = Kernel::VariableToIndex(y - j);
  const int C = Kernel::VariableToIndex(z - k);
  const int D = Kernel::VariableToIndex(t - l);

  typename CPImage::VoxelType v;
  double                      wx[2], wy[2], wz[2], wt;
  int                         ia, jb, kc, ld;

  --i, --j, --k, --l;
  for (int d = 0; d < 4; ++d) {
    ld = l + d;
    wt = Kernel::LookupTable[D][d];
    for (int c = 0; c < 4; ++c) {
      kc = k + c;
      wz[0] = Kernel::LookupTable   [C][c];
      wz[1] = Kernel::LookupTable_II[C][c];
      for (int b = 0; b < 4; ++b) {
        jb = j + b;
        wy[0] = Kernel::LookupTable   [B][b];
        wy[1] = Kernel::LookupTable_II[B][b];
        for (int a = 0; a < 4; ++a) {
          ia = i + a;
          wx[0] = Kernel::LookupTable   [A][a];
          wx[1] = Kernel::LookupTable_II[A][a];
          v += (wx[1] * wy[0] * wz[0] * wt) * coeff->Get(ia, jb, kc, ld);
          v += (wx[0] * wy[1] * wz[0] * wt) * coeff->Get(ia, jb, kc, ld);
          v += (wx[0] * wy[0] * wz[1] * wt) * coeff->Get(ia, jb, kc, ld);
        }
      }
    }
  }

  x = v._x, y = v._y, z = v._z;
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformation4D
::EvaluateLaplacian(double laplacian[3], double x, double y, double z, double t) const
{
  if (_FFD.IsInside(x, y, z, t)) mirtk::EvaluateLaplacian(&_CPImage, x, y, z, t);
  else                           mirtk::EvaluateLaplacian( _CPValue, x, y, z, t);
  laplacian[0] = x, laplacian[1] = y, laplacian[2] = z;
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformation4D
::EvaluateLaplacian(double &x, double &y, double &z, double t) const
{
  if (_FFD.IsInside(x, y, z, t)) mirtk::EvaluateLaplacian(&_CPImage, x, y, z, t);
  else                           mirtk::EvaluateLaplacian( _CPValue, x, y, z, t);
}

// =============================================================================
// Derivatives
// =============================================================================

// -----------------------------------------------------------------------------
class BSplineFreeFormTransformation4DParametricGradientBody
{
  typedef BSplineFreeFormTransformation4D::Kernel Kernel;

public:
  const BSplineFreeFormTransformation4D *_FFD;
  const GenericImage<double>            *_Input;
  const WorldCoordsImage                *_Image2World;
  const WorldCoordsImage                *_WorldCoords;
  double                                *_Output;
  double                                 _Weight;

  int _X, _Y, _N;
  double _t;  ///< Time corrresponding to input gradient image (in ms)
  int    _cl; ///< Current temporal control point index
  double _BL; ///< Pre-computed 1D B-Spline kernel value in temporal domain

  // ---------------------------------------------------------------------------
  void operator ()(const blocked_range3d<int> &re) const
  {
    int           cp, xdof, ydof, zdof;   // indices of DoFs corresponding to control point
    const double *gx, *gy, *gz;           // pointers to input gradient data
    double        x, y, z, jac;           // lattice coordinates and B-spline kernel
    FreeFormTransformation4D::CPStatus status;

    // With transformed world coordinates
    if (_WorldCoords) {
      const double *wx, *wy, *wz;
      double        x1, y1, z1, x2, y2, z2;
      // Loop over control points
      for (int ck = re.pages().begin(); ck != re.pages().end(); ++ck)
      for (int cj = re.rows ().begin(); cj != re.rows ().end(); ++cj)
      for (int ci = re.cols ().begin(); ci != re.cols ().end(); ++ci) {
        // Get index of the control point
        cp = _FFD->LatticeToIndex(ci, cj, ck, _cl);
        // Check if any DoF corresponding to the control point is active
        _FFD->GetStatus(cp, status);
        if (status._x == Passive && status._y == Passive && status._z == Passive) continue;
        // Calculate bounding box of control point in world coordinates
        _FFD->BoundingBox(cp, x1, y1, z1, x2, y2, z2, 1.0 / _FFD->SpeedupFactor());
        // Get indices of DoFs corresponding to the control point
        _FFD->IndexToDOFs(cp, xdof, ydof, zdof);
        // Loop over target voxels
        wx = _WorldCoords->Data(), wy = wx + _N, wz = wy + _N;
        gx = _Input      ->Data(), gy = gx + _N, gz = gy + _N;
        for (int i = 0; i < _N; ++i, ++wx, ++wy, ++wz, ++gx, ++gy, ++gz) {
          // Check whether reference point is valid
          if (*gx == .0 && *gy == .0 && *gz == .0) continue;
          // Check if world coordinate is in support region of control point
          if (x1 <= *wx && *wx <= x2 &&
              y1 <= *wy && *wy <= y2 &&
              z1 <= *wz && *wz <= z2) {
            // Convert point to lattice coordinates
            x = *wx, y = *wy, z = *wz;
            _FFD->WorldToLattice(x, y, z);
            // Convert non-parametric gradient into parametric gradient
            jac = Kernel::B(x - ci) * Kernel::B(y - cj) * Kernel::B(z - ck) * _BL;
            if (status._x == Active) _Output[xdof] += _Weight * jac * (*gx);
            if (status._y == Active) _Output[ydof] += _Weight * jac * (*gy);
            if (status._z == Active) _Output[zdof] += _Weight * jac * (*gz);
          }
        }
      }
    }
    // With pre-computed world coordinates
    else if (_Image2World) {
      const double *wx, *wy, *wz;
      int           i1, i2, j1, j2, k1, k2; // bounding box of target voxels affected by control point
      int           s2, s3;                 // stride of data/increment of data pointers
      // Loop over control points
      for (int ck = re.pages().begin(); ck != re.pages().end(); ++ck)
      for (int cj = re.rows ().begin(); cj != re.rows ().end(); ++cj)
      for (int ci = re.cols ().begin(); ci != re.cols ().end(); ++ci) {
        // Compute DoFs corresponding to the control point
        cp = _FFD->LatticeToIndex(ci, cj, ck, _cl);
        // Check if any DoF corresponding to the control point is active
        _FFD->GetStatus(cp, status);
        if (status._x == Passive && status._y == Passive && status._z == Passive) continue;
        // Calculate spatial bounding box of control point in image coordinates
        if (!_FFD->BoundingBox(_Input, cp, i1, j1, k1, i2, j2, k2, 1.0 / _FFD->SpeedupFactor())) continue;
        _FFD->IndexToDOFs(cp, xdof, ydof, zdof);
        // Loop over all voxels in the target (reference) volume
        s2 =  _X - (i2 - i1 + 1);
        s3 = (_Y - (j2 - j1 + 1)) * _X;
        wx = _Image2World->Data(i1, j1, k1), wy = wx + _N, wz = wy + _N;
        gx = _Input      ->Data(i1, j1, k1), gy = gx + _N, gz = gy + _N;
        for (int k = k1; k <= k2; k++, wx += s3, wy += s3, wz += s3, gx += s3, gy += s3, gz += s3)
        for (int j = j1; j <= j2; j++, wx += s2, wy += s2, wz += s2, gx += s2, gy += s2, gz += s2)
        for (int i = i1; i <= i2; i++, wx +=  1, wy +=  1, wz +=  1, gx +=  1, gy +=  1, gz +=  1) {
          // Check whether reference point is valid
          if (*gx == .0 && *gy == .0 && *gz == .0) continue;
          // Convert point to lattice coordinates
          x = *wx, y = *wy, z = *wz;
          _FFD->WorldToLattice(x, y, z);
          // Convert non-parametric gradient into parametric gradient
          jac = Kernel::B(x - ci) * Kernel::B(y - cj) * Kernel::B(z - ck) * _BL;
          if (status._x == Active) _Output[xdof] += _Weight * jac * (*gx);
          if (status._y == Active) _Output[ydof] += _Weight * jac * (*gy);
          if (status._z == Active) _Output[zdof] += _Weight * jac * (*gz);
        }
      }
    }
    // Without pre-computed world coordinates
    else {
      int i1, i2, j1, j2, k1, k2; // bounding box of target voxels affected by control point
      int s2, s3;                 // stride of data/increment of data pointers
      // Loop over control points
      for (int ck = re.pages().begin(); ck != re.pages().end(); ++ck)
      for (int cj = re.rows ().begin(); cj != re.rows ().end(); ++cj)
      for (int ci = re.cols ().begin(); ci != re.cols ().end(); ++ci) {
        // Compute DoFs corresponding to the control point
        cp = _FFD->LatticeToIndex(ci, cj, ck, _cl);
        // Check if any DoF corresponding to the control point is active
        _FFD->GetStatus(cp, status);
        if (status._x == Passive && status._y == Passive && status._z == Passive) continue;
        // Calculate spatial bounding box of control point in image coordinates
        if (!_FFD->BoundingBox(_Input, cp, i1, j1, k1, i2, j2, k2, 1.0 / _FFD->SpeedupFactor())) continue;
        // Get indices of DoFs corresponding to the control point
        _FFD->IndexToDOFs(cp, xdof, ydof, zdof);
        // Loop over all voxels in the target (reference) volume
        s2 =  _X - (i2 - i1 + 1);
        s3 = (_Y - (j2 - j1 + 1)) * _X;
        gx = _Input->Data(i1, j1, k1), gy = gx + _N, gz = gy + _N;
        for (int k = k1; k <= k2; k++, gx += s3, gy += s3, gz += s3)
        for (int j = j1; j <= j2; j++, gx += s2, gy += s2, gz += s2)
        for (int i = i1; i <= i2; i++, gx +=  1, gy +=  1, gz +=  1) {
          // Check whether reference point is valid
          if (*gx == .0 && *gy == .0 && *gz == .0) continue;
          // Convert voxel coordinates to world coordinates
          x = i, y = j, z = k;
          _Input->ImageToWorld(x, y, z);
          // Convert point to lattice coordinates
          _FFD->WorldToLattice(x, y, z);
          // Convert non-parametric gradient into parametric gradient
          jac = Kernel::B(x - ci) * Kernel::B(y - cj) * Kernel::B(z - ck) * _BL;
          if (status._x == Active) _Output[xdof] += _Weight * jac * (*gx);
          if (status._y == Active) _Output[ydof] += _Weight * jac * (*gy);
          if (status._z == Active) _Output[zdof] += _Weight * jac * (*gz);
        }
      }
    }
  }

  // ---------------------------------------------------------------------------
  void operator ()()
  {
    // Initialize members to often accessed data
    _X          = _Input->X();
    _Y          = _Input->Y();
    const int Z = _Input->Z();
    _N          = _X * _Y * Z;

    // Check input
    if (_Image2World && (_Image2World->X() != _X || _Image2World->Y() != _Y ||
                         _Image2World->Z() !=  Z || _Image2World->T() !=  3)) {
      cerr << "BSplineFreeFormTransformation4D::ParametricGradient: Invalid voxel coordinates map" << endl;
      exit(1);
    }
    if (_WorldCoords && (_WorldCoords->X() != _X || _WorldCoords->Y() != _Y ||
                         _WorldCoords->Z() !=  Z || _WorldCoords->T() !=  3)) {
      cerr << "BSplineFreeFormTransformation4D::ParametricGradient: Invalid world coordinates map" << endl;
      exit(1);
    }

    // Convert time to lattice units
    const double t = _FFD->TimeToLattice(_t);

    // Check for periodicity in temporal dimension
    ExtrapolationMode m = _FFD->ExtrapolationMode();
    const bool periodic = (m == ExtrapolationWithPeriodicTime(m));

    // Calculate parametric gradient
    for (_cl = 0; _cl < _FFD->T(); _cl++) {
      // Pre-compute B-spline kernel weight in t dimension
      _BL = Kernel::B(t - _cl);
      // Add weights corresponding to periodic repetition of control points
      if (periodic) {
        int    cl;
        double wt;
        // Periodic repetition in the past (i.e., cl < _cl)
        cl = _cl;
        do {
          cl  -= _FFD->T();
          _BL += (wt = Kernel::B(t - cl));
        } while (wt > .0);
        // Periodic repetition in the future (i.e., cl > _cl)
        cl = _cl;
        do {
          cl  += _FFD->T();
          _BL += (wt = Kernel::B(t - cl));
        } while (wt > .0);
      }
      // Skip if outside of support region - no temporal bounding box check required later
      if (_BL == .0) continue;
      // Otherwise, add derivatives w.r.t DoFs of all control points with temporal index _cl
      blocked_range3d<int> cps(0, _FFD->Z(), 0, _FFD->Y(), 0, _FFD->X());
      parallel_for(cps, BSplineFreeFormTransformation4DParametricGradientBody(*this));
    }
  }

}; // BSplineFreeFormTransformation4DParametricGradientBody

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformation4D
::ParametricGradient(const GenericImage<double> *in, double *out,
                     const WorldCoordsImage *i2w, const WorldCoordsImage *wc,
                     double, double w) const
{
  MIRTK_START_TIMING();
  BSplineFreeFormTransformation4DParametricGradientBody body;
  body._FFD         = this;
  body._Input       = in;
  body._Output      = out;
  body._Weight      = w;
  body._Image2World = i2w;
  body._WorldCoords = wc;
  body._t           = in->ImageToTime(.0);
  body();
  MIRTK_DEBUG_TIMING(2, "parametric gradient computation (4D B-spline FFD)");
}

// =============================================================================
// Properties
// =============================================================================

// -----------------------------------------------------------------------------
int BSplineFreeFormTransformation4D::KernelSize() const
{
  return 4;
}

// -----------------------------------------------------------------------------
double BSplineFreeFormTransformation4D::BendingEnergy(double x, double y, double z, double t, double, bool wrt_world) const
{
  // Convert to lattice coordinates
  this->WorldToLattice(x, y, z);
  t = this->TimeToLattice(t);
  // Calculate 2nd order derivatives
  Matrix hessian[3];
  EvaluateHessian(hessian, x, y, z, t);
  // Convert derivatives to world coordinates
  if (wrt_world) this->HessianToWorld(hessian);
  // Calculate bending energy
  return Bending3D(hessian);
}

// -----------------------------------------------------------------------------
double BSplineFreeFormTransformation4D::BendingEnergy(bool incl_passive, bool wrt_world) const
{
  int    nactive = 0;
  double bending = .0;
  Matrix hessian[3];

  for (int l = 0; l < _t; ++l)
  for (int k = 0; k < _z; ++k)
  for (int j = 0; j < _y; ++j)
  for (int i = 0; i < _x; ++i) {
    if (incl_passive || IsActive(i, j, k, l)) {
      EvaluateHessian(hessian, i, j, k, l);
      if (wrt_world) HessianToWorld(hessian);
      bending += Bending3D(hessian);
      ++nactive;
    }
  }

  if (nactive > 0) bending /= nactive;
  return bending;
}

// -----------------------------------------------------------------------------
double BSplineFreeFormTransformation4D::BendingEnergy(const ImageAttributes &attr, double, bool wrt_world) const
{
  const int N = attr.NumberOfLatticePoints();
  if (N == 0) return .0;

  double bending = .0;
  double x, y, z, t;
  Matrix hessian[3];

  for (int l = 0; l < attr._t; ++l) {
    t = this->TimeToLattice(attr.LatticeToTime(l));
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i) {
      x = i, y = j, z = k;
      attr .LatticeToWorld(x, y, z);
      this->WorldToLattice(x, y, z);
      EvaluateHessian(hessian, x, y, z, t);
      if (wrt_world) HessianToWorld(hessian);
      bending += Bending3D(hessian);
    }
  }

  return bending / N;
}


namespace BSplineFreeFormTransformation4DUtils {

// -----------------------------------------------------------------------------
/// Voxel function which evaluates the 2nd order derivatives of one component
/// of a B-spline FFD at control point lattice coordinates
struct Evaluate2ndOrderBSplineFFDDerivatives : public VoxelFunction
{
  typedef BSplineFreeFormTransformation4D::CPExtrapolator Extrapolator;
  typedef BSplineFreeFormTransformation4D::Vector         Vector;
  typedef BSplineFreeFormTransformation4D::Kernel         Kernel;

  const BSplineFreeFormTransformation4D *_Transformation; ///< B-spline free-form deformation
  const Extrapolator                    *_CPValue;        ///< Coefficients of B-spline FFD
  bool                                   _WrtWorld;       ///< Whether to compute derivatives
                                                              ///< w.r.t world or lattice coordinates

  /// Derivatives of 2nd order derivatives w.r.t control point parameters
  static double LookupTable[81][6];

  /// Initialize 4D lookup table of 2nd order derivatives w.r.t control point
  /// parameter for each voxel in the support region of a basis kernel
  static void InitializeLookupTable();

  /// Constructor
  Evaluate2ndOrderBSplineFFDDerivatives(const BSplineFreeFormTransformation4D *ffd,
                                        bool wrt_world = false)
  :
    _Transformation(ffd), _CPValue(ffd->Extrapolator()), _WrtWorld(wrt_world)
  {}

  /// Evaluate 2nd order derivatives of 4D FFD at given lattice coordinate
  void operator()(int i, int j, int k, int l,
                  Vector *dxx, Vector *dxy, Vector *dxz, Vector *dyy, Vector *dyz, Vector *dzz)
  {
    // Note: Derivatives are evaluated on a lattice that has an
    //       additional boundary margin of one voxel. Therefore,
    //       CP indices I, J, K, and L are shifted by an offset of -1.
    int n = 0;
    for (int L = l-2; L <= l; ++L)
    for (int K = k-2; K <= k; ++K)
    for (int J = j-2; J <= j; ++J)
    for (int I = i-2; I <= i; ++I, ++n) {
      *dxx += _CPValue->Get(I, J, K, L) * LookupTable[n][0];
      *dxy += _CPValue->Get(I, J, K, L) * LookupTable[n][1];
      *dxz += _CPValue->Get(I, J, K, L) * LookupTable[n][2];
      *dyy += _CPValue->Get(I, J, K, L) * LookupTable[n][3];
      *dyz += _CPValue->Get(I, J, K, L) * LookupTable[n][4];
      *dzz += _CPValue->Get(I, J, K, L) * LookupTable[n][5];
    }
    // Apply product and chain rule to convert derivatives to ones w.r.t the world
    if (_WrtWorld) {
      _Transformation->HessianToWorld(dxx->_x, dxy->_x, dxz->_x, dyy->_x, dyz->_x, dzz->_x);
      _Transformation->HessianToWorld(dxx->_y, dxy->_y, dxz->_y, dyy->_y, dyz->_y, dzz->_y);
      _Transformation->HessianToWorld(dxx->_z, dxy->_z, dxz->_z, dyy->_z, dyz->_z, dzz->_z);
    }
  }
};

// -----------------------------------------------------------------------------
double Evaluate2ndOrderBSplineFFDDerivatives::LookupTable[81][6] = {{.0}};
void Evaluate2ndOrderBSplineFFDDerivatives::InitializeLookupTable()
{
  static bool initialized = false;
  if (initialized) return;

  const double *w[3] = {
    Kernel::LatticeWeights,
    Kernel::LatticeWeights_I,
    Kernel::LatticeWeights_II
  };

  int n = 0;
  for (int d = 0; d < 3; ++d)
  for (int c = 0; c < 3; ++c)
  for (int b = 0; b < 3; ++b)
  for (int a = 0; a < 3; ++a, ++n) {
    LookupTable[n][0] = w[2][a] * w[0][b] * w[0][c] * w[0][d];
    LookupTable[n][1] = w[1][a] * w[1][b] * w[0][c] * w[0][d];
    LookupTable[n][2] = w[1][a] * w[0][b] * w[1][c] * w[0][d];
    LookupTable[n][3] = w[0][a] * w[2][b] * w[0][c] * w[0][d];
    LookupTable[n][4] = w[0][a] * w[1][b] * w[1][c] * w[0][d];
    LookupTable[n][5] = w[0][a] * w[0][b] * w[2][c] * w[0][d];
  }

  initialized = true;
}


} // namespace BSplineFreeFormTransformation4DUtils

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformation4D::BendingEnergyGradient(double *gradient, double weight, bool incl_passive, bool wrt_world) const
{
  using namespace BSplineFreeFormTransformation4DUtils;

  Vector sum;
  int    m, n;

  MIRTK_START_TIMING();

  // Pre-multiply weight by derivative of square function (2) and normalization factor
  const int ncps = this->NumberOfActiveCPs();
  if (ncps == 0) return;
  weight *= 2.0 / ncps;

  // Initialize static lookup table
  Evaluate2ndOrderBSplineFFDDerivatives::InitializeLookupTable();

  // Add a layer of boundary voxels to avoid additional boundary conditions
  ImageAttributes attr = this->Attributes();
  attr._x += 2, attr._y += 2, attr._z += 2, attr._t += 2;

  // Compute 2nd order derivatives w.r.t. lattice or world coordinates,
  // respectively, evaluated at each control point of the lattice
  GenericImage<Vector> dxx(attr);
  GenericImage<Vector> dxy(attr);
  GenericImage<Vector> dxz(attr);
  GenericImage<Vector> dyy(attr);
  GenericImage<Vector> dyz(attr);
  GenericImage<Vector> dzz(attr);

  Evaluate2ndOrderBSplineFFDDerivatives eval(this, wrt_world);
  ParallelForEachVoxel(attr, dxx, dxy, dxz, dyy, dyz, dzz, eval);

  // Compute 3rd order derivatives, twice w.r.t. lattice or world coordinate
  // and once w.r.t. transformation parameters of control point
  double w[81][9];
  if (wrt_world) {
    // Loop over support region (3x3x3x3) of a control point
    //
    // Note that the following terms are independent of the transformation
    // parameters and therefore can be pre-computed here as they are identical for
    // all control point positions at which the bending gradient is evaluated.
    for (n = 0; n < 81; ++n) {
      m = 0; // index of 2nd order derivative ( d^2 T(x[n]) / dx_i dx_j )
             // which is multiplied by weight w[n][m] according to the chain rule
      for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j, ++m) {
        w[n][m]  = Evaluate2ndOrderBSplineFFDDerivatives::LookupTable[n][0] * _matW2L(0, i) * _matW2L(0, j); // ( d^3 T(x[n]) / dudu dPhi ) * du/dx_i * du/dx_j
        w[n][m] += Evaluate2ndOrderBSplineFFDDerivatives::LookupTable[n][1] * _matW2L(0, i) * _matW2L(1, j); // ( d^3 T(x[n]) / dudv dPhi ) * du/dx_i * dv/dx_j
        w[n][m] += Evaluate2ndOrderBSplineFFDDerivatives::LookupTable[n][2] * _matW2L(0, i) * _matW2L(2, j); // ( d^3 T(x[n]) / dudw dPhi ) * du/dx_i * dw/dx_j
        w[n][m] += Evaluate2ndOrderBSplineFFDDerivatives::LookupTable[n][1] * _matW2L(1, i) * _matW2L(0, j); // ( d^3 T(x[n]) / dvdu dPhi ) * dv/dx_i * du/dx_j
        w[n][m] += Evaluate2ndOrderBSplineFFDDerivatives::LookupTable[n][3] * _matW2L(1, i) * _matW2L(1, j); // ( d^3 T(x[n]) / dvdv dPhi ) * dv/dx_i * dv/dx_j
        w[n][m] += Evaluate2ndOrderBSplineFFDDerivatives::LookupTable[n][4] * _matW2L(1, i) * _matW2L(2, j); // ( d^3 T(x[n]) / dvdw dPhi ) * dv/dx_i * dw/dx_j
        w[n][m] += Evaluate2ndOrderBSplineFFDDerivatives::LookupTable[n][2] * _matW2L(2, i) * _matW2L(0, j); // ( d^3 T(x[n]) / dwdu dPhi ) * dw/dx_i * du/dx_j
        w[n][m] += Evaluate2ndOrderBSplineFFDDerivatives::LookupTable[n][4] * _matW2L(2, i) * _matW2L(1, j); // ( d^3 T(x[n]) / dwdv dPhi ) * dw/dx_i * dv/dx_j
        w[n][m] += Evaluate2ndOrderBSplineFFDDerivatives::LookupTable[n][5] * _matW2L(2, i) * _matW2L(2, j); // ( d^3 T(x[n]) / dwdw dPhi ) * dw/dx_i * dw/dx_j
      }
      // Add weights of ( d^2 T(x[n]) / dx_i dx_j ) and ( d^2 T(x[n]) / dx_j dx_i )
      // for i != j as these derivatives are identical and re-order remaining weights.
      w[n][1] = w[n][1] + w[n][3];
      w[n][2] = w[n][2] + w[n][6];
      w[n][3] = w[n][4];
      w[n][4] = w[n][5] + w[n][7];
      w[n][5] = w[n][8];
    }
  } else {
    for (n = 0; n < 81; ++n) {
      w[n][0] =       Evaluate2ndOrderBSplineFFDDerivatives::LookupTable[n][0];
      w[n][1] = 2.0 * Evaluate2ndOrderBSplineFFDDerivatives::LookupTable[n][1];
      w[n][2] = 2.0 * Evaluate2ndOrderBSplineFFDDerivatives::LookupTable[n][2];
      w[n][3] =       Evaluate2ndOrderBSplineFFDDerivatives::LookupTable[n][3];
      w[n][4] = 2.0 * Evaluate2ndOrderBSplineFFDDerivatives::LookupTable[n][4];
      w[n][5] =       Evaluate2ndOrderBSplineFFDDerivatives::LookupTable[n][5];
    }
  }

  // Compute derivative of bending energy w.r.t each control point
  int xdof, ydof, zdof;
  for (int cl = 0; cl < _t; ++cl)
  for (int ck = 0; ck < _z; ++ck)
  for (int cj = 0; cj < _y; ++cj)
  for (int ci = 0; ci < _x; ++ci) {
    if (incl_passive || IsActive(ci, cj, ck, cl)) {
      sum = .0;
      // Loop over support region (3x3x3x3) of control point
      //
      // Note: Derivatives were evaluated on a lattice that has an
      //       additional boundary margin of one voxel. Therefore,
      //       indices i, j, k, and l are shifted by an offset of +1.
      n = 0;
      for (int l = cl; l <= cl+2; ++l)
      for (int k = ck; k <= ck+2; ++k)
      for (int j = cj; j <= cj+2; ++j)
      for (int i = ci; i <= ci+2; ++i, ++n) {
        sum += dxx(i, j, k, l) * w[n][0];
        sum += dxy(i, j, k, l) * w[n][1];
        sum += dxz(i, j, k, l) * w[n][2];
        sum += dyy(i, j, k, l) * w[n][3];
        sum += dyz(i, j, k, l) * w[n][4];
        sum += dzz(i, j, k, l) * w[n][5];
      }
      this->IndexToDOFs(this->LatticeToIndex(ci, cj, ck, cl), xdof, ydof, zdof);
      gradient[xdof] += weight * sum._x;
      gradient[ydof] += weight * sum._y;
      gradient[zdof] += weight * sum._z;
    }
  }

  MIRTK_DEBUG_TIMING(2, "bending gradient computation");
}

// =============================================================================
// I/O
// =============================================================================

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformation4D::Print(Indent indent) const
{
  cout << indent << "4D B-spline FFD:" << endl;
  FreeFormTransformation4D::Print(indent + 1);
}

// -----------------------------------------------------------------------------
bool BSplineFreeFormTransformation4D::CanRead(TransformationType format) const
{
  switch (format) {
    case TRANSFORMATION_BSPLINE_FFD_4D_v1:
    case TRANSFORMATION_BSPLINE_FFD_4D_v2:
    case TRANSFORMATION_BSPLINE_FFD_4D_v3:
      return true;
    default:
      return false;
  }
}


} // namespace mirtk
