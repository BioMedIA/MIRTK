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

#include "mirtk/ImageAttributes.h"

#include "mirtk/Math.h"
#include "mirtk/Array.h"
#include "mirtk/Memory.h"
#include "mirtk/Stream.h"
#include "mirtk/Indent.h"
#include "mirtk/Vector.h"
#include "mirtk/Vector3.h"
#include "mirtk/Matrix.h"
#include "mirtk/Matrix3x3.h"


namespace mirtk {


// -----------------------------------------------------------------------------
ImageAttributes::ImageAttributes()
{
  // Default image size
  _x  = 0;
  _y  = 0;
  _z  = 1;
  _t  = 1;

  // Default voxel size
  _dx = 1;
  _dy = 1;
  _dz = 1;
  _dt = 1;

  // Default origin
  _xorigin = 0;
  _yorigin = 0;
  _zorigin = 0;
  _torigin = 0;

  // Default x-axis
  _xaxis[0] = 1;
  _xaxis[1] = 0;
  _xaxis[2] = 0;

  // Default y-axis
  _yaxis[0] = 0;
  _yaxis[1] = 1;
  _yaxis[2] = 0;

  // Default z-axis
  _zaxis[0] = 0;
  _zaxis[1] = 0;
  _zaxis[2] = 1;

  // Default affine transformation
  _smat.Initialize(4, 4);
  _smat.Ident();

  _i2w = nullptr;
  _w2i = nullptr;
}

// -----------------------------------------------------------------------------
ImageAttributes::ImageAttributes(int x, int y, double dx, double dy)
{
  // Default image size
  _x  = x;
  _y  = y;
  _z  = 1;
  _t  = 1;

  // Default voxel size
  _dx = dx;
  _dy = dy;
  _dz = 1;
  _dt = 1;

  // Default origin
  _xorigin = 0;
  _yorigin = 0;
  _zorigin = 0;
  _torigin = 0;

  // Default x-axis
  _xaxis[0] = 1;
  _xaxis[1] = 0;
  _xaxis[2] = 0;

  // Default y-axis
  _yaxis[0] = 0;
  _yaxis[1] = 1;
  _yaxis[2] = 0;

  // Default z-axis
  _zaxis[0] = 0;
  _zaxis[1] = 0;
  _zaxis[2] = 1;

  // Default affine transformation
  _smat.Initialize(4, 4);
  _smat.Ident();

  _i2w = nullptr;
  _w2i = nullptr;
}

// -----------------------------------------------------------------------------
ImageAttributes::ImageAttributes(int x, int y, int z, double dx, double dy, double dz)
{
  _x  = x;
  _y  = y;
  _z  = z;
  _t  = 1;

  _dx = dx;
  _dy = dy;
  _dz = dz;
  _dt = 1;

  _xorigin = 0;
  _yorigin = 0;
  _zorigin = 0;
  _torigin = 0;

  _xaxis[0] = 1;
  _xaxis[1] = 0;
  _xaxis[2] = 0;

  _yaxis[0] = 0;
  _yaxis[1] = 1;
  _yaxis[2] = 0;

  _zaxis[0] = 0;
  _zaxis[1] = 0;
  _zaxis[2] = 1;

  _smat.Initialize(4, 4);
  _smat.Ident();

  _i2w = nullptr;
  _w2i = nullptr;
}

// -----------------------------------------------------------------------------
ImageAttributes::ImageAttributes(int x, int y, int z, int t, double dx, double dy, double dz, double dt)
{
  _x  = x;
  _y  = y;
  _z  = z;
  _t  = t;

  _dx = dx;
  _dy = dy;
  _dz = dz;
  _dt = dt;

  _xorigin = 0;
  _yorigin = 0;
  _zorigin = 0;
  _torigin = 0;

  _xaxis[0] = 1;
  _xaxis[1] = 0;
  _xaxis[2] = 0;

  _yaxis[0] = 0;
  _yaxis[1] = 1;
  _yaxis[2] = 0;

  _zaxis[0] = 0;
  _zaxis[1] = 0;
  _zaxis[2] = 1;

  _smat.Initialize(4, 4);
  _smat.Ident();

  _i2w = nullptr;
  _w2i = nullptr;
}

// -----------------------------------------------------------------------------
ImageAttributes::ImageAttributes(const ImageAttributes &attr)
{
  _x = attr._x;
  _y = attr._y;
  _z = attr._z;
  _t = attr._t;

  _dx = attr._dx;
  _dy = attr._dy;
  _dz = attr._dz;
  _dt = attr._dt;

  _xorigin = attr._xorigin;
  _yorigin = attr._yorigin;
  _zorigin = attr._zorigin;
  _torigin = attr._torigin;

  _xaxis[0] = attr._xaxis[0];
  _xaxis[1] = attr._xaxis[1];
  _xaxis[2] = attr._xaxis[2];

  _yaxis[0] = attr._yaxis[0];
  _yaxis[1] = attr._yaxis[1];
  _yaxis[2] = attr._yaxis[2];

  _zaxis[0] = attr._zaxis[0];
  _zaxis[1] = attr._zaxis[1];
  _zaxis[2] = attr._zaxis[2];

  _smat = attr._smat;

  // Do not copy pointers to pre-computed matrices as we don't know how
  // long the objects pointed to by the other attributes will be valid
  _i2w = nullptr;
  _w2i = nullptr;
}

// -----------------------------------------------------------------------------
ImageAttributes &ImageAttributes::operator =(const ImageAttributes &attr)
{
  if (this != &attr) {
    _x = attr._x;
    _y = attr._y;
    _z = attr._z;
    _t = attr._t;

    _dx = attr._dx;
    _dy = attr._dy;
    _dz = attr._dz;
    _dt = attr._dt;

    _xorigin = attr._xorigin;
    _yorigin = attr._yorigin;
    _zorigin = attr._zorigin;
    _torigin = attr._torigin;

    _xaxis[0] = attr._xaxis[0];
    _xaxis[1] = attr._xaxis[1];
    _xaxis[2] = attr._xaxis[2];

    _yaxis[0] = attr._yaxis[0];
    _yaxis[1] = attr._yaxis[1];
    _yaxis[2] = attr._yaxis[2];

    _zaxis[0] = attr._zaxis[0];
    _zaxis[1] = attr._zaxis[1];
    _zaxis[2] = attr._zaxis[2];

    _smat = attr._smat;
  }

  // Do not copy pointers to pre-computed matrices as we don't know how
  // long the objects pointed to by the other attributes will be valid
  //
  // Reset pointers even when assigning to itself such that behavior is
  // the same as when assigning from another instance
  _i2w = nullptr;
  _w2i = nullptr;

  return *this;
}

// -----------------------------------------------------------------------------
void ImageAttributes::LatticeToWorld(double *x, double *y, double *z) const
{
  int idx = 0;
  for (int k = 0; k < _z; ++k)
  for (int j = 0; j < _y; ++j)
  for (int i = 0; i < _x; ++i, ++idx) {
    x[idx] = i, y[idx] = j, z[idx] = k;
    LatticeToWorld(x[idx], y[idx], z[idx]);
  }
}

// -----------------------------------------------------------------------------
void ImageAttributes::LatticeToWorld(double *x, double *y, double *z, double *t) const
{
  const int xyz = _x * _y * _z;

  LatticeToWorld(x, y, z);

  if (_dt == .0) {
    double tl = LatticeToTime(0);
    for (int idx = 0; idx < xyz; ++idx) {
      (*t++) = tl;
    }
  } else {
    int idx = xyz;
    for (int l = 1; l < _t; ++l) {
      memcpy(x + idx, x, xyz * sizeof(double));
      memcpy(y + idx, y, xyz * sizeof(double));
      memcpy(z + idx, z, xyz * sizeof(double));
      idx += xyz;
    }
    for (int l = 0; l < _t; ++l) {
      double tl = LatticeToTime(l);
      for (int idx = 0; idx < xyz; ++idx) {
        (*t++) = tl;
      }
    }
  }
}

// -----------------------------------------------------------------------------
bool ImageAttributes::ContainsInSpace(const ImageAttributes &attr) const
{
  double x = 0, y = 0, z = 0;
  attr.LatticeToWorld(x, y, z);
  this->WorldToLattice(x, y, z);
  if (x <= -.5 || x >= _x - .5 || y <= -.5 || y >= _y - .5 || z <= -.5 || z >= _z - .5) return false;
  x = attr._x - 1, y = attr._y - 1, z = attr._z - 1;
  attr.LatticeToWorld(x, y, z);
  this->WorldToLattice(x, y, z);
  if (x <= -.5 || x >= _x - .5 || y <= -.5 || y >= _y - .5 || z <= -.5 || z >= _z - .5) return false;
  return true;
}

// -----------------------------------------------------------------------------
bool ImageAttributes::EqualInSpace(const ImageAttributes &attr) const
{
  // Spatial dimensions
  if (_x != attr._x || _y != attr._y || _z != attr._z) return false;
  // Spatial resolution
  if (!fequal(_dx, attr._dx) || !fequal(_dy, attr._dy) || !fequal(_dz, attr._dz)) return false;
  // Spatial orientation
  for (int i = 0; i < 3; ++i) if (!fequal(_xaxis[i], attr._xaxis[i])) return false;
  for (int i = 0; i < 3; ++i) if (!fequal(_yaxis[i], attr._yaxis[i])) return false;
  for (int i = 0; i < 3; ++i) if (!fequal(_zaxis[i], attr._zaxis[i])) return false;
  // Spatial origin
  if (!fequal(_xorigin, attr._xorigin) || !fequal(_yorigin, attr._yorigin) || !fequal(_zorigin, attr._zorigin)) return false;
  // Spatial transformation
  for (int r = 0; r < 4; ++r) {
    for (int c = 0; c < 4; ++c) {
      if (!fequal(_smat(r, c), attr._smat(r, c))) return false;
    }
  }
  return true;
}

// -----------------------------------------------------------------------------
bool ImageAttributes::EqualInTime(const ImageAttributes &attr) const
{
  return (_t == attr._t) && fequal(_dt, attr._dt) && fequal(_torigin, attr._torigin);
}

// -----------------------------------------------------------------------------
double ImageAttributes::Area() const
{
  return _x * _dx * _y * _dy;
}

// -----------------------------------------------------------------------------
double ImageAttributes::Volume() const
{
  return _x * _dx * _y * _dy * _z * _dz;
}

// -----------------------------------------------------------------------------
double ImageAttributes::Space() const
{
  if (static_cast<bool>(*this)) {
    if (_dz > 0.) return Volume();
    else          return Area();
  } else {
    return 0.;
  }
}

// -----------------------------------------------------------------------------
void ImageAttributes::Print(ostream &os, Indent indent) const
{
  // Change output stream settings
  const streamsize    w = os.width    (0);
  const streamsize    p = os.precision(5);
  const ios::fmtflags f = os.flags    ();
  // Print attributes of lattice
  bool sz = (_z > 1 && _dz);
  bool st = (_t > 1 && _dt);
  bool dz = (_dz && (_dz != 1.0 || _z > 1));
  bool dt = (_dt && (_dt != 1.0 || _t > 1));
  if (_t > 1 && !_dt) os << indent << "Channels: " << setw(10) << _t << "\n";
  os << indent << "Size:     "     << setw(10) << _x
                          << " x " << setw(10) << _y;
  if (sz || st) os        << " x " << setw(10) << _z;
  if (      st) os        << " x " << setw(10) << _t;
  os << "\n";
  os << indent << "Spacing:  "       << fixed << setw(10) << _dx
                            << " x " << fixed << setw(10) << _dy;
  if (dz || dt) os     << " x " << fixed << setw(10) << _dz;
  if (      dt) os     << " x " << fixed << setw(10) << _dt;
  os << "\n";
  os << indent << "Origin:  [" << fixed << setw(10) << _xorigin
                      << " , " << fixed << setw(10) << _yorigin
                      << " , " << fixed << setw(10) << _zorigin
                      << " , " << fixed << setw(10) << _torigin
                      <<   "]\n";
  os << indent << "X-axis:  [" << fixed << setw(10) << _xaxis[0]
                      << " , " << fixed << setw(10) << _xaxis[1]
                      << " , " << fixed << setw(10) << _xaxis[2]
                      <<   "]\n";
  os << indent << "Y-axis:  [" << fixed << setw(10) << _yaxis[0]
                      << " , " << fixed << setw(10) << _yaxis[1]
                      << " , " << fixed << setw(10) << _yaxis[2]
                      <<   "]\n";
  os << indent << "Z-axis:  [" << fixed << setw(10) << _zaxis[0]
                      << " , " << fixed << setw(10) << _zaxis[1]
                      << " , " << fixed << setw(10) << _zaxis[2]
                      <<   "]\n";
  if (!_smat.IsIdentity()) {
    os << indent << "Affine";
    _smat.Print(os, indent);
  }
  // Restore output stream settings
  os.width    (w);
  os.precision(p);
  os.flags    (f);
}

// -----------------------------------------------------------------------------
void ImageAttributes::Print(Indent indent) const
{
  Print(cout, indent);
}

// -----------------------------------------------------------------------------
Matrix ImageAttributes::GetLatticeToWorldMatrix() const
{
  // Final mat = A * T * R * S * T0
  Matrix mat(4, 4), tmp(4, 4);

  // T0: Translate image origin
  mat.Ident();
  mat(0, 3) = - (_x - 1) / 2.0;
  mat(1, 3) = - (_y - 1) / 2.0;
  mat(2, 3) = - (_z - 1) / 2.0;

  // S: Convert to world units
  tmp.Ident();
  tmp(0, 0) = _dx;
  tmp(1, 1) = _dy;
  tmp(2, 2) = _dz;
  mat = tmp * mat;

  // R: Orientation
  tmp(0, 0) = _xaxis[0];
  tmp(1, 0) = _xaxis[1];
  tmp(2, 0) = _xaxis[2];
  tmp(0, 1) = _yaxis[0];
  tmp(1, 1) = _yaxis[1];
  tmp(2, 1) = _yaxis[2];
  tmp(0, 2) = _zaxis[0];
  tmp(1, 2) = _zaxis[1];
  tmp(2, 2) = _zaxis[2];
  mat = tmp * mat;

  // T: Translate world origin
  tmp.Ident();
  tmp(0, 3) = _xorigin;
  tmp(1, 3) = _yorigin;
  tmp(2, 3) = _zorigin;
  mat = tmp * mat;

  // A: Transform world coordinate by optional affine transformation
  //    to align (scanner) coordinates with another anatomy
  mat = _smat * mat;

  return mat;
}

// -----------------------------------------------------------------------------
Matrix ImageAttributes::GetWorldToLatticeMatrix() const
{
  // Final mat = inv(A * T * R * S * T0)
  //           = inv(T0) * inv(S) * inv(R) * inv(T) * inv(A),
  Matrix mat(4, 4), tmp(4, 4);

  // inv(A): Transform world coordinate by optional inverse of the affine
  //         transformation which aligns (scanner) coordinates with another anatomy
  if (_smat.IsIdentity()) {
    mat.Ident(); // avoid numerical imprecisions caused by Matrix::Invert
  } else {
    mat = _smat.Inverse();
  }

  // inv(T): Translate world origin
  tmp.Ident();
  tmp(0, 3) = - _xorigin;
  tmp(1, 3) = - _yorigin;
  tmp(2, 3) = - _zorigin;
  mat = tmp * mat;

  // inv(R): Orientation
  tmp(0, 0) = _xaxis[0];
  tmp(0, 1) = _xaxis[1];
  tmp(0, 2) = _xaxis[2];
  tmp(0, 3) = .0;
  tmp(1, 0) = _yaxis[0];
  tmp(1, 1) = _yaxis[1];
  tmp(1, 2) = _yaxis[2];
  tmp(1, 3) = .0;
  tmp(2, 0) = _zaxis[0];
  tmp(2, 1) = _zaxis[1];
  tmp(2, 2) = _zaxis[2];
  tmp(2, 3) = .0;
  mat = tmp * mat;

  // inv(S): Convert to voxel units
  tmp.Ident();
  if (_dx) tmp(0, 0) = 1.0 / _dx;
  if (_dy) tmp(1, 1) = 1.0 / _dy;
  if (_dz) tmp(2, 2) = 1.0 / _dz;
  mat = tmp * mat;

  // inv(T0): Translate image origin
  tmp(0, 0) = tmp(1, 1) = tmp(2, 2) = 1.0;
  tmp(0, 3) = (_x - 1) / 2.0;
  tmp(1, 3) = (_y - 1) / 2.0;
  tmp(2, 3) = (_z - 1) / 2.0;
  mat = tmp * mat;

  return mat;
}

// -----------------------------------------------------------------------------
void ImageAttributes::PutAffineMatrix(const Matrix &m, bool apply)
{
  if (m.IsIdentity()) {
    _smat.Ident();
  } else if (apply) {

    // Lattice to world matrix is: A * T * R * S * T0
    const Matrix B = GetLatticeToWorldMatrix();

    // T0: Translate image origin
    Matrix T0(4, 4);
    T0.Ident();
    T0(0, 3) = - (_x - 1) / 2.;
    T0(1, 3) = - (_y - 1) / 2.;
    T0(2, 3) = - (_z - 1) / 2.;

    // Get new lattice to world matrix, excluding T0
    Matrix C = m * B * T0.Inverse();

    // Decompose C into T * R * S * K components
    // (K: shear, S: scale, R: rotation, T: translation)
    double tx, ty, tz, rx, ry, rz, sx, sy, sz, sxy, sxz, syz;
    MatrixToAffineParameters(C, tx, ty, tz, rx, ry, rz, sx, sy, sz, sxy, sxz, syz);

    // Remove shearing part
    const Matrix K = AffineParametersToMatrix(0., 0., 0., 0., 0., 0., 1., 1., 1, sxy, sxz, syz);
    const bool has_shearing = !K.IsIdentity();
    if (has_shearing) {
      C = C * K.Inverse();
    }

    // Decompose again, this time without shearing
    MatrixToAffineParameters(C, tx, ty, tz, rx, ry, rz, sx, sy, sz, sxy, sxz, syz);

    // Voxel size must be positive
    Matrix flip(4, 4);
    flip.Ident();
    if (sx < 0.) {
      flip(0, 0) = -1.;
      sx = -sx;
    }
    if (sy < 0.) {
      flip(1, 1) = -1.;
      sy = -sy;
    }
    if (sz < 0.) {
      flip(2, 2) = -1.;
      sz = -sz;
    }

    // Set new lattice attributes
    _dx = sx;
    _dy = sy;
    _dz = sz;

    // R: Orientation
    const Matrix R = RigidParametersToMatrix(0., 0., 0., rx, ry, rz) * flip;
    _xaxis[0] = R(0, 0);
    _xaxis[1] = R(1, 0);
    _xaxis[2] = R(2, 0);
    _yaxis[0] = R(0, 1);
    _yaxis[1] = R(1, 1);
    _yaxis[2] = R(2, 1);
    _zaxis[0] = R(0, 2);
    _zaxis[1] = R(1, 2);
    _zaxis[2] = R(2, 2);

    // T: Translate world origin
    _xorigin = tx;
    _yorigin = ty;
    _zorigin = tz;

    // Set remaining difference due to shearing as post-image to world mapping
    if (has_shearing) {
      _smat = C * K * T0 * GetWorldToLatticeMatrix();
    }
  } else {
    _smat = m;
  }
}

// =============================================================================
// Image domain helpers
// =============================================================================

// -----------------------------------------------------------------------------
/// Adjust size of output domain if necessary such that all corners of the
/// input image domain are within the bounds of the output image domain
///
/// \param[in]     in  Image grid that should fit into \p out domain.
/// \param[in,out] out Image grid which should fully contain \p in without
///                    changing its orientation nor voxel size.
void AdjustFieldOfView(const ImageAttributes &in, ImageAttributes &out)
{
  const Matrix i2w = in .GetImageToWorldMatrix();
  Matrix       w2i = out.GetWorldToImageMatrix();
  const int di = (in._x > 1 ? in._x - 1 : 1);
  const int dj = (in._y > 1 ? in._y - 1 : 1);
  const int dk = (in._z > 1 ? in._z - 1 : 1);
  for (int k1 = 0; k1 < in._z; k1 += dk)
  for (int j1 = 0; j1 < in._y; j1 += dj)
  for (int i1 = 0; i1 < in._x; i1 += di) {
    // Convert corner voxel indices of input domain to world coordinates
    double x  = i2w(0, 0) * i1 + i2w(0, 1) * j1 + i2w(0, 2) * k1 + i2w(0, 3);
    double y  = i2w(1, 0) * i1 + i2w(1, 1) * j1 + i2w(1, 2) * k1 + i2w(1, 3);
    double z  = i2w(2, 0) * i1 + i2w(2, 1) * j1 + i2w(2, 2) * k1 + i2w(2, 3);
    // Convert world coordinates to voxel index w.r.t output domain
    double i2 = w2i(0, 0) * x  + w2i(0, 1) * y  + w2i(0, 2) * z  + w2i(0, 3);
    double j2 = w2i(1, 0) * x  + w2i(1, 1) * y  + w2i(1, 2) * z  + w2i(1, 3);
    double k2 = w2i(2, 0) * x  + w2i(2, 1) * y  + w2i(2, 2) * z  + w2i(2, 3);
    // If point is outside the output domain...
    if (i2 < 0 || i2 > out._x-1 || j2 < 0 || j2 > out._y-1 || k2 < 0 || k2 > out._z-1) {
      // Increase output domain by difference along the axes directions
      if (i2 < 0) {
        int dv = iceil(abs(i2));
        out._x       += dv;
        out._xorigin -= out._xaxis[0] * out._dx * dv / 2.0;
        out._yorigin -= out._xaxis[1] * out._dx * dv / 2.0;
        out._zorigin -= out._xaxis[2] * out._dx * dv / 2.0;
      } else if (i2 > out._x - 1) {
        int dv = iceil(i2 - out._x + 1);
        out._x       += dv;
        out._xorigin += out._xaxis[0] * out._dx * dv / 2.0;
        out._yorigin += out._xaxis[1] * out._dx * dv / 2.0;
        out._zorigin += out._xaxis[2] * out._dx * dv / 2.0;
      }
      if (j2 < 0) {
        int dv = iceil(abs(j2));
        out._y       += dv;
        out._xorigin -= out._yaxis[0] * out._dy * dv / 2.0;
        out._yorigin -= out._yaxis[1] * out._dy * dv / 2.0;
        out._zorigin -= out._yaxis[2] * out._dy * dv / 2.0;
      } else if (j2 > out._y - 1) {
        int dv = iceil(j2 - out._y + 1);
        out._y       += dv;
        out._xorigin += out._yaxis[0] * out._dy * dv / 2.0;
        out._yorigin += out._yaxis[1] * out._dy * dv / 2.0;
        out._zorigin += out._yaxis[2] * out._dy * dv / 2.0;
      }
      if (k2 < 0) {
        int dv = iceil(abs(k2));
        out._z       += dv;
        out._xorigin -= out._zaxis[0] * out._dz * dv / 2.0;
        out._yorigin -= out._zaxis[1] * out._dz * dv / 2.0;
        out._zorigin -= out._zaxis[2] * out._dz * dv / 2.0;
      } else if (k2 > out._z - 1) {
        int dv = iceil(k2 - out._z + 1);
        out._z       += dv;
        out._xorigin += out._zaxis[0] * out._dz * dv / 2.0;
        out._yorigin += out._zaxis[1] * out._dz * dv / 2.0;
        out._zorigin += out._zaxis[2] * out._dz * dv / 2.0;
      }
      // Update transformation matrix
      w2i = out.GetWorldToImageMatrix();
    }
  }
}

// -----------------------------------------------------------------------------
void GetCorners(const ImageAttributes &attr, double *corners[3])
{
  int idx = 0;
  for (int k = 0; k <= 1; ++k)
  for (int j = 0; j <= 1; ++j)
  for (int i = 0; i <= 1; ++i, ++idx) {
    corners[idx][0] = i * (attr._x - 1);
    corners[idx][1] = j * (attr._y - 1);
    corners[idx][2] = k * (attr._z - 1);
    attr.LatticeToWorld(corners[idx][0], corners[idx][1], corners[idx][2]);
  }
}

// -----------------------------------------------------------------------------
void GetCorners(const ImageAttributes &attr, double corners[][3])
{
  int idx = 0;
  for (int k = 0; k <= 1; ++k)
  for (int j = 0; j <= 1; ++j)
  for (int i = 0; i <= 1; ++i, ++idx) {
    corners[idx][0] = i * (attr._x - 1);
    corners[idx][1] = j * (attr._y - 1);
    corners[idx][2] = k * (attr._z - 1);
    attr.LatticeToWorld(corners[idx][0], corners[idx][1], corners[idx][2]);
  }
}

// -----------------------------------------------------------------------------
ImageAttributes OrthogonalFieldOfView(const ImageAttributes &attr)
{
  ImageAttributes out = attr;
  // Apply any additional affine transformation to the attributes of the image domain
  if (!attr._smat.IsIdentity()) {
    // Reset affine transformation of output attributes
    out._smat.Ident();
    // Get corners of transformed image lattice
    double corners[8][3];
    GetCorners(attr, corners);
    // Set output origin
    out._xorigin = out._yorigin = out._zorigin = .0;
    for (int i = 0; i < 8; ++i) {
      out._xorigin += corners[i][0];
      out._yorigin += corners[i][1];
      out._zorigin += corners[i][2];
    }
    out._xorigin /= 8, out._yorigin /= 8, out._zorigin /= 8;
    // Set output orientation
    Matrix3x3 covar(.0);
    for (int i = 0; i < 8; ++i) {
      corners[i][0] -= out._xorigin;
      corners[i][1] -= out._yorigin;
      corners[i][2] -= out._zorigin;
      for (int r = 0; r < 3; ++r)
      for (int c = 0; c < 3; ++c) {
        covar[r][c] += corners[i][r] * corners[i][c];
      }
    }
    double  eigen[3];
    Vector3 axis [3];
    covar.EigenSolveSymmetric(eigen, axis);
    Vector3::GenerateOrthonormalBasis(axis[0], axis[1], axis[2]);
    for (int d = 0; d < 3; ++d) {
      out._xaxis[d] = axis[0][d];
      out._yaxis[d] = axis[1][d];
      out._zaxis[d] = axis[2][d];
    }
    // Set output size
    out._x = out._y = out._z = 1;
    AdjustFieldOfView(attr, out);
  }
  return out;
}

// -----------------------------------------------------------------------------
ImageAttributes OverallFieldOfView(const Array<ImageAttributes> &attr)
{
  ImageAttributes out;
  const int N = static_cast<int>(attr.size());
  if (N == 0) return out;
  if (N == 1) return attr[0];
  // Set output resolution
  out._dx = out._dy = out._dz = .0;
  for (int n = 0; n < N; ++n) {
    out._dx += attr[n]._dx;
    out._dy += attr[n]._dy;
    out._dz += attr[n]._dz;
  }
  out._dx /= N, out._dy /= N, out._dz /= N;
  // Set output origin
  out._xorigin = out._yorigin = out._zorigin = .0;
  for (int n = 0; n < N; ++n) {
    out._xorigin += attr[n]._xorigin;
    out._yorigin += attr[n]._yorigin;
    out._zorigin += attr[n]._zorigin;
  }
  out._xorigin /= N, out._yorigin /= N, out._zorigin /= N;
  // Set output orientation
  bool have_same_orientation = true;
  const Matrix R = attr[0].GetLatticeToWorldOrientation();
  for (int n = 1; n < N; ++n) {
    if (attr[n].GetLatticeToWorldOrientation() != R) {
      have_same_orientation = false;
      break;
    }
  }
  if (have_same_orientation) {
    for (int d = 0; d < 3; ++d) {
        out._xaxis[d] = attr[0]._xaxis[d];
        out._yaxis[d] = attr[0]._yaxis[d];
        out._zaxis[d] = attr[0]._zaxis[d];
    }
  } else {
    const int ncorners = 8 * N;
    double   **corners = Allocate<double>(3, ncorners);
    for (int n = 0; n < N; ++n) {
      GetCorners(attr[n], corners + 8 * n);
    }
    Matrix3x3 covar(.0);
    for (int i = 0; i < ncorners; ++i) {
      corners[i][0] -= out._xorigin;
      corners[i][1] -= out._yorigin;
      corners[i][2] -= out._zorigin;
      for (int r = 0; r < 3; ++r)
      for (int c = 0; c < 3; ++c) {
        covar[r][c] += corners[i][r] * corners[i][c];
      }
    }
    Deallocate(corners);
    double  eigen[3];
    Vector3 axis [3];
    covar.EigenSolveSymmetric(eigen, axis);
    Vector3::GenerateOrthonormalBasis(axis[0], axis[1], axis[2]);
    for (int d = 0; d < 3; ++d) {
      out._xaxis[d] = axis[0][d];
      out._yaxis[d] = axis[1][d];
      out._zaxis[d] = axis[2][d];
    }
  }
  // Set output size
  out._x = out._y = out._z = 1;
  for (int n = 0; n < N; ++n) {
    AdjustFieldOfView(attr[n], out);
  }
  return out;
}


} // namespace mirtk
