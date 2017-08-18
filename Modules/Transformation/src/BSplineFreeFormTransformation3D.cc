/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2017 Imperial College London
 * Copyright 2008-2013 Daniel Rueckert, Julia Schnabel
 * Copyright 2013-2017 Andreas Schuh
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

#include "mirtk/BSplineFreeFormTransformation3D.h"

#include "mirtk/Math.h"
#include "mirtk/Memory.h"
#include "mirtk/Profiling.h"
#include "mirtk/LinearInterpolateImageFunction3D.h"
#include "mirtk/ImageToInterpolationCoefficients.h"
#include "mirtk/CubicBSplineConvolution.h"


// The approximation of B-spline coefficients using parallel_reduce gives slightly
// different results depending on the number of threads; also, it did not demonstrate
// to be faster than a single-threaded execution. Thus, this code is disabled here,
// but kept for future reference.
#define _MIRTK_BSPLINEFFD_PARALLEL_APPROXIMATION 0


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
BSplineFreeFormTransformation3D::BSplineFreeFormTransformation3D()
:
  FreeFormTransformation3D(_FFD, &_FFD2D),
  _ParametricGradientCalculation(PG_Default)
{
}

// -----------------------------------------------------------------------------
BSplineFreeFormTransformation3D
::BSplineFreeFormTransformation3D(double x1, double y1, double z1,
                                  double x2, double y2, double z2,
                                  double dx, double dy, double dz,
                                  double *xaxis, double *yaxis, double *zaxis)
:
  FreeFormTransformation3D(_FFD, &_FFD2D),
  _ParametricGradientCalculation(PG_Default)
{
  Initialize(DefaultAttributes(x1, y1, z1, x2, y2, z2, dx, dy, dz, xaxis, yaxis, zaxis));
}

// -----------------------------------------------------------------------------
BSplineFreeFormTransformation3D
::BSplineFreeFormTransformation3D(const ImageAttributes &attr, double dx, double dy, double dz)
:
  FreeFormTransformation3D(_FFD, &_FFD2D),
  _ParametricGradientCalculation(PG_Default)
{
  Initialize(attr, dx, dy, dz);
}

// -----------------------------------------------------------------------------
BSplineFreeFormTransformation3D
::BSplineFreeFormTransformation3D(const BaseImage &target, double dx, double dy, double dz)
:
  FreeFormTransformation3D(_FFD, &_FFD2D),
  _ParametricGradientCalculation(PG_Default)
{
  Initialize(target.Attributes(), dx, dy, dz);
}

// -----------------------------------------------------------------------------
BSplineFreeFormTransformation3D
::BSplineFreeFormTransformation3D(const GenericImage<double> &image, bool disp)
:
  FreeFormTransformation3D(_FFD, &_FFD2D),
  _ParametricGradientCalculation(PG_Default)
{
  Initialize(image, disp);
}

// -----------------------------------------------------------------------------
BSplineFreeFormTransformation3D
::BSplineFreeFormTransformation3D(const BSplineFreeFormTransformation3D &ffd)
:
  FreeFormTransformation3D(ffd, _FFD, &_FFD2D),
  _ParametricGradientCalculation(PG_Default)
{
  if (_NumberOfDOFs > 0) InitializeInterpolator();
}

// -----------------------------------------------------------------------------
BSplineFreeFormTransformation3D::~BSplineFreeFormTransformation3D()
{
}

// =============================================================================
// Parameters (non-DoFs)
// =============================================================================

// -----------------------------------------------------------------------------
bool BSplineFreeFormTransformation3D::Set(const char *name, const char *value)
{
  if (strcmp(name, "BSplineFFD gradient calculation")   == 0 ||
      strcmp(name, "BSpline FFD gradient calculation")  == 0 ||
      strcmp(name, "B-spline FFD gradient calculation") == 0) {
    return FromString(value, _ParametricGradientCalculation);
  }
  return FreeFormTransformation3D::Set(name, value);
}

// -----------------------------------------------------------------------------
ParameterList BSplineFreeFormTransformation3D::Parameter() const
{
  ParameterList params = FreeFormTransformation3D::Parameter();
  Insert(params, "B-spline FFD gradient calculationn", _ParametricGradientCalculation);
  return params;
}

// =============================================================================
// Approximation/Interpolation
// =============================================================================

#if _MIRTK_BSPLINEFFD_PARALLEL_APPROXIMATION
namespace {

/// 2D cubic B-spline scattered data approximation
template <class Kernel>
class ApproximateSplineCoefficients2D
{
  typedef BSplineFreeFormTransformation3D::Vector Vector;

  // Scattered data to approximate
  int _NumberOfPoints;
  const double *_PointX;
  const double *_PointY;
  const double *_PointZ;
  const double *_ValueX;
  const double *_ValueY;
  const double *_ValueZ;

  // Cubic B-spline function
  const ImageAttributes  *_Attr;
  Vector                **_Data;
  double                **_Norm;

public:

  /// Default constructor
  ApproximateSplineCoefficients2D(const ImageAttributes *attr,
                                  const double *wx, const double *wy, const double *wz,
                                  const double *dx, const double *dy, const double *dz, int no)
  :
    _NumberOfPoints(no),
    _PointX(wx), _PointY(wy), _PointZ(wz),
    _ValueX(dx), _ValueY(dy), _ValueZ(dz),
    _Attr(attr)
  {
    _Data = CAllocate<Vector>(_Attr->X(), _Attr->Y());
    _Norm = CAllocate<double>(_Attr->X(), _Attr->Y());
  }

  /// Split constructor
  ApproximateSplineCoefficients2D(const ApproximateSplineCoefficients2D &rhs, split)
  :
    _NumberOfPoints(rhs._NumberOfPoints),
    _PointX(rhs._PointX), _PointY(rhs._PointY), _PointZ(rhs._PointZ),
    _ValueX(rhs._ValueX), _ValueY(rhs._ValueY), _ValueZ(rhs._ValueZ),
    _Attr(rhs._Attr)
  {
    _Data = CAllocate<Vector>(_Attr->X(), _Attr->Y());
    _Norm = CAllocate<double>(_Attr->X(), _Attr->Y());
  }

  /// Destructor
  ~ApproximateSplineCoefficients2D()
  {
    Deallocate(_Data);
    Deallocate(_Norm);
  }

  /// Joint results of parallel reduction
  void join(const ApproximateSplineCoefficients2D &rhs)
  {
    for (int j = 0; j < _Attr->Y(); ++j)
    for (int i = 0; i < _Attr->X(); ++i) {
      _Data[j][i] += rhs._Data[j][i];
      _Norm[j][i] += rhs._Norm[j][i];
    }
  }

  /// Add coefficients approximating subset of values
  void operator ()(const blocked_range<int> &re)
  {
    int    i, j, ci, cj, A, B;
    double x, y, z, w[2], basis, sum;

    const double *wx = _PointX + re.begin();
    const double *wy = _PointY + re.begin();
    const double *wz = _PointZ + re.begin();

    const double *dx = _ValueX + re.begin();
    const double *dy = _ValueY + re.begin();
    const double *dz = _ValueZ + re.begin();

    for (int n = re.begin(); n != re.end(); ++n, ++wx, ++wy, ++wz, ++dx, ++dy, ++dz) {
      if (AreEqual(*dx, 0.) && AreEqual(*dy, 0.) && AreEqual(*dz, 0.)) continue;

      x = *wx, y = *wy, z = *wz;
      _Attr->WorldToLattice(x, y, z);

      i = ifloor(x);
      j = ifloor(y);

      A = Kernel::VariableToIndex(x - i);
      B = Kernel::VariableToIndex(y - j);
      --i, --j;

      sum = 0.;
      for (int b = 0; b <= 3; ++b) {
        w[1] = Kernel::LookupTable[B][b];
        for (int a = 0; a <= 3; ++a) {
          w[0] = Kernel::LookupTable[A][a] * w[1];
          sum += w[0] * w[0];
        }
      }

      for (int b = 0; b <= 3; ++b) {
        cj = j + b;
        if (cj < 0 || cj >= _Attr->Y()) continue;
        w[1] = Kernel::LookupTable[B][b];
        for (int a = 0; a <= 3; ++a) {
          ci = i + a;
          if (ci < 0 || ci >= _Attr->X()) continue;
          w[0]  = Kernel::LookupTable[A][a] * w[1];
          basis = w[0] * w[0];
          _Norm[cj][ci] += basis;
          basis *= w[0] / sum;
          _Data[cj][ci]._x += basis * (*dx);
          _Data[cj][ci]._y += basis * (*dy);
          _Data[cj][ci]._z += basis * (*dz);
        }
      }
    }
  }

  void Approximate()
  {
    parallel_reduce(blocked_range<int>(0, _NumberOfPoints, _NumberOfPoints / 4), *this);
    for (int j = 0; j < _Attr->Y(); ++j)
    for (int i = 0; i < _Attr->X(); ++i) {
      if (!AreEqual(_Norm[j][i], 0.)) {
        _Data[j][i] /= _Norm[j][i];
      }
    }
  }

  const Vector &Get(int cp) const
  {
    return *(_Data[0] + cp);
  }

  const Vector &Get(int i, int j) const
  {
    return *(_Data[j][i]);
  }
};

/// 3D cubic B-spline scattered data approximation
template <class Kernel>
class ApproximateSplineCoefficients3D
{
  typedef BSplineFreeFormTransformation3D::Vector Vector;

  // Scattered data to approximate
  int _NumberOfPoints;
  const double *_PointX;
  const double *_PointY;
  const double *_PointZ;
  const double *_ValueX;
  const double *_ValueY;
  const double *_ValueZ;

  // Cubic B-spline function
  const ImageAttributes   *_Attr;
  Vector                ***_Data;
  double                ***_Norm;

public:

  /// Default constructor
  ApproximateSplineCoefficients3D(const ImageAttributes *attr,
                                  const double *wx, const double *wy, const double *wz,
                                  const double *dx, const double *dy, const double *dz, int no)
  :
    _NumberOfPoints(no),
    _PointX(wx), _PointY(wy), _PointZ(wz),
    _ValueX(dx), _ValueY(dy), _ValueZ(dz),
    _Attr(attr)
  {
    _Data = CAllocate<Vector>(_Attr->X(), _Attr->Y(), _Attr->Z());
    _Norm = CAllocate<double>(_Attr->X(), _Attr->Y(), _Attr->Z());
  }

  /// Split constructor
  ApproximateSplineCoefficients3D(const ApproximateSplineCoefficients3D &rhs, split)
  :
    _NumberOfPoints(rhs._NumberOfPoints),
    _PointX(rhs._PointX), _PointY(rhs._PointY), _PointZ(rhs._PointZ),
    _ValueX(rhs._ValueX), _ValueY(rhs._ValueY), _ValueZ(rhs._ValueZ),
    _Attr(rhs._Attr)
  {
    _Data = CAllocate<Vector>(_Attr->X(), _Attr->Y(), _Attr->Z());
    _Norm = CAllocate<double>(_Attr->X(), _Attr->Y(), _Attr->Z());
  }

  /// Destructor
  ~ApproximateSplineCoefficients3D()
  {
    Deallocate(_Data);
    Deallocate(_Norm);
  }

  /// Joint results of parallel reduction
  void join(const ApproximateSplineCoefficients3D &rhs)
  {
    for (int k = 0; k < _Attr->Z(); ++k)
    for (int j = 0; j < _Attr->Y(); ++j)
    for (int i = 0; i < _Attr->X(); ++i) {
      _Data[k][j][i] += rhs._Data[k][j][i];
      _Norm[k][j][i] += rhs._Norm[k][j][i];
    }
  }

  /// Add coefficients approximating subset of values
  void operator ()(const blocked_range<int> &re)
  {
    int    i, j, k, ci, cj, ck, A, B, C;
    double x, y, z, w[3], basis, sum;

    const double *wx = _PointX + re.begin();
    const double *wy = _PointY + re.begin();
    const double *wz = _PointZ + re.begin();

    const double *dx = _ValueX + re.begin();
    const double *dy = _ValueY + re.begin();
    const double *dz = _ValueZ + re.begin();

    for (int n = re.begin(); n != re.end(); ++n, ++wx, ++wy, ++wz, ++dx, ++dy, ++dz) {
      if (AreEqual(*dx, 0.) && AreEqual(*dy, 0.) && AreEqual(*dz, 0.)) continue;

      x = *wx, y = *wy, z = *wz;
      _Attr->WorldToLattice(x, y, z);

      i = ifloor(x);
      j = ifloor(y);
      k = ifloor(z);

      A = Kernel::VariableToIndex(x - i);
      B = Kernel::VariableToIndex(y - j);
      C = Kernel::VariableToIndex(z - k);
      --i, --j, --k;

      sum = 0.;
      for (int c = 0; c <= 3; ++c) {
        w[2] = Kernel::LookupTable[C][c];
        for (int b = 0; b <= 3; ++b) {
          w[1] = Kernel::LookupTable[B][b] * w[2];
          for (int a = 0; a <= 3; ++a) {
            w[0] = Kernel::LookupTable[A][a] * w[1];
            sum += w[0] * w[0];
          }
        }
      }

      for (int c = 0; c <= 3; ++c) {
        ck = k + c;
        if (ck < 0 || ck >= _Attr->Z()) continue;
        w[2] = Kernel::LookupTable[C][c];
        for (int b = 0; b <= 3; ++b) {
          cj = j + b;
          if (cj < 0 || cj >= _Attr->Y()) continue;
          w[1] = Kernel::LookupTable[B][b] * w[2];
          for (int a = 0; a <= 3; ++a) {
            ci = i + a;
            if (ci < 0 || ci >= _Attr->X()) continue;
            w[0]  = Kernel::LookupTable[A][a] * w[1];
            basis = w[0] * w[0];
            _Norm[ck][cj][ci] += basis;
            basis *= w[0] / sum;
            _Data[ck][cj][ci]._x += basis * (*dx);
            _Data[ck][cj][ci]._y += basis * (*dy);
            _Data[ck][cj][ci]._z += basis * (*dz);
          }
        }
      }
    }
  }

  void Approximate()
  {
    parallel_reduce(blocked_range<int>(0, _NumberOfPoints, _NumberOfPoints / 4), *this);
    for (int k = 0; k < _Attr->Z(); ++k)
    for (int j = 0; j < _Attr->Y(); ++j)
    for (int i = 0; i < _Attr->X(); ++i) {
      if (!AreEqual(_Norm[k][j][i], 0.)) {
        _Data[k][j][i] /= _Norm[k][j][i];
      }
    }
  }

  const Vector &Get(int cp) const
  {
    return *(_Data[0][0] + cp);
  }

  const Vector &Get(int i, int j, int k) const
  {
    return *(_Data[k][j][i]);
  }
};

} // namespace

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformation3D
::AddApproximateSplineCoefficients(const double *wx, const double *wy, const double *wz,
                                   const double *dx, const double *dy, const double *dz,
                                   int no, double *gradient, double weight, bool incl_passive) const
{
  if (AreEqual(weight, 0.)) return;
  int xdof, ydof, zdof;
  if (_z == 1) {
    ApproximateSplineCoefficients2D<Kernel> spline(&_attr, wx, wy, wz, dx, dy, dz, no);
    spline.Approximate();
    for (int cp = 0; cp < NumberOfCPs(); ++cp) {
      if (incl_passive || IsActive(cp)) {
        auto &coeff = spline.Get(cp);
        this->IndexToDOFs(cp, xdof, ydof, zdof);
        gradient[xdof] += weight * coeff._x;
        gradient[ydof] += weight * coeff._y;
        gradient[zdof] += weight * coeff._z;
      }
    }
  } else {
    ApproximateSplineCoefficients3D<Kernel> spline(&_attr, wx, wy, wz, dx, dy, dz, no);
    spline.Approximate();
    for (int cp = 0; cp < NumberOfCPs(); ++cp) {
      if (incl_passive || IsActive(cp)) {
        auto &coeff = spline.Get(cp);
        this->IndexToDOFs(cp, xdof, ydof, zdof);
        gradient[xdof] += weight * coeff._x;
        gradient[ydof] += weight * coeff._y;
        gradient[zdof] += weight * coeff._z;
      }
    }
  }
}
#else // _MIRTK_BSPLINEFFD_PARALLEL_APPROXIMATION
// -----------------------------------------------------------------------------
void BSplineFreeFormTransformation3D
::AddApproximateSplineCoefficients(const double *wx, const double *wy, const double *wz,
                                   const double *dx, const double *dy, const double *dz,
                                   int no, double *gradient, double weight, bool incl_passive) const
{
  if (AreEqual(weight, 0.)) return;

  int    i, j, k, ci, cj, ck, A, B, C;
  double x, y, z, w[3], basis, sum;

  // Allocate memory
  Vector ***data = CAllocate<Vector>(_x, _y, _z);
  double ***norm = CAllocate<double>(_x, _y, _z);

  // Initial loop: Calculate change of control points
  for (int idx = 0; idx < no; ++idx) {
    if (AreEqual(dx[idx], 0.) && AreEqual(dy[idx], 0.) && AreEqual(dz[idx], 0.)) continue;

    x = wx[idx], y = wy[idx], z = wz[idx];
    this->WorldToLattice(x, y, z);

    i = ifloor(x);
    j = ifloor(y);

    A = Kernel::VariableToIndex(x - i);
    B = Kernel::VariableToIndex(y - j);
    --i, --j;

    // 2D
    if (_z == 1) {

      sum = .0;
      for (int b = 0; b <= 3; ++b) {
        w[1] = Kernel::LookupTable[B][b];
        for (int a = 0; a <= 3; ++a) {
          w[0] = Kernel::LookupTable[A][a] * w[1];
          sum += w[0] * w[0];
        }
      }

      for (int b = 0; b <= 3; ++b) {
        cj = j + b;
        if (cj < 0 || cj >= _y) continue;
        w[1] = Kernel::LookupTable[B][b];
        for (int a = 0; a <= 3; ++a) {
          ci = i + a;
          if (ci < 0 || ci >= _x) continue;
          w[0]  = Kernel::LookupTable[A][a] * w[1];
          basis = w[0] * w[0];
          norm[0][cj][ci] += basis;
          basis *= w[0] / sum;
          data[0][cj][ci]._x += basis * dx[idx];
          data[0][cj][ci]._y += basis * dy[idx];
          data[0][cj][ci]._z += basis * dz[idx];
        }
      }

    // 3D
    } else {

      k = ifloor(z);
      C = Kernel::VariableToIndex(z - k);
      --k;

      sum = .0;
      for (int c = 0; c <= 3; ++c) {
        w[2] = Kernel::LookupTable[C][c];
        for (int b = 0; b <= 3; ++b) {
          w[1] = Kernel::LookupTable[B][b] * w[2];
          for (int a = 0; a <= 3; ++a) {
            w[0] = Kernel::LookupTable[A][a] * w[1];
            sum += w[0] * w[0];
          }
        }
      }

      for (int c = 0; c <= 3; ++c) {
        ck = k + c;
        if (ck < 0 || ck >= _z) continue;
        w[2] = Kernel::LookupTable[C][c];
        for (int b = 0; b <= 3; ++b) {
          cj = j + b;
          if (cj < 0 || cj >= _y) continue;
          w[1] = Kernel::LookupTable[B][b] * w[2];
          for (int a = 0; a <= 3; ++a) {
            ci = i + a;
            if (ci < 0 || ci >= _x) continue;
            w[0]  = Kernel::LookupTable[A][a] * w[1];
            basis = w[0] * w[0];
            norm[ck][cj][ci] += basis;
            basis *= w[0] / sum;
            data[ck][cj][ci]._x += basis * dx[idx];
            data[ck][cj][ci]._y += basis * dy[idx];
            data[ck][cj][ci]._z += basis * dz[idx];
          }
        }
      }

    }
  }

  // Final loop: Calculate new control points
  int xdof, ydof, zdof;

  Vector *in  = data[0][0];
  double *div = norm[0][0];

  for (int cp = 0; cp < NumberOfCPs(); ++cp, ++in, ++div) {
    if ((incl_passive || IsActive(cp)) && !AreEqual(*div, 0.)) {
      this->IndexToDOFs(cp, xdof, ydof, zdof);
      gradient[xdof] += weight * in->_x / (*div);
      gradient[ydof] += weight * in->_y / (*div);
      gradient[zdof] += weight * in->_z / (*div);
    }
  }

  // Deallocate memory
  Deallocate(data);
  Deallocate(norm);
}
#endif // _MIRTK_BSPLINEFFD_PARALLEL_APPROXIMATION

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformation3D
::ApproximateDOFs(const double *wx, const double *wy, const double *wz, const double *,
                  const double *dx, const double *dy, const double *dz, int no)
{
  memset(_Param, 0, _NumberOfDOFs * sizeof(double));
  AddApproximateSplineCoefficients(wx, wy, wz, dx, dy, dz, no, _Param);
  this->Changed(true);
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformation3D
::ApproximateDOFsGradient(const double *wx, const double *wy, const double *wz, const double *,
                          const double *dx, const double *dy, const double *dz,
                          int no, double *gradient, double weight) const
{
  int    i, j, k, ci, cj, ck, A, B, C;
  double x, y, z, w[3];

  // Allocate memory
  Vector ***data = CAllocate<Vector>(_x, _y, _z);

  // Initial loop: Calculate change of control points
  for (int idx = 0; idx < no; ++idx) {
    if (dx[idx] == .0 && dy[idx] == .0 && dz[idx] == .0) continue;

    x = wx[idx], y = wy[idx], z = wz[idx];
    this->WorldToLattice(x, y, z);

    i = ifloor(x);
    j = ifloor(y);
    A = Kernel::VariableToIndex(x - i);
    B = Kernel::VariableToIndex(y - j);
    --i, --j;

    // 2D
    if (_z == 1) {

      for (int b = 0; b <= 3; ++b) {
        cj = j + b;
        if (cj < 0 || cj >= _y) continue;
        w[1] = Kernel::LookupTable[B][b];
        for (int a = 0; a <= 3; ++a) {
          ci = i + a;
          if (ci < 0 || ci >= _x) continue;
          w[0] = Kernel::LookupTable[A][a] * w[1];
          data[0][cj][ci]._x += w[0] * dx[idx];
          data[0][cj][ci]._y += w[0] * dy[idx];
          data[0][cj][ci]._z += w[0] * dz[idx];
        }
      }

    // 3D
    } else {

      k = ifloor(z);
      C = Kernel::VariableToIndex(z - k);
      --k;

      for (int c = 0; c <= 3; ++c) {
        ck = k + c;
        if (ck < 0 || ck >= _z) continue;
        w[2] = Kernel::LookupTable[C][c];
        for (int b = 0; b <= 3; ++b) {
          cj = j + b;
          if (cj < 0 || cj >= _y) continue;
          w[1] = Kernel::LookupTable[B][b] * w[2];
          for (int a = 0; a <= 3; ++a) {
            ci = i + a;
            if (ci < 0 || ci >= _x) continue;
            w[0] = Kernel::LookupTable[A][a] * w[1];
            data[ck][cj][ci]._x += w[0] * dx[idx];
            data[ck][cj][ci]._y += w[0] * dy[idx];
            data[ck][cj][ci]._z += w[0] * dz[idx];
          }
        }
      }

    }
  }

  // Final loop
  int           xdof, ydof, zdof;
  const Vector *grad = data[0][0];

  for (int cp = 0; cp < NumberOfCPs(); ++cp, ++grad) {
    if (!IsActive(cp)) continue;
    this->IndexToDOFs(cp, xdof, ydof, zdof);
    gradient[xdof] += weight * grad->_x;
    gradient[ydof] += weight * grad->_y;
    gradient[zdof] += weight * grad->_z;
  }

  // Deallocate memory
  Deallocate(data);
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformation3D
::Interpolate(const double *dx, const double *dy, const double *dz)
{
  if (dz) {
    for (int idx = 0; idx < NumberOfCPs(); ++idx) {
      _CPImage(idx) = Vector(dx[idx], dy[idx], dz[idx]);
    }
  } else {
    for (int idx = 0; idx < NumberOfCPs(); ++idx) {
      _CPImage(idx) = Vector(dx[idx], dy[idx], .0);
    }
  }
  if (_z == 1) ConvertToSplineCoefficients(3, _CPImage, 0, 0);
  else         ConvertToSplineCoefficients(3, _CPImage,    0);
  this->Changed(true);
}

// =============================================================================
// Lattice
// =============================================================================

// -----------------------------------------------------------------------------
int BSplineFreeFormTransformation3D::GetXAfterSubdivision() const
{
  return (_x == 1 ? _x : 2 * _x - 1);
}

// -----------------------------------------------------------------------------
int BSplineFreeFormTransformation3D::GetYAfterSubdivision() const
{
  return (_y == 1 ? _y : 2 * _y - 1);
}

// -----------------------------------------------------------------------------
int BSplineFreeFormTransformation3D::GetZAfterSubdivision() const
{
  return (_z == 1 ? _z : 2 * _z - 1);
}

// -----------------------------------------------------------------------------
int BSplineFreeFormTransformation3D::GetTAfterSubdivision() const
{
  return _t;
}

// -----------------------------------------------------------------------------
double BSplineFreeFormTransformation3D::GetXSpacingAfterSubdivision() const
{
  return (_x == 1 ? _dx : 0.5 * _dx);
}

// -----------------------------------------------------------------------------
double BSplineFreeFormTransformation3D::GetYSpacingAfterSubdivision() const
{
  return (_y == 1 ? _dy : 0.5 * _dy);
}

// -----------------------------------------------------------------------------
double BSplineFreeFormTransformation3D::GetZSpacingAfterSubdivision() const
{
  return (_z == 1 ? _dz : 0.5 * _dz);
}

// -----------------------------------------------------------------------------
double BSplineFreeFormTransformation3D::GetTSpacingAfterSubdivision() const
{
  return _dt;
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformation3D
::Subdivide(bool subdivide_x, bool subdivide_y, bool subdivide_z, bool)
{
  if (!subdivide_x && !subdivide_y && !subdivide_z) return;

  if (_x == 1) subdivide_x = false;
  if (_y == 1) subdivide_y = false;
  if (_z == 1) subdivide_z = false;

  // Weights for subdivision
  const double w1[2][3] = {{0.125, 0.75, 0.125}, {0.0, 0.5, 0.5}};
  const double w2[2][3] = {{0.0,   1.0,  0.0},   {0.0, 0.0, 0.0}};

  const double (*wx)[3] = (subdivide_x ? w1 : w2);
  const double (*wy)[3] = (subdivide_y ? w1 : w2);
  const double (*wz)[3] = (subdivide_z ? w1 : w2);

  // Size of new control point grid
  const int X = (subdivide_x ? (2 * _x - 1) : _x);
  const int Y = (subdivide_y ? (2 * _y - 1) : _y);
  const int Z = (subdivide_z ? (2 * _z - 1) : _z);

  // Allocate memory for new control points
  CPValue ***data = Allocate<CPValue>(X, Y, Z);

  // Limits for inner loops
  const int I2 = (subdivide_x ? 2 : 1); // s.t. i+i2-1 == i if no subdivision
  const int J2 = (subdivide_y ? 2 : 1);
  const int K2 = (subdivide_z ? 2 : 1);

  // Compute control point values on subdivided grid
  int si, sj, sk, I1, J1, K1;
  for (int k = 0; k < _z; ++k) {
    if (subdivide_z) sk = 2 * k, K1 = int(k < _z-1);
    else             sk =     k, K1 = 0;
    for (int j = 0; j < _y; ++j) {
      if (subdivide_y) sj = 2 * j, J1 = int(j < _y-1);
      else             sj =     j, J1 = 0;
      for (int i = 0; i < _x; ++i) {
        if (subdivide_x) si = 2 * i, I1 = int(i < _x-1);
        else             si =     i, I1 = 0;
        for (int k1 = 0; k1 <= K1; ++k1)
        for (int j1 = 0; j1 <= J1; ++j1)
        for (int i1 = 0; i1 <= I1; ++i1) {
          for (int k2 = (subdivide_z ? 0 : K2); k2 <= K2; ++k2)
          for (int j2 = (subdivide_y ? 0 : J2); j2 <= J2; ++j2)
          for (int i2 = (subdivide_x ? 0 : I2); i2 <= I2; ++i2) {
            // If a dimension is not subdivided,
            // 1) si + i1    == i
            // 2) i + i2 - 1 == i
            // 3) w[i1][i2]  == 1.0
            data[sk+k1][sj+j1][si+i1] += wx[i1][i2] * wy[j1][j2] * wz[k1][k2]
                                       * _CPValue->Get(i+i2-1, j+j2-1, k+k2-1);
          }
        }
      }
    }
  }

  // Initialize subdivided free-form transformation
  ImageAttributes attr = this->Attributes();
  attr._x = X;
  attr._y = Y;
  attr._z = Z;
  if (subdivide_x) attr._dx *= 0.5;
  if (subdivide_y) attr._dy *= 0.5;
  if (subdivide_z) attr._dz *= 0.5;
  this->Initialize(attr);

  // Copy subdivided control point data
  memcpy(_CPImage.Data(), data[0][0], X * Y * Z * sizeof(CPValue));

  // Deallocate temporary memory
  Deallocate(data);
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformation3D
::BoundingBox(int cp, double &x1, double &y1, double &z1,
                      double &x2, double &y2, double &z2, double fraction) const
{
  int i, j, k;
  IndexToLattice(cp, i, j, k);

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

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
template <class CPImage>
void Evaluate(const CPImage *coeff, double &dx, double &dy, double &dz, int i, int j)
{
  typedef BSplineFreeFormTransformation3D::Kernel Kernel;

  const double *w = Kernel::LatticeWeights;

  typename CPImage::VoxelType d;
  int                         jb;

  --i, --j;
  for (int b = 0; b < 3; ++b) {
    jb = j + b;
    for (int a = 0; a < 3; ++a) {
      d += (w[a] * w[b]) * coeff->Get(i + a, jb);
    }
  }

  dx = d._x, dy = d._y, dz = d._z;
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformation3D
::Evaluate(double &dx, double &dy, double &dz, int i, int j) const
{
  if (_FFD.IsInside(i, j)) mirtk::Evaluate(&_CPImage, dx, dy, dz, i, j);
  else                     mirtk::Evaluate( _CPValue, dx, dy, dz, i, j);
}

// -----------------------------------------------------------------------------
template <class CPImage>
void Evaluate(const CPImage *coeff, double &dx, double &dy, double &dz, int i, int j, int k)
{
  typedef BSplineFreeFormTransformation3D::Kernel Kernel;

  const double *w = Kernel::LatticeWeights;

  typename CPImage::VoxelType d;
  int                         jb, kc;

  --i, --j, --k;
  for (int c = 0; c < 3; ++c) {
    kc = k + c;
    for (int b = 0; b < 3; ++b) {
      jb = j + b;
      for (int a = 0; a < 3; ++a) {
        d += (w[a] * w[b] * w[c]) * coeff->Get(i + a, jb, kc);
      }
    }
  }

  dx = d._x, dy = d._y, dz = d._z;
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformation3D
::Evaluate(double &dx, double &dy, double &dz, int i, int j, int k) const
{
  if (_FFD.IsInside(i, j, k)) mirtk::Evaluate(&_CPImage, dx, dy, dz, i, j, k);
  else                        mirtk::Evaluate( _CPValue, dx, dy, dz, i, j, k);
}

// -----------------------------------------------------------------------------
template <class CPImage>
void EvaluateJacobian(const CPImage *coeff, Matrix &jac, int i, int j)
{
  typedef BSplineFreeFormTransformation3D::Kernel Kernel;

  const double *w[2] = {
    Kernel::LatticeWeights,
    Kernel::LatticeWeights_I
  };

  typename CPImage::VoxelType dx, dy;
  int                         ia, jb;

  --i, --j;
  for (int b = 0; b < 3; ++b) {
    jb = j + b;
    for (int a = 0; a < 3; ++a) {
      ia = i + a;
      dx += (w[1][a] * w[0][b]) * coeff->Get(ia, jb);
      dy += (w[0][a] * w[1][b]) * coeff->Get(ia, jb);
    }
  }

  jac.Initialize(3, 3);
  jac(0, 0) = dx._x; jac(0, 1) = dy._x;
  jac(1, 0) = dx._y; jac(1, 1) = dy._y;
  jac(2, 0) = dx._z; jac(2, 1) = dy._z;
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformation3D
::EvaluateJacobian(Matrix &jac, int i, int j) const
{
  if (_FFD.IsInside(i, j)) mirtk::EvaluateJacobian(&_CPImage, jac, i, j);
  else                     mirtk::EvaluateJacobian( _CPValue, jac, i, j);
}

// -----------------------------------------------------------------------------
template <class CPImage>
void EvaluateJacobian(const CPImage *coeff, Matrix &jac, int i, int j, int k)
{
  typedef BSplineFreeFormTransformation3D::Kernel Kernel;

  const double *w[2] = {
    Kernel::LatticeWeights,
    Kernel::LatticeWeights_I
  };

  typename CPImage::VoxelType dx, dy, dz;
  int                         ia, jb, kc;

  --i, --j, --k;
  for (int c = 0; c < 3; ++c) {
    kc = k + c;
    for (int b = 0; b < 3; ++b) {
      jb = j + b;
      for (int a = 0; a < 3; ++a) {
        ia = i + a;
        dx += (w[1][a] * w[0][b] * w[0][c]) * coeff->Get(ia, jb, kc);
        dy += (w[0][a] * w[1][b] * w[0][c]) * coeff->Get(ia, jb, kc);
        dz += (w[0][a] * w[0][b] * w[1][c]) * coeff->Get(ia, jb, kc);
      }
    }
  }

  jac.Initialize(3, 3);
  jac(0, 0) = dx._x; jac(0, 1) = dy._x; jac(0, 2) = dz._x;
  jac(1, 0) = dx._y; jac(1, 1) = dy._y; jac(1, 2) = dz._y;
  jac(2, 0) = dx._z; jac(2, 1) = dy._z; jac(2, 2) = dz._z;
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformation3D
::EvaluateJacobian(Matrix &jac, int i, int j, int k) const
{
  if (_FFD.IsInside(i, j, k)) mirtk::EvaluateJacobian(&_CPImage, jac, i, j, k);
  else                        mirtk::EvaluateJacobian( _CPValue, jac, i, j, k);
}

// -----------------------------------------------------------------------------
template <class CPImage>
void EvaluateJacobian(const CPImage *coeff, Matrix &jac, double x, double y)
{
  typedef BSplineFreeFormTransformation3D::Kernel Kernel;

  int i = ifloor(x);
  int j = ifloor(y);

  const int A = Kernel::VariableToIndex(x - i);
  const int B = Kernel::VariableToIndex(y - j);

  typename CPImage::VoxelType dx, dy;
  double                      wx[2], wy[2];
  int                         ia, jb;

  --i, --j;
  for (int b = 0; b < 4; ++b) {
    jb = j + b;
    wy[0] = Kernel::LookupTable  [B][b];
    wy[1] = Kernel::LookupTable_I[B][b];
    for (int a = 0; a < 4; ++a) {
      ia = i + a;
      wx[0] = Kernel::LookupTable  [A][a];
      wx[1] = Kernel::LookupTable_I[A][a];
      dx += (wx[1] * wy[0]) * coeff->Get(ia, jb);
      dy += (wx[0] * wy[1]) * coeff->Get(ia, jb);
    }
  }

  jac.Initialize(3, 3);
  jac(0, 0) = dx._x; jac(0, 1) = dy._x;
  jac(1, 0) = dx._y; jac(1, 1) = dy._y;
  jac(2, 0) = dx._z; jac(2, 1) = dy._z;
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformation3D
::EvaluateJacobian(Matrix &jac, double x, double y) const
{
  if (_FFD.IsInside(x, y)) mirtk::EvaluateJacobian(&_CPImage, jac, x, y);
  else                     mirtk::EvaluateJacobian( _CPValue, jac, x, y);
}

// -----------------------------------------------------------------------------
template <class CPImage>
void EvaluateJacobian(const CPImage *coeff, Matrix &jac, double x, double y, double z)
{
  typedef BSplineFreeFormTransformation3D::Kernel Kernel;

  int i = ifloor(x);
  int j = ifloor(y);
  int k = ifloor(z);

  const int A = Kernel::VariableToIndex(x - i);
  const int B = Kernel::VariableToIndex(y - j);
  const int C = Kernel::VariableToIndex(z - k);

  typename CPImage::VoxelType dx, dy, dz;
  double                      wx[2], wy[2], wz[2];
  int                         ia, jb, kc;

  --i, --j, --k;
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
        dx += (wx[1] * wy[0] * wz[0]) * coeff->Get(ia, jb, kc);
        dy += (wx[0] * wy[1] * wz[0]) * coeff->Get(ia, jb, kc);
        dz += (wx[0] * wy[0] * wz[1]) * coeff->Get(ia, jb, kc);
      }
    }
  }

  jac.Initialize(3, 3);
  jac(0, 0) = dx._x; jac(0, 1) = dy._x; jac(0, 2) = dz._x;
  jac(1, 0) = dx._y; jac(1, 1) = dy._y; jac(1, 2) = dz._y;
  jac(2, 0) = dx._z; jac(2, 1) = dy._z; jac(2, 2) = dz._z;
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformation3D
::EvaluateJacobian(Matrix &jac, double x, double y, double z) const
{
  if (_FFD.IsInside(x, y, z)) mirtk::EvaluateJacobian(&_CPImage, jac, x, y, z);
  else                        mirtk::EvaluateJacobian( _CPValue, jac, x, y, z);
}

// -----------------------------------------------------------------------------
template <class CPImage>
void EvaluateHessian(const CPImage *coeff, Matrix hessian[3], int i, int j)
{
  typedef BSplineFreeFormTransformation3D::Kernel Kernel;

  const double *w[3] = {
    Kernel::LatticeWeights,
    Kernel::LatticeWeights_I,
    Kernel::LatticeWeights_II
  };

  typename CPImage::VoxelType dxx, dxy, dyy;
  int                         ia, jb;

  --i, --j;
  for (int b = 0; b < 3; ++b) {
    jb = j + b;
    for (int a = 0; a < 3; ++a) {
      ia = i + a;
      dxx += (w[2][a] * w[0][b]) * coeff->Get(ia, jb);
      dxy += (w[1][a] * w[1][b]) * coeff->Get(ia, jb);
      dyy += (w[0][a] * w[2][b]) * coeff->Get(ia, jb);
    }
  }

  Matrix &hx = hessian[0];
  hx.Initialize(3, 3);
  hx(0, 0) = dxx._x; hx(0, 1) = dxy._x;
  hx(1, 0) = dxy._x; hx(1, 1) = dyy._x;

  Matrix &hy = hessian[1];
  hy.Initialize(3, 3);
  hy(0, 0) = dxx._y; hy(0, 1) = dxy._y;
  hy(1, 0) = dxy._y; hy(1, 1) = dyy._y;

  Matrix &hz = hessian[2];
  hz.Initialize(3, 3);
  hz(0, 0) = dxx._z; hz(0, 1) = dxy._z;
  hz(1, 0) = dxy._z; hz(1, 1) = dyy._z;
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformation3D
::EvaluateHessian(Matrix hessian[3], int i, int j) const
{
  if (_FFD.IsInside(i, j)) mirtk::EvaluateHessian(&_CPImage, hessian, i, j);
  else                     mirtk::EvaluateHessian( _CPValue, hessian, i, j);
}

// -----------------------------------------------------------------------------
template <class CPImage>
void EvaluateHessian(const CPImage *coeff, Matrix hessian[3], int i, int j, int k)
{
  typedef BSplineFreeFormTransformation3D::Kernel Kernel;

  const double *w[3] = {
    Kernel::LatticeWeights,
    Kernel::LatticeWeights_I,
    Kernel::LatticeWeights_II
  };

  typename CPImage::VoxelType dxx, dxy, dxz, dyy, dyz, dzz;
  int                         ia, jb, kc;

  --i, --j, --k;
  for (int c = 0; c < 3; ++c) {
    kc = k + c;
    for (int b = 0; b < 3; ++b) {
      jb = j + b;
      for (int a = 0; a < 3; ++a) {
        ia = i + a;
        dxx += (w[2][a] * w[0][b] * w[0][c]) * coeff->Get(ia, jb, kc);
        dxy += (w[1][a] * w[1][b] * w[0][c]) * coeff->Get(ia, jb, kc);
        dxz += (w[1][a] * w[0][b] * w[1][c]) * coeff->Get(ia, jb, kc);
        dyy += (w[0][a] * w[2][b] * w[0][c]) * coeff->Get(ia, jb, kc);
        dyz += (w[0][a] * w[1][b] * w[1][c]) * coeff->Get(ia, jb, kc);
        dzz += (w[0][a] * w[0][b] * w[2][c]) * coeff->Get(ia, jb, kc);
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
void BSplineFreeFormTransformation3D
::EvaluateHessian(Matrix hessian[3], int i, int j, int k) const
{
  if (_FFD.IsInside(i, j, k)) mirtk::EvaluateHessian(&_CPImage, hessian, i, j, k);
  else                        mirtk::EvaluateHessian( _CPValue, hessian, i, j, k);
}

// -----------------------------------------------------------------------------
template <class CPImage>
void EvaluateHessian(const CPImage *coeff, Matrix hessian[3], double x, double y)
{
  typedef BSplineFreeFormTransformation3D::Kernel Kernel;

  int i = ifloor(x);
  int j = ifloor(y);

  const int A = Kernel::VariableToIndex(x - i);
  const int B = Kernel::VariableToIndex(y - j);

  typename CPImage::VoxelType dxx, dxy, dyy;
  double                      wx[3], wy[3];
  int                         ia, jb;

  --i, --j;
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
      dxx += (wx[2] * wy[0]) * coeff->Get(ia, jb);
      dxy += (wx[1] * wy[1]) * coeff->Get(ia, jb);
      dyy += (wx[0] * wy[2]) * coeff->Get(ia, jb);
    }
  }

  Matrix &hx = hessian[0];
  hx.Initialize(3, 3);
  hx(0, 0) = dxx._x; hx(0, 1) = dxy._x;
  hx(1, 0) = dxy._x; hx(1, 1) = dyy._x;

  Matrix &hy = hessian[1];
  hy.Initialize(3, 3);
  hy(0, 0) = dxx._y; hy(0, 1) = dxy._y;
  hy(1, 0) = dxy._y; hy(1, 1) = dyy._y;

  Matrix &hz = hessian[2];
  hz.Initialize(3, 3);
  hz(0, 0) = dxx._z; hz(0, 1) = dxy._z;
  hz(1, 0) = dxy._z; hz(1, 1) = dyy._z;
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformation3D
::EvaluateHessian(Matrix hessian[3], double x, double y) const
{
  if (_FFD.IsInside(x, y)) mirtk::EvaluateHessian(&_CPImage, hessian, x, y);
  else                     mirtk::EvaluateHessian( _CPValue, hessian, x, y);
}

// -----------------------------------------------------------------------------
template <class CPImage>
void EvaluateHessian(const CPImage *coeff, Matrix hessian[3], double x, double y, double z)
{
  typedef BSplineFreeFormTransformation3D::Kernel Kernel;

  int i = ifloor(x);
  int j = ifloor(y);
  int k = ifloor(z);

  const int A = Kernel::VariableToIndex(x - i);
  const int B = Kernel::VariableToIndex(y - j);
  const int C = Kernel::VariableToIndex(z - k);

  typename CPImage::VoxelType dxx, dxy, dxz, dyy, dyz, dzz;
  double                      wx[3], wy[3], wz[3];
  int                         ia, jb, kc;

  --i, --j, --k;
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
        dxx += (wx[2] * wy[0] * wz[0]) * coeff->Get(ia, jb, kc);
        dxy += (wx[1] * wy[1] * wz[0]) * coeff->Get(ia, jb, kc);
        dxz += (wx[1] * wy[0] * wz[1]) * coeff->Get(ia, jb, kc);
        dyy += (wx[0] * wy[2] * wz[0]) * coeff->Get(ia, jb, kc);
        dyz += (wx[0] * wy[1] * wz[1]) * coeff->Get(ia, jb, kc);
        dzz += (wx[0] * wy[0] * wz[2]) * coeff->Get(ia, jb, kc);
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
void BSplineFreeFormTransformation3D
::EvaluateHessian(Matrix hessian[3], double x, double y, double z) const
{
  if (_FFD.IsInside(x, y, z)) mirtk::EvaluateHessian(&_CPImage, hessian, x, y, z);
  else                        mirtk::EvaluateHessian( _CPValue, hessian, x, y, z);
}

// -----------------------------------------------------------------------------
template <class CPImage>
void EvaluateLaplacian(const CPImage *coeff, double laplacian[3], int i, int j, int k)
{
  typedef BSplineFreeFormTransformation3D::Kernel Kernel;

  const double *w[3] = {
    Kernel::LatticeWeights,
    Kernel::LatticeWeights_I,
    Kernel::LatticeWeights_II
  };

  typename CPImage::VoxelType v = .0;
  int                         ia, jb, kc;

  --i, --j, --k;
  for (int c = 0; c < 3; ++c) {
    kc = k + c;
    for (int b = 0; b < 3; ++b) {
      jb = j + b;
      for (int a = 0; a < 3; ++a) {
        ia = i + a;
        v += (w[2][a] * w[0][b] * w[0][c]) * coeff->Get(ia, jb, kc);
        v += (w[0][a] * w[2][b] * w[0][c]) * coeff->Get(ia, jb, kc);
        v += (w[0][a] * w[0][b] * w[2][c]) * coeff->Get(ia, jb, kc);
      }
    }
  }

  laplacian[0] = v._x;
  laplacian[1] = v._y;
  laplacian[2] = v._z;
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformation3D
::EvaluateLaplacian(double laplacian[3], int i, int j, int k) const
{
  if (_FFD.IsInside(i, j, k)) mirtk::EvaluateLaplacian(&_CPImage, laplacian, i, j, k);
  else                        mirtk::EvaluateLaplacian( _CPValue, laplacian, i, j, k);
}

// -----------------------------------------------------------------------------
template <class CPImage>
void EvaluateLaplacian(const CPImage *coeff, double &x, double &y, double &z)
{
  typedef BSplineFreeFormTransformation3D::Kernel Kernel;

  int i = ifloor(x);
  int j = ifloor(y);
  int k = ifloor(z);

  const int A = Kernel::VariableToIndex(x - i);
  const int B = Kernel::VariableToIndex(y - j);
  const int C = Kernel::VariableToIndex(z - k);

  typename CPImage::VoxelType v;
  double                      wx[2], wy[2], wz[2];
  int                         ia, jb, kc;

  --i, --j, --k;
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
        v += (wx[1] * wy[0] * wz[0]) * coeff->Get(ia, jb, kc);
        v += (wx[0] * wy[1] * wz[0]) * coeff->Get(ia, jb, kc);
        v += (wx[0] * wy[0] * wz[1]) * coeff->Get(ia, jb, kc);
      }
    }
  }

  x = v._x, y = v._y, z = v._z;
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformation3D
::EvaluateLaplacian(double laplacian[3], double x, double y, double z) const
{
  if (_FFD.IsInside(x, y, z)) mirtk::EvaluateLaplacian(&_CPImage, x, y, z);
  else                        mirtk::EvaluateLaplacian( _CPValue, x, y, z);
  laplacian[0] = x, laplacian[1] = y, laplacian[2] = z;
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformation3D
::EvaluateLaplacian(double &x, double &y, double &z) const
{
  if (_FFD.IsInside(x, y, z)) mirtk::EvaluateLaplacian(&_CPImage, x, y, z);
  else                        mirtk::EvaluateLaplacian( _CPValue, x, y, z);
}

// =============================================================================
// Point transformation
// =============================================================================

// -----------------------------------------------------------------------------
bool BSplineFreeFormTransformation3D::CanModifyDisplacement(int) const
{
  return true;
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformation3D
::DisplacementAfterDOFChange(int dof, double dv, GenericImage<double> &disp,
           /* unused --> */ double, double, const WorldCoordsImage *) const
{
  // Corresponding control point index and dimension of displacement field
  const int cp  = this->DOFToIndex    (dof);
  const int dim = this->DOFToDimension(dof);

  // Bounding box of control point in image coordinates
  double                x1, y1, z1, x2, y2, z2;
  this->BoundingBox(cp, x1, y1, z1, x2, y2, z2);
  disp.WorldToImage(x1, y1, z1);
  disp.WorldToImage(x2, y2, z2);

  if (x1 > x2) swap(x1, x2);
  if (y1 > y2) swap(y1, y2);
  if (z1 > z2) swap(z1, z2);

  // Calculate increment for voxel offset to kernel function weight lookup table
  const double dx = (Kernel::LookupTableSize - 1) / (x2 - x1);
  const double dy = (Kernel::LookupTableSize - 1) / (y2 - y1);
  const double dz = (Kernel::LookupTableSize - 1) / (z2 - z1);

  // Bounding box of control point in voxel indices
  int                          i1, j1, k1, i2, j2, k2;
  this->BoundingBox(&disp, cp, i1, j1, k1, i2, j2, k2, 1.0 / _SpeedupFactor);

  // Get pointer to displacement of lower-left voxel in bounding box and
  // increments for update of pointer while looping over the voxels
  const int s1 = 1;
  const int s2 = (disp.X() - (i2 - i1 + 1));
  const int s3 = (disp.Y() - (j2 - j1 + 1)) * disp.X();
  double    *d = disp.Data(i1, j1, k1, dim);

  // Loop over voxels in bounding box of control point and add additional
  // displacement along dim induced by change of control point parameter
  double di, dj, dk; // displacement change = (B_i * (B_j * (B_k * dv)))
  for (int k = k1; k <= k2; ++k, d += s3) {
    dk = Kernel::WeightLookupTable[iround((k - z1) * dz)] * dv;
    for (int j = j1; j <= j2; ++j, d += s2) {
      dj = Kernel::WeightLookupTable[iround((j - y1) * dy)] * dk;
      for (int i = i1; i <= i2; ++i, d += s1) {
        di = Kernel::WeightLookupTable[iround((i - x1) * dx)] * dj;
        (*d) += di;
      }
    }
  }
}

// =============================================================================
// Derivatives
// =============================================================================

namespace {


template <class TReal>
class SampleGradientImage
{
public:

  typedef GenericImage<TReal>                                      GradientImage;
  typedef GenericLinearInterpolateImageFunction3D<GradientImage>   GradientFunction;

  const BSplineFreeFormTransformation3D *_FFD;
  const GradientFunction                *_Gradient;
  double                                *_Output;

  void operator ()(const blocked_range3d<int> &re) const
  {
    double g[3], x, y, z;
    int cp, xdof, ydof, zdof;
    BSplineFreeFormTransformation3D::DOFStatus status_x, status_y, status_z;
    for (int ck = re.pages().begin(); ck != re.pages().end(); ++ck)
    for (int cj = re.rows ().begin(); cj != re.rows ().end(); ++cj)
    for (int ci = re.cols ().begin(); ci != re.cols ().end(); ++ci) {
      _FFD->GetStatus(ci, cj, ck, status_x, status_y, status_z);
      if (status_x == Active || status_y == Active || status_z == Active) {
        cp = _FFD->LatticeToIndex(ci, cj, ck);
        x = ci, y = cj, z = ck;
        _FFD->LatticeToWorld(x, y, z);
        _Gradient->WorldToImage(x, y, z);
        _Gradient->Evaluate(g, x, y, z);
        _FFD->IndexToDOFs(cp, xdof, ydof, zdof);
        if (status_x == Active) _Output[xdof] = g[0];
        if (status_y == Active) _Output[ydof] = g[1];
        if (status_z == Active) _Output[zdof] = g[2];
      }
    }
  }
};


}


// -----------------------------------------------------------------------------
void BSplineFreeFormTransformation3D
::ParametricGradient(const GenericImage<double> *in, double *out,
                     const WorldCoordsImage *i2w, const WorldCoordsImage *wc,
                     double t0, double w) const
{
  switch (_ParametricGradientCalculation) {
    case PG_Default:
    case PG_Analytic: {
      FreeFormTransformation3D::ParametricGradient(in, out, i2w, wc, t0, w);
    } break;

    case PG_Convolution: {
      if (wc != nullptr) {
        Throw(ERR_NotImplemented, __FUNCTION__, "Convolution-based B-spline FFD gradient calculation not implemented for wc != nullptr, i.e., 'Fluid' MFFD!");
      }

      MIRTK_START_TIMING();

      typedef SampleGradientImage<double>        GradientSampler;
      typedef GradientSampler::GradientFunction  GradientFunction;

      // Note: in may have more components that are unused here, see SVFFD subclass
      GenericImage<double> tmp;
      CubicBSplineConvolution<double> conv(_dx, _dy, _dz);
      conv.Components(3);
      conv.Input(in);
      conv.Output(&tmp);
      conv.Run();

      GradientFunction grad;
      grad.Input(conv.Output());
      grad.Initialize();

      GradientSampler sample;
      sample._FFD      = this;
      sample._Gradient = &grad;
      sample._Output   = out;
      parallel_for(blocked_range3d<int>(0, _z, 0, _y, 0, _x), sample);

      MIRTK_DEBUG_TIMING(2, "parametric gradient computation (3D B-spline convolution)");
    } break;

    case PG_Approximation: {
      MIRTK_START_TIMING();

      UniquePtr<WorldCoordsImage> _wc;
      if (wc == nullptr) {
        if (i2w == nullptr) {
          _wc.reset(new WorldCoordsImage(in->Attributes(), 3));
          in->ImageToWorld(*_wc.get());
          wc = _wc.get();
        } else {
          wc = i2w;
        }
      }
      AddApproximateSplineCoefficients(wc->Data(0, 0, 0, 0), wc->Data(0, 0, 0, 1), wc->Data(0, 0, 0, 2),
                                       in->Data(0, 0, 0, 0), in->Data(0, 0, 0, 1), in->Data(0, 0, 0, 2),
                                       in->NumberOfSpatialVoxels(), out, w);

      MIRTK_DEBUG_TIMING(2, "parametric gradient computation (3D B-spline DMFFD)");
    } break;
  }
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformation3D
::EvaluateJacobianDetDerivative(double dJ[3], const Matrix &adj, double a, double b, double c) const
{
  // Values of the B-spline basis functions and its 1st derivatives
  // Note: Calling Kernel::B faster/not slower than Kernel::VariableToIndex + Kernel:LookupTable.
  double wa = Kernel::B(a), da = Kernel::B_I(a);
  double wb = Kernel::B(b), db = Kernel::B_I(b);
  double wc = Kernel::B(c), dc = Kernel::B_I(c);

  // Derivatives of 3D cubic B-spline kernel w.r.t. lattice distance from control point
  double du = da * wb * wc;
  double dv = wa * db * wc;
  double dw = wa * wb * dc;

  // Apply chain rule for da/dx, db/dx, dc/dx and sum terms (same for y and z)
  // to get the derivative of the 3D cubic B-spline kernel centered at this
  // control point w.r.t. each spatial world coordinate
  JacobianToWorld(du, dv, dw);

  // Apply Jacobi's formula to get derivatives of Jacobian determinant
  dJ[0] = adj(0, 0) * du + adj(1, 0) * dv + adj(2, 0) * dw;
  dJ[1] = adj(0, 1) * du + adj(1, 1) * dv + adj(2, 1) * dw;
  dJ[2] = adj(0, 2) * du + adj(1, 2) * dv + adj(2, 2) * dw;
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformation3D
::EvaluateJacobianDetDerivative(double dJ[3], const Matrix &adj, int a, int b, int c) const
{
  // Values of the B-spline basis functions and its 1st derivatives at lattice points
  // Note: The order of the cubic B-spline pieces is *not* B0, B1, B2, B3! Hence, the minus signs.
  double wa, wb, wc, da, db, dc;
  if (-1 <= a && a <= 1) {
    wa = Kernel::LatticeWeights[a + 1];
    da = - Kernel::LatticeWeights_I[a + 1];
  }
  if (-1 <= b && b <= 1) {
    wb = Kernel::LatticeWeights[b + 1];
    db = - Kernel::LatticeWeights_I[b + 1];
  }
  if (-1 <= c && c <= 1) {
    wc = Kernel::LatticeWeights[c + 1];
    dc = - Kernel::LatticeWeights_I[c + 1];
  }

  // Product terms of 3D cubic B-spline kernel required for total derivative
  double du = da * wb * wc;
  double dv = wa * db * wc;
  double dw = wa * wb * dc;

  // Apply chain rule for da/dx, db/dx, dc/dx and sum terms (same for y and z)
  // to get the derivative of the 3D cubic B-spline kernel centered at this
  // control point w.r.t. each spatial coordinate
  JacobianToWorld(du, dv, dw);

  // Apply Jacobi's formula to get derivatives of Jacobian determinant
  dJ[0] = adj(0, 0) * du + adj(1, 0) * dv + adj(2, 0) * dw;
  dJ[1] = adj(0, 1) * du + adj(1, 1) * dv + adj(2, 1) * dw;
  dJ[2] = adj(0, 2) * du + adj(1, 2) * dv + adj(2, 2) * dw;
}

// =============================================================================
// Properties
// =============================================================================

// -----------------------------------------------------------------------------
int BSplineFreeFormTransformation3D::KernelSize() const
{
  return 4;
}

// -----------------------------------------------------------------------------
double BSplineFreeFormTransformation3D
::BendingEnergy(double x, double y, double z, double, double, bool wrt_world) const
{
  // Convert to lattice coordinates
  this->WorldToLattice(x, y, z);
  // Calculate 2nd order derivatives
  Matrix hessian[3];
  if (_z == 1) EvaluateHessian(hessian, x, y);
  else         EvaluateHessian(hessian, x, y, z);
  // Convert derivatives to world coordinates
  if (wrt_world) HessianToWorld(hessian);
  // Calculate bending energy
  return Bending3D(hessian);
}

// -----------------------------------------------------------------------------
double BSplineFreeFormTransformation3D::BendingEnergy(bool incl_passive, bool wrt_world) const
{
  Matrix hessian[3];
  double bending = .0;
  int    nactive = 0;

  for (int k = 0; k < _z; ++k)
  for (int j = 0; j < _y; ++j)
  for (int i = 0; i < _x; ++i) {
    if (incl_passive || IsActive(i, j, k)) {
      if (_z == 1) EvaluateHessian(hessian, i, j);
      else         EvaluateHessian(hessian, i, j, k);
      if (wrt_world) HessianToWorld(hessian);
      bending += Bending3D(hessian);
      ++nactive;
    }
  }

  if (nactive) bending /= nactive;
  return bending;
}

// -----------------------------------------------------------------------------
double BSplineFreeFormTransformation3D
::BendingEnergy(const ImageAttributes &attr, double, bool wrt_world) const
{
  const int nvox = attr.NumberOfSpatialPoints();
  if (nvox == 0) return .0;

  double bending = .0;
  double x, y, z;
  Matrix hessian[3];

  for (int k = 0; k < attr._z; ++k)
  for (int j = 0; j < attr._y; ++j)
  for (int i = 0; i < attr._x; ++i) {
    x = i, y = j, z = k;
    attr .LatticeToWorld(x, y, z);
    this->WorldToLattice(x, y, z);
    if (_z == 1) EvaluateHessian(hessian, x, y);
    else         EvaluateHessian(hessian, x, y, z);
    if (wrt_world) HessianToWorld(hessian);
    bending += Bending3D(hessian);
  }

  return bending / nvox;
}


namespace BSplineFreeFormTransformation3DUtils {

// -----------------------------------------------------------------------------
/// Voxel function which evaluates the 2nd order derivatives of one component
/// of a B-spline FFD at control point lattice coordinates
struct Evaluate2ndOrderBSplineFFDDerivatives : public VoxelFunction
{
  typedef BSplineFreeFormTransformation3D::CPExtrapolator Extrapolator;
  typedef BSplineFreeFormTransformation3D::Vector         Vector;
  typedef BSplineFreeFormTransformation3D::Kernel         Kernel;

  const BSplineFreeFormTransformation3D *_Transformation; ///< B-spline free-form deformation
  const Extrapolator                    *_CPValue;        ///< Coefficients of B-spline FFD
  bool                                   _WrtWorld;       ///< Whether to compute derivatives
                                                          ///< w.r.t world or lattice coordinates

  /// Derivatives of 2nd order derivatives w.r.t control point parameters
  static double LookupTable_2D[9][3];

  /// Derivatives of 2nd order derivatives w.r.t control point parameters
  static double LookupTable_3D[27][6];

  /// Initialize static lookup tables of 2nd order derivatives w.r.t control
  /// point parameters evaluated for each voxel in the 2D kernel support region
  static void InitializeLookupTable2D();

  /// Initialize static lookup tables of 2nd order derivatives w.r.t control
  /// point parameters evaluated for each voxel in the 3D kernel support region
  static void InitializeLookupTable3D();

  /// Constructor
  Evaluate2ndOrderBSplineFFDDerivatives(const BSplineFreeFormTransformation3D *ffd,
                                        bool wrt_world = false)
  :
    _Transformation(ffd), _CPValue(ffd->Extrapolator()), _WrtWorld(wrt_world)
  {}

  /// Evaluate 2nd order derivatives of 2D FFD at given lattice coordinate
  void operator()(int i, int j, int k, int, Vector *dxx, Vector *dxy, Vector *dyy)
  {
    // Note: Derivatives are evaluated on a lattice that has an
    //       additional boundary margin of one voxel. Therefore,
    //       CP indices I and J are shifted by an offset of -1.
    int n = 0;
    for (int J = j-2; J <= j; ++J)
    for (int I = i-2; I <= i; ++I, ++n) {
      *dxx += _CPValue->Get(I, J) * LookupTable_2D[n][0];
      *dxy += _CPValue->Get(I, J) * LookupTable_2D[n][1];
      *dyy += _CPValue->Get(I, J) * LookupTable_2D[n][2];
    }
    // Apply product and chain rule to convert derivatives to ones w.r.t the world
    if (_WrtWorld) {
      _Transformation->HessianToWorld(dxx->_x, dxy->_x, dyy->_x);
      _Transformation->HessianToWorld(dxx->_y, dxy->_y, dyy->_y);
      _Transformation->HessianToWorld(dxx->_z, dxy->_z, dyy->_z);
    }
  }

  /// Evaluate 2nd order derivatives of 3D FFD at given lattice coordinate
  void operator()(int i, int j, int k, int,
                  Vector *dxx, Vector *dxy, Vector *dxz, Vector *dyy, Vector *dyz, Vector *dzz)
  {
    // Note: Derivatives are evaluated on a lattice that has an
    //       additional boundary margin of one voxel. Therefore,
    //       CP indices I, J and K are shifted by an offset of -1.
    int n = 0;
    for (int K = k-2; K <= k; ++K)
    for (int J = j-2; J <= j; ++J)
    for (int I = i-2; I <= i; ++I, ++n) {
      *dxx += _CPValue->Get(I, J, K) * LookupTable_3D[n][0];
      *dxy += _CPValue->Get(I, J, K) * LookupTable_3D[n][1];
      *dxz += _CPValue->Get(I, J, K) * LookupTable_3D[n][2];
      *dyy += _CPValue->Get(I, J, K) * LookupTable_3D[n][3];
      *dyz += _CPValue->Get(I, J, K) * LookupTable_3D[n][4];
      *dzz += _CPValue->Get(I, J, K) * LookupTable_3D[n][5];
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
double Evaluate2ndOrderBSplineFFDDerivatives::LookupTable_2D[9][3] = {{.0}};
void Evaluate2ndOrderBSplineFFDDerivatives::InitializeLookupTable2D()
{
  static bool initialized = false;
  if (initialized) return;

  const double *w[3] = {
    Kernel::LatticeWeights,
    Kernel::LatticeWeights_I,
    Kernel::LatticeWeights_II
  };

  int n = 0;
  for (int b = 0; b < 3; ++b)
  for (int a = 0; a < 3; ++a, ++n) {
    LookupTable_2D[n][0] = w[2][a] * w[0][b];
    LookupTable_2D[n][1] = w[1][a] * w[1][b];
    LookupTable_2D[n][2] = w[0][a] * w[2][b];
  }

  initialized = true;
}

// -----------------------------------------------------------------------------
double Evaluate2ndOrderBSplineFFDDerivatives::LookupTable_3D[27][6] = {{.0}};
void Evaluate2ndOrderBSplineFFDDerivatives::InitializeLookupTable3D()
{
  static bool initialized = false;
  if (initialized) return;

  const double *w[3] = {
    Kernel::LatticeWeights,
    Kernel::LatticeWeights_I,
    Kernel::LatticeWeights_II
  };

  int n = 0;
  for (int c = 0; c < 3; ++c)
  for (int b = 0; b < 3; ++b)
  for (int a = 0; a < 3; ++a, ++n) {
    LookupTable_3D[n][0] = w[2][a] * w[0][b] * w[0][c];
    LookupTable_3D[n][1] = w[1][a] * w[1][b] * w[0][c];
    LookupTable_3D[n][2] = w[1][a] * w[0][b] * w[1][c];
    LookupTable_3D[n][3] = w[0][a] * w[2][b] * w[0][c];
    LookupTable_3D[n][4] = w[0][a] * w[1][b] * w[1][c];
    LookupTable_3D[n][5] = w[0][a] * w[0][b] * w[2][c];
  }

  initialized = true;
}

} // namespace BSplineFreeFormTransformation3DUtils

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformation3D::BendingEnergyGradient(double *gradient, double weight, bool incl_passive, bool wrt_world) const
{
  using namespace BSplineFreeFormTransformation3DUtils;
  
  Vector sum;
  int    m, n;

  MIRTK_START_TIMING();

  // Pre-multiply weight by derivative of square function (2) and normalization factor
  const int ncps = this->NumberOfActiveCPs();
  if (ncps == 0) return;
  weight *= 2.0 / ncps;

  // ---------------------------------------------------------------------------
  // Bending energy of 2D FFD
  if (_z == 1) {

    // Initialize static lookup table
    Evaluate2ndOrderBSplineFFDDerivatives::InitializeLookupTable2D();

    // Add a layer of boundary voxels to avoid additional boundary conditions
    ImageAttributes attr = this->Attributes();
    attr._x += 2, attr._y += 2;

    // Compute 2nd order derivatives w.r.t control point lattice coordinates,
    // evaluated at each control point of the lattice
    GenericImage<Vector> dxx(attr);
    GenericImage<Vector> dxy(attr);
    GenericImage<Vector> dyy(attr);

    Evaluate2ndOrderBSplineFFDDerivatives eval(this, wrt_world);
    ParallelForEachVoxel(attr, dxx, dxy, dyy, eval);

    // Compute 3rd order derivatives, twice w.r.t. lattice or world coordinate
    // and once w.r.t. transformation parameters of control point
    double w[9][4];
    if (wrt_world) {
      // Loop over support region (3x3) of a control point
      //
      // Note that the following terms are independent of the transformation
      // parameters and therefore can be pre-computed here as they are identical for
      // all control point positions at which the bending gradient is evaluated.
      for (n = 0; n < 9; ++n) {
        m = 0; // index of 2nd order derivative ( d^2 T(x[n]) / dx_i dx_j )
               // which is multiplied by weight w[n][m] according to the chain rule
        for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j, ++m) {
            w[n][m]  = Evaluate2ndOrderBSplineFFDDerivatives::LookupTable_2D[n][0] * _matW2L(0, i) * _matW2L(0, j); // ( d^3 T(x[n]) / dudu dPhi ) * du/dx_i * du/dx_j
            w[n][m] += Evaluate2ndOrderBSplineFFDDerivatives::LookupTable_2D[n][1] * _matW2L(0, i) * _matW2L(1, j); // ( d^3 T(x[n]) / dudv dPhi ) * du/dx_i * dv/dx_j
            w[n][m] += Evaluate2ndOrderBSplineFFDDerivatives::LookupTable_2D[n][1] * _matW2L(1, i) * _matW2L(0, j); // ( d^3 T(x[n]) / dvdu dPhi ) * dv/dx_i * du/dx_j
            w[n][m] += Evaluate2ndOrderBSplineFFDDerivatives::LookupTable_2D[n][2] * _matW2L(1, i) * _matW2L(1, j); // ( d^3 T(x[n]) / dvdv dPhi ) * dv/dx_i * dv/dx_j
        }
        // Add weights of ( d^2 T(x[n]) / dx_i dx_j ) and ( d^2 T(x[n]) / dx_j dx_i )
        // for i != j as these derivatives are identical and re-order remaining weights.
        w[n][1] = w[n][1] + w[n][2];
        w[n][2] = w[n][3];
      }
    } else {
      for (n = 0; n < 9; ++n) {
        w[n][0] =       Evaluate2ndOrderBSplineFFDDerivatives::LookupTable_2D[n][0];
        w[n][1] = 2.0 * Evaluate2ndOrderBSplineFFDDerivatives::LookupTable_2D[n][1];
        w[n][2] =       Evaluate2ndOrderBSplineFFDDerivatives::LookupTable_2D[n][2];
      }
    }

    // Compute derivative of bending energy w.r.t each control point
    int xdof, ydof, zdof;
    for (int cj = 0; cj < _y; ++cj)
    for (int ci = 0; ci < _x; ++ci) {
      if (incl_passive || IsActive(ci, cj)) {
        sum = .0;
        // Loop over support region (3x3) of control point
        //
        // Note: Derivatives were evaluated on a lattice that has an
        //       additional boundary margin of one voxel. Therefore,
        //       indices i and j are shifted by an offset of +1.
        n = 0;
        for (int j = cj; j <= cj+2; ++j)
        for (int i = ci; i <= ci+2; ++i, ++n) {
            sum += dxx(i, j) * w[n][0];
            sum += dxy(i, j) * w[n][1];
            sum += dyy(i, j) * w[n][2];
        }
        this->IndexToDOFs(this->LatticeToIndex(ci, cj), xdof, ydof, zdof);
        gradient[xdof] += weight * sum._x;
        gradient[ydof] += weight * sum._y;
        gradient[zdof] += weight * sum._z;
      }
    }

  // ---------------------------------------------------------------------------
  // Bending energy of 3D FFD
  } else {

    // Initialize static lookup table
    Evaluate2ndOrderBSplineFFDDerivatives::InitializeLookupTable3D();

    // Add a layer of boundary voxels to avoid additional boundary conditions
    ImageAttributes attr = this->Attributes();
    attr._x += 2, attr._y += 2, attr._z += 2;

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
    double w[27][9];
    if (wrt_world) {
      // Loop over support region (3x3x3) of a control point
      //
      // Note that the following terms are independent of the transformation
      // parameters and therefore can be pre-computed here as they are identical for
      // all control point positions at which the bending gradient is evaluated.
      for (n = 0; n < 27; ++n) {
        m = 0; // index of 2nd order derivative ( d^2 T(x[n]) / dx_i dx_j )
               // which is multiplied by weight w[n][m] according to the chain rule
        for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j, ++m) {
          w[n][m]  = Evaluate2ndOrderBSplineFFDDerivatives::LookupTable_3D[n][0] * _matW2L(0, i) * _matW2L(0, j); // ( d^3 T(x[n]) / dudu dPhi ) * du/dx_i * du/dx_j
          w[n][m] += Evaluate2ndOrderBSplineFFDDerivatives::LookupTable_3D[n][1] * _matW2L(0, i) * _matW2L(1, j); // ( d^3 T(x[n]) / dudv dPhi ) * du/dx_i * dv/dx_j
          w[n][m] += Evaluate2ndOrderBSplineFFDDerivatives::LookupTable_3D[n][2] * _matW2L(0, i) * _matW2L(2, j); // ( d^3 T(x[n]) / dudw dPhi ) * du/dx_i * dw/dx_j
          w[n][m] += Evaluate2ndOrderBSplineFFDDerivatives::LookupTable_3D[n][1] * _matW2L(1, i) * _matW2L(0, j); // ( d^3 T(x[n]) / dvdu dPhi ) * dv/dx_i * du/dx_j
          w[n][m] += Evaluate2ndOrderBSplineFFDDerivatives::LookupTable_3D[n][3] * _matW2L(1, i) * _matW2L(1, j); // ( d^3 T(x[n]) / dvdv dPhi ) * dv/dx_i * dv/dx_j
          w[n][m] += Evaluate2ndOrderBSplineFFDDerivatives::LookupTable_3D[n][4] * _matW2L(1, i) * _matW2L(2, j); // ( d^3 T(x[n]) / dvdw dPhi ) * dv/dx_i * dw/dx_j
          w[n][m] += Evaluate2ndOrderBSplineFFDDerivatives::LookupTable_3D[n][2] * _matW2L(2, i) * _matW2L(0, j); // ( d^3 T(x[n]) / dwdu dPhi ) * dw/dx_i * du/dx_j
          w[n][m] += Evaluate2ndOrderBSplineFFDDerivatives::LookupTable_3D[n][4] * _matW2L(2, i) * _matW2L(1, j); // ( d^3 T(x[n]) / dwdv dPhi ) * dw/dx_i * dv/dx_j
          w[n][m] += Evaluate2ndOrderBSplineFFDDerivatives::LookupTable_3D[n][5] * _matW2L(2, i) * _matW2L(2, j); // ( d^3 T(x[n]) / dwdw dPhi ) * dw/dx_i * dw/dx_j
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
      for (n = 0; n < 27; ++n) {
        w[n][0] =       Evaluate2ndOrderBSplineFFDDerivatives::LookupTable_3D[n][0];
        w[n][1] = 2.0 * Evaluate2ndOrderBSplineFFDDerivatives::LookupTable_3D[n][1];
        w[n][2] = 2.0 * Evaluate2ndOrderBSplineFFDDerivatives::LookupTable_3D[n][2];
        w[n][3] =       Evaluate2ndOrderBSplineFFDDerivatives::LookupTable_3D[n][3];
        w[n][4] = 2.0 * Evaluate2ndOrderBSplineFFDDerivatives::LookupTable_3D[n][4];
        w[n][5] =       Evaluate2ndOrderBSplineFFDDerivatives::LookupTable_3D[n][5];
      }
    }

    // Compute derivative of bending energy w.r.t each control point
    int xdof, ydof, zdof;
    for (int ck = 0; ck < _z; ++ck)
    for (int cj = 0; cj < _y; ++cj)
    for (int ci = 0; ci < _x; ++ci) {
      if (incl_passive || IsActive(ci, cj, ck)) {
        sum = .0;
        // Loop over support region (3x3x3) of control point
        //
        // Note: Derivatives were evaluated on a lattice that has an
        //       additional boundary margin of one voxel. Therefore,
        //       indices i, j and k are shifted by an offset of +1.
        n = 0;
        for (int k = ck; k <= ck+2; ++k)
        for (int j = cj; j <= cj+2; ++j)
        for (int i = ci; i <= ci+2; ++i, ++n) {
          sum += dxx(i, j, k) * w[n][0];
          sum += dxy(i, j, k) * w[n][1];
          sum += dxz(i, j, k) * w[n][2];
          sum += dyy(i, j, k) * w[n][3];
          sum += dyz(i, j, k) * w[n][4];
          sum += dzz(i, j, k) * w[n][5];
        }
        this->IndexToDOFs(this->LatticeToIndex(ci, cj, ck), xdof, ydof, zdof);
        gradient[xdof] += weight * sum._x;
        gradient[ydof] += weight * sum._y;
        gradient[zdof] += weight * sum._z;
      }
    }

  }
  MIRTK_DEBUG_TIMING(2, "bending gradient computation");
}

// =============================================================================
// I/O
// =============================================================================

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformation3D::Print(ostream &os, Indent indent) const
{
  os << indent << "3D B-spline FFD:" << endl;
  FreeFormTransformation3D::Print(os, indent + 1);
}

// -----------------------------------------------------------------------------
bool BSplineFreeFormTransformation3D::CanRead(TransformationType format) const
{
  switch (format) {
    case TRANSFORMATION_BSPLINE_FFD_2D_v1:
    case TRANSFORMATION_BSPLINE_FFD_3D_v1:
    case TRANSFORMATION_BSPLINE_FFD_3D_v2:
    case TRANSFORMATION_BSPLINE_FFD_3D_v3:
    case TRANSFORMATION_BSPLINE_FFD_3D_v4:
      return true;
    default:
      return false;
  }
}


} // namespace mirtk
