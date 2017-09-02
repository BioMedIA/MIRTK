/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2017 Imperial College London
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

#include "mirtk/ScalingAndSquaring.h"

#include "mirtk/Math.h"
#include "mirtk/VoxelFunction.h"
#include "mirtk/InterpolateImageFunction.hxx"
#include "mirtk/LinearInterpolateImageFunction.hxx"
#include "mirtk/GaussianBlurring.h"
#include "mirtk/Matrix.h"
#include "mirtk/Parallel.h"
#include "mirtk/Profiling.h"
#include "mirtk/Deallocate.h"

#include <algorithm>


namespace mirtk {


// =============================================================================
// Auxiliary voxel functions for parallel execution
// =============================================================================

namespace {


// -----------------------------------------------------------------------------
/// Multiply vector from left with matrix
class InPlaceMatrixVectorProduct : public VoxelFunction
{
  const Matrix *_Matrix;

public:

  InPlaceMatrixVectorProduct(const Matrix *m) : _Matrix(m) {}

  template <class TImage, class TReal>
  void operator ()(const TImage &im, int, TReal *x)
  {
    const Matrix &m = *_Matrix;
    const int n = im.NumberOfSpatialVoxels();
    auto y = x + n, z = y + n;
    double u = m(0, 0) * (*x) + m(0, 1) * (*y) + m(0, 2) * (*z);
    double v = m(1, 0) * (*x) + m(1, 1) * (*y) + m(1, 2) * (*z);
    double w = m(2, 0) * (*x) + m(2, 1) * (*y) + m(2, 2) * (*z);
    (*x) = static_cast<TReal>(u);
    (*y) = static_cast<TReal>(v);
    (*z) = static_cast<TReal>(w);
  }
};

// -----------------------------------------------------------------------------
template <class TReal>
struct ConvertToDisplacement3D : public VoxelFunction
{
  const int              _x, _y, _z;
  const ImageAttributes &_Domain;

  ConvertToDisplacement3D(const ImageAttributes &domain)
  :
    _x(0), _y(domain.NumberOfSpatialPoints()), _z(2 * _y), _Domain(domain)
  {}

  void operator()(int i, int j, int k, int, const TReal *d, TReal *u)
  {
    double x = i, y = j, z = k;
    _Domain.LatticeToWorld(x, y, z);
    u[_x] = d[_x] - static_cast<TReal>(x);
    u[_y] = d[_y] - static_cast<TReal>(y);
    u[_z] = d[_z] - static_cast<TReal>(z);
  }
};

// -----------------------------------------------------------------------------
template <class TReal>
struct ConvertToDeformation3D : public VoxelFunction
{
  const int              _x, _y, _z;
  const ImageAttributes &_Domain;

  ConvertToDeformation3D(const ImageAttributes &domain)
  :
    _x(0), _y(domain.NumberOfSpatialPoints()), _z(2 * _y), _Domain(domain)
  {}

  void operator()(int i, int j, int k, int, const TReal *u, TReal *d)
  {
    double x = i, y = j, z = k;
    _Domain.LatticeToWorld(x, y, z);
    d[_x] = static_cast<TReal>(x) + u[_x];
    d[_y] = static_cast<TReal>(y) + u[_y];
    d[_z] = static_cast<TReal>(z) + u[_z];
  }
};

// -----------------------------------------------------------------------------
/// Voxel function for the evaluation of the Jacobian w.r.t. x during scaling step
template <class TReal>
struct EvaluateJacobianBase : public VoxelFunction
{
  typedef typename ScalingAndSquaring<TReal>::VelocityField VelocityField;

  static const int     _xx = 0;
  int                  _xy, _xz, _yx, _yy, _yz, _zx, _zy, _zz;
  ImageAttributes      _Domain;
  Matrix               _MatW2L;
  const VelocityField *_VelocityField;
  double               _Scale;

  // ---------------------------------------------------------------------------
  inline void Initialize(const ImageAttributes &domain,
                         const VelocityField   *v,
                         double                 s = 1.0)
  {
    const int n = domain.NumberOfSpatialPoints();
                 _xy = 1 * n; _xz = 2 * n;
    _yx = 3 * n; _yy = 4 * n; _yz = 5 * n;
    _zx = 6 * n; _zy = 7 * n; _zz = 8 * n;
    _Domain        = domain;
    _MatW2L        = domain.GetWorldToLatticeMatrix();
    _VelocityField = v;
    _Scale         = s;
  }

  // ---------------------------------------------------------------------------
  inline void JacobianToWorld(double &du, double &dv, double &dw) const
  {
    double dx = du * _MatW2L(0, 0) + dv * _MatW2L(1, 0) + dw * _MatW2L(2, 0);
    double dy = du * _MatW2L(0, 1) + dv * _MatW2L(1, 1) + dw * _MatW2L(2, 1);
    double dz = du * _MatW2L(0, 2) + dv * _MatW2L(1, 2) + dw * _MatW2L(2, 2);
    du = dx, dv = dy, dw = dz;
  }

  // ---------------------------------------------------------------------------
  inline Matrix Jacobian(int i, int j, int k)
  {
    // Convert output voxel indices to velocity field voxel coordinates
    double x = i, y = j, z = k;
    _Domain.LatticeToWorld(x, y, z);
    _VelocityField->WorldToImage(x, y, z);
    // Evaluate Jacobian of velocity field
    Matrix dvx, dvy, dvz;
    _VelocityField->EvaluateJacobian(dvx, x, y, z, 0);
    _VelocityField->EvaluateJacobian(dvy, x, y, z, 1);
    _VelocityField->EvaluateJacobian(dvz, x, y, z, 2);
    JacobianToWorld(dvx(0, 0), dvx(0, 1), dvx(0, 2));
    JacobianToWorld(dvy(0, 0), dvy(0, 1), dvy(0, 2));
    JacobianToWorld(dvz(0, 0), dvz(0, 1), dvz(0, 2));
    // Compute Jacobian of (scaled) transformation
    Matrix jac(3, 3);
    jac(0, 0) = _Scale * dvx(0, 0) + 1.0;
    jac(0, 1) = _Scale * dvx(0, 1);
    jac(0, 2) = _Scale * dvx(0, 2);
    jac(1, 0) = _Scale * dvy(0, 0);
    jac(1, 1) = _Scale * dvy(0, 1) + 1.0;
    jac(1, 2) = _Scale * dvy(0, 2);
    jac(2, 0) = _Scale * dvz(0, 0);
    jac(2, 1) = _Scale * dvz(0, 1);
    jac(2, 2) = _Scale * dvz(0, 2) + 1.0;
    return jac;
  }

  // ---------------------------------------------------------------------------
  inline void PutJacobian(const Matrix &jac, TReal *out)
  {
    out[_xx] = static_cast<TReal>(jac(0, 0));
    out[_xy] = static_cast<TReal>(jac(0, 1));
    out[_xz] = static_cast<TReal>(jac(0, 2));
    out[_yx] = static_cast<TReal>(jac(1, 0));
    out[_yy] = static_cast<TReal>(jac(1, 1));
    out[_yz] = static_cast<TReal>(jac(1, 2));
    out[_zx] = static_cast<TReal>(jac(2, 0));
    out[_zy] = static_cast<TReal>(jac(2, 1));
    out[_zz] = static_cast<TReal>(jac(2, 2));
  }

  // ---------------------------------------------------------------------------
  inline void PutDetJacobian(const Matrix &jac, TReal *dj)
  {
    *dj = static_cast<TReal>(max(.0001, jac.Det3x3()));
  }

  // ---------------------------------------------------------------------------
  inline void PutLogJacobian(const Matrix &jac, TReal *lj)
  {
    double dj = jac.Det3x3();
    *lj = static_cast<TReal>(max(.0001, log(dj)));
  }

  // ---------------------------------------------------------------------------
  inline void PutDetAndLogJacobian(const Matrix &jac, TReal *dj, TReal *lj)
  {
    *dj = static_cast<TReal>(jac.Det3x3());
    if (*dj < TReal(.0001)) *dj = TReal(.0001);
    *lj = log(*dj);
  }
};

// -----------------------------------------------------------------------------
template <class TReal>
struct EvaluateJacobian : public EvaluateJacobianBase<TReal>
{
  void operator()(int i, int j, int k, int, TReal *out)
  {
    Matrix jac = this->Jacobian(i, j, k);
    this->PutJacobian(jac, out);
  }
};

// -----------------------------------------------------------------------------
template <class TReal>
struct EvaluateDetJacobian : public EvaluateJacobianBase<TReal>
{
  void operator()(int i, int j, int k, int, TReal *dj)
  {
    Matrix jac = this->Jacobian(i, j, k);
    this->PutDetJacobian(jac, dj);
  }
};

// -----------------------------------------------------------------------------
template <class TReal>
struct EvaluateLogJacobian : public EvaluateJacobianBase<TReal>
{
  void operator()(int i, int j, int k, int, TReal *lj)
  {
    Matrix jac = this->Jacobian(i, j, k);
    this->PutLogJacobian(jac, lj);
  }
};

// -----------------------------------------------------------------------------
template <class TReal>
struct EvaluateJacobianAndDet : public EvaluateJacobianBase<TReal>
{
  void operator()(int i, int j, int k, int, TReal *out, TReal *dj)
  {
    Matrix jac = this->Jacobian(i, j, k);
    this->PutJacobian   (jac, out);
    this->PutDetJacobian(jac, dj);
  }
};

// -----------------------------------------------------------------------------
template <class TReal>
struct EvaluateJacobianAndLog : public EvaluateJacobianBase<TReal>
{
  void operator()(int i, int j, int k, int, TReal *out, TReal *lj)
  {
    Matrix jac = this->Jacobian(i, j, k);
    this->PutJacobian   (jac, out);
    this->PutLogJacobian(jac, lj);
  }
};

// -----------------------------------------------------------------------------
template <class TReal>
struct EvaluateDetJacobianAndLog : public EvaluateJacobianBase<TReal>
{
  void operator()(int i, int j, int k, int, TReal *dj, TReal *lj)
  {
    Matrix jac = this->Jacobian(i, j, k);
    this->PutDetAndLogJacobian(jac, dj, lj);
  }
};

// -----------------------------------------------------------------------------
template <class TReal>
struct EvaluateJacobianAndDetAndLog : public EvaluateJacobianBase<TReal>
{
  void operator()(int i, int j, int k, int, TReal *out, TReal *dj, TReal *lj)
  {
    Matrix jac = this->Jacobian(i, j, k);
    this->PutJacobian         (jac, out);
    this->PutDetAndLogJacobian(jac, dj, lj);
  }
};

// -----------------------------------------------------------------------------
/// Voxel function for update of displacement field at each squaring step
template <class TInterpolator>
struct UpdateDisplacement : public VoxelFunction
{
  typedef typename TInterpolator::VoxelType TReal;

  const int            _x, _y, _z;
  const TInterpolator *_Displacement;

  UpdateDisplacement(const TInterpolator *f)
  :
    _x(0), _y(f->Input()->NumberOfSpatialVoxels()), _z(2 * _y),
    _Displacement(f)
  {}

  void operator ()(int i, int j, int k, int, const TReal *in, TReal *out)
  {
    double u[3] = {.0, .0, .0};
    _Displacement->Evaluate(u, static_cast<double>(i + in[_x]),
                               static_cast<double>(j + in[_y]),
                               static_cast<double>(k + in[_z]));
    out[_x] = in[_x] + static_cast<TReal>(u[0]);
    out[_y] = in[_y] + static_cast<TReal>(u[1]);
    out[_z] = in[_z] + static_cast<TReal>(u[2]);
  }
};

// -----------------------------------------------------------------------------
/// Voxel function for update of gradient image at each squaring step
template <class TDispInterpolator, class TGradInterpolator>
class UpdateGradientVoxelFunction : public VoxelFunction
{
protected:

  typedef typename TGradInterpolator::VoxelType TGradient;

  int _x, _y, _z;
  const TDispInterpolator *_Displacement;
  const TGradInterpolator *_Gradient;

  /// Initialize voxel function attributes
  void Initialize(const TDispInterpolator *u, const TGradInterpolator *g)
  {
    _x = 0;
    _y = g->Attributes().NumberOfSpatialPoints();
    _z = 2 * _y;
    _Displacement = u;
    _Gradient = g;
  }

  /// Transform displacement field lattice point using interpolated displacement field
  inline void Transform(double &u, double &v, double &w) const
  {
    double d[3];
    _Displacement->Evaluate(d, u, v, w);
    u += d[0], v += d[1], w += d[2];
  }

  /// Evaluate first order derivatives of deformation field at given lattice point
  inline void EvaluateJacobian(Matrix &jac, double x, double y, double z) const
  {
    _Displacement->Jacobian3D(jac, x, y, z);
  }

  /// Update gradient vector at current voxel given interpolated values
  inline void EvaluateGradient(TGradient *out, const TGradient *in, const Matrix &jac, double x, double y, double z) const
  {
    double g[3];
    _Gradient->Evaluate(g, x, y, z);
    out[_x] = in[_x] + jac(0, 0) * g[0] + jac(0, 1) * g[1] + jac(0, 2) * g[2];
    out[_y] = in[_y] + jac(1, 0) * g[0] + jac(1, 1) * g[1] + jac(1, 2) * g[2];
    out[_z] = in[_z] + jac(2, 0) * g[0] + jac(2, 1) * g[1] + jac(2, 2) * g[2];
  }
};

// -----------------------------------------------------------------------------
/// Voxel function for update of gradient image at each squaring step
template <class TDispInterpolator, class TGradInterpolator>
struct UpdateGradientWithEqualOutputLattice
: public UpdateGradientVoxelFunction<TDispInterpolator, TGradInterpolator>
{
  typedef typename TGradInterpolator::VoxelType TGradient;

  // ---------------------------------------------------------------------------
  /// Initialize voxel function attributes
  UpdateGradientWithEqualOutputLattice(const TDispInterpolator *u, const TGradInterpolator *g)
  {
    if (!u->Attributes().EqualInSpace(g->Attributes())) {
      Throw(ERR_InvalidArgument, __FUNCTION__, "Lattice attributes of displacement field and gradient image must be identical");
    }
    this->Initialize(u, g);
  }

  // ---------------------------------------------------------------------------
  void operator ()(int i, int j, int k, int, const TGradient *in, TGradient *out)
  {
    Matrix jac;
    double x = i, y = j, z = k;
    this->EvaluateJacobian(jac, x, y, z);
    this->Transform(x, y, z);
    this->EvaluateGradient(out, in, jac, x, y, z);
  }
};

// -----------------------------------------------------------------------------
/// Voxel function for update of gradient image at each squaring step
template <class TDispInterpolator, class TGradInterpolator>
struct UpdateGradientWithDifferentOutputLattice
: public UpdateGradientVoxelFunction<TDispInterpolator, TGradInterpolator>
{
  typedef typename TGradInterpolator::VoxelType TGradient;

  // ---------------------------------------------------------------------------
  /// Initialize voxel function attributes
  UpdateGradientWithDifferentOutputLattice(const TDispInterpolator *u, const TGradInterpolator *g)
  {
    this->Initialize(u, g);
  }

  /// Map coordinates from gradient image voxel to displacement field lattice
  inline void MapFromGradientImageToDisplacementField(double &u, double &v, double &w) const
  {
    this->_Gradient->ImageToWorld(u, v, w);
    this->_Displacement->WorldToImage(u, v, w);
  }

  /// Map coordinates from displacement field lattice to gradient image voxel
  inline void MapFromDisplacementFieldToGradientImage(double &u, double &v, double &w) const
  {
    this->_Displacement->ImageToWorld(u, v, w);
    this->_Gradient->WorldToImage(u, v, w);
  }

  // ---------------------------------------------------------------------------
  void operator ()(int i, int j, int k, int, const TGradient *in, TGradient *out)
  {
    Matrix jac;
    double x = i, y = j, z = k;
    this->MapFromGradientImageToDisplacementField(x, y, z);
    this->EvaluateJacobian(jac, x, y, z);
    this->Transform(x, y, z);
    this->MapFromDisplacementFieldToGradientImage(x, y, z);
    this->EvaluateGradient(out, in, jac, x, y, z);
  }
};

// -----------------------------------------------------------------------------
/// Voxel function for update of Jacobian w.r.t. x at each squaring step
template <class TReal>
struct UpdateJacobianBase : public VoxelFunction
{
  static const int     _x = 0, _xx = 0;
  int                  _y, _z, _xy, _xz, _yx, _yy, _yz, _zx, _zy, _zz;

  inline void Initialize(const ImageAttributes &attr)
  {
    // number of voxels
    const int n = attr.NumberOfSpatialPoints();
    // vector element offsets
    /* _x  = 0 */_y  = 1 * n; _z  = 2 * n;
    // matrix element offsets
    /* _xx = 0 */_xy = 1 * n; _xz = 2 * n;
    _yx = 3 * n; _yy = 4 * n; _yz = 5 * n;
    _zx = 6 * n; _zy = 7 * n; _zz = 8 * n;
  }

  // ---------------------------------------------------------------------------
  /// Evaluate Jacobian of current displacement field using linear interpolation
  inline void EvaluateJac(Matrix &jac, int i, int j, int k, const TReal *u)
  {
    jac.Initialize(3, 3);

    double x = static_cast<double>(i) + static_cast<double>(u[_x]);
    double y = static_cast<double>(j) + static_cast<double>(u[_y]);
    double z = static_cast<double>(k) + static_cast<double>(u[_z]);
    u -= _Domain->LatticeToIndex(i, j, k);

    i = ifloor(x);
    j = ifloor(y);
    k = ifloor(z);
    if (0 <= i && i < _Domain->X() - 1 &&
        0 <= j && j < _Domain->Y() - 1 &&
        0 <= k && k < _Domain->Z() - 1) {
      double w[3], wx[2], wy[2], wz[2], wd[2] = {-1., +1.};
      wx[1] = x - static_cast<double>(i), wx[0] = 1. - wx[1];
      wy[1] = y - static_cast<double>(j), wy[0] = 1. - wy[1];
      wz[1] = z - static_cast<double>(k), wz[0] = 1. - wz[1];
      for (int c = 0; c < 2; ++c) {
        int K = k + c;
        if (0 <= K && K < _Domain->_z) {
          for (int b = 0; b < 2; ++b) {
            int J = j + b;
            if (0 <= J && J < _Domain->_y) {
              for (int a = 0; a < 2; ++a) {
                int I = i + a;
                if (0 <= I && I < _Domain->_x) {
                  const TReal *d = u + _Domain->LatticeToIndex(I, J, K);
                  w[0] = wd[a] * wy[b] * wz[c];
                  w[1] = wx[a] * wd[b] * wz[c];
                  w[2] = wx[a] * wy[b] * wd[c];
                  for (int c = 0; c < 3; ++c) {
                    jac(0, c) += w[c] * d[_x];
                    jac(1, c) += w[c] * d[_y];
                    jac(2, c) += w[c] * d[_z];
                  }
                }
              }
            }
          }
        }
      }
    }

    jac(0, 0) += 1.;
    jac(1, 1) += 1.;
    jac(2, 2) += 1.;
  }

  // ---------------------------------------------------------------------------
  inline void UpdateJac(const Matrix &jac, const TReal *in, TReal *out)
  {
    out[_xx] = static_cast<TReal>(jac(0, 0)) * in[_xx]
             + static_cast<TReal>(jac(0, 1)) * in[_yx]
             + static_cast<TReal>(jac(0, 2)) * in[_zx];
    out[_xy] = static_cast<TReal>(jac(0, 0)) * in[_xy]
             + static_cast<TReal>(jac(0, 1)) * in[_yy]
             + static_cast<TReal>(jac(0, 2)) * in[_zy];
    out[_xz] = static_cast<TReal>(jac(0, 0)) * in[_xz]
             + static_cast<TReal>(jac(0, 1)) * in[_yz]
             + static_cast<TReal>(jac(0, 2)) * in[_zz];
    out[_yx] = static_cast<TReal>(jac(1, 0)) * in[_xx]
             + static_cast<TReal>(jac(1, 1)) * in[_yx]
             + static_cast<TReal>(jac(1, 2)) * in[_zx];
    out[_yy] = static_cast<TReal>(jac(1, 0)) * in[_xy]
             + static_cast<TReal>(jac(1, 1)) * in[_yy]
             + static_cast<TReal>(jac(1, 2)) * in[_zy];
    out[_yz] = static_cast<TReal>(jac(1, 0)) * in[_xz]
             + static_cast<TReal>(jac(1, 1)) * in[_yz]
             + static_cast<TReal>(jac(1, 2)) * in[_zz];
    out[_zx] = static_cast<TReal>(jac(2, 0)) * in[_xx]
             + static_cast<TReal>(jac(2, 1)) * in[_yx]
             + static_cast<TReal>(jac(2, 2)) * in[_zx];
    out[_zy] = static_cast<TReal>(jac(2, 0)) * in[_xy]
             + static_cast<TReal>(jac(2, 1)) * in[_yy]
             + static_cast<TReal>(jac(2, 2)) * in[_zy];
    out[_zz] = static_cast<TReal>(jac(2, 0)) * in[_xz]
             + static_cast<TReal>(jac(2, 1)) * in[_yz]
             + static_cast<TReal>(jac(2, 2)) * in[_zz];
  }

  // ---------------------------------------------------------------------------
  inline void UpdateDet(const Matrix &jac, const TReal *dj_in, TReal *dj_out)
  {
    (*dj_out) = (*dj_in) * static_cast<TReal>(max(.0001, jac.Det3x3()));
  }

  // ---------------------------------------------------------------------------
  inline void UpdateLog(const Matrix &jac, const TReal *lj_in, TReal *lj_out)
  {
    (*lj_out) = (*lj_in) + static_cast<TReal>(max(/*log(.0001)=*/-4.0, log(jac.Det3x3())));
  }

  // ---------------------------------------------------------------------------
  // Lorenzi, M., Ayache, N., Frisoni, G. B., & Pennec, X. (2013).
  // LCC-Demons: a robust and accurate symmetric diffeomorphic registration algorithm.
  // NeuroImage, 81, 470–83. doi:10.1016/j.neuroimage.2013.04.114
  inline void UpdateDetAndLog(const Matrix &jac, const TReal *, const TReal *lj_in, TReal *dj_out, TReal *lj_out)
  {
    (*lj_out) = (*lj_in) + static_cast<TReal>(log(max(.0001, jac.Det3x3())));
    (*dj_out) = exp(*lj_out);
  }

  // ---------------------------------------------------------------------------
  inline void Transform(double &x, double &y, double &z, const TReal *u) const
  {
    x += static_cast<double>(u[_x]);
    y += static_cast<double>(u[_y]);
    z += static_cast<double>(u[_z]);
  }
};

// -----------------------------------------------------------------------------
template <class TReal>
struct UpdateJacobian : public UpdateJacobianBase<TReal>
{
  void operator()(int i, int j, int k, int, const TReal *u, const TReal *jac_in, TReal *jac_out)
  {
    Matrix jac;
    this->EvaluateJac(jac, i, j, k, u);
    this->UpdateJac(jac, jac_in, jac_out);
  }
};

// -----------------------------------------------------------------------------
template <class TReal>
struct UpdateDetJacobian : public UpdateJacobianBase<TReal>
{
  void operator()(int i, int j, int k, int, const TReal *u, const TReal *dj_in, TReal *dj_out)
  {
    Matrix jac;
    this->EvaluateJac(jac, i, j, k, u);
    this->UpdateDet(jac, dj_in, dj_out);
  }
};

// -----------------------------------------------------------------------------
template <class TReal>
struct UpdateLogJacobian : public UpdateJacobianBase<TReal>
{
  void operator()(int i, int j, int k, int, const TReal *u, const TReal *lj_in, TReal *lj_out)
  {
    Matrix jac;
    this->EvaluateJac(jac, i, j, k, u);
    this->UpdateLog(jac, lj_in, lj_out);
  }
};

// -----------------------------------------------------------------------------
template <class TReal>
struct UpdateDetJacobianAndLog : public UpdateJacobianBase<TReal>
{
  void operator()(int i, int j, int k, int, const TReal *u,
                  const TReal *dj_in,  const TReal *lj_in,
                  const TReal *dj_out, TReal       *lj_out)
  {
    Matrix jac;
    this->EvaluateJac(jac, i, j, k, u);
    this->UpdateDetAndLog(jac, dj_in, lj_in, const_cast<TReal *>(dj_out), lj_out);
  }
};

// -----------------------------------------------------------------------------
template <class TReal>
struct UpdateJacobianAndDet : public UpdateJacobianBase<TReal>
{
  void operator()(int i, int j, int k, int, const TReal *u,
                  const TReal *jac_in,  const TReal *dj_in,
                  const TReal *jac_out, TReal       *dj_out)
  {
    Matrix jac;
    this->EvaluateJac(jac, i, j, k, u);
    this->UpdateJac(jac, jac_in, const_cast<TReal *>(jac_out));
    this->UpdateDet(jac, dj_in,  dj_out);
  }
};

// -----------------------------------------------------------------------------
template <class TReal>
struct UpdateJacobianAndLog : public UpdateJacobianBase<TReal>
{
  void operator()(int i, int j, int k, int, const TReal *u,
                  const TReal *jac_in,  const TReal *lj_in,
                  const TReal *jac_out, TReal       *lj_out)
  {
    Matrix jac;
    this->EvaluateJac(jac, i, j, k, u);
    this->UpdateJac(jac, jac_in, const_cast<TReal *>(jac_out));
    this->UpdateLog(jac, lj_in,  lj_out);
  }
};

// -----------------------------------------------------------------------------
template <class TReal>
struct UpdateJacobianAndDetAndLog : public UpdateJacobianBase<TReal>
{
  void operator()(int i, int j, int k, int, const TReal *u,
                  const TReal *jac_in,  const TReal *dj_in,  const TReal *lj_in,
                  const TReal *jac_out, const TReal *dj_out, TReal       *lj_out)
  {
    Matrix jac;
    this->EvaluateJac(jac, i, j, k, u);
    this->UpdateJac(jac, jac_in, const_cast<TReal *>(jac_out));
    this->UpdateDetAndLog(jac, dj_in, lj_in, const_cast<TReal *>(dj_out), lj_out);
  }
};

// -----------------------------------------------------------------------------
/// Voxel function for composition of output displacement with input displacement
template <class TInterpolator>
class ApplyInputDisplacement : public VoxelFunction
{
  typedef typename TInterpolator::VoxelType TReal;
  const TInterpolator *_OutputDisplacement;

public:

  ApplyInputDisplacement(const TInterpolator *disp)
  :
    _OutputDisplacement(disp)
  {}

  void operator ()(int i, int j, int k, int, const TReal *din, TReal *out)
  {
    // Get world coordinates of this voxel
    double x = i, y = j, z = k;
    _Domain->LatticeToWorld(x, y, z);
    // Add input displacements
    const int n = _Domain->NumberOfSpatialPoints();
    x += (*din), din += n;
    y += (*din), din += n;
    z += (*din);
    // Evaluate output displacement at mapped point
    double v[3];
    _OutputDisplacement->WorldToImage(x, y, z);
    _OutputDisplacement->Evaluate(v, x, y, z);
    (*out) = static_cast<TReal>(v[0]), out += n;
    (*out) = static_cast<TReal>(v[1]), out += n;
    (*out) = static_cast<TReal>(v[2]);
  }
};

// -----------------------------------------------------------------------------
/// Voxel function for composition of output displacement with input deformation
template <class TInterpolator>
class ApplyInputDeformation : public VoxelFunction
{
  typedef typename TInterpolator::VoxelType TReal;
  const TInterpolator *_OutputDisplacement;

  void Apply(const TReal *def, TReal *out) const
  {
    const int n = _Domain->NumberOfSpatialPoints();
    // Get mapped world coordinates of this voxel
    double x, y, z;
    x = static_cast<double>(*def), def += n;
    y = static_cast<double>(*def), def += n;
    z = static_cast<double>(*def);
    // Evaluate output displacement at mapped point
    double v[3];
    _OutputDisplacement->WorldToImage(x, y, z);
    _OutputDisplacement->Evaluate(v, x, y, z);
    (*out) = static_cast<TReal>(v[0]), out += n;
    (*out) = static_cast<TReal>(v[1]), out += n;
    (*out) = static_cast<TReal>(v[2]);
  }

public:

  ApplyInputDeformation(const TInterpolator *disp)
  :
    _OutputDisplacement(disp)
  {}

  template <class TImage>
  void operator ()(const TImage &, int, const TReal *def, TReal *out)
  {
    Apply(def, out);
  }

  void operator ()(int, int, int, int, const TReal *def, TReal *out)
  {
    Apply(def, out);
  }
};

// -----------------------------------------------------------------------------
/// Voxel function for resampling of output images
template <class TInterpolator>
struct ResampleOutput : public VoxelFunction
{
  typedef typename TInterpolator::VoxelType TReal;

  const ImageAttributes &_Domain;
  const TInterpolator   *_Image;
  const int              _NumberOfVoxels;
  const int              _NumberOfComponents;

  ResampleOutput(const ImageAttributes &domain, const TInterpolator *f)
  :
    _Domain(domain), _Image(f),
    _NumberOfVoxels(domain.NumberOfSpatialPoints()),
    _NumberOfComponents(_Image->Input()->T())
  {}

  void operator()(int i, int j, int k, int, TReal *value)
  {
    double x = i, y = j, z = k;
    _Domain.LatticeToWorld(x, y, z);
    _Image->Input()->WorldToImage(x, y, z);
    for (int l = 0; l < _NumberOfComponents; ++l, value += _NumberOfVoxels) {
      (*value) = static_cast<TReal>(_Image->Evaluate(x, y, z, l));
    }
  }
};

// -----------------------------------------------------------------------------
ExtrapolationMode GetExtrapolationMode(InterpolationMode imode)
{
  return (imode == Interpolation_BSpline ||
          imode == Interpolation_CubicBSpline ||
          imode == Interpolation_FastCubicBSpline)
         ? Extrapolation_Mirror : Extrapolation_NN;
}


} // anonymous namespace

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
template <class TReal>
ScalingAndSquaring<TReal>::ScalingAndSquaring()
:
  _InputVelocity(nullptr),
  _InputDisplacement(nullptr),
  _InputDeformation(nullptr),
  _InputGradient(nullptr),
  _OutputDisplacement(nullptr),
  _OutputDeformation(nullptr),
  _OutputJacobian(nullptr),
  _OutputDetJacobian(nullptr),
  _OutputLogJacobian(nullptr),
  _OutputGradient(nullptr),
  _Interpolation(Interpolation_Linear),
  _ComputeInterpolationCoefficients(true),
  _ComputeInverse(false),
  _UpperIntegrationLimit(1.0),
  _NumberOfSteps(0),
  _NumberOfSquaringSteps(0),
  _MaxScaledVelocity(.0),
  _Upsample(false),
  _SmoothBeforeDownsampling(false)
{
}

// -----------------------------------------------------------------------------
template <class TReal>
void ScalingAndSquaring<TReal>::Clear()
{
  _InterimDisplacement.reset();
  _InterimJacobian.reset();
  _InterimDetJacobian.reset();
  _InterimLogJacobian.reset();
}

// -----------------------------------------------------------------------------
template <class TReal>
ScalingAndSquaring<TReal>::~ScalingAndSquaring()
{
  Clear();
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
template <class TReal>
void ScalingAndSquaring<TReal>::Initialize()
{
  MIRTK_START_TIMING();

  // Check input/output
  if (!_InputVelocity) {
    cerr << "ScalingAndSquaring::Initialize: Missing input velocity field!" << endl;
    exit(1);
  }
  if (!_OutputAttributes) {
    if      (_InputDisplacement) _OutputAttributes = _InputDisplacement->Attributes();
    else if (_InputDeformation ) _OutputAttributes = _InputDeformation ->Attributes();
    else                         _OutputAttributes = _InputVelocity    ->Attributes();
  }
  _OutputAttributes._t = 1, _OutputAttributes._dt = .0;
  if (_InputDisplacement && !_InputDisplacement->Attributes().EqualInSpace(_OutputAttributes)) {
    cerr << "ScalingAndSquaring::Initialize: Output attributes do not match those of input displacement field" << endl;
    exit(1);
  }
  if (_InputDeformation && !_InputDeformation->Attributes().EqualInSpace(_OutputAttributes)) {
    cerr << "ScalingAndSquaring::Initialize: Output attributes do not match those of input deformation field" << endl;
    exit(1);
  }
  if (_InputGradient) {
    if (_InputGradient->T() != 3) {
      cerr << "ScalingAndSquaring::Initialize: Input gradient field must have 3 vector components (_t)" << endl;
      exit(1);
    }
  } else if (_OutputGradient) {
    cerr << "ScalingAndSquaring::Initialize: Input gradient field not specified, but output gradient image is set" << endl;
    exit(1);
  }
  // Check settings
  if (_NumberOfSteps <= 0 && _NumberOfSquaringSteps <= 0 && _MaxScaledVelocity <= .0) {
    cerr << "ScalingAndSquaring::Initialize: Either number of integration or squaring steps or maximum scaled velocity must be positive!" << endl;
    exit(1);
  }

  // Initialize input interpolator
  VelocityField velocity;
  velocity.Input     (_InputVelocity);
  velocity.Initialize(!_ComputeInterpolationCoefficients);

  // Attributes of intermediate and output images
  //
  // According to
  //   Bossa, M., Zacur, E., & Olmos, S. (2008). Algorithms for
  //   computing the group exponential of diffeomorphisms: Performance evaluation.
  //   In 2008 IEEE Computer Society Conference on Computer Vision and Pattern
  //   Recognition Workshops (pp. 1–8). IEEE. doi:10.1109/CVPRW.2008.4563005
  // an upsampling of the input velocity field by a factor of r * 2^d, where
  // r is the interpolation radius (i.e., 1 in case of linear interpolation)
  // and d is the dimension of the vector field (i.e., 3 for 3D) would reduce
  // the error made by using the approximation that the exponential of the
  // scaled velocity field is equal to the scaled velocity itself and especially
  // the larger sampling error at later stages of the squaring steps.
  if (!_InterimAttributes) _InterimAttributes = _OutputAttributes;
  _InterimAttributes._t = 1, _InterimAttributes._dt = .0;
  ImageAttributes attr = _InterimAttributes;
  if (_Upsample) {
    if (attr._x > 1) attr._x *= 2, attr._dx /= 2.0;
    if (attr._y > 1) attr._y *= 2, attr._dy /= 2.0;
    if (attr._z > 1) attr._z *= 2, attr._dz /= 2.0;
  }

  // Number of squaring steps
  if (_NumberOfSquaringSteps <= 0) {
    _NumberOfSquaringSteps = iceil(log(static_cast<double>(_NumberOfSteps)) / log(2.0));
  }
  if (_NumberOfSquaringSteps < 0) {
    _NumberOfSquaringSteps = 0; // i.e., 1 integration step only
  }

  // Initialize deformation field and increase number of squaring steps if needed
  // Note that input image may contain precomputed interpolation coefficients!
  {
    _InterimDisplacement.reset(new ImageType(attr, 3));
    velocity.Evaluate(*_InterimDisplacement);
  }

  TReal  vmax(0);
  TReal  scale = static_cast<TReal>(_UpperIntegrationLimit / pow(2.0, _NumberOfSquaringSteps));
  TReal *v     = _InterimDisplacement->Data();
  const int n  = 3 * attr.NumberOfSpatialPoints();
  for (int idx = 0; idx < n; ++idx, ++v) {
    (*v) *= scale;
    if (abs(*v) > vmax) vmax = abs(*v);
  }

  // Continue halfing input velocities as long as maximum absolute velocity
  // exceeds the specified maximum; skip if fixed number of steps
  if (_MaxScaledVelocity > .0) {
    TReal s(1);
    while ((vmax * s) > _MaxScaledVelocity) {
      s *= TReal(.5);
      _NumberOfSquaringSteps++;
    }
    if (s != TReal(1)) {
      (*_InterimDisplacement) *= s;
      vmax                    *= s;
      scale                   *= s;
    }
  }

  // Update number of steps
  _NumberOfSteps = static_cast<int>(pow(2, _NumberOfSquaringSteps));

  // Convert scaled displacements to voxel units
  {
    const Matrix m = _InterimAttributes.GetWorldToImageMatrix();
    ParallelForEachVoxel(InPlaceMatrixVectorProduct(&m), _InterimDisplacement.get());
  }

  // Compute derivatives of initial deformation w.r.t. x and/or its (log) determinant
  int jac_mode = 0;
  if (_OutputJacobian)    jac_mode += 1;
  if (_OutputDetJacobian) jac_mode += 2;
  if (_OutputLogJacobian) jac_mode += 4;

  if (jac_mode > 0) {
    if (jac_mode & 1) {
      _InterimJacobian.reset(new ImageType(attr, 9));
    } else {
      _InterimJacobian.reset();
    }
    if (jac_mode & 2) {
      _InterimDetJacobian.reset(new ImageType(attr));
    } else {
      _InterimDetJacobian.reset();
    }
    if (jac_mode & 4) {
      _InterimLogJacobian.reset(new ImageType(attr));
    } else {
      _InterimLogJacobian.reset();
    }
    switch (jac_mode) {
      case 1: {
        EvaluateJacobian<TReal> eval;
        eval.Initialize(attr, &velocity, scale);
        ParallelForEachVoxel(attr, _InterimJacobian.get(), eval);
      } break;
      case 2: {
        EvaluateDetJacobian<TReal> eval;
        eval.Initialize(attr, &velocity, scale);
        ParallelForEachVoxel(attr, _InterimDetJacobian.get(), eval);
      } break;
      case 3: {
        EvaluateJacobianAndDet<TReal> eval;
        eval.Initialize(attr, &velocity, scale);
        ParallelForEachVoxel(attr, _InterimJacobian.get(), _InterimDetJacobian.get(), eval);
      } break;
      case 4: {
        EvaluateLogJacobian<TReal> eval;
        eval.Initialize(attr, &velocity, scale);
        ParallelForEachVoxel(attr, _InterimLogJacobian.get(), eval);
      } break;
      case 5: {
        EvaluateJacobianAndLog<TReal> eval;
        eval.Initialize(attr, &velocity, scale);
        ParallelForEachVoxel(attr, _InterimJacobian.get(), _InterimLogJacobian.get(), eval);
      } break;
      case 6: {
        EvaluateDetJacobianAndLog<TReal> eval;
        eval.Initialize(attr, &velocity, scale);
        ParallelForEachVoxel(attr, _InterimDetJacobian.get(), _InterimLogJacobian.get(), eval);
      } break;
      case 7: {
        EvaluateJacobianAndDetAndLog<TReal> eval;
        eval.Initialize(attr, &velocity, scale);
        ParallelForEachVoxel(attr, _InterimJacobian.get(), _InterimDetJacobian.get(), _InterimLogJacobian.get(), eval);
      } break;
    }
  } else {
    _InterimJacobian.reset();
    _InterimDetJacobian.reset();
    _InterimLogJacobian.reset();
  }

  // Initialize exponentiated output gradient
  if (_OutputGradient) {
    _OutputGradient->Initialize(_InputGradient->Attributes(), _InputGradient->T());
    _OutputGradient->CopyFrom(_InputGradient->Data());
    // Convert input gradient to derivatives w.r.t. displacement field lattice
    Matrix m = _InterimAttributes.GetWorldToImageMatrix();
    m *= scale; // pre-multiply matrix elements with scaling factor
    ParallelForEachVoxel(InPlaceMatrixVectorProduct(&m), _OutputGradient);
  }

  // Use output deformation field as temporary output displacement field
  if (!_OutputDisplacement) _OutputDisplacement = _OutputDeformation;

  MIRTK_DEBUG_TIMING(5, "scaling step");
}

// -----------------------------------------------------------------------------
template <class TReal>
void ScalingAndSquaring<TReal>
::Resample(const ImageType *interim, ImageType *output)
{
  if (output) {
    // Attributes of final output images
    const ImageAttributes &attr = _OutputAttributes;
    if (_InputDisplacement || _InputDeformation || attr != _InterimAttributes) {
      // Initialize interim vector field interpolator
      UniquePtr<VectorField> interp;
      interp.reset(VectorField::New(_Interpolation, GetExtrapolationMode(_Interpolation), interim));
      // Either apply input deformation field while resampling output...
      if (_InputDeformation) {
        interp->Initialize();
        output->Initialize(attr, interim->T());
        ApplyInputDeformation<VectorField> warp(interp.get());
        ParallelForEachVoxel(_InputDeformation, output, warp);
      // ...or apply input displacement field while resampling output
      } else if (_InputDisplacement) {
        interp->Initialize();
        output->Initialize(attr, interim->T());
        ApplyInputDisplacement<VectorField> warp(interp.get());
        ParallelForEachVoxel(attr, _InputDisplacement, output, warp);
      // ...or resample intermediate output if necessary
      } else {
        if (_SmoothBeforeDownsampling) {
          const double sx = (attr._dx > interim->XSize()) ? attr._dx / 2.0 : .0;
          const double sy = (attr._dy > interim->YSize()) ? attr._dy / 2.0 : .0;
          const double sz = (attr._dz > interim->ZSize()) ? attr._dz / 2.0 : .0;
          if (sx || sy || sz) {
            GaussianBlurring<TReal> blurring(sx, sy, sz);
            blurring.Input (interim);
            blurring.Output(const_cast<ImageType *>(interim));
            blurring.Run();
          }
        }
        interp->Initialize();
        output->Initialize(attr, interim->T());
        ResampleOutput<VectorField> resample(attr, interp.get());
        ParallelForEachVoxel(attr, output, resample);
        TReal min_val, max_val;
        interim->GetMinMax(min_val, max_val);
        for (int vox = 0; vox < output->NumberOfVoxels(); ++vox) {
          output->Put(vox, clamp(output->Get(vox), min_val, max_val));
        }
      }
    } else {
      // Otherwise, just copy intermediate output
      (*output) = (*interim);
    }
  }
}

// -----------------------------------------------------------------------------
template <class TReal>
void ScalingAndSquaring<TReal>::Finalize()
{
  MIRTK_START_TIMING();

  if (_OutputDisplacement) {
    // Convert final displacements back to world units if output requested
    const Matrix m = _InterimAttributes.GetImageToWorldMatrix();
    ParallelForEachVoxel(InPlaceMatrixVectorProduct(&m), _InterimDisplacement.get());
    // Resample final displacement field to requested output size
    Resample(_InterimDisplacement.get(), _OutputDisplacement);
  }

  // Resample final output fields to requested output size
  Resample(_InterimJacobian.get(),     _OutputJacobian);
  Resample(_InterimDetJacobian.get(),  _OutputDetJacobian);
  Resample(_InterimLogJacobian.get(),  _OutputLogJacobian);

  // Compute output deformation field
  if (_OutputDeformation) {
    if (_OutputDeformation != _OutputDisplacement) {
      _OutputDeformation->Initialize(_OutputAttributes, 3);
    }
    ConvertToDeformation3D<TReal> plusx(_OutputAttributes);
    ParallelForEachVoxel(_OutputAttributes, _OutputDisplacement, _OutputDeformation, plusx);
    if (_OutputDeformation == _OutputDisplacement) _OutputDisplacement = nullptr;
  }

  // Multiply output Jacobian by Jacobian of input deformation
  if ((_OutputJacobian || _OutputDetJacobian || _OutputLogJacobian) && (_InputDisplacement || _InputDeformation)) {
    // TODO: Multiply output Jacobian by Jacobian of input deformation
    // TODO: Multiply output determinant of Jacobian by determinant of Jacobian of input deformation
    cerr << "WARNING: ScalingAndSquaring::Finalize: Output Jacobian does not include the input deformation!" << endl;
    cerr << "         Be aware that this may be fixed in the future and thus change the output." << endl;
    exit(1);
  }

  // Convert output gradients back to world vectors
  if (_OutputGradient) {
    Matrix A = _InterimAttributes.GetImageToWorldMatrix();
    ParallelForEachVoxel(InPlaceMatrixVectorProduct(&A), _OutputGradient);
  }

  // Free allocated memory
  Clear();

  MIRTK_DEBUG_TIMING(5, "finalization of scaling and squaring");
}

// -----------------------------------------------------------------------------
template <class TReal>
void ScalingAndSquaring<TReal>::Run()
{
  // Do the initial set up and scaling
  this->Initialize();

  MIRTK_START_TIMING();

  // Get common attributes of intermediate images
  const ImageAttributes &attr = _InterimDisplacement->Attributes();

  int jac_mode = 0;
  if (_InterimJacobian   ) jac_mode += 1;
  if (_InterimDetJacobian) jac_mode += 2;
  if (_InterimLogJacobian) jac_mode += 4;

  // Either allocate required temporary images or use provided output
  ImageType *disp   = _OutputDisplacement;
  ImageType *jac3x3 = nullptr;
  ImageType *detjac = nullptr;
  ImageType *logjac = nullptr;

  UniquePtr<ImageType> t_disp;
  if (!disp) {
    disp = new ImageType;
    t_disp.reset(disp);
  }
  disp->Initialize(attr, 3);
  UniquePtr<VectorField> f_disp(VectorField::New(_Interpolation, GetExtrapolationMode(_Interpolation), _InterimDisplacement.get()));

  UniquePtr<ImageType> t_jac3x3;
  if (_InterimJacobian) {
    if (_OutputJacobian) {
      jac3x3 = _OutputJacobian;
    } else {
      jac3x3 = new ImageType;
      t_jac3x3.reset(jac3x3);
    }
    jac3x3->Initialize(attr, 9);
  }

  UniquePtr<ImageType> t_detjac;
  if (_InterimDetJacobian) {
    if (_OutputDetJacobian) {
      detjac = _OutputDetJacobian;
    } else {
      detjac = new ImageType;
      t_detjac.reset(detjac);
    }
    detjac->Initialize(attr, 1);
  }

  UniquePtr<ImageType> t_logjac;
  if (_InterimLogJacobian) {
    if (_OutputLogJacobian) {
      logjac = _OutputLogJacobian;
    } else {
      logjac = new ImageType;
      t_logjac.reset(logjac);
    }
    logjac->Initialize(attr, 1);
  }

  UniquePtr<ImageType> grad;
  UniquePtr<VectorField> f_grad;
  bool grad_attr_equal_interim_attr = false;
  if (_OutputGradient) {
    grad.reset(new ImageType(_OutputGradient->Attributes()));
    f_grad.reset(VectorField::New(_Interpolation, Extrapolation_Const, grad.get()));
    grad_attr_equal_interim_attr = _OutputGradient->Attributes().EqualInSpace(attr);
  }

  // Do the squaring steps
  int n = _NumberOfSquaringSteps;
  while (n--) {
    // (Re-)initialize interpolators
    f_disp->Initialize();
    // Compose displacement field with itself
    {
      UpdateDisplacement<VectorField> update(f_disp.get());
      ParallelForEachVoxel(attr, _InterimDisplacement.get(), disp, update);
    }
    // Update Jacobian and/or determinants
    switch (jac_mode) {
      case 1: {
        UpdateJacobian<TReal> update;
        update.Initialize(attr);
        ParallelForEachVoxel(attr, _InterimDisplacement.get(), _InterimJacobian.get(), jac3x3, update);
      } break;
      case 2: {
        UpdateDetJacobian<TReal> update;
        update.Initialize(attr);
        ParallelForEachVoxel(attr, _InterimDisplacement.get(), _InterimDetJacobian.get(), detjac, update);
      } break;
      case 3: {
        UpdateJacobianAndDet<TReal> update;
        update.Initialize(attr);
        ParallelForEachVoxel(attr, _InterimDisplacement.get(), _InterimJacobian.get(), _InterimDetJacobian.get(), jac3x3, detjac, update);
      } break;
      case 4: {
        UpdateLogJacobian<TReal> update;
        update.Initialize(attr);
        ParallelForEachVoxel(attr, _InterimDisplacement.get(), _InterimLogJacobian.get(), logjac, update);
      } break;
      case 5: {
        UpdateJacobianAndLog<TReal> update;
        update.Initialize(attr);
        ParallelForEachVoxel(attr, _InterimDisplacement.get(), _InterimJacobian.get(), _InterimLogJacobian.get(), jac3x3, logjac, update);
      } break;
      case 6: {
        UpdateDetJacobianAndLog<TReal> update;
        update.Initialize(attr);
        ParallelForEachVoxel(attr, _InterimDisplacement.get(), _InterimDetJacobian.get(), _InterimLogJacobian.get(), detjac, logjac, update);
      } break;
      case 7: {
        UpdateJacobianAndDetAndLog<TReal> update;
        update.Initialize(attr);
        ParallelForEachVoxel(attr, _InterimDisplacement.get(), _InterimJacobian.get(), _InterimDetJacobian.get(), _InterimLogJacobian.get(), jac3x3, detjac, logjac, update);
      } break;
    }
    // Propage input gradient using current deformation
    if (_OutputGradient && n >= 0) {
      typedef GenericLinearInterpolateImageFunction<ImageType> DispInterpolator;
      UniquePtr<DispInterpolator> f_disp(new DispInterpolator());
      f_disp->Input(_InterimDisplacement.get());
      f_disp->Initialize();
      typedef GenericLinearInterpolateImageFunction<ImageType> GradInterpolator;
      UniquePtr<GradInterpolator> f_grad(new GradInterpolator());
      f_grad->Input(_OutputGradient);
      f_grad->Initialize();
      if (grad_attr_equal_interim_attr) {
        UpdateGradientWithEqualOutputLattice<DispInterpolator, GradInterpolator> update(f_disp.get(), f_grad.get());
        ParallelForEachVoxel(attr, _OutputGradient, grad.get(), update);
      } else {
        UpdateGradientWithDifferentOutputLattice<DispInterpolator, GradInterpolator> update(f_disp.get(), f_grad.get());
        ParallelForEachVoxel(_OutputGradient->Attributes(), _OutputGradient, grad.get(), update);
      }
      _OutputGradient->CopyFrom(grad->Data());
    }
    // Update intermediate output images
    _InterimDisplacement                        ->CopyFrom(disp  ->Data());
    if (_InterimJacobian   ) _InterimJacobian   ->CopyFrom(jac3x3->Data());
    if (_InterimDetJacobian) _InterimDetJacobian->CopyFrom(detjac->Data());
    if (_InterimLogJacobian) _InterimLogJacobian->CopyFrom(logjac->Data());
  }

  MIRTK_DEBUG_TIMING(5, "squaring steps"
                           " (d="    << (_OutputDisplacement ? "on" : "off")
                        << ", J="    << (_OutputJacobian     ? "on" : "off")
                        << ", detJ=" << (_OutputDetJacobian  ? "on" : "off")
                        << ", logJ=" << (_OutputLogJacobian  ? "on" : "off")
                        << ", grad=" << (_OutputGradient     ? "on" : "off") << ")");

  this->Finalize();
}

// =============================================================================
// Explicit template instantiations
// =============================================================================

template class ScalingAndSquaring<float>;
template class ScalingAndSquaring<double>;


} // namespace mirtk
