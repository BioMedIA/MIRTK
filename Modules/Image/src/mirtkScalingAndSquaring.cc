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

#include <mirtkScalingAndSquaring.h>

#include <mirtkMath.h>
#include <mirtkVoxelFunction.h>
#include <mirtkInterpolateImageFunction.hxx>
#include <mirtkGaussianBlurring.h>
#include <mirtkMatrix.h>
#include <mirtkParallel.h>
#include <mirtkProfiling.h>
#include <mirtkDeallocate.h>

#include <algorithm>


namespace mirtk {


// =============================================================================
// Auxiliary voxel functions for parallel execution
// =============================================================================

namespace ScalingAndSquaringUtils {


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
template <class TReal>
struct ConvertToVoxelUnits3D : public VoxelFunction
{
  const int              _x, _y, _z;
  const ImageAttributes &_Domain;

  ConvertToVoxelUnits3D(const ImageAttributes &domain)
  :
    _x(0), _y(domain.NumberOfSpatialPoints()), _z(2 * _y), _Domain(domain)
  {}

  void operator()(int i, int j, int k, int, const TReal *d, TReal *u)
  {
    double x = i, y = j, z = k;
    _Domain.LatticeToWorld(x, y, z);
    x += static_cast<double>(d[_x]);
    y += static_cast<double>(d[_y]);
    z += static_cast<double>(d[_z]);
    _Domain.WorldToLattice(x, y, z);
    u[_x] = static_cast<TReal>(x - i);
    u[_y] = static_cast<TReal>(y - j);
    u[_z] = static_cast<TReal>(z - k);
  }
};

// -----------------------------------------------------------------------------
template <class TReal>
struct ConvertToWorldUnits3D : public VoxelFunction
{
  const int              _x, _y, _z;
  const ImageAttributes &_Domain;

  ConvertToWorldUnits3D(const ImageAttributes &domain)
  :
    _x(0), _y(domain.NumberOfSpatialPoints()), _z(2 * _y), _Domain(domain)
  {}

  void operator()(int i, int j, int k, int, const TReal *u, TReal *d)
  {
    double x1 = i, y1 = j, z1 = k;
    _Domain.LatticeToWorld(x1, y1, z1);
    double x2 = i + u[_x];
    double y2 = j + u[_y];
    double z2 = k + u[_z];
    _Domain.LatticeToWorld(x2, y2, z2);
    d[_x] = static_cast<TReal>(x2 - x1);
    d[_y] = static_cast<TReal>(y2 - y1);
    d[_z] = static_cast<TReal>(z2 - z1);
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
                         const VelocityField       *v,
                         double                     s = 1.0)
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
    _VelocityField->Input()->WorldToImage(x, y, z);
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
    *dj = static_cast<TReal>(jac.Det3x3());
  }

  // ---------------------------------------------------------------------------
  inline void PutLogJacobian(const Matrix &jac, TReal *lj)
  {
    double dj = jac.Det3x3();
    if (dj < .0001) dj = .0001;
    *lj = static_cast<TReal>(log(dj));
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

  void operator()(int i, int j, int k, int, const TReal *in, TReal *out)
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
/// Voxel function for update of Jacobian w.r.t. x at each squaring step
template <class TInterpolator>
struct UpdateJacobianBase : public VoxelFunction
{
  typedef typename TInterpolator::VoxelType TReal;

  static const int     _x = 0, _xx = 0;
  int                  _y, _z, _xy, _xz, _yx, _yy, _yz, _zx, _zy, _zz;
  const TInterpolator *_Jacobian;
  const TInterpolator *_DetJacobian;
  const TInterpolator *_LogJacobian;

  inline void Initialize(const TInterpolator *j,
                         const TInterpolator *dj,
                         const TInterpolator *lj)
  {
    // number of voxels
    int n = 0;
    if (j ) n = j ->Input()->NumberOfSpatialVoxels();
    if (dj) n = dj->Input()->NumberOfSpatialVoxels();
    if (lj) n = lj->Input()->NumberOfSpatialVoxels();
    // vector element offsets
    /* _x  = 0 */_y  = 1 * n; _z  = 2 * n;
    // matrix element offsets
    /* _xx = 0 */_xy = 1 * n; _xz = 2 * n;
    _yx = 3 * n; _yy = 4 * n; _yz = 5 * n;
    _zx = 6 * n; _zy = 7 * n; _zz = 8 * n;
    // interpolators
    _Jacobian    = j;
    _DetJacobian = dj;
    _LogJacobian = lj;
  }

  // ---------------------------------------------------------------------------
  inline void Transform(double &x, double &y, double &z, const TReal *u)
  {
    x += static_cast<double>(u[_x]);
    y += static_cast<double>(u[_y]);
    z += static_cast<double>(u[_z]);
  }

  // ---------------------------------------------------------------------------
  inline void UpdateJac(double x, double y, double z, const TReal *in, TReal *out)
  {
    if (_Jacobian->IsInside(x, y, z)) {
      double jac[9];
      _Jacobian->EvaluateInside(jac, x, y, z);
      out[_xx] = static_cast<TReal>(jac[0]) * in[_xx]
               + static_cast<TReal>(jac[1]) * in[_yx]
               + static_cast<TReal>(jac[2]) * in[_zx];
      out[_xy] = static_cast<TReal>(jac[0]) * in[_xy]
               + static_cast<TReal>(jac[1]) * in[_yy]
               + static_cast<TReal>(jac[2]) * in[_zy];
      out[_xz] = static_cast<TReal>(jac[0]) * in[_xz]
               + static_cast<TReal>(jac[1]) * in[_yz]
               + static_cast<TReal>(jac[2]) * in[_zz];
      out[_yx] = static_cast<TReal>(jac[3]) * in[_xx]
               + static_cast<TReal>(jac[4]) * in[_yx]
               + static_cast<TReal>(jac[5]) * in[_zx];
      out[_yy] = static_cast<TReal>(jac[3]) * in[_xy]
               + static_cast<TReal>(jac[4]) * in[_yy]
               + static_cast<TReal>(jac[5]) * in[_zy];
      out[_yz] = static_cast<TReal>(jac[3]) * in[_xz]
               + static_cast<TReal>(jac[4]) * in[_yz]
               + static_cast<TReal>(jac[5]) * in[_zz];
      out[_zx] = static_cast<TReal>(jac[6]) * in[_xx]
               + static_cast<TReal>(jac[7]) * in[_yx]
               + static_cast<TReal>(jac[8]) * in[_zx];
      out[_zy] = static_cast<TReal>(jac[6]) * in[_xy]
               + static_cast<TReal>(jac[7]) * in[_yy]
               + static_cast<TReal>(jac[8]) * in[_zy];
      out[_zz] = static_cast<TReal>(jac[6]) * in[_xz]
               + static_cast<TReal>(jac[7]) * in[_yz]
               + static_cast<TReal>(jac[8]) * in[_zz];
    } else {
      out[_xx] = in[_xx]; out[_xy] = in[_xy]; out[_xz] = in[_xz];
      out[_yx] = in[_yx]; out[_yy] = in[_yy]; out[_yz] = in[_yz];
      out[_zx] = in[_zx]; out[_zy] = in[_zy]; out[_zz] = in[_zz];
    }
  }

  // ---------------------------------------------------------------------------
  inline void UpdateDet(double x, double y, double z, const TReal *dj_in, TReal *dj_out)
  {
    if (_DetJacobian->IsInside(x, y, z)) {
      (*dj_out) = (*dj_in) * static_cast<TReal>(max(.0001, _DetJacobian->EvaluateInside(x, y, z)));
    } else {
      (*dj_out) = (*dj_in);
    }
  }

  // ---------------------------------------------------------------------------
  inline void UpdateLog(double x, double y, double z, const TReal *lj_in, TReal *lj_out)
  {
    if (_LogJacobian->IsInside(x, y, z)) {
      (*lj_out) = (*lj_in) + static_cast<TReal>(max(/*log(.0001)=*/-4.0, _LogJacobian->EvaluateInside(x, y, z)));
    } else {
      (*lj_out) = (*lj_in);
    }
  }

  // ---------------------------------------------------------------------------
  // Lorenzi, M., Ayache, N., Frisoni, G. B., & Pennec, X. (2013).
  // LCC-Demons: a robust and accurate symmetric diffeomorphic registration algorithm.
  // NeuroImage, 81, 470–83. doi:10.1016/j.neuroimage.2013.04.114
  inline void UpdateDetAndLog(double x, double y, double z,
                              const TReal *dj_in, const TReal *lj_in,
                              TReal *dj_out, TReal *lj_out)
  {
    if (_DetJacobian->IsInside(x, y, z)) {
      (*lj_out) = (*lj_in) + static_cast<TReal>(log(max(.0001, _DetJacobian->EvaluateInside(x, y, z))));
      (*dj_out) = exp(*lj_out);
    } else {
      (*lj_out) = (*lj_in);
      (*dj_out) = (*dj_in);
    }
  }
};

// -----------------------------------------------------------------------------
template <class TInterpolator>
struct UpdateJacobian : public UpdateJacobianBase<TInterpolator>
{
  typedef typename TInterpolator::VoxelType TReal;
  void operator()(int i, int j, int k, int,
                  const TReal *u, const TReal *in, TReal *out)
  {
    double x = i, y = j, z = k;
    this->Transform(x, y, z, u);
    this->UpdateJac(x, y, z, in, out);
  }
};

// -----------------------------------------------------------------------------
template <class TInterpolator>
struct UpdateDetJacobian : public UpdateJacobianBase<TInterpolator>
{
  typedef typename TInterpolator::VoxelType TReal;
  void operator()(int i, int j, int k, int,
                  const TReal *u, const TReal *dj_in, TReal *dj_out)
  {
    double x = i, y = j, z = k;
    this->Transform(x, y, z, u);
    this->UpdateDet(x, y, z, dj_in, dj_out);
  }
};

// -----------------------------------------------------------------------------
template <class TInterpolator>
struct UpdateLogJacobian : public UpdateJacobianBase<TInterpolator>
{
  typedef typename TInterpolator::VoxelType TReal;
  void operator()(int i, int j, int k, int,
                  const TReal *u, const TReal *lj_in, TReal *lj_out)
  {
    double x = i, y = j, z = k;
    this->Transform(x, y, z, u);
    this->UpdateLog(x, y, z, lj_in, lj_out);
  }
};

// -----------------------------------------------------------------------------
template <class TInterpolator>
struct UpdateDetJacobianAndLog : public UpdateJacobianBase<TInterpolator>
{
  typedef typename TInterpolator::VoxelType TReal;
  void operator()(int i, int j, int k, int, const TReal *u,
                  const TReal *dj_in,  const TReal *lj_in,
                  const TReal *dj_out, TReal       *lj_out)
  {
    double x = i, y = j, z = k;
    this->Transform(x, y, z, u);
    this->UpdateDetAndLog(x, y, z, dj_in, lj_in, const_cast<TReal *>(dj_out), lj_out);
  }
};

// -----------------------------------------------------------------------------
template <class TInterpolator>
struct UpdateJacobianAndDet : public UpdateJacobianBase<TInterpolator>
{
  typedef typename TInterpolator::VoxelType TReal;
  void operator()(int i, int j, int k, int, const TReal *u,
                  const TReal *in,  const TReal *dj_in,
                  const TReal *out, TReal       *dj_out)
  {
    double x = i, y = j, z = k;
    this->Transform(x, y, z, u);
    this->UpdateJac(x, y, z, in, const_cast<TReal *>(out));
    this->UpdateDet(x, y, z, dj_in, dj_out);
  }
};

// -----------------------------------------------------------------------------
template <class TInterpolator>
struct UpdateJacobianAndLog : public UpdateJacobianBase<TInterpolator>
{
  typedef typename TInterpolator::VoxelType TReal;
  void operator()(int i, int j, int k, int, const TReal *u,
                  const TReal *in,  const TReal *lj_in,
                  const TReal *out, TReal       *lj_out)
  {
    double x = i, y = j, z = k;
    this->Transform(x, y, z, u);
    this->UpdateJac(x, y, z, in, const_cast<TReal *>(out));
    this->UpdateLog(x, y, z, lj_in, lj_out);
  }
};

// -----------------------------------------------------------------------------
template <class TInterpolator>
struct UpdateJacobianAndDetAndLog : public UpdateJacobianBase<TInterpolator>
{
  typedef typename TInterpolator::VoxelType TReal;
  void operator()(int i, int j, int k, int, const TReal *u,
                  const TReal *in,  const TReal *dj_in,  const TReal *lj_in,
                  const TReal *out, const TReal *dj_out, TReal       *lj_out)
  {
    double x = i, y = j, z = k;
    this->Transform(x, y, z, u);
    this->UpdateJac(x, y, z, in, const_cast<TReal *>(out));
    this->UpdateDetAndLog(x, y, z, dj_in, lj_in, const_cast<TReal *>(dj_out), lj_out);
  }
};

// -----------------------------------------------------------------------------
/// Voxel function for update of Jacobian w.r.t. v at each squaring step
template <class TInterpolator>
struct UpdateJacobianDOFs : public VoxelFunction
{
  typedef typename TInterpolator::VoxelType TReal;

  const int            _x, _y, _z;
  const int            _xx, _xy, _xz, _yx, _yy, _yz, _zx, _zy, _zz;
  const TInterpolator *_Jacobian;
  const TInterpolator *_JacobianDOFs;

  UpdateJacobianDOFs(const TInterpolator *f, const TInterpolator *g)
  :
    _x(0), _y(f->Input()->NumberOfSpatialVoxels()), _z(2 * _y),
    _xx(0      ), _xy(     _y), _xz(2 * _xy),
    _yx(3 * _xy), _yy(4 * _xy), _yz(5 * _xy),
    _zx(6 * _xy), _zy(7 * _xy), _zz(8 * _xy),
    _Jacobian(f), _JacobianDOFs(g)
  {}

  void operator()(int i, int j, int k, int, const TReal *u, const TReal *in, TReal *out)
  {
    double x = i + u[_x], y = j + u[_y], z = k + u[_z], jac[9], jacdof[9];
    _Jacobian    ->Evaluate(jac,    x, y, z);
    _JacobianDOFs->Evaluate(jacdof, x, y, z);
    out[_xx] = static_cast<TReal>(jac[0]) * in[_xx]
             + static_cast<TReal>(jac[1]) * in[_yx]
             + static_cast<TReal>(jac[2]) * in[_zx]
             + static_cast<TReal>(jacdof[0]);
    out[_xy] = static_cast<TReal>(jac[0]) * in[_xy]
             + static_cast<TReal>(jac[1]) * in[_yy]
             + static_cast<TReal>(jac[2]) * in[_zy]
             + static_cast<TReal>(jacdof[1]);
    out[_xz] = static_cast<TReal>(jac[0]) * in[_xz]
             + static_cast<TReal>(jac[1]) * in[_yz]
             + static_cast<TReal>(jac[2]) * in[_zz]
             + static_cast<TReal>(jacdof[2]);
    out[_yx] = static_cast<TReal>(jac[3]) * in[_xx]
             + static_cast<TReal>(jac[4]) * in[_yx]
             + static_cast<TReal>(jac[5]) * in[_zx]
             + static_cast<TReal>(jacdof[3]);
    out[_yy] = static_cast<TReal>(jac[3]) * in[_xy]
             + static_cast<TReal>(jac[4]) * in[_yy]
             + static_cast<TReal>(jac[5]) * in[_zy]
             + static_cast<TReal>(jacdof[4]);
    out[_yz] = static_cast<TReal>(jac[3]) * in[_xz]
             + static_cast<TReal>(jac[4]) * in[_yz]
             + static_cast<TReal>(jac[5]) * in[_zz]
             + static_cast<TReal>(jacdof[5]);
    out[_zx] = static_cast<TReal>(jac[6]) * in[_xx]
             + static_cast<TReal>(jac[7]) * in[_yx]
             + static_cast<TReal>(jac[8]) * in[_zx]
             + static_cast<TReal>(jacdof[6]);
    out[_zy] = static_cast<TReal>(jac[6]) * in[_xy]
             + static_cast<TReal>(jac[7]) * in[_yy]
             + static_cast<TReal>(jac[8]) * in[_zy]
             + static_cast<TReal>(jacdof[7]);
    out[_zz] = static_cast<TReal>(jac[6]) * in[_xz]
             + static_cast<TReal>(jac[7]) * in[_yz]
             + static_cast<TReal>(jac[8]) * in[_zz]
             + static_cast<TReal>(jacdof[8]);
  }
};

// -----------------------------------------------------------------------------
/// Voxel function for composition of output displacement with input displacement
template <class TInterpolator>
struct ApplyInputDisplacement : public VoxelFunction
{
  typedef typename TInterpolator::VoxelType TReal;

  const int              _x, _y, _z;
  const ImageAttributes &_Domain;
  const TInterpolator   *_Image;
  const int              _NumberOfComponents;

  ApplyInputDisplacement(const ImageAttributes &domain, const TInterpolator *f)
  :
    _x(0), _y(domain.NumberOfSpatialPoints()), _z(2 * _y),
    _Domain(domain), _Image(f), _NumberOfComponents(_Image->Input()->T())
  {}

  void operator()(int i, int j, int k, int, const TReal *d, TReal *out)
  {
    double x = i, y = j, z = k;
    _Domain.LatticeToWorld(x, y, z);
    x += d[_x], y += d[_y], z += d[_z];
    _Image->Input()->WorldToImage(x, y, z);
    for (int l = 0; l < _NumberOfComponents; ++l, out += _y /* =X*Y*Z */) {
      (*out) = static_cast<TReal>(_Image->Evaluate(x, y, z, l));
    }
  }
};

// -----------------------------------------------------------------------------
/// Voxel function for composition of output displacement with input deformation
template <class TInterpolator>
struct ApplyInputDeformation : public VoxelFunction
{
  typedef typename TInterpolator::VoxelType TReal;

  const int            _x, _y, _z;
  const TInterpolator *_Image;
  const int            _NumberOfComponents;

  ApplyInputDeformation(const ImageAttributes &domain, const TInterpolator *f)
  :
    _x(0), _y(domain.NumberOfSpatialPoints()), _z(2 * _y),
    _Image(f), _NumberOfComponents(_Image->Input()->T())
  {}

  void operator()(int, int, int, int, const TReal *pos, TReal *out)
  {
    double x = pos[_x], y = pos[_y], z = pos[_z];
    _Image->Input()->WorldToImage(x, y, z);
    for (int l = 0; l < _NumberOfComponents; ++l, out += _y /* =X*Y*Z */) {
      (*out) = static_cast<TReal>(_Image->Evaluate(x, y, z, l));
    }
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


} // namespace ScalingAndSquaringUtils
using namespace ScalingAndSquaringUtils;

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
template <class TReal>
ScalingAndSquaring<TReal>::ScalingAndSquaring()
:
  _InputVelocity(NULL),
  _InputDisplacement(NULL),
  _InputDeformation(NULL),
  _InterimDisplacement(NULL),
  _InterimJacobian(NULL),
  _InterimDetJacobian(NULL),
  _InterimLogJacobian(NULL),
  _InterimJacobianDOFs(NULL),
  _Displacement(NULL),
  _Jacobian(NULL),
  _DetJacobian(NULL),
  _LogJacobian(NULL),
  _JacobianDOFs(NULL),
  _OutputDisplacement(NULL),
  _OutputDeformation(NULL),
  _OutputJacobian(NULL),
  _OutputDetJacobian(NULL),
  _OutputLogJacobian(NULL),
  _OutputJacobianDOFs(NULL),
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
  Delete(_Displacement);
  Delete(_Jacobian);
  Delete(_DetJacobian);
  Delete(_LogJacobian);
  Delete(_JacobianDOFs);
  Delete(_InterimDisplacement);
  Delete(_InterimJacobian);
  Delete(_InterimDetJacobian);
  Delete(_InterimLogJacobian);
  Delete(_InterimJacobianDOFs);
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
    cerr << "ScalingAndSquaring::Initialize: Output attributes due not match those of input displacement field" << endl;
    exit(1);
  }
  if (_InputDeformation && !_InputDeformation->Attributes().EqualInSpace(_OutputAttributes)) {
    cerr << "ScalingAndSquaring::Initialize: Output attributes due not match those of input deformation field" << endl;
    exit(1);
  }
  // Check settings
  if (_NumberOfSteps <= 0 && _NumberOfSquaringSteps <= 0 && _MaxScaledVelocity <= .0) {
    cerr << "ScalingAndSquaring::Initialize: Either number of integration or squaring steps or maximum scaled velocity must be positive!" << endl;
    exit(1);
  }

  // Leave no potential memory leaks
  Clear();

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
  _InterimDisplacement = new ImageType(attr, 3);
  velocity.Evaluate(*_InterimDisplacement);

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

  // Compute derivatives of initial deformation w.r.t. x and/or its (log) determinant
  int jac_mode = 0;
  if (_OutputJacobian || _OutputJacobianDOFs) jac_mode += 1;
  if (_OutputDetJacobian                    ) jac_mode += 2;
  if (_OutputLogJacobian                    ) jac_mode += 4;

  if (jac_mode & 1) _InterimJacobian    = new ImageType(attr, 9);
  if (jac_mode & 2) _InterimDetJacobian = new ImageType(attr);
  if (jac_mode & 4) _InterimLogJacobian = new ImageType(attr);

  switch (jac_mode) {
    case 1: {
      EvaluateJacobian<TReal> eval;
      eval.Initialize(attr, &velocity, scale);
      ParallelForEachVoxel(attr, _InterimJacobian, eval);
    } break;
    case 2: {
      EvaluateDetJacobian<TReal> eval;
      eval.Initialize(attr, &velocity, scale);
      ParallelForEachVoxel(attr, _InterimDetJacobian, eval);
    } break;
    case 3: {
      EvaluateJacobianAndDet<TReal> eval;
      eval.Initialize(attr, &velocity, scale);
      ParallelForEachVoxel(attr, _InterimJacobian, _InterimDetJacobian, eval);
    } break;
    case 4: {
      EvaluateLogJacobian<TReal> eval;
      eval.Initialize(attr, &velocity, scale);
      ParallelForEachVoxel(attr, _InterimLogJacobian, eval);
    } break;
    case 5: {
      EvaluateJacobianAndLog<TReal> eval;
      eval.Initialize(attr, &velocity, scale);
      ParallelForEachVoxel(attr, _InterimJacobian, _InterimLogJacobian, eval);
    } break;
    case 6: {
      EvaluateDetJacobianAndLog<TReal> eval;
      eval.Initialize(attr, &velocity, scale);
      ParallelForEachVoxel(attr, _InterimDetJacobian, _InterimLogJacobian, eval);
    } break;
    case 7: {
      EvaluateJacobianAndDetAndLog<TReal> eval;
      eval.Initialize(attr, &velocity, scale);
      ParallelForEachVoxel(attr, _InterimJacobian, _InterimDetJacobian, _InterimLogJacobian, eval);
    } break;
  }

  // Compute derivatives of initial deformation w.r.t. v
  if (_OutputJacobianDOFs) {
    _InterimJacobianDOFs = new ImageType;
    _InterimJacobianDOFs->Initialize(attr, 9);
    TReal *dxx = _InterimJacobianDOFs->Data(0, 0, 0, 0);
    TReal *dyy = _InterimJacobianDOFs->Data(0, 0, 0, 4);
    TReal *dzz = _InterimJacobianDOFs->Data(0, 0, 0, 8);
    const int n = attr.NumberOfSpatialPoints();
    for (int i = 0; i < n; ++i) {
      dxx[i] = dyy[i] = dzz[i] = scale;
    }
  }

  // Initialize interpolators used during squaring steps
  InterpolationMode imode = _Interpolation;
  ExtrapolationMode emode = (imode == Interpolation_BSpline ||
                             imode == Interpolation_CubicBSpline ||
                             imode == Interpolation_FastCubicBSpline)
                            ? Extrapolation_Mirror : Extrapolation_NN;
  _Displacement = DisplacementField::New(imode, emode, _InterimDisplacement);
  if (_InterimJacobian    ) _Jacobian     = JacobianField::New(imode, emode, _InterimJacobian);
  if (_InterimDetJacobian ) _DetJacobian  = JacobianField::New(imode, emode, _InterimDetJacobian);
  if (_InterimLogJacobian ) _LogJacobian  = JacobianField::New(imode, emode, _InterimLogJacobian);
  if (_InterimJacobianDOFs) _JacobianDOFs = JacobianField::New(imode, emode, _InterimJacobianDOFs);

  // Use output deformation field as temporary output displacement field
  if (!_OutputDisplacement) _OutputDisplacement = _OutputDeformation;

  MIRTK_DEBUG_TIMING(5, "scaling step");
}

// -----------------------------------------------------------------------------
template <class TReal> template <class TInterpolator>
void ScalingAndSquaring<TReal>
::Resample(ImageType *output, TInterpolator *interp)
{
  // Do nothing if output is not requested
  if (!output) return;

  // Attributes of final output images
  const ImageAttributes &attr    = _OutputAttributes;
  const ImageType       *interim = interp->Input();

  // Either apply input deformation field while resampling output...
  if (_InputDeformation) {
    interp->Initialize();
    output->Initialize(attr, interim->T());
    ApplyInputDeformation<TInterpolator> warp(attr, interp);
    ParallelForEachVoxel(attr, _InputDeformation, output, warp);
  // ...or apply input displacement field while resampling output
  } else if (_InputDisplacement) {
    interp->Initialize();
    output->Initialize(attr, interim->T());
    ApplyInputDisplacement<TInterpolator> warp(attr, interp);
    ParallelForEachVoxel(attr, _InputDisplacement, output, warp);
  // ...or resample intermediate output if necessary
  } else if (attr != _InterimAttributes) {
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
    ResampleOutput<TInterpolator> resample(attr, interp);
    ParallelForEachVoxel(attr, output, resample);
  // Otherwise, just copy intermediate output
  } else {
    (*output) = (*interim);
  }
}

// -----------------------------------------------------------------------------
template <class TReal>
void ScalingAndSquaring<TReal>::Finalize()
{
  MIRTK_START_TIMING();

  // Finalize output displacement field and its derivatives
  Resample(_OutputDisplacement, _Displacement);
  Resample(_OutputJacobian,     _Jacobian);
  Resample(_OutputDetJacobian,  _DetJacobian);
  Resample(_OutputLogJacobian,  _LogJacobian);
  Resample(_OutputJacobianDOFs, _JacobianDOFs);

  // Compute output deformation field
  if (_OutputDeformation) {
    if (_OutputDeformation != _OutputDisplacement) {
      _OutputDeformation->Initialize(_OutputAttributes, 3);
    }
    ConvertToDeformation3D<TReal> plusx(_OutputAttributes);
    ParallelForEachVoxel(_OutputAttributes, _OutputDisplacement, _OutputDeformation, plusx);
    if (_OutputDeformation == _OutputDisplacement) _OutputDisplacement = NULL;
  }

  // Multiply output Jacobian by Jacobian of input deformation
  if ((_OutputJacobian || _OutputDetJacobian || _OutputLogJacobian) && (_InputDeformation || _InputDeformation)) {
    // TODO: Multiply output Jacobian by Jacobian of input deformation
    // TODO: Multiply output determinant of Jacobian by determinant of Jacobian of input deformation
    cerr << "WARNING: ScalingAndSquaring::Finalize: Output Jacobian does not include the input deformation!" << endl;
    cerr << "         Be aware that this will be fixed in the future and thus change the output." << endl;
    exit(1);
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
  ImageType *jac3x3 = NULL;
  ImageType *detjac = NULL;
  ImageType *logjac = NULL;
  ImageType *dofjac = NULL;

  if (!disp) disp = new ImageType;
  disp->Initialize(attr, 3);

  if (_InterimJacobian) {
    if (_OutputJacobian) jac3x3 = _OutputJacobian;
    else                 jac3x3 = new ImageType;
    jac3x3->Initialize(attr, 9);
  }
  if (_InterimDetJacobian) {
    if (_OutputDetJacobian) detjac = _OutputDetJacobian;
    else                    detjac = new ImageType;
    detjac->Initialize(attr, 1);
  }
  if (_InterimLogJacobian) {
    if (_OutputLogJacobian) logjac = _OutputLogJacobian;
    else                    logjac = new ImageType;
    logjac->Initialize(attr, 1);
  }
  if (_InterimJacobianDOFs) {
    if (_OutputJacobianDOFs) dofjac = _OutputJacobianDOFs;
    else                     dofjac = new ImageType;
    dofjac->Initialize(attr, 9);
  }

  // Convert scaled displacements to voxel units
  ConvertToVoxelUnits3D<TReal> w2i(attr);
  ParallelForEachVoxel(attr, _InterimDisplacement, _InterimDisplacement, w2i);

  // Do the squaring steps
  int n = _NumberOfSquaringSteps;
  while (n--) {
    // (Re-)initialize interpolators
    _Displacement                   ->Initialize();
    if (_Jacobian    ) _Jacobian    ->Initialize();
    if (_DetJacobian ) _DetJacobian ->Initialize();
    if (_LogJacobian ) _LogJacobian ->Initialize();
    if (_JacobianDOFs) _JacobianDOFs->Initialize();
    // Compute updates
    {
      UpdateDisplacement<DisplacementField> update(_Displacement);
      ParallelForEachVoxel(attr, _InterimDisplacement, disp, update);
    }
    switch (jac_mode) {
      case 1: {
        UpdateJacobian<JacobianField> update;
        update.Initialize(_Jacobian, NULL, NULL);
        ParallelForEachVoxel(attr, _InterimDisplacement, _InterimJacobian, jac3x3, update);
      } break;
      case 2: {
        UpdateDetJacobian<JacobianField> update;
        update.Initialize(NULL, _DetJacobian, NULL);
        ParallelForEachVoxel(attr, _InterimDisplacement, _InterimDetJacobian, detjac, update);
      } break;
      case 3: {
        UpdateJacobianAndDet<JacobianField> update;
        update.Initialize(_Jacobian, _DetJacobian, NULL);
        ParallelForEachVoxel(attr, _InterimDisplacement, _InterimJacobian, _InterimDetJacobian, jac3x3, detjac, update);
      } break;
      case 4: {
        UpdateLogJacobian<JacobianField> update;
        update.Initialize(NULL, NULL, _LogJacobian);
        ParallelForEachVoxel(attr, _InterimDisplacement, _InterimLogJacobian, logjac, update);
      } break;
      case 5: {
        UpdateJacobianAndLog<JacobianField> update;
        update.Initialize(_Jacobian, NULL, _LogJacobian);
        ParallelForEachVoxel(attr, _InterimDisplacement, _InterimJacobian, _InterimLogJacobian, jac3x3, logjac, update);
      } break;
      case 6: {
        UpdateDetJacobianAndLog<JacobianField> update;
        update.Initialize(NULL, _DetJacobian, _LogJacobian);
        ParallelForEachVoxel(attr, _InterimDisplacement, _InterimDetJacobian, _InterimLogJacobian, detjac, logjac, update);
      } break;
      case 7: {
        UpdateJacobianAndDetAndLog<JacobianField> update;
        update.Initialize(_Jacobian, _DetJacobian, _LogJacobian);
        ParallelForEachVoxel(attr, _InterimDisplacement, _InterimJacobian, _InterimDetJacobian, _InterimLogJacobian, jac3x3, detjac, logjac, update);
      } break;
    }
    if (_InterimJacobianDOFs) {
      UpdateJacobianDOFs<JacobianField> update(_Jacobian, _JacobianDOFs);
      ParallelForEachVoxel(attr, _InterimDisplacement, _InterimJacobianDOFs, dofjac, update);
    }
    // Update intermediate output images
    _InterimDisplacement                          ->CopyFrom(disp  ->Data());
    if (_InterimJacobian    ) _InterimJacobian    ->CopyFrom(jac3x3->Data());
    if (_InterimDetJacobian ) _InterimDetJacobian ->CopyFrom(detjac->Data());
    if (_InterimLogJacobian ) _InterimLogJacobian ->CopyFrom(logjac->Data());
    if (_InterimJacobianDOFs) _InterimJacobianDOFs->CopyFrom(dofjac->Data());
  }

  // Convert final displacements back to world units if output requested
  if (_OutputDisplacement) {
    ConvertToWorldUnits3D<TReal> i2w(attr);
    ParallelForEachVoxel(attr, _InterimDisplacement, _InterimDisplacement, i2w);
  }

  // Do the final cleaning up
  if (disp   != _OutputDisplacement) Delete(disp);
  if (jac3x3 != _OutputJacobian    ) Delete(jac3x3);
  if (detjac != _OutputDetJacobian ) Delete(detjac);
  if (logjac != _OutputLogJacobian ) Delete(logjac);
  if (dofjac != _OutputJacobianDOFs) Delete(dofjac);

  MIRTK_DEBUG_TIMING(5, "squaring steps"
                           " (d="    << (_OutputDisplacement ? "on" : "off")
                        << ", J="    << (_OutputJacobian     ? "on" : "off")
                        << ", detJ=" << (_OutputDetJacobian  ? "on" : "off")
                        << ", logJ=" << (_OutputLogJacobian  ? "on" : "off")
                        << ", dv="   << (_OutputJacobianDOFs ? "on" : "off") << ")");

  this->Finalize();
}

// =============================================================================
// Explicit template instantiations
// =============================================================================

template class ScalingAndSquaring<float>;
template class ScalingAndSquaring<double>;


} // namespace mirtk
