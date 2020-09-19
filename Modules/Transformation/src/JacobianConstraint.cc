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

#include "mirtk/JacobianConstraint.h"

#include "mirtk/Memory.h"
#include "mirtk/Matrix.h"
#include "mirtk/Parallel.h"
#include "mirtk/FreeFormTransformation.h"
#include "mirtk/BSplineFreeFormTransformation3D.h"
#include "mirtk/BSplineFreeFormTransformationSV.h"
#include "mirtk/MultiLevelTransformation.h"
#include "mirtk/GaussianBlurring.h"
#include "mirtk/Resampling.h"
#include "mirtk/LinearInterpolateImageFunction.h"
#include "mirtk/Profiling.h"

#include "mirtk/CommonExport.h" 

// The sum of the individual contributions to the gradient w.r.t. a control point
// has been normalized by dividing by the number of terms. It is disabled with this flag.
#define _USE_INNER_GRADIENT_SUM_NORM 0

// Whether to divide gradient vectors by number of points in (sub-)domain.
// This was done in the IRTK nreg2 implementation, but the gradient magnitude is too
// small in this case and would need a high weight (>10). The volume preservation and
// topology preservation penalties for the classic FFD model seemed to work better
// without dividing the gradient vectors by the number of sub-domain points.
// Instead, the sum of gradient contributions is divided by the number of sub-domain
// points within the local support of the respective control point.
#define _DIVIDE_GRADIENT_BY_NPOINTS 1


namespace mirtk {


// global "debug" flag (cf. mirtk/Options.h)
MIRTK_Common_EXPORT extern int debug;


// =============================================================================
// Auxiliaries
// =============================================================================

namespace {


using WeightImage = JacobianConstraint::WeightImage;


// -----------------------------------------------------------------------------
/// Subdivide lattice by inserting points in-between previous lattice points
void SubdivDomain(ImageAttributes &domain)
{
  if (domain._x > 1) {
    domain._x  *= 2;
    domain._dx /= 2.;
  }
  if (domain._y > 1) {
    domain._y  *= 2;
    domain._dy /= 2.;
  }
  if (domain._z > 1) {
    domain._z  *= 2;
    domain._dz /= 2.;
  }
  if (domain._t > 1) {
    domain._t  *= 2;
    domain._dt /= 2.;
  }
}

// -----------------------------------------------------------------------------
/// Add one voxel to each dimension, while center remains fixed.
/// This results in an evaluation of the penalty in-between lattice points
void ShiftDomain(ImageAttributes &domain)
{
  if (domain._x > 1) domain._x += 1;
  if (domain._y > 1) domain._y += 1;
  if (domain._z > 1) domain._z += 1;
  if (domain._t > 1) domain._t += 1;
}

// -----------------------------------------------------------------------------
/// Check if world point is within local support of an active control point
inline bool IsActiveWorld(const FreeFormTransformation *ffd, double x, double y, double z, double t)
{
  ffd->WorldToLattice(x, y, z);
  t = ffd->TimeToLattice(t);

  const double r = static_cast<double>(ffd->KernelRadius());

  int i1 = max(0, ifloor(x - r));
  int j1 = max(0, ifloor(y - r));
  int k1 = max(0, ifloor(z - r));
  int l1 = max(0, ifloor(t - r));

  int i2 = min(iceil(x + r), ffd->X() - 1);
  int j2 = min(iceil(y + r), ffd->Y() - 1);
  int k2 = min(iceil(z + r), ffd->Z() - 1);
  int l2 = min(iceil(t + r), ffd->T() - 1);

  for (int l = l1; l <= l2; ++l)
  for (int k = k1; k <= k2; ++k)
  for (int j = j1; j <= j2; ++j)
  for (int i = i1; i <= i2; ++i) {
    if (ffd->IsActive(i, j, k, l)) return true;
  }
  return false;
}

// -----------------------------------------------------------------------------
/// Check if lattice point is within local support of an active control point
inline bool IsActiveLattice(const FreeFormTransformation3D *ffd, double x, double y, double z)
{
  const double r = static_cast<double>(ffd->KernelRadius());

  int i1 = max(0, ifloor(x - r));
  int j1 = max(0, ifloor(y - r));
  int k1 = max(0, ifloor(z - r));

  int i2 = min(iceil(x + r), ffd->X() - 1);
  int j2 = min(iceil(y + r), ffd->Y() - 1);
  int k2 = min(iceil(z + r), ffd->Z() - 1);

  for (int k = k1; k <= k2; ++k)
  for (int j = j1; j <= j2; ++j)
  for (int i = i1; i <= i2; ++i) {
    if (ffd->IsActive(i, j, k)) return true;
  }
  return false;
}

// -----------------------------------------------------------------------------
/// Check if lattice point is within local support of an active control point
inline bool IsActiveBSplineLatticePoint(const BSplineFreeFormTransformation3D *ffd, int ci, int cj, int ck)
{
  int i1 = max(0, ci - 1);
  int j1 = max(0, cj - 1);
  int k1 = max(0, ck - 1);

  int i2 = min(ci + 1, ffd->X() - 1);
  int j2 = min(cj + 1, ffd->Y() - 1);
  int k2 = min(ck + 1, ffd->Z() - 1);

  for (int k = k1; k <= k2; ++k)
  for (int j = j1; j <= j2; ++j)
  for (int i = i1; i <= i2; ++i) {
    if (ffd->IsActive(i, j, k)) return true;
  }
  return false;
}

// -----------------------------------------------------------------------------
/// Multiply 2D row vector from right with 2x2 matrix
inline void MultiplyRight2x2(double &du, double &dv, const Matrix &m)
{
  double dx = du * m(0, 0) + dv * m(1, 0);
  double dy = du * m(0, 1) + dv * m(1, 1);
  du = dx, dv = dy;
}

// -----------------------------------------------------------------------------
/// Multiply 2x2 Jacobian matrix from right with reorientation (and scaling) matrix
inline void MultiplyRight2x2(Matrix &jac, const Matrix &m)
{
  MultiplyRight2x2(jac(0, 0), jac(0, 1), m);
  MultiplyRight2x2(jac(1, 0), jac(1, 1), m);
}

// -----------------------------------------------------------------------------
/// Multiply 3D row vector from right with 3x3 matrix
inline void MultiplyRight3x3(double &du, double &dv, double &dw, const Matrix &m)
{
  double dx = du * m(0, 0) + dv * m(1, 0) + dw * m(2, 0);
  double dy = du * m(0, 1) + dv * m(1, 1) + dw * m(2, 1);
  double dz = du * m(0, 2) + dv * m(1, 2) + dw * m(2, 2);
  du = dx, dv = dy, dw = dz;
}

// -----------------------------------------------------------------------------
/// Multiply 3x3 Jacobian matrix from right with reorientation (and scaling) matrix
inline void MultiplyRight3x3(Matrix &jac, const Matrix &m)
{
  MultiplyRight3x3(jac(0, 0), jac(0, 1), jac(0, 2), m);
  MultiplyRight3x3(jac(1, 0), jac(1, 1), jac(1, 2), m);
  MultiplyRight3x3(jac(2, 0), jac(2, 1), jac(2, 2), m);
}


// -----------------------------------------------------------------------------
/// Body of generic Update function
struct EvaluateJacobian_AnyFFD
{
private:

  const JacobianConstraint     *_Constraint;
  const FreeFormTransformation *_FFD;
  const ImageAttributes        *_Domain;
  const WeightImage            *_Mask;
  Matrix                       *_Jacobian;

public:

  void operator ()(const blocked_range3d<int> &re) const
  {
    int idx;
    double x, y, z, t;
    if (_Constraint->ConstrainPassiveDoFs()) {
      if (_Constraint->ConstrainParameterization()) {
        for (int l = 0; l < _Domain->T(); ++l) {
          t = _Domain->LatticeToTime(l);
          for (int k = re.pages().begin(); k != re.pages().end(); ++k)
          for (int j = re.rows ().begin(); j != re.rows ().end(); ++j)
          for (int i = re.cols ().begin(); i != re.cols ().end(); ++i) {
            idx = _Domain->LatticeToIndex(i, j, k, l);
            if (_Mask && IsZero(_Mask->Get(i, j, k))) {
              _Jacobian[idx].Clear();
              continue;
            }
            x = i, y = j, z = k;
            _Domain->LatticeToWorld(x, y, z);
            _FFD->FFDJacobianWorld(_Jacobian[idx], x, y, z, t);
          }
        }
      } else {
        for (int l = 0; l < _Domain->T(); ++l) {
          t = _Domain->LatticeToTime(l);
          for (int k = re.pages().begin(); k != re.pages().end(); ++k)
          for (int j = re.rows ().begin(); j != re.rows ().end(); ++j)
          for (int i = re.cols ().begin(); i != re.cols ().end(); ++i) {
            idx = _Domain->LatticeToIndex(i, j, k, l);
            if (_Mask && IsZero(_Mask->Get(i, j, k))) {
              _Jacobian[idx].Clear();
              continue;
            }
            x = i, y = j, z = k;
            _Domain->LatticeToWorld(x, y, z);
            _FFD->LocalJacobian(_Jacobian[idx], x, y, z, t);
          }
        }
      }
    } else {
      if (_Constraint->ConstrainParameterization()) {
        for (int l = 0; l < _Domain->T(); ++l) {
          t = _Domain->LatticeToTime(l);
          for (int k = re.pages().begin(); k != re.pages().end(); ++k)
          for (int j = re.rows ().begin(); j != re.rows ().end(); ++j)
          for (int i = re.cols ().begin(); i != re.cols ().end(); ++i) {
            idx = _Domain->LatticeToIndex(i, j, k, l);
            if (_Mask && IsZero(_Mask->Get(i, j, k))) {
              _Jacobian[idx].Clear();
              continue;
            }
            x = i, y = j, z = k;
            _Domain->LatticeToWorld(x, y, z);
            if (IsActiveWorld(_FFD, x, y, z, t)) {
              _FFD->FFDJacobianWorld(_Jacobian[idx], x, y, z, t);
            } else {
              _Jacobian[idx].Clear();
            }
          }
        }
      } else {
        for (int l = 0; l < _Domain->T(); ++l) {
          t = _Domain->LatticeToTime(l);
          for (int k = re.pages().begin(); k != re.pages().end(); ++k)
          for (int j = re.rows ().begin(); j != re.rows ().end(); ++j)
          for (int i = re.cols ().begin(); i != re.cols ().end(); ++i) {
            idx = _Domain->LatticeToIndex(i, j, k, l);
            if (_Mask && IsZero(_Mask->Get(i, j, k))) {
              _Jacobian[idx].Clear();
              continue;
            }
            x = i, y = j, z = k;
            _Domain->LatticeToWorld(x, y, z);
            if (IsActiveWorld(_FFD, x, y, z, t)) {
              _FFD->LocalJacobian(_Jacobian[idx], x, y, z, t);
            } else {
              _Jacobian[idx].Clear();
            }
          }
        }
      }
    }
  }

  static void Run(const JacobianConstraint     *constraint,
                  const FreeFormTransformation *ffd,
                  const ImageAttributes        &domain,
                  const WeightImage            &mask,
                  Matrix                       *jac)
  {
    EvaluateJacobian_AnyFFD body;
    body._Constraint = constraint;
    body._FFD        = ffd;
    body._Domain     = &domain;
    body._Mask       = mask.IsEmpty() ? nullptr : &mask;
    body._Jacobian   = jac;
    blocked_range3d<int> range(0, domain.Z(), 0, domain.Y(), 0, domain.X());
    parallel_for(range, body);
  }
};


// -----------------------------------------------------------------------------
/// Reorient 2x2 Jacobian matrices computed by EvaluateJacobian_AnyFFD
struct ReorientJacobian2D
{
private:

  Matrix       *_Jacobian;
  const Matrix *_Orient;

public:

  void operator ()(const blocked_range<int> &re) const
  {
    for (int idx = re.begin(); idx != re.end(); ++idx) {
      Matrix &jac = _Jacobian[idx];
      jac(0, 0) -= 1.;
      jac(1, 1) -= 1.;
      MultiplyRight2x2(jac, *_Orient);
      jac(0, 0) += 1.;
      jac(1, 1) += 1.;
    }
  }

  static void Run(int n, Matrix *jac, const Matrix *orient)
  {
    ReorientJacobian2D body;
    body._Jacobian = jac;
    body._Orient   = orient;
    parallel_for(blocked_range<int>(0, n), body);
  }
};


// -----------------------------------------------------------------------------
/// Reorient 3x3 Jacobian matrices computed by EvaluateJacobian_AnyFFD
struct ReorientJacobian3D
{
private:

  Matrix       *_Jacobian;
  const Matrix *_Orient;

public:

  void operator ()(const blocked_range<int> &re) const
  {
    for (int idx = re.begin(); idx != re.end(); ++idx) {
      Matrix &jac = _Jacobian[idx];
      jac(0, 0) -= 1.;
      jac(1, 1) -= 1.;
      jac(2, 2) -= 1.;
      MultiplyRight3x3(jac, *_Orient);
      jac(0, 0) += 1.;
      jac(1, 1) += 1.;
      jac(2, 2) += 1.;
    }
  }

  static void Run(int n, Matrix *jac, const Matrix *orient)
  {
    ReorientJacobian3D body;
    body._Jacobian = jac;
    body._Orient   = orient;
    parallel_for(blocked_range<int>(0, n), body);
  }
};


// -----------------------------------------------------------------------------
/// Body of Update function specialized for BSplineFreeFormTransformation3D
/// for evaluation of Jacobian determinants of spline function at lattice points
/// of aribtrary image domain in the vicinity of active control points.
struct EvaluateJacobian_Domain_ActiveCPs_BSplineFFD_3D
{
private:

  const BSplineFreeFormTransformation3D *_FFD;
  const ImageAttributes                 *_Domain;
  const WeightImage                     *_Mask;
  const Matrix                          *_Orient;
  Matrix                                *_Jacobian;

public:

  void operator ()(const blocked_range3d<int> &re) const
  {
    int idx;
    double x, y, z;
    for (int k = re.pages().begin(); k != re.pages().end(); ++k)
    for (int j = re.rows ().begin(); j != re.rows ().end(); ++j)
    for (int i = re.cols ().begin(); i != re.cols ().end(); ++i) {
      idx = _Domain->LatticeToIndex(i, j, k);
      if (_Mask && IsZero(_Mask->Get(i, j, k))) {
        _Jacobian[idx].Clear();
        continue;
      }
      x = i, y = j, z = k;
      _Domain->LatticeToWorld(x, y, z);
      _FFD->WorldToLattice(x, y, z);
      if (IsActiveLattice(_FFD, x, y, z)) {
        Matrix &jac = _Jacobian[idx];
        _FFD->EvaluateJacobian(jac, x, y, z);
        if (_Orient) {
          MultiplyRight3x3(jac, *_Orient);
        }
        jac(0, 0) += 1.;
        jac(1, 1) += 1.;
        jac(2, 2) += 1.;
      } else {
        _Jacobian[idx].Clear();
      }
    }
  }

  static void Run(const BSplineFreeFormTransformation3D *ffd,
                  const ImageAttributes &domain,
                  const WeightImage &mask,
                  Matrix *jac,
                  const Matrix *orient = nullptr)
  {
    EvaluateJacobian_Domain_ActiveCPs_BSplineFFD_3D body;
    body._FFD      = ffd;
    body._Domain   = &domain;
    body._Mask     = mask.IsEmpty() ? nullptr : &mask;
    body._Orient   = orient;
    body._Jacobian = jac;
    blocked_range3d<int> range(0, domain.Z(), 0, domain.Y(), 0, domain.X());
    parallel_for(range, body);
  }
};


// -----------------------------------------------------------------------------
/// Body of Update function specialized for BSplineFreeFormTransformationSV
/// for evaluation of Jacobian determinants of spline function at lattice points
/// of aribtrary image domain in the vicinity of active control points.
struct EvaluateJacobian_Domain_ActiveCPs_BSplineFFD_SV
{
private:

  const BSplineFreeFormTransformationSV *_FFD;
  const ImageAttributes                 *_Domain;
  const WeightImage                     *_Mask;
  const Matrix                          *_Orient;
  Matrix                                *_Jacobian;

public:

  void operator ()(const blocked_range3d<int> &re) const
  {
    int idx1, idx2;
    double x, y, z;
    for (int k = re.pages().begin(); k != re.pages().end(); ++k)
    for (int j = re.rows ().begin(); j != re.rows ().end(); ++j)
    for (int i = re.cols ().begin(); i != re.cols ().end(); ++i) {
      idx1 = _Domain->LatticeToIndex(i, j, k, 0);
      idx2 = _Domain->LatticeToIndex(i, j, k, 1);
      if (_Mask && IsZero(_Mask->Get(i, j, k))) {
        _Jacobian[idx1].Clear();
        _Jacobian[idx2].Clear();
        continue;
      }
      x = i, y = j, z = k;
      _Domain->LatticeToWorld(x, y, z);
      _FFD->WorldToLattice(x, y, z);
      if (IsActiveLattice(_FFD, x, y, z)) {
        Matrix &jac1 = _Jacobian[idx1];
        Matrix &jac2 = _Jacobian[idx2];
        _FFD->EvaluateJacobian(jac1, x, y, z);
        if (_Orient) {
          MultiplyRight3x3(jac1, *_Orient);
        }
        jac2 = jac1, jac2 *= -1.;
        jac1(0, 0) += 1.;
        jac1(1, 1) += 1.;
        jac1(2, 2) += 1.;
        jac2(0, 0) += 1.;
        jac2(1, 1) += 1.;
        jac2(2, 2) += 1.;
      } else {
        _Jacobian[idx1].Clear();
        _Jacobian[idx2].Clear();
      }
    }
  }

  static void Run(const BSplineFreeFormTransformationSV *ffd,
                  const ImageAttributes &domain,
                  const WeightImage &mask,
                  Matrix *jac,
                  const Matrix *orient = nullptr)
  {
    EvaluateJacobian_Domain_ActiveCPs_BSplineFFD_SV body;
    body._FFD      = ffd;
    body._Domain   = &domain;
    body._Mask     = mask.IsEmpty() ? nullptr : &mask;
    body._Orient   = orient;
    body._Jacobian = jac;
    blocked_range3d<int> range(0, domain.Z(), 0, domain.Y(), 0, domain.X());
    parallel_for(range, body);
  }
};


// -----------------------------------------------------------------------------
/// Body of Update function specialized for BSplineFreeFormTransformation3D
/// for evaluation of Jacobian determinants of spline function at lattice points
/// of control point grid in the vicinity of active control points.
struct EvaluateJacobian_Lattice_ActiveCPs_BSplineFFD_3D
{
private:

  const BSplineFreeFormTransformation3D *_FFD;
  const WeightImage                     *_Mask;
  const Matrix                          *_Orient;
  Matrix                                *_Jacobian;

public:

  void operator ()(const blocked_range3d<int> &re) const
  {
    int idx;
    for (int k = re.pages().begin(); k != re.pages().end(); ++k)
    for (int j = re.rows ().begin(); j != re.rows ().end(); ++j)
    for (int i = re.cols ().begin(); i != re.cols ().end(); ++i) {
      idx = _FFD->LatticeToIndex(i, j, k);
      bool skip = _Mask && IsZero(_Mask->Get(i, j, k));
      if (!skip && IsActiveBSplineLatticePoint(_FFD, i, j, k)) {
        Matrix &jac = _Jacobian[idx];
        _FFD->EvaluateJacobian(jac, i, j, k);
        if (_Orient) {
          MultiplyRight3x3(jac, *_Orient);
        }
        jac(0, 0) += 1.;
        jac(1, 1) += 1.;
        jac(2, 2) += 1.;
      } else {
        _Jacobian[idx].Clear();
      }
    }
  }

  static void Run(const BSplineFreeFormTransformation3D *ffd, const WeightImage &mask,
                  Matrix *jac, const Matrix *orient = nullptr)
  {
    EvaluateJacobian_Lattice_ActiveCPs_BSplineFFD_3D body;
    body._FFD      = ffd;
    body._Mask     = mask.IsEmpty() ? nullptr : &mask;
    body._Orient   = orient;
    body._Jacobian = jac;
    blocked_range3d<int> range(0, ffd->Z(), 0, ffd->Y(), 0, ffd->X());
    parallel_for(range, body);
  }
};


// -----------------------------------------------------------------------------
/// Body of Update function specialized for BSplineFreeFormTransformationSV
/// for evaluation of Jacobian determinants of spline function at lattice points
/// of control point grid in the vicinity of active control points.
struct EvaluateJacobian_Lattice_ActiveCPs_BSplineFFD_SV
{
private:

  const BSplineFreeFormTransformationSV *_FFD;
  const WeightImage                     *_Mask;
  const Matrix                          *_Orient;
  Matrix                                *_Jacobian;

public:

  void operator ()(const blocked_range3d<int> &re) const
  {
    int idx1, idx2;
    for (int k = re.pages().begin(); k != re.pages().end(); ++k)
    for (int j = re.rows ().begin(); j != re.rows ().end(); ++j)
    for (int i = re.cols ().begin(); i != re.cols ().end(); ++i) {
      idx1 = _FFD->LatticeToIndex(i, j, k, 0);
      idx2 = _FFD->LatticeToIndex(i, j, k, 1);
      bool skip = _Mask && IsZero(_Mask->Get(i, j, k));
      if (!skip && IsActiveBSplineLatticePoint(_FFD, i, j, k)) {
        Matrix &jac1 = _Jacobian[idx1];
        Matrix &jac2 = _Jacobian[idx2];
        _FFD->EvaluateJacobian(jac1, i, j, k);
        if (_Orient) {
          MultiplyRight3x3(jac1, *_Orient);
        }
        jac2 = jac1, jac2 *= -1.;
        jac1(0, 0) += 1.;
        jac1(1, 1) += 1.;
        jac1(2, 2) += 1.;
        jac2(0, 0) += 1.;
        jac2(1, 1) += 1.;
        jac2(2, 2) += 1.;
      } else {
        _Jacobian[idx1].Clear();
        _Jacobian[idx2].Clear();
      }
    }
  }

  static void Run(const BSplineFreeFormTransformationSV *ffd, const WeightImage &mask,
                  Matrix *jac, const Matrix *orient = nullptr)
  {
    EvaluateJacobian_Lattice_ActiveCPs_BSplineFFD_SV body;
    body._FFD      = ffd;
    body._Mask     = mask.IsEmpty() ? nullptr : &mask;
    body._Orient   = orient;
    body._Jacobian = jac;
    blocked_range3d<int> range(0, ffd->Z(), 0, ffd->Y(), 0, ffd->X());
    parallel_for(range, body);
  }
};


// -----------------------------------------------------------------------------
/// Body of Update function specialized for BSplineFreeFormTransformation3D
/// for evaluation of Jacobian determinants of spline function at lattice points
/// of aribtrary image domain for both active and passive control points.
struct EvaluateJacobian_Domain_InclPassiveCPs_BSplineFFD_3D
{
private:

  const BSplineFreeFormTransformation3D *_FFD;
  const ImageAttributes                 *_Domain;
  const WeightImage                     *_Mask;
  const Matrix                          *_Orient;
  Matrix                                *_Jacobian;

public:

  void operator ()(const blocked_range3d<int> &re) const
  {
    int idx;
    double x, y, z;
    for (int k = re.pages().begin(); k != re.pages().end(); ++k)
    for (int j = re.rows ().begin(); j != re.rows ().end(); ++j)
    for (int i = re.cols ().begin(); i != re.cols ().end(); ++i) {
      idx = _Domain->LatticeToIndex(i, j, k);
      Matrix &jac = _Jacobian[idx];
      if (_Mask && IsZero(_Mask->Get(i, j, k))) {
        jac.Clear();
        continue;
      }
      x = i, y = j, z = k;
      _Domain->LatticeToWorld(x, y, z);
      _FFD->WorldToLattice(x, y, z);
      _FFD->EvaluateJacobian(jac, x, y, z);
      if (_Orient) {
        MultiplyRight3x3(jac, *_Orient);
      }
      jac(0, 0) += 1.;
      jac(1, 1) += 1.;
      jac(2, 2) += 1.;
    }
  }

  static void Run(const BSplineFreeFormTransformation3D *ffd,
                  const ImageAttributes &domain, const WeightImage &mask,
                  Matrix *jac, const Matrix *orient = nullptr)
  {
    EvaluateJacobian_Domain_InclPassiveCPs_BSplineFFD_3D body;
    body._FFD      = ffd;
    body._Domain   = &domain;
    body._Mask     = mask.IsEmpty() ? nullptr : &mask;
    body._Orient   = orient;
    body._Jacobian = jac;
    blocked_range3d<int> range(0, domain.Z(), 0, domain.Y(), 0, domain.X());
    parallel_for(range, body);
  }
};


// -----------------------------------------------------------------------------
/// Body of Update function specialized for BSplineFreeFormTransformationSV
/// for evaluation of Jacobian determinants of spline function at lattice points
/// of aribtrary image domain for both active and passive control points.
struct EvaluateJacobian_Domain_InclPassiveCPs_BSplineFFD_SV
{
private:

  const BSplineFreeFormTransformationSV *_FFD;
  const ImageAttributes                 *_Domain;
  const WeightImage                     *_Mask;
  const Matrix                          *_Orient;
  Matrix                                *_Jacobian;

public:

  void operator ()(const blocked_range3d<int> &re) const
  {
    int idx1, idx2;
    double x, y, z;
    for (int k = re.pages().begin(); k != re.pages().end(); ++k)
    for (int j = re.rows ().begin(); j != re.rows ().end(); ++j)
    for (int i = re.cols ().begin(); i != re.cols ().end(); ++i) {
      idx1 = _Domain->LatticeToIndex(i, j, k, 0);
      idx2 = _Domain->LatticeToIndex(i, j, k, 1);
      Matrix &jac1 = _Jacobian[idx1];
      Matrix &jac2 = _Jacobian[idx2];
      if (_Mask && IsZero(_Mask->Get(i, j, k))) {
        jac1.Clear();
        jac2.Clear();
        continue;
      }
      x = i, y = j, z = k;
      _Domain->LatticeToWorld(x, y, z);
      _FFD->WorldToLattice(x, y, z);
      _FFD->EvaluateJacobian(jac1, x, y, z);
      if (_Orient) {
        MultiplyRight3x3(jac1, *_Orient);
      }
      jac2 = jac1, jac2 *= -1.;
      jac1(0, 0) += 1.;
      jac1(1, 1) += 1.;
      jac1(2, 2) += 1.;
      jac2(0, 0) += 1.;
      jac2(1, 1) += 1.;
      jac2(2, 2) += 1.;
    }
  }

  static void Run(const BSplineFreeFormTransformationSV *svffd,
                  const ImageAttributes &domain, const WeightImage &mask,
                  Matrix *jac, const Matrix *orient = nullptr)
  {
    EvaluateJacobian_Domain_InclPassiveCPs_BSplineFFD_SV body;
    body._FFD      = svffd;
    body._Domain   = &domain;
    body._Mask     = mask.IsEmpty() ? nullptr : &mask;
    body._Orient   = orient;
    body._Jacobian = jac;
    blocked_range3d<int> range(0, domain.Z(), 0, domain.Y(), 0, domain.X());
    parallel_for(range, body);
  }
};


// -----------------------------------------------------------------------------
/// Body of Update function specialized for BSplineFreeFormTransformation3D
/// for evaluation of Jacobian determinants of spline function at lattice points
/// of control point grid for both passive and active control points.
struct EvaluateJacobian_Lattice_InclPassiveCPs_BSplineFFD_3D
{
private:

  const BSplineFreeFormTransformation3D *_FFD;
  const WeightImage                      *_Mask;
  const Matrix                          *_Orient;
  Matrix                                *_Jacobian;

public:

  void operator ()(const blocked_range3d<int> &re) const
  {
    int idx;
    for (int k = re.pages().begin(); k != re.pages().end(); ++k)
    for (int j = re.rows ().begin(); j != re.rows ().end(); ++j)
    for (int i = re.cols ().begin(); i != re.cols ().end(); ++i) {
      idx = _FFD->LatticeToIndex(i, j, k);
      Matrix &jac = _Jacobian[idx];
      if (_Mask && IsZero(_Mask->Get(i, j, k))) {
        jac.Clear();
        continue;
      }
      _FFD->EvaluateJacobian(jac, i, j, k);
      if (_Orient) {
        MultiplyRight3x3(jac, *_Orient);
      }
      jac(0, 0) += 1.;
      jac(1, 1) += 1.;
      jac(2, 2) += 1.;
    }
  }

  static void Run(const BSplineFreeFormTransformation3D *ffd, const WeightImage &mask,
                  Matrix *jac, const Matrix *orient = nullptr)
  {
    EvaluateJacobian_Lattice_InclPassiveCPs_BSplineFFD_3D body;
    body._FFD      = ffd;
    body._Mask     = mask.IsEmpty() ? nullptr : &mask;
    body._Orient   = orient;
    body._Jacobian = jac;
    blocked_range3d<int> range(0, ffd->Z(), 0, ffd->Y(), 0, ffd->X());
    parallel_for(range, body);
  }
};


// -----------------------------------------------------------------------------
/// Body of Update function specialized for BSplineFreeFormTransformationSV
/// for evaluation of Jacobian determinants of spline function at lattice points
/// of control point grid for both passive and active control points.
struct EvaluateJacobian_Lattice_InclPassiveCPs_BSplineFFD_SV
{
private:

  const BSplineFreeFormTransformationSV *_FFD;
  const WeightImage                     *_Mask;
  const Matrix                          *_Orient;
  Matrix                                *_Jacobian;

public:

  void operator ()(const blocked_range3d<int> &re) const
  {
    int idx1, idx2;
    for (int k = re.pages().begin(); k != re.pages().end(); ++k)
    for (int j = re.rows ().begin(); j != re.rows ().end(); ++j)
    for (int i = re.cols ().begin(); i != re.cols ().end(); ++i) {
      idx1 = _FFD->LatticeToIndex(i, j, k, 0);
      idx2 = _FFD->LatticeToIndex(i, j, k, 1);
      Matrix &jac1 = _Jacobian[idx1];
      Matrix &jac2 = _Jacobian[idx2];
      if (_Mask && IsZero(_Mask->Get(i, j, k))) {
        jac1.Clear();
        jac2.Clear();
        continue;
      }
      _FFD->EvaluateJacobian(jac1, i, j, k);
      if (_Orient) {
        MultiplyRight3x3(jac1, *_Orient);
      }
      jac2 = jac1, jac2 *= -1.;
      jac1(0, 0) += 1.;
      jac1(1, 1) += 1.;
      jac1(2, 2) += 1.;
      jac2(0, 0) += 1.;
      jac2(1, 1) += 1.;
      jac2(2, 2) += 1.;
    }
  }

  static void Run(const BSplineFreeFormTransformationSV *ffd, const WeightImage &mask,
                  Matrix *jac, const Matrix *orient = nullptr)
  {
    EvaluateJacobian_Lattice_InclPassiveCPs_BSplineFFD_SV body;
    body._FFD      = ffd;
    body._Mask     = mask.IsEmpty() ? nullptr : &mask;
    body._Orient   = orient;
    body._Jacobian = jac;
    blocked_range3d<int> range(0, ffd->Z(), 0, ffd->Y(), 0, ffd->X());
    parallel_for(range, body);
  }
};


// -----------------------------------------------------------------------------
/// Evaluate Jacobian matrices for given FFD
int EvaluateJacobian(const JacobianConstraint *constraint,
                     const FreeFormTransformation *ffd,
                     const ImageAttributes &domain,
                     const WeightImage &mask,
                     Matrix *jac)
{
  const int npts = domain.NumberOfPoints();
  auto bffd  = dynamic_cast<const BSplineFreeFormTransformation3D *>(ffd);
  auto svffd = dynamic_cast<const BSplineFreeFormTransformationSV *>(ffd);
  if (!constraint->ConstrainParameterization()) {
    if (svffd) {
      bffd  = nullptr;
      svffd = nullptr;
    }
  }
  if (!constraint->Symmetric()) {
    svffd = nullptr;
  }
  if (svffd) {
    UniquePtr<Matrix> orient;
    if (constraint->WithRespectToWorld()) {
      if (constraint->UseLatticeSpacing()) {
        orient.reset(new Matrix(svffd->Attributes().GetWorldToImageMatrix()));
      } else {
        orient.reset(new Matrix(svffd->Attributes().GetWorldToImageOrientation()));
      }
    }
    if (constraint->SubDomain() == JacobianConstraint::SD_Lattice) {
      if (constraint->ConstrainPassiveDoFs()) {
        typedef EvaluateJacobian_Lattice_InclPassiveCPs_BSplineFFD_SV Update;
        Update::Run(svffd, mask, jac, orient.get());
      } else {
        typedef EvaluateJacobian_Lattice_ActiveCPs_BSplineFFD_SV Update;
        Update::Run(svffd, mask, jac, orient.get());
      }
    } else {
      if (constraint->ConstrainPassiveDoFs()) {
        typedef EvaluateJacobian_Domain_InclPassiveCPs_BSplineFFD_SV Update;
        Update::Run(svffd, domain, mask, jac, orient.get());
      } else {
        typedef EvaluateJacobian_Domain_ActiveCPs_BSplineFFD_SV Update;
        Update::Run(svffd, domain, mask, jac, orient.get());
      }
    }
    return 2 * npts;
  }
  if (bffd && bffd->Z() > 1) {
    UniquePtr<Matrix> orient;
    if (constraint->WithRespectToWorld()) {
      if (constraint->UseLatticeSpacing()) {
        orient.reset(new Matrix(bffd->Attributes().GetWorldToImageMatrix()));
      } else {
        orient.reset(new Matrix(bffd->Attributes().GetWorldToImageOrientation()));
      }
    }
    if (constraint->SubDomain() == JacobianConstraint::SD_Lattice) {
      if (constraint->ConstrainPassiveDoFs()) {
        typedef EvaluateJacobian_Lattice_InclPassiveCPs_BSplineFFD_3D Update;
        Update::Run(bffd, mask, jac, orient.get());
      } else {
        typedef EvaluateJacobian_Lattice_ActiveCPs_BSplineFFD_3D Update;
        Update::Run(bffd, mask, jac, orient.get());
      }
    } else {
      if (constraint->ConstrainPassiveDoFs()) {
        typedef EvaluateJacobian_Domain_InclPassiveCPs_BSplineFFD_3D Update;
        Update::Run(bffd, domain, mask, jac, orient.get());
      } else {
        typedef EvaluateJacobian_Domain_ActiveCPs_BSplineFFD_3D Update;
        Update::Run(bffd, domain, mask, jac, orient.get());
      }
    }
    return npts;
  }
  EvaluateJacobian_AnyFFD::Run(constraint, ffd, domain, mask, jac);
  if (!constraint->WithRespectToWorld() || !constraint->UseLatticeSpacing()) {
    // Jacobian matrix w.r.t. FFD lattice is right multiplied by world to image matrix;
    // this is reverted here if needed by multiplying with the inverse of it from right
    Matrix orient;
    if (!constraint->WithRespectToWorld()) {
      orient = ffd->Attributes().GetImageToWorldMatrix();
    } else {
      orient = ffd->Attributes().GetImageToWorldMatrix() * ffd->Attributes().GetWorldToImageOrientation();
    }
    if (ffd->Z() == 1) {
      ReorientJacobian2D::Run(npts, jac, &orient);
    } else {
      ReorientJacobian3D::Run(npts, jac, &orient);
    }
  }
  return npts;
}


// -----------------------------------------------------------------------------
/// Compute determinants of 3x3 Jacobian matrices
struct ComputeDeterminant
{
private:

  const Matrix *_Jacobian;
  double       *_Determinant;

public:

  void operator ()(const blocked_range<int> &re) const
  {
    for (int idx = re.begin(); idx != re.end(); ++idx) {
      if (_Jacobian[idx].Rows() == 3) {
        _Determinant[idx] = _Jacobian[idx].Det3x3();
      } else {
        _Determinant[idx] = NaN;
      }
    }
  }

  static void Run(int n, const Matrix *jac, double *det)
  {
    ComputeDeterminant body;
    body._Jacobian    = jac;
    body._Determinant = det;
    parallel_for(blocked_range<int>(0, n), body);
  }
};


// -----------------------------------------------------------------------------
/// Compute adjugate of 3x3 Jacobian matrices
struct ComputeAdjugate
{
private:

  Matrix *_Jacobian;

  /// Compute determinant of 2x2 cofactor matrix
  inline double Det(double a11, double a12, double a21, double a22) const
  {
    return a11 * a22 - a12 * a21;
  }

public:

  void operator ()(const blocked_range<int> &re) const
  {
    Matrix adj(3, 3);
    for (int idx = re.begin(); idx != re.end(); ++idx) {
      const auto &jac = _Jacobian[idx];
      if (jac.Rows() == 3) {
        adj(0, 0) = + Det(jac(1, 1), jac(1, 2), jac(2, 1), jac(2, 2));
        adj(0, 1) = - Det(jac(0, 1), jac(0, 2), jac(2, 1), jac(2, 2));
        adj(0, 2) = + Det(jac(0, 1), jac(0, 2), jac(1, 1), jac(1, 2));
        adj(1, 0) = - Det(jac(1, 0), jac(1, 2), jac(2, 0), jac(2, 2));
        adj(1, 1) = + Det(jac(0, 0), jac(0, 2), jac(2, 0), jac(2, 2));
        adj(1, 2) = - Det(jac(0, 0), jac(0, 2), jac(1, 0), jac(1, 2));
        adj(2, 0) = + Det(jac(1, 0), jac(1, 1), jac(2, 0), jac(2, 1));
        adj(2, 1) = - Det(jac(0, 0), jac(0, 1), jac(2, 0), jac(2, 1));
        adj(2, 2) = + Det(jac(0, 0), jac(0, 1), jac(1, 0), jac(1, 1));
        _Jacobian[idx] = adj;
      }
    }
  }

  static void Run(int n, Matrix *jac)
  {
    ComputeAdjugate body;
    body._Jacobian = jac;
    parallel_for(blocked_range<int>(0, n), body);
  }
};


// -----------------------------------------------------------------------------
/// Evaluate value of Jacobian based penalty term
struct EvaluateConstraint
{
private:

  const JacobianConstraint     *_Constraint;
  const FreeFormTransformation *_FFD;
  const ImageAttributes        *_Domain;
  const WeightImage            *_Mask;
  const double                 *_DetJacobian;
  double                        _Penalty;
  double                        _Norm;

public:

  EvaluateConstraint() {}

  EvaluateConstraint(const EvaluateConstraint &lhs, split)
  :
    _Constraint(lhs._Constraint),
    _FFD(lhs._FFD),
    _Domain(lhs._Domain),
    _Mask(lhs._Mask),
    _DetJacobian(lhs._DetJacobian),
    _Penalty(0.),
    _Norm(0)
  {}

  void join(const EvaluateConstraint &rhs)
  {
    _Penalty += rhs._Penalty;
    _Norm    += rhs._Norm;
  }

  void operator ()(const blocked_range3d<int> &re)
  {
    int idx;
    double w = 1.;
    for (int l = 0; l < _Domain->T(); ++l) {
      for (int k = re.pages().begin(); k != re.pages().end(); ++k)
      for (int j = re.rows ().begin(); j != re.rows ().end(); ++j)
      for (int i = re.cols ().begin(); i != re.cols ().end(); ++i) {
        idx = _Domain->LatticeToIndex(i, j, k, l);
        if (!IsNaN(_DetJacobian[idx])) {
          if (_Mask) w = _Mask->Get(i, j, k);
          _Penalty += w * _Constraint->Penalty(_DetJacobian[idx]);
          _Norm += w;
        }
      }
    }
  }

  static void Run(const JacobianConstraint     *constraint,
                  const FreeFormTransformation *ffd,
                  const ImageAttributes        &domain,
                  const WeightImage            &mask,
                  const double                 *det,
                  double                       &penalty)
  {
    EvaluateConstraint body;
    body._Constraint    = constraint;
    body._FFD           = ffd;
    body._Domain        = &domain;
    body._Mask          = mask.IsEmpty() ? nullptr : &mask;
    body._DetJacobian   = det;
    body._Penalty       = 0.;
    body._Norm          = 0.;
    blocked_range3d<int> range(0, domain.Z(), 0, domain.Y(), 0, domain.X());
    parallel_reduce(range, body);
    if (body._Norm > 0.) {
      penalty += body._Penalty / body._Norm;
    }
  }
};


// -----------------------------------------------------------------------------
/// Evaluate value of symmetric Jacobian based penalty term for SVFFD
struct EvaluateSymmetricConstraint
{
private:

  const JacobianConstraint     *_Constraint;
  const FreeFormTransformation *_FFD;
  const ImageAttributes        *_Domain;
  const WeightImage            *_Mask;
  const double                 *_DetJacobian;
  double                        _Penalty;
  double                        _Norm;

public:

  EvaluateSymmetricConstraint() {}

  EvaluateSymmetricConstraint(const EvaluateSymmetricConstraint &lhs, split)
  :
    _Constraint(lhs._Constraint),
    _FFD(lhs._FFD),
    _Domain(lhs._Domain),
    _Mask(lhs._Mask),
    _DetJacobian(lhs._DetJacobian),
    _Penalty(0.),
    _Norm(0.)
  {}

  void join(const EvaluateSymmetricConstraint &rhs)
  {
    _Penalty += rhs._Penalty;
    _Norm    += rhs._Norm;
  }

  void operator ()(const blocked_range3d<int> &re)
  {
    int idx1, idx2;
    double w = 1.;
    for (int k = re.pages().begin(); k != re.pages().end(); ++k)
    for (int j = re.rows ().begin(); j != re.rows ().end(); ++j)
    for (int i = re.cols ().begin(); i != re.cols ().end(); ++i) {
      idx1 = _Domain->LatticeToIndex(i, j, k, 0);
      idx2 = _Domain->LatticeToIndex(i, j, k, 1);
      if (!IsNaN(_DetJacobian[idx1])) {
        if (_Mask) w = _Mask->Get(i, j, k);
        _Penalty += w * _Constraint->Penalty(_DetJacobian[idx1]);
        _Penalty += w * _Constraint->Penalty(_DetJacobian[idx2]);
        _Norm += 2. * w;
      }
    }
  }

  static void Run(const JacobianConstraint     *constraint,
                  const FreeFormTransformation *ffd,
                  const ImageAttributes        &domain,
                  const WeightImage            &mask,
                  const double                 *det,
                  double                       &penalty)
  {
    EvaluateSymmetricConstraint body;
    body._Constraint    = constraint;
    body._FFD           = ffd;
    body._Domain        = &domain;
    body._Mask          = mask.IsEmpty() ? nullptr : &mask;
    body._DetJacobian   = det;
    body._Penalty       = 0.;
    body._Norm          = 0;
    blocked_range3d<int> range(0, domain.Z(), 0, domain.Y(), 0, domain.X());
    parallel_reduce(range, body);
    if (body._Norm > 0) penalty += body._Penalty / body._Norm;
  }
};


// -----------------------------------------------------------------------------
/// Generic function for evaluation of Jacobian constraint gradient
struct EvaluateGradient_AnyFFD
{
  const JacobianConstraint     *_Constraint;
  const FreeFormTransformation *_FFD;
  const ImageAttributes        *_Domain;
  const double                 *_DetJacobian;
  const Matrix                 *_AdjJacobian;
  double                       *_Gradient;
  const WeightImage            *_Mask;
  double                        _Weight;

  void operator ()(const blocked_range3d<int> &re) const
  {
    int    idx, xdof, ydof, zdof, cp, i1, j1, k1, l1, i2, j2, k2, l2;
    double x, y, z, t, df, dp[3], gradient[3], w = 1.;

    #if _USE_INNER_GRADIENT_SUM_NORM
      double norm;
    #endif

    // Loop over control points
    for (int cl = 0; cl < _FFD->T(); ++cl) {
      for (int ck = re.pages().begin(); ck != re.pages().end(); ++ck)
      for (int cj = re.rows ().begin(); cj != re.rows ().end(); ++cj)
      for (int ci = re.cols ().begin(); ci != re.cols ().end(); ++ci) {
        cp = _FFD->LatticeToIndex(ci, cj, ck, cl);
        if (_Constraint->ConstrainPassiveDoFs() || _FFD->IsActive(cp)) {
          if (_FFD->BoundingBox(*_Domain, cp, i1, j1, k1, l1, i2, j2, k2, l2)) {

            gradient[0] = gradient[1] = gradient[2] = 0.;

            #if _USE_INNER_GRADIENT_SUM_NORM
              norm = 0.;
            #endif

            // Loop over image domain
            for (int l = l1; l <= l2; ++l)
            for (int k = k1; k <= k2; ++k)
            for (int j = j1; j <= j2; ++j)
            for (int i = i1; i <= i2; ++i) {

              // Lattice point index
              idx = _Domain->LatticeToIndex(i, j, k, l);
              if (IsNaN(_DetJacobian[idx])) {
                continue;
              }

              // Derivative of penalty term w.r.t. Jacobian determinant
              df = _Constraint->DerivativeWrtJacobianDet(_DetJacobian[idx]);
              if (!IsZero(df)) {

                x = i, y = j, z = k;
                _Domain->LatticeToWorld(x, y, z);
                t = _Domain->LatticeToTime(l);

                // Spatially varying weight
                if (_Mask) {
                  w = _Mask->Get(i, j, k);
                }

                // Derivatives of Jacobian determinant w.r.t. DoFs of control point
                // (https://en.wikipedia.org/wiki/Jacobi's_formula)
                if (_Constraint->ConstrainParameterization()) {
                  _FFD->FFDJacobianDetDerivative(dp, _AdjJacobian[idx], cp, x, y, z, t,
                                                 _Constraint->WithRespectToWorld(),
                                                 _Constraint->UseLatticeSpacing());
                } else {
                  _FFD->JacobianDetDerivative(dp, _AdjJacobian[idx], cp, x, y, z, t,
                                              _Constraint->WithRespectToWorld(),
                                              _Constraint->UseLatticeSpacing());
                }

                // Apply chain rule and add to sum of gradient vectors
                gradient[0] += w * df * dp[0];
                gradient[1] += w * df * dp[1];
                gradient[2] += w * df * dp[2];

                #if _USE_INNER_GRADIENT_SUM_NORM
                  norm += w;
                #endif
              }
            }

            // Divide by size of support region
            #if _USE_INNER_GRADIENT_SUM_NORM
              if (norm > 0.) {
                gradient[0] /= norm;
                gradient[1] /= norm;
                gradient[2] /= norm;
              }
            #endif

            // Add to total gradient
            _FFD->IndexToDOFs(cp, xdof, ydof, zdof);
            _Gradient[xdof] += _Weight * gradient[0];
            _Gradient[ydof] += _Weight * gradient[1];
            _Gradient[zdof] += _Weight * gradient[2];
          }
        }
      }
    }
  }

  static void Run(const JacobianConstraint     *constraint,
                  const FreeFormTransformation *ffd,
                  const ImageAttributes        &domain,
                  const double                 *det,
                  const Matrix                 *adj,
                  double                       *gradient,
                  const WeightImage            &mask,
                  double                        weight)
  {
    EvaluateGradient_AnyFFD body;
    body._Constraint    = constraint;
    body._FFD           = ffd;
    body._Domain        = &domain;
    body._DetJacobian   = det;
    body._AdjJacobian   = adj;
    body._Gradient      = gradient;
    body._Mask          = mask.IsEmpty() ? nullptr : &mask;
    body._Weight        = weight;
    parallel_for(blocked_range3d<int>(0, ffd->Z(), 0, ffd->Y(), 0, ffd->X()), body);
  }
};


// -----------------------------------------------------------------------------
/// Evaluate gradient of Jacobian constraint evaluated at lattice points of
/// arbitrary domain w.r.t. coefficients of 3D cubic B-spline FFD.
struct EvaluateGradient_Domain_BSplineFFD_3D
{
  const JacobianConstraint              *_Constraint;
  const BSplineFreeFormTransformation3D *_FFD;
  const ImageAttributes                 *_Domain;
  const double                          *_DetJacobian;
  const Matrix                          *_AdjJacobian;
  double                                *_Gradient;
  const WeightImage                     *_Mask;
  double                                 _Weight;

  void operator ()(const blocked_range3d<int> &re) const
  {
    int    idx, xdof, ydof, zdof, cp, i1, j1, k1, i2, j2, k2;
    double x, y, z, df, dp[3], gradient[3];
    double w = 1.;

    #if _USE_INNER_GRADIENT_SUM_NORM
      double norm;
    #endif

    // Loop over control points
    for (int ck = re.pages().begin(); ck != re.pages().end(); ++ck)
    for (int cj = re.rows ().begin(); cj != re.rows ().end(); ++cj)
    for (int ci = re.cols ().begin(); ci != re.cols ().end(); ++ci) {
      cp = _FFD->LatticeToIndex(ci, cj, ck);
      if (_Constraint->ConstrainPassiveDoFs() || _FFD->IsActive(cp)) {
        if (_FFD->BoundingBox(*_Domain, cp, i1, j1, k1, i2, j2, k2)) {

          gradient[0] = gradient[1] = gradient[2] = 0.;

          #if _USE_INNER_GRADIENT_SUM_NORM
            norm = 0.;
          #endif

          // Loop over image domain
          for (int k = k1; k <= k2; ++k)
          for (int j = j1; j <= j2; ++j)
          for (int i = i1; i <= i2; ++i) {

            // Lattice point index
            idx = _Domain->LatticeToIndex(i, j, k);
            if (IsNaN(_DetJacobian[idx])) {
              continue;
            }

            // Derivative of penalty term w.r.t. Jacobian determinant
            df = _Constraint->DerivativeWrtJacobianDet(_DetJacobian[idx]);
            if (IsZero(df)) {
              continue;
            }

            // Spatially varying weight
            if (_Mask) {
              w = _Mask->Get(i, j, k);
            }

            // FFD lattice coordinates
            x = i, y = j, z = k;
            _Domain->LatticeToWorld(x, y, z);
            _FFD->WorldToLattice(x, y, z);

            // Derivatives of Jacobian determinant w.r.t. DoFs of control point
            // (https://en.wikipedia.org/wiki/Jacobi's_formula)
            _FFD->EvaluateJacobianDetDerivative(dp, _AdjJacobian[idx],
                                                x - ci, y - cj, z - ck,
                                                _Constraint->WithRespectToWorld(),
                                                _Constraint->UseLatticeSpacing());

            // Apply chain rule and add to sum of gradient vectors
            gradient[0] += w * df * dp[0];
            gradient[1] += w * df * dp[1];
            gradient[2] += w * df * dp[2];

            #if _USE_INNER_GRADIENT_SUM_NORM
              norm += w;
            #endif
          }

          // Divide by size of support region
          #if _USE_INNER_GRADIENT_SUM_NORM
            if (norm > 0.) {
              gradient[0] /= norm;
              gradient[1] /= norm;
              gradient[2] /= norm;
            }
          #endif

          // Add to total gradient
          _FFD->IndexToDOFs(cp, xdof, ydof, zdof);
          _Gradient[xdof] += _Weight * gradient[0];
          _Gradient[ydof] += _Weight * gradient[1];
          _Gradient[zdof] += _Weight * gradient[2];
        }
      }
    }
  }

  static void Run(const JacobianConstraint              *constraint,
                  const BSplineFreeFormTransformation3D *ffd,
                  const ImageAttributes                 &domain,
                  const double                          *det,
                  const Matrix                          *adj,
                  double                                *gradient,
                  const WeightImage                     &mask,
                  double                                 weight)
  {
    EvaluateGradient_Domain_BSplineFFD_3D body;
    body._Constraint    = constraint;
    body._FFD           = ffd;
    body._Domain        = &domain;
    body._DetJacobian   = det;
    body._AdjJacobian   = adj;
    body._Gradient      = gradient;
    body._Mask          = mask.IsEmpty() ? nullptr : &mask;
    body._Weight        = weight;
    parallel_for(blocked_range3d<int>(0, ffd->Z(), 0, ffd->Y(), 0, ffd->X()), body);
  }
};


// -----------------------------------------------------------------------------
/// Evaluate gradient of Jacobian constraint evaluated at lattice points of
/// arbitrary domain w.r.t. coefficients of 3D cubic B-spline SVF FFD.
struct EvaluateGradient_Domain_BSplineFFD_SV
{
  const JacobianConstraint              *_Constraint;
  const BSplineFreeFormTransformationSV *_FFD;
  const ImageAttributes                 *_Domain;
  const double                          *_DetJacobian;
  const Matrix                          *_AdjJacobian;
  double                                *_Gradient;
  const WeightImage                     *_Mask;
  double                                 _Weight;

  void operator ()(const blocked_range3d<int> &re) const
  {
    Vector3D<double> grad1, grad2;
    int    idx1, idx2, xdof, ydof, zdof, cp, i1, j1, k1, i2, j2, k2;
    double x, y, z, df, a, b, c, wa, wb, wc, da, db, dc, du, dv, dw;
    double w = 1.;

    #if _USE_INNER_GRADIENT_SUM_NORM
      double norm;
    #endif

    // Loop over control points
    for (int ck = re.pages().begin(); ck != re.pages().end(); ++ck)
    for (int cj = re.rows ().begin(); cj != re.rows ().end(); ++cj)
    for (int ci = re.cols ().begin(); ci != re.cols ().end(); ++ci) {
      cp = _FFD->LatticeToIndex(ci, cj, ck);
      if (_Constraint->ConstrainPassiveDoFs() || _FFD->IsActive(cp)) {
        if (_FFD->BoundingBox(*_Domain, cp, i1, j1, k1, i2, j2, k2)) {

          grad1 = 0.;
          grad2 = 0.;

          #if _USE_INNER_GRADIENT_SUM_NORM
            norm = 0.;
          #endif

          // Loop over image domain
          for (int k = k1; k <= k2; ++k)
          for (int j = j1; j <= j2; ++j)
          for (int i = i1; i <= i2; ++i) {

            // Lattice point index
            idx1 = _Domain->LatticeToIndex(i, j, k, 0);
            idx2 = _Domain->LatticeToIndex(i, j, k, 1);
            if (IsNaN(_DetJacobian[idx1])) {
              continue;
            }

            // Spatially varying weight
            if (_Mask) {
              w = _Mask->Get(i, j, k);
            }

            // FFD lattice coordinates
            x = i, y = j, z = k;
            _Domain->LatticeToWorld(x, y, z);
            _FFD->WorldToLattice(x, y, z);

            // Values of the B-spline basis functions and its 1st derivatives
            // Note: Calling Kernel::B faster/not slower than Kernel::VariableToIndex + Kernel:LookupTable.
            a = x - ci;
            wa = BSpline<double>::B  (a);
            da = BSpline<double>::B_I(a);

            b = y - cj;
            wb = BSpline<double>::B  (b);
            db = BSpline<double>::B_I(b);

            c = z - ck;
            wc = BSpline<double>::B  (c);
            dc = BSpline<double>::B_I(c);

            // Derivatives of 3D cubic B-spline kernel w.r.t. lattice distance from control point
            du = da * wb * wc;
            dv = wa * db * wc;
            dw = wa * wb * dc;

            // Apply chain rule for da/dx, db/dx, dc/dx and sum terms (same for y and z)
            // to get the derivative of the 3D cubic B-spline kernel centered at this
            // control point w.r.t. each spatial world coordinate
            if (_Constraint->WithRespectToWorld()) {
              if (_Constraint->UseLatticeSpacing()) {
                _FFD->JacobianToWorld(du, dv, dw);
              } else {
                _FFD->JacobianToWorldOrientation(du, dv, dw);
              }
            }

            // Derivative of penalty term w.r.t. Jacobian determinant
            df = _Constraint->DerivativeWrtJacobianDet(_DetJacobian[idx1]);

            // Apply chain rule and Jacobi's formula to get derivatives of Jacobian determinant
            // (https://en.wikipedia.org/wiki/Jacobi's_formula)
            if (!IsZero(df)) {
              const auto &adj = _AdjJacobian[idx1];
              grad1._x += w * df * (adj(0, 0) * du + adj(1, 0) * dv + adj(2, 0) * dw);
              grad1._y += w * df * (adj(0, 1) * du + adj(1, 1) * dv + adj(2, 1) * dw);
              grad1._z += w * df * (adj(0, 2) * du + adj(1, 2) * dv + adj(2, 2) * dw);
              #if _USE_INNER_GRADIENT_SUM_NORM
                norm += w;
              #endif
            }

            // Derivative of penalty term w.r.t. Jacobian determinant of inverse mapping
            df = _Constraint->DerivativeWrtJacobianDet(_DetJacobian[idx2]);

            // Apply chain rule and Jacobi's formula to get derivatives of Jacobian determinant
            // (https://en.wikipedia.org/wiki/Jacobi's_formula)
            if (!IsZero(df)) {
              const auto &adj = _AdjJacobian[idx2];
              grad2._x -= w * df * (adj(0, 0) * du + adj(1, 0) * dv + adj(2, 0) * dw);
              grad2._y -= w * df * (adj(0, 1) * du + adj(1, 1) * dv + adj(2, 1) * dw);
              grad2._z -= w * df * (adj(0, 2) * du + adj(1, 2) * dv + adj(2, 2) * dw);
              #if _USE_INNER_GRADIENT_SUM_NORM
                norm += w;
              #endif
            }
          }

          // Divide by size of support region
          #if _USE_INNER_GRADIENT_SUM_NORM
            if (norm > 0.) {
              grad1 /= norm;
              grad2 /= norm;
            }
          #endif

          // Add to total gradient
          _FFD->IndexToDOFs(cp, xdof, ydof, zdof);
          _Gradient[xdof] += _Weight * (grad1._x + grad2._x);
          _Gradient[ydof] += _Weight * (grad1._y + grad2._y);
          _Gradient[zdof] += _Weight * (grad1._z + grad2._z);
        }
      }
    }
  }

  static void Run(const JacobianConstraint              *constraint,
                  const BSplineFreeFormTransformationSV *ffd,
                  const ImageAttributes                 &domain,
                  const double                          *det,
                  const Matrix                          *adj,
                  double                                *gradient,
                  const WeightImage                     &mask,
                  double                                 weight)
  {
    EvaluateGradient_Domain_BSplineFFD_SV body;
    body._Constraint    = constraint;
    body._FFD           = ffd;
    body._Domain        = &domain;
    body._DetJacobian   = det;
    body._AdjJacobian   = adj;
    body._Gradient      = gradient;
    body._Mask          = mask.IsEmpty() ? nullptr : &mask;
    body._Weight        = .5 * weight;
    parallel_for(blocked_range3d<int>(0, ffd->Z(), 0, ffd->Y(), 0, ffd->X()), body);
  }
};


// -----------------------------------------------------------------------------
/// Evaluate gradient of Jacobian constraint evaluated at lattice points of
/// control point grid w.r.t. coefficients of 3D cubic B-spline FFD.
struct EvaluateGradient_Lattice_BSplineFFD_3D
{
  const JacobianConstraint              *_Constraint;
  const BSplineFreeFormTransformation3D *_FFD;
  const double                          *_DetJacobian;
  const Matrix                          *_AdjJacobian;
  double                                *_Gradient;
  const WeightImage                     *_Mask;
  double                                 _Weight;

  void operator ()(const blocked_range3d<int> &re) const
  {
    int    idx, xdof, ydof, zdof, cp, i1, j1, k1, i2, j2, k2;
    double df, dp[3], gradient[3], w = 1.;

    #if _USE_INNER_GRADIENT_SUM_NORM
      double norm;
    #endif

    // Loop over control points
    for (int ck = re.pages().begin(); ck != re.pages().end(); ++ck)
    for (int cj = re.rows ().begin(); cj != re.rows ().end(); ++cj)
    for (int ci = re.cols ().begin(); ci != re.cols ().end(); ++ci) {
      cp = _FFD->LatticeToIndex(ci, cj, ck);
      if (_Constraint->ConstrainPassiveDoFs() || _FFD->IsActive(cp)) {

        gradient[0] = gradient[1] = gradient[2] = 0.;

        i1 = max(ci - 1, 0);
        j1 = max(cj - 1, 0);
        k1 = max(ck - 1, 0);

        i2 = min(ci + 1, _FFD->X() - 1);
        j2 = min(cj + 1, _FFD->Y() - 1);
        k2 = min(ck + 1, _FFD->Z() - 1);

        #if _USE_INNER_GRADIENT_SUM_NORM
          norm = 0.;
        #endif

        // Loop over lattice domain
        for (int k = k1; k <= k2; ++k)
        for (int j = j1; j <= j2; ++j)
        for (int i = i1; i <= i2; ++i) {

          // Lattice point index
          idx = _FFD->LatticeToIndex(i, j, k);
          if (IsNaN(_DetJacobian[idx])) {
            continue;
          }

          // Derivative of penalty term w.r.t. Jacobian determinant
          df = _Constraint->DerivativeWrtJacobianDet(_DetJacobian[idx]);
          if (!IsZero(df)) {

            // Derivatives of Jacobian determinant w.r.t. DoFs of control point
            // (https://en.wikipedia.org/wiki/Jacobi's_formula)
            _FFD->EvaluateJacobianDetDerivative(dp, _AdjJacobian[idx],
                                                i - ci, j - cj, k - ck,
                                                _Constraint->WithRespectToWorld(),
                                                _Constraint->UseLatticeSpacing());

            // Spatially varying weight
            if (_Mask) {
              w = _Mask->Get(i, j, k);
            }

            // Apply chain rule and add to sum of gradient vectors
            gradient[0] += w * df * dp[0];
            gradient[1] += w * df * dp[1];
            gradient[2] += w * df * dp[2];

            #if _USE_INNER_GRADIENT_SUM_NORM
              norm += w;
            #endif
          }
        }

        // Divide by sum of spatial weights
        #if _USE_INNER_GRADIENT_SUM_NORM
          if (norm > 0.) {
            gradient[0] /= norm;
            gradient[1] /= norm;
            gradient[2] /= norm;
          }
        #endif

        // Add to total gradient
        _FFD->IndexToDOFs(cp, xdof, ydof, zdof);
        _Gradient[xdof] += _Weight * gradient[0];
        _Gradient[ydof] += _Weight * gradient[1];
        _Gradient[zdof] += _Weight * gradient[2];
      }
    }
  }

  static void Run(const JacobianConstraint              *constraint,
                  const BSplineFreeFormTransformation3D *ffd,
                  const double                          *det,
                  const Matrix                          *adj,
                  double                                *gradient,
                  const WeightImage                     &mask,
                  double                                 weight)
  {
    EvaluateGradient_Lattice_BSplineFFD_3D body;
    body._Constraint    = constraint;
    body._FFD           = ffd;
    body._DetJacobian   = det;
    body._AdjJacobian   = adj;
    body._Gradient      = gradient;
    body._Mask          = mask.IsEmpty() ? nullptr : &mask;
    body._Weight        = weight;
    parallel_for(blocked_range3d<int>(0, ffd->Z(), 0, ffd->Y(), 0, ffd->X()), body);
  }
};


// -----------------------------------------------------------------------------
/// Evaluate gradient of Jacobian constraint evaluated at lattice points of
/// control point grid w.r.t. coefficients of 3D cubic B-spline SVF FFD.
struct EvaluateGradient_Lattice_BSplineFFD_SV
{
  const JacobianConstraint              *_Constraint;
  const BSplineFreeFormTransformationSV *_FFD;
  const double                          *_DetJacobian;
  const Matrix                          *_AdjJacobian;
  double                                *_Gradient;
  const WeightImage                     *_Mask;
  double                                 _Weight;

  void operator ()(const blocked_range3d<int> &re) const
  {
    Vector3D<double> grad1, grad2;
    int    idx1, idx2, xdof, ydof, zdof, cp, i1, j1, k1, i2, j2, k2, a, b, c;
    double df, wa, wb, wc, da, db, dc, du, dv, dw, w = 1.;

    #if _USE_INNER_GRADIENT_SUM_NORM
      double norm;
    #endif

    // Loop over control points
    for (int ck = re.pages().begin(); ck != re.pages().end(); ++ck)
    for (int cj = re.rows ().begin(); cj != re.rows ().end(); ++cj)
    for (int ci = re.cols ().begin(); ci != re.cols ().end(); ++ci) {
      cp = _FFD->LatticeToIndex(ci, cj, ck);
      if (_Constraint->ConstrainPassiveDoFs() || _FFD->IsActive(cp)) {

        grad1 = 0.;
        grad2 = 0.;

        i1 = max(ci - 1, 0);
        j1 = max(cj - 1, 0);
        k1 = max(ck - 1, 0);

        i2 = min(ci + 1, _FFD->X() - 1);
        j2 = min(cj + 1, _FFD->Y() - 1);
        k2 = min(ck + 1, _FFD->Z() - 1);

        #if _USE_INNER_GRADIENT_SUM_NORM
          norm = 0.;
        #endif

        // Loop over lattice domain
        for (int k = k1; k <= k2; ++k)
        for (int j = j1; j <= j2; ++j)
        for (int i = i1; i <= i2; ++i) {

          // Lattice point index
          idx1 = _FFD->LatticeToIndex(i, j, k, 0);
          idx2 = _FFD->LatticeToIndex(i, j, k, 1);
          if (IsNaN(_DetJacobian[idx1])) {
            continue;
          }

          // Spatially varying weight
          if (_Mask) {
            w = _Mask->Get(i, j, k);
          }

          // Values of the B-spline basis functions and its 1st derivatives at lattice points
          // Note: The order of the cubic B-spline pieces is *not* B0, B1, B2, B3! Hence, the minus signs.
          a = i - ci;
          if (-1 <= a && a <= 1) {
            wa = BSpline<double>::LatticeWeights[a + 1];
            da = - BSpline<double>::LatticeWeights_I[a + 1];
          }

          b = j - cj;
          if (-1 <= b && b <= 1) {
            wb = BSpline<double>::LatticeWeights[b + 1];
            db = - BSpline<double>::LatticeWeights_I[b + 1];
          }

          c = k - ck;
          if (-1 <= c && c <= 1) {
            wc = BSpline<double>::LatticeWeights[c + 1];
            dc = - BSpline<double>::LatticeWeights_I[c + 1];
          }

          // Derivatives of 3D cubic B-spline kernel w.r.t. lattice distance from control point
          du = da * wb * wc;
          dv = wa * db * wc;
          dw = wa * wb * dc;

          // Apply chain rule for da/dx, db/dx, dc/dx and sum terms (same for y and z)
          // to get the derivative of the 3D cubic B-spline kernel centered at this
          // control point w.r.t. each spatial world coordinate
          if (_Constraint->WithRespectToWorld()) {
            if (_Constraint->UseLatticeSpacing()) {
              _FFD->JacobianToWorld(du, dv, dw);
            } else {
              _FFD->JacobianToWorldOrientation(du, dv, dw);
            }
          }

          // Derivative of penalty term w.r.t. Jacobian determinant
          df = _Constraint->DerivativeWrtJacobianDet(_DetJacobian[idx1]);

          // Apply chain rule and Jacobi's formula to get derivatives of Jacobian determinant
          // (https://en.wikipedia.org/wiki/Jacobi's_formula)
          if (!IsZero(df)) {
            const auto &adj = _AdjJacobian[idx1];
            grad1._x += w * df * (adj(0, 0) * du + adj(1, 0) * dv + adj(2, 0) * dw);
            grad1._y += w * df * (adj(0, 1) * du + adj(1, 1) * dv + adj(2, 1) * dw);
            grad1._z += w * df * (adj(0, 2) * du + adj(1, 2) * dv + adj(2, 2) * dw);
            #if _USE_INNER_GRADIENT_SUM_NORM
              norm += w;
            #endif
          }

          // Derivative of penalty term w.r.t. Jacobian determinant of inverse mapping
          df = _Constraint->DerivativeWrtJacobianDet(_DetJacobian[idx2]);

          // Apply chain rule and Jacobi's formula to get derivatives of Jacobian determinant
          // (https://en.wikipedia.org/wiki/Jacobi's_formula)
          if (!IsZero(df)) {
            const auto &adj = _AdjJacobian[idx2];
            grad2._x -= w * df * (adj(0, 0) * du + adj(1, 0) * dv + adj(2, 0) * dw);
            grad2._y -= w * df * (adj(0, 1) * du + adj(1, 1) * dv + adj(2, 1) * dw);
            grad2._z -= w * df * (adj(0, 2) * du + adj(1, 2) * dv + adj(2, 2) * dw);
            #if _USE_INNER_GRADIENT_SUM_NORM
              norm += w;
            #endif
          }
        }

        // Divide by size of support region
        #if _USE_INNER_GRADIENT_SUM_NORM
          if (norm > 0.) {
            grad1 /= norm;
            grad2 /= norm;
          }
        #endif

        // Add to total gradient
        _FFD->IndexToDOFs(cp, xdof, ydof, zdof);
        _Gradient[xdof] += _Weight * (grad1._x + grad2._x);
        _Gradient[ydof] += _Weight * (grad1._y + grad2._y);
        _Gradient[zdof] += _Weight * (grad1._z + grad2._z);
      }
    }
  }

  static void Run(const JacobianConstraint              *constraint,
                  const BSplineFreeFormTransformationSV *ffd,
                  const double                          *det,
                  const Matrix                          *adj,
                  double                                *gradient,
                  const WeightImage                     &mask,
                  double                                 weight)
  {
    EvaluateGradient_Lattice_BSplineFFD_SV body;
    body._Constraint    = constraint;
    body._FFD           = ffd;
    body._DetJacobian   = det;
    body._AdjJacobian   = adj;
    body._Gradient      = gradient;
    body._Mask          = mask.IsEmpty() ? nullptr : &mask;
    body._Weight        = .5 * weight;
    parallel_for(blocked_range3d<int>(0, ffd->Z(), 0, ffd->Y(), 0, ffd->X()), body);
  }
};


// -----------------------------------------------------------------------------
/// Evaluate gradient of Jacobian constraint for given FFD
int EvaluateConstraintGradient(const JacobianConstraint *constraint,
                               const FreeFormTransformation *ffd,
                               const ImageAttributes &domain,
                               const double *det,
                               const Matrix *adj,
                               double *gradient,
                               const WeightImage &mask,
                               double weight)
{
  const int npts = domain.NumberOfPoints();
  auto bffd  = dynamic_cast<const BSplineFreeFormTransformation3D *>(ffd);
  auto svffd = dynamic_cast<const BSplineFreeFormTransformationSV *>(ffd);
  if (!constraint->ConstrainParameterization()) {
    if (svffd) {
      bffd  = nullptr;
      svffd = nullptr;
    }
  }
  if (!constraint->Symmetric()) {
    svffd = nullptr;
  }
  if (svffd) {
    if (constraint->SubDomain() == JacobianConstraint::SD_Lattice) {
      typedef EvaluateGradient_Lattice_BSplineFFD_SV EvaluateGradient;
      EvaluateGradient::Run(constraint, svffd, det, adj, gradient, mask, weight);
    } else {
      typedef EvaluateGradient_Domain_BSplineFFD_SV EvaluateGradient;
      EvaluateGradient::Run(constraint, svffd, domain, det, adj, gradient, mask, weight);
    }
    return 2 * npts;
  }
  if (bffd && bffd->Z() > 1) {
    if (constraint->SubDomain() == JacobianConstraint::SD_Lattice) {
      typedef EvaluateGradient_Lattice_BSplineFFD_3D EvaluateGradient;
      EvaluateGradient::Run(constraint, bffd, det, adj, gradient, mask, weight);
    } else {
      typedef EvaluateGradient_Domain_BSplineFFD_3D EvaluateGradient;
      EvaluateGradient::Run(constraint, bffd, domain, det, adj, gradient, mask, weight);
    }
    return npts;
  }
  typedef EvaluateGradient_AnyFFD EvaluateGradient;
  EvaluateGradient::Run(constraint, ffd, domain, det, adj, gradient, mask, weight);
  return npts;
}


// -----------------------------------------------------------------------------
/// Write image of Jacobian determinants to specified file
void WriteJacobianImage(const char *fname, const ImageAttributes &domain, const double *det)
{
  GenericImage<double> jac(domain, const_cast<double *>(det));
  jac.Write(fname);
}


} // namespace

// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
JacobianConstraint::JacobianConstraint(const char *name, bool constrain_spline)
:
  TransformationConstraint(name),
  _SubDomain(SD_Lattice),
  _MaskBlurring(0.2),
  _ConstrainParameterization(constrain_spline),
  _WithRespectToWorld(true),
  _UseLatticeSpacing(false),
  _Symmetric(true),
  _DetJacobian(nullptr),
  _AdjJacobian(nullptr)
{
}

// -----------------------------------------------------------------------------
JacobianConstraint::~JacobianConstraint()
{
  Deallocate(_DetJacobian);
  Deallocate(_AdjJacobian);
}

// =============================================================================
// Parameters
// =============================================================================

// -----------------------------------------------------------------------------
bool JacobianConstraint::SetWithoutPrefix(const char *param, const char *value)
{
  if (strcmp(param, "Domain") == 0) {
    return FromString(value, _SubDomain);
  }
  if (strstr(param, "W.r.t world"          ) == param ||
      strstr(param, "W.r.t. world"         ) == param ||
      strstr(param, "Wrt world"            ) == param ||
      strstr(param, "With respect to world") == param) {
    return FromString(value, _WithRespectToWorld);
  }
  if (strcmp(param, "Use lattice spacing") == 0 ||
      strcmp(param, "Use grid spacing")    == 0 ||
      strcmp(param, "Use spacing")         == 0) {
    return FromString(value, _UseLatticeSpacing);
  }
  if (strcmp(param, "Symmetric") == 0) {
    return FromString(value, _Symmetric);
  }
  if (strcmp(param, "Mask") == 0) {
    _Mask.Clear();
    if (strcmp(value, "none") == 0 || strcmp(value, "None") == 0) {
      _MaskFileName = "";
    } else {
      _MaskFileName = value;
    }
    return true;
  }
  if (strcmp(param, "Mask blurring") == 0) {
    double sigma;
    if (FromString(value, sigma) && sigma >= 0.) {
      _MaskBlurring = sigma;
      return true;
    }
    return false;
  }
  return TransformationConstraint::SetWithoutPrefix(param, value);
}

// -----------------------------------------------------------------------------
ParameterList JacobianConstraint::Parameter() const
{
  ParameterList params = TransformationConstraint::Parameter();
  InsertWithPrefix(params, "Domain", _SubDomain);
  InsertWithPrefix(params, "W.r.t. world", _WithRespectToWorld);
  InsertWithPrefix(params, "Use lattice spacing", _UseLatticeSpacing);
  InsertWithPrefix(params, "Symmetric", _Symmetric);
  if (!_MaskFileName.empty()) {
    InsertWithPrefix(params, "Mask", _MaskFileName);
  }
  InsertWithPrefix(params, "Mask blurring", _MaskBlurring);
  return params;
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
void JacobianConstraint::Initialize()
{
  ImageAttributes domain;
  const FreeFormTransformation   *ffd  = FFD();
  const MultiLevelTransformation *mffd = MFFD();
  bool symmetric = false;

  _SubDomains.clear();
  _MatW2L.clear();
  _MatL2W.clear();
  Deallocate(_DetJacobian);
  Deallocate(_AdjJacobian);
  _NumJacobian = 0;

  if (ffd) {
    _SubDomains.reserve(1);
    if (_SubDomain == SD_Image) {
      domain = _Domain;
    } else {
      domain = ffd->Attributes();
      if (_SubDomain == SD_Shifted) {
        ShiftDomain(domain);
      } else if (_SubDomain == SD_SubDiv) {
        SubdivDomain(domain);
      }
    }
    _SubDomains.push_back(domain);
    if (dynamic_cast<const BSplineFreeFormTransformationSV *>(ffd) != nullptr) {
      symmetric = _Symmetric;
    }
  } else if (mffd) {
    _SubDomains.reserve(mffd->NumberOfActiveLevels());
    for (int n = 0; n < mffd->NumberOfLevels(); ++n) {
      if (mffd->LocalTransformationIsActive(n)) {
        if (_SubDomain == SD_Image) {
          domain = _Domain;
        } else {
          ffd = mffd->GetLocalTransformation(n);
          domain = ffd->Attributes();
          if (_SubDomain == SD_Shifted) {
            ShiftDomain(domain);
          } else if (_SubDomain == SD_SubDiv) {
            SubdivDomain(domain);
          }
        }
        _SubDomains.push_back(domain);
        if (dynamic_cast<const BSplineFreeFormTransformationSV *>(ffd) != nullptr) {
          symmetric = _Symmetric;
        }
      }
    }
  }

  if (!_SubDomains.empty()) {
    _MatW2L.resize(_SubDomains.size());
    _MatL2W.resize(_SubDomains.size());
    for (size_t i = 0; i < _SubDomains.size(); ++i) {
      _NumJacobian += _SubDomains[i].NumberOfPoints();
      _MatW2L[i] = _SubDomains[i].GetWorldToLatticeMatrix();
      _MatL2W[i] = _SubDomains[i].GetLatticeToWorldMatrix();
      _SubDomains[i]._w2i = &_MatW2L[i];
      _SubDomains[i]._i2w = &_MatL2W[i];
    }
    if (_NumJacobian > 0) {
      if (symmetric) _NumJacobian *= 2;
      _DetJacobian = Allocate<double>(_NumJacobian);
      _AdjJacobian = Allocate<Matrix>(_NumJacobian);
    }
  }

  if (!_MaskFileName.empty()) {
    _Mask.Read(_MaskFileName.c_str());
  }
  _Masks.resize(_SubDomains.size());
  if (!_Mask.IsEmpty()) {
    const double min_scale = 1.2;
    for (size_t i = 0; i < _SubDomains.size(); ++i) {
      const auto domain = _SubDomains[i];
      double sigma_x = 0.;
      double sigma_y = 0.;
      double sigma_z = 0.;
      if (domain._dx > min_scale * _Mask.XSize()) {
        sigma_x = _MaskBlurring * domain._dx;
      }
      if (domain._dy > min_scale * _Mask.YSize()) {
        sigma_y = _MaskBlurring * domain._dy;
      }
      if (domain._dz > min_scale * _Mask.ZSize()) {
        sigma_z = _MaskBlurring * domain._dz;
      }
      using VoxelType = WeightImage::VoxelType;
      GaussianBlurring<VoxelType> blur(sigma_x, sigma_y, sigma_z);
      WeightImage mask;
      blur.Input(&_Mask);
      blur.Output(&mask);
      blur.Run();
      Resampling<VoxelType> resampler(
        domain._x,  domain._y,  domain._z,
        domain._dx, domain._dy, domain._dz
      );
      GenericLinearInterpolateImageFunction<WeightImage> interp;
      resampler.Interpolator(&interp);
      resampler.Input(&mask);
      resampler.Output(&_Masks[i]);
      resampler.Run();
    }
  }
}

// -----------------------------------------------------------------------------
void JacobianConstraint::Update(bool gradient)
{
  if (_SubDomains.empty()) {
    return;
  }

  // Evaluate Jacobian matrices
  const FreeFormTransformation   *ffd  = FFD();
  const MultiLevelTransformation *mffd = MFFD();

  double *det = _DetJacobian;
  Matrix *jac = _AdjJacobian;

  if (mffd) {
    int i = 0;
    for (int n = 0; n < mffd->NumberOfLevels(); ++n) {
      if (mffd->LocalTransformationIsActive(n)) {
        ffd = mffd->GetLocalTransformation(n);
        auto npts = EvaluateJacobian(this, ffd, _SubDomains[i], _Masks[i], jac);
        det += npts;
        jac += npts;
        ++i;
      }
    }
  } else if (ffd) {
    EvaluateJacobian(this, ffd, _SubDomains[0], _Masks[0], jac);
  }

  // Evaluate determinants
  ComputeDeterminant::Run(_NumJacobian, _AdjJacobian, _DetJacobian);

  // Compute adjugate matrices required for gradient calculation
  if (gradient) {
    ComputeAdjugate::Run(_NumJacobian, _AdjJacobian);
  }
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
double JacobianConstraint::Evaluate()
{
  if (_SubDomains.empty()) {
    return 0.;
  }

  MIRTK_START_TIMING();

  const FreeFormTransformation   *ffd  = FFD();
  const MultiLevelTransformation *mffd = MFFD();

  double penalty = 0.;
  double *det = _DetJacobian;

  if (mffd) {
    int i = 0;
    for (int n = 0; n < mffd->NumberOfLevels(); ++n) {
      if (mffd->LocalTransformationIsActive(n)) {
        auto npts = _SubDomains[i].NumberOfPoints();
        ffd = mffd->GetLocalTransformation(n);
        auto svffd = dynamic_cast<const BSplineFreeFormTransformationSV *>(ffd);
        if (svffd) {
          EvaluateSymmetricConstraint::Run(this, ffd, _SubDomains[i], _Masks[i], det, penalty);
          det += 2 * npts;
        } else {
          EvaluateConstraint::Run(this, ffd, _SubDomains[i], _Masks[i], det, penalty);
          det += npts;
        }
        ++i;
      }
    }
  } else if (ffd) {
    auto svffd = dynamic_cast<const BSplineFreeFormTransformationSV *>(ffd);
    if (svffd) {
      EvaluateSymmetricConstraint::Run(this, ffd, _SubDomains[0], _Masks[0], det, penalty);
    } else {
      EvaluateConstraint::Run(this, ffd, _SubDomains[0], _Masks[0], det, penalty);
    }
  }

  MIRTK_DEBUG_TIMING(2, (this->HasName() ? this->Name() : this->DefaultName()) << " evaluation");
  return penalty / _SubDomains.size();
}

// -----------------------------------------------------------------------------
void JacobianConstraint::EvaluateGradient(double *gradient, double, double weight)
{
  if (_SubDomains.empty()) {
    return;
  }

  MIRTK_START_TIMING();

  const FreeFormTransformation   *ffd  = FFD();
  const MultiLevelTransformation *mffd = MFFD();

  const double *det = _DetJacobian;
  const Matrix *adj = _AdjJacobian;

  if (mffd) {
    int i = 0;
    for (int n = 0; n < mffd->NumberOfLevels(); ++n) {
      if (mffd->LocalTransformationIsActive(n)) {
        double w = weight;
        #if _DIVIDE_GRADIENT_BY_NPOINTS
          w /= _SubDomains[i].NumberOfPoints();
        #endif
        ffd = mffd->GetLocalTransformation(n);
        auto npts = EvaluateConstraintGradient(this, ffd, _SubDomains[i], det, adj, gradient, _Masks[i], w);
        gradient += ffd->NumberOfDOFs();
        det += npts;
        adj += npts;
        ++i;
      }
    }
  } else if (ffd) {
    double w = weight;
    #if _DIVIDE_GRADIENT_BY_NPOINTS
      w /= _SubDomains[0].NumberOfPoints();
    #endif
    EvaluateConstraintGradient(this, ffd, _SubDomains[0], det, adj, gradient, _Masks[0], w);
  }

  MIRTK_DEBUG_TIMING(2, (this->HasName() ? this->Name() : this->DefaultName()) << " gradient computation");
}

// =============================================================================
// Debugging
// =============================================================================

// -----------------------------------------------------------------------------
void JacobianConstraint::WriteDataSets(const char *p, const char *suffix, bool) const
{
  if (debug < 3 || _SubDomains.empty()) return;

  const int   sz = 1024;
  char        fname[sz];
  string _prefix = Prefix(p);
  const char  *prefix = _prefix.c_str();

  const FreeFormTransformation   *ffd  = FFD();
  const MultiLevelTransformation *mffd = MFFD();

  if (mffd) {
    int i = 0;
    double *det = _DetJacobian;
    for (int n = 0; n < mffd->NumberOfLevels(); ++n) {
      if (mffd->LocalTransformationIsActive(n)) {
        if (!_Masks[i].IsEmpty()) {
          if (mffd->NumberOfActiveLevels() == 1) {
            snprintf(fname, sz, "%smask%s", prefix, suffix);
          } else {
            snprintf(fname, sz, "%smask_for_ffd_at_level_%d%s", prefix, n+1, suffix);
          }
          _Masks[i].Write(fname);
        }
        if (mffd->NumberOfActiveLevels() == 1) {
          snprintf(fname, sz, "%sjacobian_determinant%s", prefix, suffix);
        } else {
          snprintf(fname, sz, "%sjacobian_determinant_of_ffd_at_level_%d%s", prefix, n+1, suffix);
        }
        WriteJacobianImage(fname, _SubDomains[i], det);
        det += _SubDomains[i].NumberOfPoints();
        ++i;
      }
    }
  } else if (ffd) {
    if (!_Masks[0].IsEmpty()) {
      snprintf(fname, sz, "%smask%s", prefix, suffix);
      _Masks[0].Write(fname);
    }
    snprintf(fname, sz, "%sjacobian_determinant%s", prefix, suffix);
    WriteJacobianImage(fname, _SubDomains[0], _DetJacobian);
  }
}

// -----------------------------------------------------------------------------
void JacobianConstraint::WriteGradient(const char *p, const char *suffix) const
{
  if (_SubDomains.empty()) return;

  const FreeFormTransformation   *ffd  = FFD();
  const MultiLevelTransformation *mffd = MFFD();
  if (!ffd && !mffd) return;

  const int   sz = 1024;
  char        fname[sz];
  string _prefix = Prefix(p);
  const char  *prefix = _prefix.c_str();

  double *gradient = CAllocate<double>(_Transformation->NumberOfDOFs());
  const_cast<JacobianConstraint *>(this)->EvaluateGradient(gradient, 1., 1.);

  if (mffd) {
    double *grad = gradient;
    for (int n = 0; n < mffd->NumberOfLevels(); ++n) {
      if (mffd->LocalTransformationIsActive(n)) {
        ffd = mffd->GetLocalTransformation(n);
        if (mffd->NumberOfActiveLevels() == 1) {
          snprintf(fname, sz, "%sgradient%s", prefix, suffix);
        } else {
          snprintf(fname, sz, "%sgradient_of_ffd_at_level_%d%s", prefix, n+1, suffix);
        }
        WriteFFDGradient(fname, ffd, grad);
        grad += ffd->NumberOfDOFs();
      }
    }
  } else if (ffd) {
    snprintf(fname, sz, "%sgradient%s", prefix, suffix);
    WriteFFDGradient(fname, ffd, gradient);
  }

  Deallocate(gradient);
}


} // namespace mirtk
