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
inline bool IsActiveLattice(const FreeFormTransformation3D *ffd, int ci, int cj, int ck)
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
struct UpdateJacobian_AnyFFD
{
private:

  const JacobianConstraint     *_Constraint;
  const FreeFormTransformation *_FFD;
  const ImageAttributes        *_Domain;
  const Matrix                 *_Orient;
  Matrix                       *_AdjJacobian;
  double                       *_DetJacobian;

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
            x = i, y = j, z = k;
            _Domain->LatticeToWorld(x, y, z);
            idx = _Domain->LatticeToIndex(i, j, k, l);
            _FFD->FFDJacobianWorld(_AdjJacobian[idx], x, y, z, t);
            if (_Orient) {
              if (_FFD->Z() == 1) MultiplyRight2x2(_AdjJacobian[idx], *_Orient);
              else                MultiplyRight3x3(_AdjJacobian[idx], *_Orient);
            }
            _AdjJacobian[idx].Adjugate(_DetJacobian[idx]);
          }
        }
      } else {
        for (int l = 0; l < _Domain->T(); ++l) {
          t = _Domain->LatticeToTime(l);
          for (int k = re.pages().begin(); k != re.pages().end(); ++k)
          for (int j = re.rows ().begin(); j != re.rows ().end(); ++j)
          for (int i = re.cols ().begin(); i != re.cols ().end(); ++i) {
            x = i, y = j, z = k;
            _Domain->LatticeToWorld(x, y, z);
            idx = _Domain->LatticeToIndex(i, j, k, l);
            _FFD->LocalJacobian(_AdjJacobian[idx], x, y, z, t);
            if (_Orient) {
              if (_FFD->Z() == 1) MultiplyRight2x2(_AdjJacobian[idx], *_Orient);
              else                MultiplyRight3x3(_AdjJacobian[idx], *_Orient);
            }
            _AdjJacobian[idx].Adjugate(_DetJacobian[idx]);
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
            x = i, y = j, z = k;
            _Domain->LatticeToWorld(x, y, z);
            idx = _Domain->LatticeToIndex(i, j, k, l);
            if (IsActiveWorld(_FFD, x, y, z, t)) {
              _FFD->FFDJacobianWorld(_AdjJacobian[idx], x, y, z, t);
              if (_Orient) {
                if (_FFD->Z() == 1) MultiplyRight2x2(_AdjJacobian[idx], *_Orient);
                else                MultiplyRight3x3(_AdjJacobian[idx], *_Orient);
              }
              _AdjJacobian[idx].Adjugate(_DetJacobian[idx]);
            } else {
              _DetJacobian[idx] = NaN;
            }
          }
        }
      } else {
        for (int l = 0; l < _Domain->T(); ++l) {
          t = _Domain->LatticeToTime(l);
          for (int k = re.pages().begin(); k != re.pages().end(); ++k)
          for (int j = re.rows ().begin(); j != re.rows ().end(); ++j)
          for (int i = re.cols ().begin(); i != re.cols ().end(); ++i) {
            x = i, y = j, z = k;
            _Domain->LatticeToWorld(x, y, z);
            idx = _Domain->LatticeToIndex(i, j, k, l);
            if (IsActiveWorld(_FFD, x, y, z, t)) {
              _FFD->LocalJacobian(_AdjJacobian[idx], x, y, z, t);
              if (_Orient) {
                if (_FFD->Z() == 1) MultiplyRight2x2(_AdjJacobian[idx], *_Orient);
                else                MultiplyRight3x3(_AdjJacobian[idx], *_Orient);
              }
              _AdjJacobian[idx].Adjugate(_DetJacobian[idx]);
            } else {
              _DetJacobian[idx] = NaN;
            }
          }
        }
      }
    }
  }

  static void Run(const JacobianConstraint     *constraint,
                  const FreeFormTransformation *ffd,
                  const ImageAttributes        &domain,
                  double                       *det,
                  Matrix                       *adj)
  {
    // Jacobian matrix w.r.t. FFD lattice is right multiplied by world to image matrix;
    // this is reverted here if needed by multiplying with the inverse of it from right
    UniquePtr<Matrix> orient;
    if (!constraint->WithRespectToWorld()) {
      orient.reset(new Matrix(ffd->Attributes().GetImageToWorldMatrix()));
    } else if (!constraint->UseLatticeSpacing()) {
      orient.reset(new Matrix(ffd->Attributes().GetImageToWorldMatrix() * ffd->Attributes().GetWorldToImageOrientation()));
    }
    UpdateJacobian_AnyFFD body;
    body._Constraint  = constraint;
    body._FFD         = ffd;
    body._Domain      = &domain;
    body._Orient      = orient.get();
    body._AdjJacobian = adj;
    body._DetJacobian = det;
    blocked_range3d<int> range(0, domain.Z(), 0, domain.Y(), 0, domain.X());
    parallel_for(range, body);
  }
};


// -----------------------------------------------------------------------------
/// Body of Update function specialized for BSplineFreeFormTransformation3D
/// for evaluation of Jacobian determinants of spline function at lattice points
/// of aribtrary image domain in the vicinity of active control points.
struct UpdateJacobian_Domain_ActiveCPs_BSplineFFD_3D
{
private:

  const BSplineFreeFormTransformation3D *_FFD;
  const ImageAttributes                 *_Domain;
  const Matrix                          *_Orient;
  Matrix                                *_AdjJacobian;
  double                                *_DetJacobian;

public:

  void operator ()(const blocked_range3d<int> &re) const
  {
    int idx;
    double x, y, z;
    for (int k = re.pages().begin(); k != re.pages().end(); ++k)
    for (int j = re.rows ().begin(); j != re.rows ().end(); ++j)
    for (int i = re.cols ().begin(); i != re.cols ().end(); ++i) {
      idx = _Domain->LatticeToIndex(i, j, k);
      x = i, y = j, z = k;
      _Domain->LatticeToWorld(x, y, z);
      _FFD->WorldToLattice(x, y, z);
      if (IsActiveLattice(_FFD, x, y, z)) {
        _FFD->EvaluateJacobian(_AdjJacobian[idx], x, y, z);
        if (_Orient) MultiplyRight3x3(_AdjJacobian[idx], *_Orient);
        _AdjJacobian[idx].Adjugate(_DetJacobian[idx]);
      } else {
        _DetJacobian[idx] = NaN;
      }
    }
  }

  static void Run(const BSplineFreeFormTransformation3D *ffd,
                  const ImageAttributes &domain, double *det, Matrix *adj,
                  const Matrix *orient = nullptr)
  {
    UpdateJacobian_Domain_ActiveCPs_BSplineFFD_3D body;
    body._FFD         = ffd;
    body._Domain      = &domain;
    body._Orient      = orient;
    body._AdjJacobian = adj;
    body._DetJacobian = det;
    blocked_range3d<int> range(0, domain.Z(), 0, domain.Y(), 0, domain.X());
    parallel_for(range, body);
  }
};


// -----------------------------------------------------------------------------
/// Body of Update function specialized for BSplineFreeFormTransformation3D
/// for evaluation of Jacobian determinants of spline function at lattice points
/// of control point grid in the vicinity of active control points.
struct UpdateJacobian_Lattice_ActiveCPs_BSplineFFD_3D
{
private:

  const BSplineFreeFormTransformation3D *_FFD;
  const Matrix                          *_Orient;
  Matrix                                *_AdjJacobian;
  double                                *_DetJacobian;

public:

  void operator ()(const blocked_range3d<int> &re) const
  {
    int idx;
    for (int k = re.pages().begin(); k != re.pages().end(); ++k)
    for (int j = re.rows ().begin(); j != re.rows ().end(); ++j)
    for (int i = re.cols ().begin(); i != re.cols ().end(); ++i) {
      idx = _FFD->LatticeToIndex(i, j, k);
      if (IsActiveLattice(_FFD, i, j, k)) {
        _FFD->EvaluateJacobian(_AdjJacobian[idx], i, j, k);
        if (_Orient) MultiplyRight3x3(_AdjJacobian[idx], *_Orient);
        _AdjJacobian[idx].Adjugate(_DetJacobian[idx]);
      } else {
        _DetJacobian[idx] = NaN;
      }
    }
  }

  static void Run(const BSplineFreeFormTransformation3D *ffd,
                  double *det, Matrix *adj, const Matrix *orient = nullptr)
  {
    UpdateJacobian_Lattice_ActiveCPs_BSplineFFD_3D body;
    body._FFD         = ffd;
    body._Orient      = orient;
    body._AdjJacobian = adj;
    body._DetJacobian = det;
    blocked_range3d<int> range(0, ffd->Z(), 0, ffd->Y(), 0, ffd->X());
    parallel_for(range, body);
  }
};


// -----------------------------------------------------------------------------
/// Body of Update function specialized for BSplineFreeFormTransformation3D
/// for evaluation of Jacobian determinants of spline function at lattice points
/// of aribtrary image domain for both active and passive control points.
struct UpdateJacobian_Domain_InclPassiveCPs_BSplineFFD_3D
{
private:

  const BSplineFreeFormTransformation3D *_FFD;
  const ImageAttributes                 *_Domain;
  const Matrix                          *_Orient;
  Matrix                                *_AdjJacobian;
  double                                *_DetJacobian;

public:

  void operator ()(const blocked_range3d<int> &re) const
  {
    int idx;
    double x, y, z;
    for (int k = re.pages().begin(); k != re.pages().end(); ++k)
    for (int j = re.rows ().begin(); j != re.rows ().end(); ++j)
    for (int i = re.cols ().begin(); i != re.cols ().end(); ++i) {
      idx = _Domain->LatticeToIndex(i, j, k);
      x = i, y = j, z = k;
      _Domain->LatticeToWorld(x, y, z);
      _FFD->WorldToLattice(x, y, z);
      _FFD->EvaluateJacobian(_AdjJacobian[idx], x, y, z);
      if (_Orient) MultiplyRight3x3(_AdjJacobian[idx], *_Orient);
      _AdjJacobian[idx].Adjugate(_DetJacobian[idx]);
    }
  }

  static void Run(const BSplineFreeFormTransformation3D *ffd,
                  const ImageAttributes &domain, double *det, Matrix *adj,
                  const Matrix *orient = nullptr)
  {
    UpdateJacobian_Domain_InclPassiveCPs_BSplineFFD_3D body;
    body._FFD         = ffd;
    body._Domain      = &domain;
    body._Orient      = orient;
    body._AdjJacobian = adj;
    body._DetJacobian = det;
    blocked_range3d<int> range(0, domain.Z(), 0, domain.Y(), 0, domain.X());
    parallel_for(range, body);
  }
};


// -----------------------------------------------------------------------------
/// Body of Update function specialized for BSplineFreeFormTransformation3D
/// for evaluation of Jacobian determinants of spline function at lattice points
/// of control point grid for both passive and active control points.
struct UpdateJacobian_Lattice_InclPassiveCPs_BSplineFFD_3D
{
private:

  const BSplineFreeFormTransformation3D *_FFD;
  const Matrix                          *_Orient;
  Matrix                                *_AdjJacobian;
  double                                *_DetJacobian;

public:

  void operator ()(const blocked_range3d<int> &re) const
  {
    int idx;
    for (int k = re.pages().begin(); k != re.pages().end(); ++k)
    for (int j = re.rows ().begin(); j != re.rows ().end(); ++j)
    for (int i = re.cols ().begin(); i != re.cols ().end(); ++i) {
      idx = _FFD->LatticeToIndex(i, j, k);
      _FFD->EvaluateJacobian(_AdjJacobian[idx], i, j, k);
      if (_Orient) MultiplyRight3x3(_AdjJacobian[idx], *_Orient);
      _AdjJacobian[idx].Adjugate(_DetJacobian[idx]);
    }
  }

  static void Run(const BSplineFreeFormTransformation3D *ffd,
                  double *det, Matrix *adj, const Matrix *orient = nullptr)
  {
    UpdateJacobian_Lattice_InclPassiveCPs_BSplineFFD_3D body;
    body._FFD         = ffd;
    body._Orient      = orient;
    body._AdjJacobian = adj;
    body._DetJacobian = det;
    blocked_range3d<int> range(0, ffd->Z(), 0, ffd->Y(), 0, ffd->X());
    parallel_for(range, body);
  }
};


// -----------------------------------------------------------------------------
/// Update Jacobian determinants and adjugate matrices for given FFD
void UpdateJacobian(const JacobianConstraint *constraint,
                    const FreeFormTransformation *ffd,
                    const ImageAttributes &domain,
                    double *det, Matrix *adj)
{
  auto bffd = dynamic_cast<const BSplineFreeFormTransformation3D *>(ffd);
  if (!constraint->ConstrainParameterization()) {
    auto svffd = dynamic_cast<const BSplineFreeFormTransformationSV *>(ffd);
    if (svffd) bffd = nullptr;
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
        typedef UpdateJacobian_Lattice_InclPassiveCPs_BSplineFFD_3D Update;
        Update::Run(bffd, det, adj, orient.get());
      } else {
        typedef UpdateJacobian_Lattice_ActiveCPs_BSplineFFD_3D Update;
        Update::Run(bffd, det, adj, orient.get());
      }
    } else {
      if (constraint->ConstrainPassiveDoFs()) {
        typedef UpdateJacobian_Domain_InclPassiveCPs_BSplineFFD_3D Update;
        Update::Run(bffd, domain, det, adj, orient.get());
      } else {
        typedef UpdateJacobian_Domain_ActiveCPs_BSplineFFD_3D Update;
        Update::Run(bffd, domain, det, adj, orient.get());
      }
    }
  } else {
    typedef UpdateJacobian_AnyFFD Update;
    Update::Run(constraint, ffd, domain, det, adj);
  }
}


// -----------------------------------------------------------------------------
/// Evaluate value of Jacobian based penalty term
struct EvaluateConstraint
{
private:

  const JacobianConstraint     *_Constraint;
  const FreeFormTransformation *_FFD;
  const ImageAttributes        *_Domain;
  const Matrix                 *_AdjJacobian;
  const double                 *_DetJacobian;
  double                        _Penalty;
  int                           _N;

public:

  EvaluateConstraint() {}

  EvaluateConstraint(const EvaluateConstraint &lhs, split)
  :
    _Constraint(lhs._Constraint), _FFD(lhs._FFD),
    _Domain(lhs._Domain),
    _AdjJacobian(lhs._AdjJacobian),
    _DetJacobian(lhs._DetJacobian),
    _Penalty(0.), _N(0)
  {}

  void join(const EvaluateConstraint &rhs)
  {
    _Penalty += rhs._Penalty;
    _N       += rhs._N;
  }

  void operator ()(const blocked_range3d<int> &re)
  {
    int idx;
    for (int l = 0; l < _Domain->T(); ++l) {
      for (int k = re.pages().begin(); k != re.pages().end(); ++k)
      for (int j = re.rows ().begin(); j != re.rows ().end(); ++j)
      for (int i = re.cols ().begin(); i != re.cols ().end(); ++i) {
        idx = _Domain->LatticeToIndex(i, j, k, l);
        if (!IsNaN(_DetJacobian[idx])) {
          _Penalty += _Constraint->Penalty(_DetJacobian[idx]);
          ++_N;
        }
      }
    }
  }

  static void Run(const JacobianConstraint     *constraint,
                  const FreeFormTransformation *ffd,
                  const ImageAttributes        &domain,
                  const double                 *det,
                  const Matrix                 *adj,
                  double                       &penalty)
  {
    EvaluateConstraint body;
    body._Constraint  = constraint;
    body._FFD         = ffd;
    body._Domain      = &domain;
    body._AdjJacobian = adj;
    body._DetJacobian = det;
    body._Penalty     = 0.;
    body._N           = 0;
    blocked_range3d<int> range(0, domain.Z(), 0, domain.Y(), 0, domain.X());
    parallel_reduce(range, body);
    if (body._N > 0) penalty += body._Penalty / body._N;
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
  double                        _Weight;

  void operator ()(const blocked_range3d<int> &re) const
  {
    int    idx, xdof, ydof, zdof, cp, i1, j1, k1, l1, i2, j2, k2, l2;
    double x, y, z, t, df, dp[3], gradient[3];

    // Loop over control points
    for (int cl = 0; cl < _FFD->T(); ++cl) {
      for (int ck = re.pages().begin(); ck != re.pages().end(); ++ck)
      for (int cj = re.rows ().begin(); cj != re.rows ().end(); ++cj)
      for (int ci = re.cols ().begin(); ci != re.cols ().end(); ++ci) {
        cp = _FFD->LatticeToIndex(ci, cj, ck, cl);
        if (_Constraint->ConstrainPassiveDoFs() || _FFD->IsActive(cp)) {
          if (_FFD->BoundingBox(*_Domain, cp, i1, j1, k1, l1, i2, j2, k2, l2)) {

            gradient[0] = gradient[1] = gradient[2] = 0.;

            // Loop over image domain
            for (int l = l1; l <= l2; ++l)
            for (int k = k1; k <= k2; ++k)
            for (int j = j1; j <= j2; ++j)
            for (int i = i1; i <= i2; ++i) {

              // Lattice point index
              idx = _Domain->LatticeToIndex(i, j, k, l);

              // Derivative of penalty term w.r.t. Jacobian determinant
              df = _Constraint->DerivativeWrtJacobianDet(_DetJacobian[idx]);
              if (!IsZero(df)) {

                // Derivatives of Jacobian determinant w.r.t. DoFs of control point
                // (https://en.wikipedia.org/wiki/Jacobi's_formula)
                x = i, y = j, z = k;
                _Domain->LatticeToWorld(x, y, z);
                t = _Domain->LatticeToTime(l);
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
                gradient[0] += df * dp[0];
                gradient[1] += df * dp[1];
                gradient[2] += df * dp[2];
              }
            }

            // Divide by size of support region
            #if _USE_INNER_GRADIENT_SUM_NORM
              int n = (i2 - i1 + 1) * (j2 - j1 + 1) * (k2 - k1 + 1) * (l2 - l1 + 1);
              gradient[0] /= n;
              gradient[1] /= n;
              gradient[2] /= n;
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
                  double                        weight)
  {
    EvaluateGradient_AnyFFD body;
    body._Constraint  = constraint;
    body._FFD         = ffd;
    body._Domain      = &domain;
    body._DetJacobian = det;
    body._AdjJacobian = adj;
    body._Gradient    = gradient;
    body._Weight      = weight;
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
  double                                 _Weight;

  void operator ()(const blocked_range3d<int> &re) const
  {
    int    idx, xdof, ydof, zdof, cp, i1, j1, k1, i2, j2, k2;
    double x, y, z, df, dp[3], gradient[3];

    // Loop over control points
    for (int ck = re.pages().begin(); ck != re.pages().end(); ++ck)
    for (int cj = re.rows ().begin(); cj != re.rows ().end(); ++cj)
    for (int ci = re.cols ().begin(); ci != re.cols ().end(); ++ci) {
      cp = _FFD->LatticeToIndex(ci, cj, ck);
      if (_Constraint->ConstrainPassiveDoFs() || _FFD->IsActive(cp)) {
        if (_FFD->BoundingBox(*_Domain, cp, i1, j1, k1, i2, j2, k2)) {

          gradient[0] = gradient[1] = gradient[2] = 0.;

          // Loop over image domain
          for (int k = k1; k <= k2; ++k)
          for (int j = j1; j <= j2; ++j)
          for (int i = i1; i <= i2; ++i) {

            // Lattice point index
            idx = _Domain->LatticeToIndex(i, j, k);

            // Derivative of penalty term w.r.t. Jacobian determinant
            df = _Constraint->DerivativeWrtJacobianDet(_DetJacobian[idx]);
            if (!IsZero(df)) {

              // Derivatives of Jacobian determinant w.r.t. DoFs of control point
              // (https://en.wikipedia.org/wiki/Jacobi's_formula)
              x = i, y = j, z = k;
              _Domain->LatticeToWorld(x, y, z);
              _FFD->WorldToLattice(x, y, z);
              _FFD->EvaluateJacobianDetDerivative(dp, _AdjJacobian[idx], x - ci, y - cj, z - ck,
                                                  _Constraint->WithRespectToWorld(),
                                                  _Constraint->UseLatticeSpacing());

              // Apply chain rule and add to sum of gradient vectors
              gradient[0] += df * dp[0];
              gradient[1] += df * dp[1];
              gradient[2] += df * dp[2];
            }
          }

          // Divide by size of support region
          #if _USE_INNER_GRADIENT_SUM_NORM
            int n = (i2 - i1 + 1) * (j2 - j1 + 1) * (k2 - k1 + 1);
            gradient[0] /= n;
            gradient[1] /= n;
            gradient[2] /= n;
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
                  double                                 weight)
  {
    EvaluateGradient_Domain_BSplineFFD_3D body;
    body._Constraint  = constraint;
    body._FFD         = ffd;
    body._Domain      = &domain;
    body._DetJacobian = det;
    body._AdjJacobian = adj;
    body._Gradient    = gradient;
    body._Weight      = weight;
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
  double                                 _Weight;

  void operator ()(const blocked_range3d<int> &re) const
  {
    int    idx, xdof, ydof, zdof, cp, i1, j1, k1, i2, j2, k2;
    double df, dp[3], gradient[3];

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

        // Loop over lattice domain
        for (int k = k1; k <= k2; ++k)
        for (int j = j1; j <= j2; ++j)
        for (int i = i1; i <= i2; ++i) {

          // Lattice point index
          idx = _FFD->LatticeToIndex(i, j, k);

          // Derivative of penalty term w.r.t. Jacobian determinant
          df = _Constraint->DerivativeWrtJacobianDet(_DetJacobian[idx]);
          if (!IsZero(df)) {

            // Derivatives of Jacobian determinant w.r.t. DoFs of control point
            // (https://en.wikipedia.org/wiki/Jacobi's_formula)
            _FFD->EvaluateJacobianDetDerivative(dp, _AdjJacobian[idx], i - ci, j - cj, k - ck,
                                                _Constraint->WithRespectToWorld(),
                                                _Constraint->UseLatticeSpacing());

            // Apply chain rule and add to sum of gradient vectors
            gradient[0] += df * dp[0];
            gradient[1] += df * dp[1];
            gradient[2] += df * dp[2];
          }
        }

        // Divide by size of support region
        #if _USE_INNER_GRADIENT_SUM_NORM
          int n = (i2 - i1 + 1) * (j2 - j1 + 1) * (k2 - k1 + 1);
          gradient[0] /= n;
          gradient[1] /= n;
          gradient[2] /= n;
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
                  double                                 weight)
  {
    EvaluateGradient_Lattice_BSplineFFD_3D body;
    body._Constraint  = constraint;
    body._FFD         = ffd;
    body._DetJacobian = det;
    body._AdjJacobian = adj;
    body._Gradient    = gradient;
    body._Weight      = weight;
    parallel_for(blocked_range3d<int>(0, ffd->Z(), 0, ffd->Y(), 0, ffd->X()), body);
  }
};


// -----------------------------------------------------------------------------
/// Evaluate gradient of Jacobian constraint for given FFD
void EvaluateConstraintGradient(const JacobianConstraint *constraint,
                                const FreeFormTransformation *ffd,
                                const ImageAttributes &domain,
                                const double *det, const Matrix *adj,
                                double *gradient, double weight)
{
  auto bffd = dynamic_cast<const BSplineFreeFormTransformation3D *>(ffd);
  if (!constraint->ConstrainParameterization()) {
    auto svffd = dynamic_cast<const BSplineFreeFormTransformationSV *>(ffd);
    if (svffd) bffd = nullptr;
  }
  if (bffd && bffd->Z() > 1) {
    if (constraint->SubDomain() == JacobianConstraint::SD_Lattice) {
      typedef EvaluateGradient_Lattice_BSplineFFD_3D EvaluateGradient;
      EvaluateGradient::Run(constraint, bffd, det, adj, gradient, weight);
    } else {
      typedef EvaluateGradient_Domain_BSplineFFD_3D EvaluateGradient;
      EvaluateGradient::Run(constraint, bffd, domain, det, adj, gradient, weight);
    }
  } else {
    typedef EvaluateGradient_AnyFFD EvaluateGradient;
    EvaluateGradient::Run(constraint, ffd, domain, det, adj, gradient, weight);
  }
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
  _ConstrainParameterization(constrain_spline),
  _WithRespectToWorld(true),
  _UseLatticeSpacing(true),
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
  return TransformationConstraint::SetWithoutPrefix(param, value);
}

// -----------------------------------------------------------------------------
ParameterList JacobianConstraint::Parameter() const
{
  ParameterList params = TransformationConstraint::Parameter();
  InsertWithPrefix(params, "Domain", _SubDomain);
  InsertWithPrefix(params, "W.r.t. world", _WithRespectToWorld);
  InsertWithPrefix(params, "Use lattice spacing", _UseLatticeSpacing);
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

  _SubDomains.clear();
  _MatW2L.clear();
  _MatL2W.clear();
  Deallocate(_DetJacobian);
  Deallocate(_AdjJacobian);

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
  } else if (mffd) {
    _SubDomains.reserve(mffd->NumberOfActiveLevels());
    for (int n = 0; n < mffd->NumberOfLevels(); ++n) {
      if (mffd->LocalTransformationIsActive(n)) {
        if (_SubDomain == SD_Image) {
          domain = _Domain;
        } else {
          domain = mffd->GetLocalTransformation(n)->Attributes();
          if (_SubDomain == SD_Shifted) {
            ShiftDomain(domain);
          } else if (_SubDomain == SD_SubDiv) {
            SubdivDomain(domain);
          }
        }
        _SubDomains.push_back(domain);
      }
    }
  }

  if (!_SubDomains.empty()) {
    int npts = 0;
    _MatW2L.resize(_SubDomains.size());
    _MatL2W.resize(_SubDomains.size());
    for (size_t i = 0; i < _SubDomains.size(); ++i) {
      npts += _SubDomains[i].NumberOfPoints();
      _MatW2L[i] = _SubDomains[i].GetWorldToLatticeMatrix();
      _MatL2W[i] = _SubDomains[i].GetLatticeToWorldMatrix();
      _SubDomains[i]._w2i = &_MatW2L[i];
      _SubDomains[i]._i2w = &_MatL2W[i];
    }
    if (npts > 0) {
      _DetJacobian = Allocate<double>(npts);
      _AdjJacobian = Allocate<Matrix>(npts);
    }
  }
}

// -----------------------------------------------------------------------------
void JacobianConstraint::Update(bool)
{
  if (_SubDomains.empty()) {
    return;
  }

  const FreeFormTransformation   *ffd  = FFD();
  const MultiLevelTransformation *mffd = MFFD();

  double *det = _DetJacobian;
  Matrix *adj = _AdjJacobian;
  auto domain = _SubDomains.begin();

  if (mffd) {
    for (int n = 0; n < mffd->NumberOfLevels(); ++n) {
      if (mffd->LocalTransformationIsActive(n)) {
        ffd = mffd->GetLocalTransformation(n);
        UpdateJacobian(this, ffd, *domain, det, adj);
        auto npts = domain->NumberOfPoints();
        det += npts;
        adj += npts;
        ++domain;
      }
    }
  } else if (ffd) {
    UpdateJacobian(this, ffd, *domain, det, adj);
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
  Matrix *adj = _AdjJacobian;
  auto domain = _SubDomains.begin();

  if (mffd) {
    for (int n = 0; n < mffd->NumberOfLevels(); ++n) {
      if (mffd->LocalTransformationIsActive(n)) {
        ffd = mffd->GetLocalTransformation(n);
        EvaluateConstraint::Run(this, ffd, *domain, det, adj, penalty);
        auto npts = domain->NumberOfPoints();
        det += npts;
        adj += npts;
        ++domain;
      }
    }
  } else if (ffd) {
    EvaluateConstraint::Run(this, ffd, *domain, det, adj, penalty);
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
  auto domain = _SubDomains.begin();

  if (mffd) {
    for (int n = 0; n < mffd->NumberOfLevels(); ++n) {
      if (mffd->LocalTransformationIsActive(n)) {
        auto npts = domain->NumberOfPoints();
        double w = weight;
        #if _DIVIDE_GRADIENT_BY_NPOINTS
          w /= npts;
        #endif
        ffd = mffd->GetLocalTransformation(n);
        EvaluateConstraintGradient(this, ffd, *domain, det, adj, gradient, w);
        gradient += ffd->NumberOfDOFs();
        det += npts;
        adj += npts;
        ++domain;
      }
    }
  } else if (ffd) {
    double w = weight;
    #if _DIVIDE_GRADIENT_BY_NPOINTS
      w /= domain->NumberOfPoints();
    #endif
    EvaluateConstraintGradient(this, ffd, *domain, det, adj, gradient, w);
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
    double *det = _DetJacobian;
    auto domain = _SubDomains.begin();
    for (int n = 0; n < mffd->NumberOfLevels(); ++n) {
      if (mffd->LocalTransformationIsActive(n)) {
        ffd = mffd->GetLocalTransformation(n);
        if (mffd->NumberOfActiveLevels() == 1) {
          snprintf(fname, sz, "%sjacobian_determinant%s", prefix, suffix);
        } else {
          snprintf(fname, sz, "%sjacobian_determinant_of_ffd_at_level_%d%s", prefix, n+1, suffix);
        }
        WriteJacobianImage(fname, *domain, det);
        det += domain->NumberOfPoints();
        ++domain;
      }
    }
  } else if (ffd) {
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
    for (int n = 0; n < mffd->NumberOfLevels(); ++n) {
      if (mffd->LocalTransformationIsActive(n)) {
        ffd = mffd->GetLocalTransformation(n);
        if (mffd->NumberOfActiveLevels() == 1) {
          snprintf(fname, sz, "%sgradient%s", prefix, suffix);
        } else {
          snprintf(fname, sz, "%sgradient_of_ffd_at_level_%d%s", prefix, n+1, suffix);
        }
        WriteFFDGradient(fname, ffd, gradient);
        gradient += ffd->NumberOfDOFs();
      }
    }
  } else if (ffd) {
    snprintf(fname, sz, "%sgradient%s", prefix, suffix);
    WriteFFDGradient(fname, ffd, gradient);
  }

  Deallocate(gradient);
}


} // namespace mirtk
