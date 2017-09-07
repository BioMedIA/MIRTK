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
inline bool IsActive(const FreeFormTransformation *ffd, double x, double y, double z, double t)
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
/// Body of Update function for parallel execution
struct JacobianConstraintUpdate
{
private:

  const JacobianConstraint     *_This;
  const FreeFormTransformation *_FFD;
  const ImageAttributes        *_Domain;
  Matrix                       *_AdjJacobian;
  double                       *_DetJacobian;

public:

  void operator()(const blocked_range3d<int> &re) const
  {
    int idx;
    double x, y, z, t;
    if (_This->ConstrainPassiveDoFs()) {
      if (_This->ConstrainJacobianOfParameterization()) {
        for (int l = 0; l < _Domain->T(); ++l) {
          t = _Domain->LatticeToTime(l);
          for (int k = re.pages().begin(); k != re.pages().end(); ++k)
          for (int j = re.rows ().begin(); j != re.rows ().end(); ++j)
          for (int i = re.cols ().begin(); i != re.cols ().end(); ++i) {
            x = i, y = j, z = k;
            _Domain->LatticeToWorld(x, y, z);
            idx = _Domain->LatticeToIndex(i, j, k, l);
            _FFD->FFDJacobianWorld(_AdjJacobian[idx], x, y, z, t);
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
            _AdjJacobian[idx].Adjugate(_DetJacobian[idx]);
          }
        }
      }
    } else {
      if (_This->ConstrainJacobianOfParameterization()) {
        for (int l = 0; l < _Domain->T(); ++l) {
          t = _Domain->LatticeToTime(l);
          for (int k = re.pages().begin(); k != re.pages().end(); ++k)
          for (int j = re.rows ().begin(); j != re.rows ().end(); ++j)
          for (int i = re.cols ().begin(); i != re.cols ().end(); ++i) {
            x = i, y = j, z = k;
            _Domain->LatticeToWorld(x, y, z);
            idx = _Domain->LatticeToIndex(i, j, k, l);
            if (IsActive(_FFD, x, y, z, t)) {
              _FFD->FFDJacobianWorld(_AdjJacobian[idx], x, y, z, t);
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
            if (IsActive(_FFD, x, y, z, t)) {
              _FFD->LocalJacobian(_AdjJacobian[idx], x, y, z, t);
              _AdjJacobian[idx].Adjugate(_DetJacobian[idx]);
            } else {
              _DetJacobian[idx] = NaN;
            }
          }
        }
      }
    }
  }

  static void Run(const JacobianConstraint     *obj,
                  const FreeFormTransformation *ffd,
                  const ImageAttributes        &domain,
                  double                       *det,
                  Matrix                       *adj)
  {
    JacobianConstraintUpdate body;
    body._This        = obj;
    body._FFD         = ffd;
    body._Domain      = &domain;
    body._AdjJacobian = adj;
    body._DetJacobian = det;
    blocked_range3d<int> range(0, domain.Z(),
                               0, domain.Y(),
                               0, domain.X());
    parallel_for(range, body);
  }
};

// -----------------------------------------------------------------------------
struct JacobianConstraintEvaluate
{
private:

  const JacobianConstraint     *_This;
  const FreeFormTransformation *_FFD;
  const ImageAttributes        *_Domain;
  const Matrix                 *_AdjJacobian;
  const double                 *_DetJacobian;
  double                        _Penalty;
  int                           _N;

public:

  JacobianConstraintEvaluate() {}

  JacobianConstraintEvaluate(JacobianConstraintEvaluate &lhs, split)
  :
    _This(lhs._This), _FFD(lhs._FFD),
    _Domain(lhs._Domain),
    _AdjJacobian(lhs._AdjJacobian),
    _DetJacobian(lhs._DetJacobian),
    _Penalty(0.), _N(0)
  {}

  void join(JacobianConstraintEvaluate &rhs)
  {
    _Penalty += rhs._Penalty;
    _N       += rhs._N;
  }

  void operator()(const blocked_range3d<int> &re)
  {
    int idx;
    for (int l = 0; l < _Domain->T(); ++l) {
      for (int k = re.pages().begin(); k != re.pages().end(); ++k)
      for (int j = re.rows ().begin(); j != re.rows ().end(); ++j)
      for (int i = re.cols ().begin(); i != re.cols ().end(); ++i) {
        idx = _Domain->LatticeToIndex(i, j, k, l);
        if (!IsNaN(_DetJacobian[idx])) {
          _Penalty += _This->Penalty(_DetJacobian[idx]);
          ++_N;
        }
      }
    }
  }

  static void Run(const JacobianConstraint     *obj,
                  const FreeFormTransformation *ffd,
                  const ImageAttributes        &domain,
                  const double                 *det,
                  const Matrix                 *adj,
                  double                       &penalty)
  {
    JacobianConstraintEvaluate body;
    body._This        = obj;
    body._FFD         = ffd;
    body._Domain      = &domain;
    body._AdjJacobian = adj;
    body._DetJacobian = det;
    body._Penalty     = 0.;
    body._N           = 0;
    blocked_range3d<int> range(0, domain.Z(),
                               0, domain.Y(),
                               0, domain.X());
    parallel_reduce(range, body);
    if (body._N > 0) {
      penalty += body._Penalty / body._N;
    }
  }
};

// -----------------------------------------------------------------------------
struct JacobianConstraintEvaluateGradient
{
  const JacobianConstraint     *_This;
  const FreeFormTransformation *_FFD;
  const ImageAttributes        *_Domain;
  const double                 *_DetJacobian;
  const Matrix                 *_AdjJacobian;
  double                       *_Gradient;
  double                        _Weight;

  void operator()(const blocked_range3d<int> &re) const
  {
    int    idx, xdof, ydof, zdof, cp, i1, j1, k1, l1, i2, j2, k2, l2;
    double x, y, z, t, df, dp[3], gradient[3];

    // Loop over control points
    for (int cl = 0; cl < _FFD->T(); ++cl) {
      for (int ck = re.pages().begin(); ck != re.pages().end(); ++ck)
      for (int cj = re.rows ().begin(); cj != re.rows ().end(); ++cj)
      for (int ci = re.cols ().begin(); ci != re.cols ().end(); ++ci) {
        cp = _FFD->LatticeToIndex(ci, cj, ck, cl);
        if (_This->ConstrainPassiveDoFs() || _FFD->IsActive(cp)) {
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
              df = _This->DerivativeWrtJacobianDet(_DetJacobian[idx]);
              if (!AreEqual(df, 0.)) {

                // Derivatives of Jacobian determinant w.r.t. DoFs of control point
                // (https://en.wikipedia.org/wiki/Jacobi's_formula)
                x = i, y = j, z = k;
                _Domain->LatticeToWorld(x, y, z);
                t = _Domain->LatticeToTime(l);
                if (_This->ConstrainJacobianOfParameterization()) {
                  _FFD->FFDJacobianDetDerivative(dp, _AdjJacobian[idx], cp, x, y, z, t);
                } else {
                  _FFD->JacobianDetDerivative(dp, _AdjJacobian[idx], cp, x, y, z, t);
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

  static void Run(const JacobianConstraint     *obj,
                  const FreeFormTransformation *ffd,
                  const ImageAttributes        &domain,
                  const double                 *det,
                  const Matrix                 *adj,
                  double                       *gradient,
                  double                        weight)
  {
    JacobianConstraintEvaluateGradient body;
    body._This        = obj;
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
struct BSplineFFD3DJacobianConstraintEvaluateGradient
{
  const JacobianConstraint              *_This;
  const BSplineFreeFormTransformation3D *_FFD;
  const ImageAttributes                 *_Domain;
  const double                          *_DetJacobian;
  const Matrix                          *_AdjJacobian;
  double                                *_Gradient;
  double                                 _Weight;

  void operator()(const blocked_range3d<int> &re) const
  {
    int    idx, xdof, ydof, zdof, cp, i1, j1, k1, i2, j2, k2;
    double x, y, z, df, dp[3], gradient[3];

    // Loop over control points
    for (int ck = re.pages().begin(); ck != re.pages().end(); ++ck)
    for (int cj = re.rows ().begin(); cj != re.rows ().end(); ++cj)
    for (int ci = re.cols ().begin(); ci != re.cols ().end(); ++ci) {
      cp = _FFD->LatticeToIndex(ci, cj, ck);
      if (_This->ConstrainPassiveDoFs() || _FFD->IsActive(cp)) {
        if (_FFD->BoundingBox(*_Domain, cp, i1, j1, k1, i2, j2, k2)) {

          gradient[0] = gradient[1] = gradient[2] = 0.;

          // Loop over image domain
          for (int k = k1; k <= k2; ++k)
          for (int j = j1; j <= j2; ++j)
          for (int i = i1; i <= i2; ++i) {

            // Lattice point index
            idx = _Domain->LatticeToIndex(i, j, k);

            // Derivative of penalty term w.r.t. Jacobian determinant
            df = _This->DerivativeWrtJacobianDet(_DetJacobian[idx]);
            if (!AreEqual(df, 0.)) {

              // Derivatives of Jacobian determinant w.r.t. DoFs of control point
              // (https://en.wikipedia.org/wiki/Jacobi's_formula)
              x = i, y = j, z = k;
              _Domain->LatticeToWorld(x, y, z);
              _FFD->WorldToLattice(x, y, z);
              _FFD->EvaluateJacobianDetDerivative(dp, _AdjJacobian[idx], x - ci, y - cj, z - ck);

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

  static void Run(const JacobianConstraint              *obj,
                  const BSplineFreeFormTransformation3D *ffd,
                  const ImageAttributes                 &domain,
                  const double                          *det,
                  const Matrix                          *adj,
                  double                                *gradient,
                  double                                 weight)
  {
    BSplineFFD3DJacobianConstraintEvaluateGradient body;
    body._This        = obj;
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
struct BSplineFFD3DJacobianConstraintEvaluateGradientLattice
{
  const JacobianConstraint              *_This;
  const BSplineFreeFormTransformation3D *_FFD;
  const ImageAttributes                 *_Domain;
  const double                          *_DetJacobian;
  const Matrix                          *_AdjJacobian;
  double                                *_Gradient;
  double                                 _Weight;

  void operator()(const blocked_range3d<int> &re) const
  {
    int    idx, xdof, ydof, zdof, cp, i1, j1, k1, i2, j2, k2;
    double df, dp[3], gradient[3];

    // Loop over control points
    for (int ck = re.pages().begin(); ck != re.pages().end(); ++ck)
    for (int cj = re.rows ().begin(); cj != re.rows ().end(); ++cj)
    for (int ci = re.cols ().begin(); ci != re.cols ().end(); ++ci) {
      cp = _FFD->LatticeToIndex(ci, cj, ck);
      if (_This->ConstrainPassiveDoFs() || _FFD->IsActive(cp)) {

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
          idx = _Domain->LatticeToIndex(i, j, k);

          // Derivative of penalty term w.r.t. Jacobian determinant
          df = _This->DerivativeWrtJacobianDet(_DetJacobian[idx]);
          if (!AreEqual(df, 0.)) {

            // Derivatives of Jacobian determinant w.r.t. DoFs of control point
            // (https://en.wikipedia.org/wiki/Jacobi's_formula)
            _FFD->EvaluateJacobianDetDerivative(dp, _AdjJacobian[idx], i - ci, j - cj, k - ck);

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

  static void Run(const JacobianConstraint              *obj,
                  const BSplineFreeFormTransformation3D *ffd,
                  const ImageAttributes                 &domain,
                  const double                          *det,
                  const Matrix                          *adj,
                  double                                *gradient,
                  double                                 weight)
  {
    BSplineFFD3DJacobianConstraintEvaluateGradientLattice body;
    body._This        = obj;
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
  _ConstrainJacobianOfParameterization(constrain_spline),
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
  return TransformationConstraint::SetWithoutPrefix(param, value);
}

// -----------------------------------------------------------------------------
ParameterList JacobianConstraint::Parameter() const
{
  ParameterList params = TransformationConstraint::Parameter();
  InsertWithPrefix(params, "Domain", _SubDomain);
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

  typedef JacobianConstraintUpdate Body;

  const FreeFormTransformation   *ffd  = FFD();
  const MultiLevelTransformation *mffd = MFFD();

  if (mffd) {
    double *det = _DetJacobian;
    Matrix *adj = _AdjJacobian;
    auto domain = _SubDomains.begin();
    for (int n = 0; n < mffd->NumberOfLevels(); ++n) {
      if (mffd->LocalTransformationIsActive(n)) {
        ffd = mffd->GetLocalTransformation(n);
        Body::Run(this, ffd, *domain, det, adj);
        auto npts = domain->NumberOfPoints();
        det += npts;
        adj += npts;
        ++domain;
      }
    }
  } else if (ffd) {
    Body::Run(this, ffd, _SubDomains[0], _DetJacobian, _AdjJacobian);
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

  typedef JacobianConstraintEvaluate Body;

  MIRTK_START_TIMING();

  const FreeFormTransformation   *ffd  = FFD();
  const MultiLevelTransformation *mffd = MFFD();

  double penalty = .0;

  if (mffd) {
    double *det = _DetJacobian;
    Matrix *adj = _AdjJacobian;
    auto domain = _SubDomains.begin();
    for (int n = 0; n < mffd->NumberOfLevels(); ++n) {
      if (mffd->LocalTransformationIsActive(n)) {
        ffd = mffd->GetLocalTransformation(n);
        Body::Run(this, ffd, *domain, det, adj, penalty);
        auto npts = domain->NumberOfPoints();
        det += npts;
        adj += npts;
        ++domain;
      }
    }
  } else if (ffd) {
    Body::Run(this, ffd, _SubDomains[0], _DetJacobian, _AdjJacobian, penalty);
  }

  MIRTK_DEBUG_TIMING(2, (this->HasName() ? this->Name() : this->DefaultName()) << " evaluation");
  return penalty / _SubDomains.size();
}

// -----------------------------------------------------------------------------
// TODO: A more generic implementation would compute the voxel-wise gradient
//       for each point of the discrete target image domain and then use
//       the ParametricGradient function of the transformation similar to
//       the ImageSimilarity::EvaluateGradient implementation. This could be
//       implemented in the TransformationConstraint base class.
void JacobianConstraint::EvaluateGradient(double *gradient, double, double weight)
{
  typedef JacobianConstraintEvaluateGradient Body;

  if (_SubDomains.empty()) {
    return;
  }

  MIRTK_START_TIMING();

  const FreeFormTransformation   *ffd  = FFD();
  const MultiLevelTransformation *mffd = MFFD();
  const BSplineFreeFormTransformation3D *bffd3d;
  const BSplineFreeFormTransformationSV *svffd;
  double w;

  if (mffd) {
    double *det = _DetJacobian;
    Matrix *adj = _AdjJacobian;
    auto domain = _SubDomains.begin();
    for (int n = 0; n < mffd->NumberOfLevels(); ++n) {
      if (mffd->LocalTransformationIsActive(n)) {
        auto npts = domain->NumberOfPoints();
        #if _DIVIDE_GRADIENT_BY_NPOINTS
          w = weight / npts;
        #else
          w = weight;
        #endif
        ffd = mffd->GetLocalTransformation(n);
        bffd3d = dynamic_cast<const BSplineFreeFormTransformation3D *>(ffd);
        if (bffd3d) {
          if (_ConstrainJacobianOfParameterization) {
            svffd = nullptr;
          } else {
            svffd = dynamic_cast<const BSplineFreeFormTransformationSV *>(ffd);
          }
          if (svffd) {
            // Penalize Jacobian of displacement generated by SVF
            Body::Run(this, ffd, *domain, det, adj, gradient, w);
          } else if (_SubDomain == SD_Lattice) {
            typedef BSplineFFD3DJacobianConstraintEvaluateGradientLattice Body;
            Body::Run(this, bffd3d, *domain, det, adj, gradient, w);
          } else {
            typedef BSplineFFD3DJacobianConstraintEvaluateGradient Body;
            Body::Run(this, bffd3d, *domain, det, adj, gradient, w);
          }
        } else {
          Body::Run(this, ffd, *domain, det, adj, gradient, w);
        }
        gradient += ffd->NumberOfDOFs();
        det += npts;
        adj += npts;
        ++domain;
      }
    }
  } else if (ffd) {
    #if _DIVIDE_GRADIENT_BY_NPOINTS
      w = weight / _SubDomains[0].NumberOfPoints();
    #else
      w = weight;
    #endif
    bffd3d = dynamic_cast<const BSplineFreeFormTransformation3D *>(ffd);
    if (bffd3d) {
      if (_ConstrainJacobianOfParameterization) {
        svffd = nullptr;
      } else {
        svffd = dynamic_cast<const BSplineFreeFormTransformationSV *>(ffd);
      }
      if (svffd) {
        // Penalize Jacobian of displacement generated by SVF
        Body::Run(this, ffd, _SubDomains[0], _DetJacobian, _AdjJacobian, gradient, w);
      } else if (_SubDomain == SD_Lattice) {
        typedef BSplineFFD3DJacobianConstraintEvaluateGradientLattice Body;
        Body::Run(this, bffd3d, _SubDomains[0], _DetJacobian, _AdjJacobian, gradient, w);
      } else {
        typedef BSplineFFD3DJacobianConstraintEvaluateGradient Body;
        Body::Run(this, bffd3d, _SubDomains[0], _DetJacobian, _AdjJacobian, gradient, w);
      }
    } else {
      Body::Run(this, ffd, _SubDomains[0], _DetJacobian, _AdjJacobian, gradient, w);
    }
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

  const int   sz = 1024;
  char        fname[sz];
  string _prefix = Prefix(p);
  const char  *prefix = _prefix.c_str();

  const FreeFormTransformation   *ffd  = FFD();
  const MultiLevelTransformation *mffd = MFFD();
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
