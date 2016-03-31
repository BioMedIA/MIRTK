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

#include "mirtk/MinJacobianConstraint.h"

#include "mirtk/Math.h"
#include "mirtk/Matrix.h"
#include "mirtk/Parallel.h"
#include "mirtk/Profiling.h"
#include "mirtk/FreeFormTransformation.h"
#include "mirtk/BSplineFreeFormTransformation3D.h"
#include "mirtk/ObjectFactory.h"


namespace mirtk {


// Register energy term with object factory during static initialization
mirtkAutoRegisterEnergyTermMacro(MinJacobianConstraint);


// =============================================================================
// Auxiliaries
// =============================================================================

namespace MinJacobianConstraintUtils {


// -----------------------------------------------------------------------------
/// Body of Evaluate function for parallel execution
struct MinJacobianConstraintEvaluate
{
private:

  const MinJacobianConstraint  *_This;
  const FreeFormTransformation *_FFD;
  const double                 *_DetJacobian;
  double                        _Gamma;
  double                        _Penalty;
  int                           _N;

public:

  MinJacobianConstraintEvaluate() {}

  MinJacobianConstraintEvaluate(MinJacobianConstraintEvaluate &lhs, split)
  :
    _This(lhs._This),
    _FFD(lhs._FFD),
    _DetJacobian(lhs._DetJacobian),
    _Gamma(lhs._Gamma),
    _Penalty(.0), _N(0)
  {}

  void join(MinJacobianConstraintEvaluate &rhs)
  {
    _Penalty += rhs._Penalty;
    _N       += rhs._N;
  }

  void operator()(const blocked_range<int> &re)
  {
    for (int cp = re.begin(); cp != re.end(); ++cp) {
      if (_This->ConstrainPassiveDoFs() || _FFD->IsActive(cp)) {
        if (_DetJacobian[cp] <= _Gamma) {
          _Penalty += pow(_Gamma / _DetJacobian[cp], 2);// - 2.0;
        }
        ++_N;
      }
    }
  }

  static void Run(const MinJacobianConstraint  *obj,
                  const FreeFormTransformation *ffd,
                  const double                 *det,
                  double                        gamma,
                  double                       &penalty,
                  int                          &num)
  {
    MinJacobianConstraintEvaluate body;
    body._This        = obj;
    body._FFD         = ffd;
    body._DetJacobian = det;
    body._Gamma       = gamma;
    body._Penalty     = .0;
    body._N           = 0;
    parallel_reduce(blocked_range<int>(0, ffd->NumberOfCPs()), body);
    penalty += body._Penalty;
    num     += body._N;
  }
};

// -----------------------------------------------------------------------------
/// Body of EvaluateGradient function for parallel execution
struct MinJacobianConstraintEvaluateGradient
{
  const BSplineFreeFormTransformation3D *_FFD;
  const double                          *_DetJacobian;
  const Matrix                          *_AdjJacobian;
  double                                 _Gamma;
  double                                *_Gradient;
  double                                 _Weight;
  bool                                   _ConstrainPassiveDoFs;

  void operator()(const blocked_range3d<int> &re) const
  {
    int    xdof, ydof, zdof, n, cp, i1, j1, k1, i2, j2, k2;
    double pendrv, pengrad[3];
    Matrix detdrv[3];

    // Loop over control points
    for (int ck = re.pages().begin(); ck != re.pages().end(); ++ck)
    for (int cj = re.rows ().begin(); cj != re.rows ().end(); ++cj)
    for (int ci = re.cols ().begin(); ci != re.cols ().end(); ++ci) {
      cp = _FFD->LatticeToIndex(ci, cj, ck);
      if (_ConstrainPassiveDoFs || _FFD->IsActive(cp)) {
        if (_DetJacobian[cp] > _Gamma) continue;
        _FFD->IndexToDOFs(cp, xdof, ydof, zdof);

        i1 = (ci - 1 >         0 ? ci - 1 :             0);
        j1 = (cj - 1 >         0 ? cj - 1 :             0);
        k1 = (ck - 1 >         0 ? ck - 1 :             0);
        i2 = (ci + 1 < _FFD->X() ? ci + 1 : _FFD->X() - 1);
        j2 = (cj + 1 < _FFD->Y() ? cj + 1 : _FFD->Y() - 1);
        k2 = (ck + 1 < _FFD->Z() ? ck + 1 : _FFD->Z() - 1);

        pengrad[0] = pengrad[1] = pengrad[2] = .0;

        for (int nk = k1; nk <= k2; ++nk)
        for (int nj = j1; nj <= j2; ++nj)
        for (int ni = i1; ni <= i2; ++ni) {
          if (ni == ci && nj == cj && nk == ck) continue;
          cp = _FFD->LatticeToIndex(ni, nj, nk);

          // Derivative of penalty term w.r.t. Jacobian determinant
          pendrv = -2.0 * _Gamma / pow(_DetJacobian[cp], 3);

          // Derivative of Jacobian determinant w.r.t. DoFs
          _FFD->JacobianDetDerivative(detdrv, ni - ci, nj - cj, nk - ck);

          // Apply chain rule and Jacobi's formula
          // (cf. https://en.wikipedia.org/wiki/Jacobi's_formula )
          const Matrix &adj = _AdjJacobian[cp];
          pengrad[0] += pendrv * (adj(0, 0) * detdrv[0](0, 0) +
                                  adj(1, 0) * detdrv[0](0, 1) +
                                  adj(2, 0) * detdrv[0](0, 2));
          pengrad[1] += pendrv * (adj(0, 1) * detdrv[1](1, 0) +
                                  adj(1, 1) * detdrv[1](1, 1) +
                                  adj(2, 1) * detdrv[1](1, 2));
          pengrad[2] += pendrv * (adj(0, 2) * detdrv[2](2, 0) +
                                  adj(1, 2) * detdrv[2](2, 1) +
                                  adj(2, 2) * detdrv[2](2, 2));
        }

        n = (i2 - i1 + 1) * (j2 - j1 + 1) * (k2 - k1 + 1) - 1;
        _Gradient[xdof] += _Weight * pengrad[0] / n;
        _Gradient[ydof] += _Weight * pengrad[1] / n;
        _Gradient[zdof] += _Weight * pengrad[2] / n;
      }
    }
  }

  static void Run(const FreeFormTransformation *ffd,
                  const double                 *det,
                  const Matrix                 *adj,
                  double                        gamma,
                  double                       *gradient,
                  double                        weight,
                  bool                          incl_passive)
  {
    MinJacobianConstraintEvaluateGradient body;
    body._FFD = dynamic_cast<const BSplineFreeFormTransformation3D *>(ffd);
    if (body._FFD == NULL) {
      cerr << "MinJacobianConstraint::EvaluateGradient: Only implemented for 3D B-spline FFD" << endl;
      exit(1);
    }
    body._DetJacobian          = det;
    body._AdjJacobian          = adj;
    body._Gamma                = gamma;
    body._Gradient             = gradient;
    body._Weight               = weight;
    body._ConstrainPassiveDoFs = incl_passive;
    parallel_for(blocked_range3d<int>(0, ffd->Z(), 0, ffd->Y(), 0, ffd->X()), body);
  }
};


} // namespace MinJacobianConstraintUtils
using namespace MinJacobianConstraintUtils;

// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
MinJacobianConstraint::MinJacobianConstraint(const char *name)
:
  JacobianConstraint(name), _Gamma(.3)
{
}

// -----------------------------------------------------------------------------
MinJacobianConstraint::~MinJacobianConstraint()
{
}

// =============================================================================
// Parameters
// =============================================================================

// -----------------------------------------------------------------------------
bool MinJacobianConstraint::SetWithoutPrefix(const char *param, const char *value)
{
  if (strcmp(param, "Threshold") == 0 ||
      strcmp(param, "Gamma")     == 0) {
    return FromString(value, _Gamma);
  }
  return JacobianConstraint::SetWithoutPrefix(param, value);
}

// -----------------------------------------------------------------------------
ParameterList MinJacobianConstraint::Parameter() const
{
  ParameterList params = JacobianConstraint::Parameter();
  InsertWithPrefix(params, "Threshold", _Gamma);
  return params;
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
double MinJacobianConstraint::Evaluate()
{
  typedef MinJacobianConstraintEvaluate Body;

  MIRTK_START_TIMING();

  const MultiLevelTransformation *mffd = NULL;
  const FreeFormTransformation   *ffd  = NULL;

  (mffd = MFFD()) || (ffd = FFD());

  double penalty = .0;
  int    ncps    = 0;

  if (mffd) {
    double *det = _DetJacobian;
    for (int n = 0; n < mffd->NumberOfLevels(); ++n) {
      if (mffd->LocalTransformationIsActive(n)) {
        ffd = mffd->GetLocalTransformation(n);
        Body::Run(this, ffd, det, _Gamma, penalty, ncps);
        det += ffd->NumberOfCPs();
      }
    }
  } else if (ffd) {
    Body::Run(this, ffd, _DetJacobian, _Gamma, penalty, ncps);
  }

  if (ncps == 0) return .0;
  MIRTK_DEBUG_TIMING(2, "non-diffeomorphic transformation penalty evaluation");
  return penalty / ncps;
}

// -----------------------------------------------------------------------------
void MinJacobianConstraint::EvaluateGradient(double *gradient, double, double weight)
{
  typedef MinJacobianConstraintEvaluateGradient Body;

  MIRTK_START_TIMING();

  const MultiLevelTransformation *mffd = NULL;
  const FreeFormTransformation   *ffd  = NULL;

  (mffd = MFFD()) || (ffd = FFD());

  int ncps = 0;
  if (_ConstrainPassiveDoFs) {
    if      (ffd)  ncps =  ffd->NumberOfCPs();
    else if (mffd) ncps = mffd->NumberOfCPs();
  } else {
    if      (ffd)  ncps =  ffd->NumberOfActiveCPs();
    else if (mffd) ncps = mffd->NumberOfActiveCPs();
  }
  if (ncps == 0) return;

  if (mffd) {
    double     *det = _DetJacobian;
    Matrix *adj = _AdjJacobian;
    for (int n = 0; n < mffd->NumberOfLevels(); ++n) {
      if (mffd->LocalTransformationIsActive(n)) {
        ffd = mffd->GetLocalTransformation(n);
        Body::Run(ffd, det, adj, _Gamma, gradient, _Weight / ncps, _ConstrainPassiveDoFs);
        det += ffd->NumberOfCPs();
        adj += ffd->NumberOfCPs();
        gradient += ffd->NumberOfDOFs();
      }
    }
  } else if (ffd) {
    Body::Run(ffd, _DetJacobian, _AdjJacobian, _Gamma, gradient, _Weight / ncps, _ConstrainPassiveDoFs);
  }

  MIRTK_DEBUG_TIMING(2, "non-diffeomorphic transformation penalty gradient computation");
}


} // namespace mirtk
