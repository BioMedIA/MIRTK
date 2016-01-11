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

#include <mirtkLogJacobianConstraint.h>

#include <mirtkMath.h>
#include <mirtkParallel.h>
#include <mirtkProfiling.h>
#include <mirtkMatrix.h>
#include <mirtkFreeFormTransformation.h>
#include <mirtkBSplineFreeFormTransformation3D.h>
#include <mirtkObjectFactory.h>


namespace mirtk {


// Register energy term with object factory during static initialization
mirtkAutoRegisterEnergyTermMacro(LogJacobianConstraint);


// =============================================================================
// Auxiliaries
// =============================================================================

namespace LogJacobianConstraintUtils {


// -----------------------------------------------------------------------------
struct LogJacobianConstraintEvaluate
{
private:

  const LogJacobianConstraint  *_This;
  const FreeFormTransformation *_FFD;
  const Matrix                 *_AdjJacobian;
  const double                 *_DetJacobian;
  double                        _Penalty;
  int                           _N;

public:

  LogJacobianConstraintEvaluate() {}

  LogJacobianConstraintEvaluate(LogJacobianConstraintEvaluate &lhs, split)
  :
    _This(lhs._This), _FFD(lhs._FFD),
    _AdjJacobian(lhs._AdjJacobian),
    _DetJacobian(lhs._DetJacobian),
    _Penalty(.0), _N(0)
  {}

  void join(LogJacobianConstraintEvaluate &rhs)
  {
    _Penalty += rhs._Penalty;
    _N       += rhs._N;
  }

  void operator()(const blocked_range<int> &re)
  {
    double logdet;
    for (int cp = re.begin(); cp != re.end(); ++cp) {
      if (_This->ConstrainPassiveDoFs() || _FFD->IsActive(cp)) {
        logdet = log(_DetJacobian[cp]);
        _Penalty += logdet * logdet;
        ++_N;
      }
    }
  }

  static void Run(const LogJacobianConstraint  *obj,
                  const FreeFormTransformation *ffd,
                  const double                 *det,
                  const Matrix                 *adj,
                  double                       &penalty,
                  int                          &num)
  {
    LogJacobianConstraintEvaluate body;
    body._This        = obj;
    body._FFD         = ffd;
    body._AdjJacobian = adj;
    body._DetJacobian = det;
    body._Penalty     = .0;
    body._N           = 0;
    parallel_reduce(blocked_range<int>(0, ffd->NumberOfCPs()), body);
    penalty += body._Penalty;
    num     += body._N;
  }
};

// -----------------------------------------------------------------------------
// c.f. Torsten Rohlfing et al. MICCAI'01 (w/o scaling correction)
struct LogJacobianConstraintEvaluateGradient
{
  const BSplineFreeFormTransformation3D *_FFD;
  const double                          *_DetJacobian;
  const Matrix                          *_AdjJacobian;
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
          pendrv = -2.0 * log(_DetJacobian[cp]) / _DetJacobian[cp];

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
                  double                       *gradient,
                  double                        weight,
                  bool                          incl_passive)
  {
    LogJacobianConstraintEvaluateGradient body;
    body._FFD = dynamic_cast<const BSplineFreeFormTransformation3D *>(ffd);
    if (body._FFD == NULL) {
      cerr << "LogJacobianConstraint::EvaluateGradient: Only implemented for 3D B-spline (SV) FFD" << endl;
      exit(1);
    }
    body._DetJacobian          = det;
    body._AdjJacobian          = adj;
    body._Gradient             = gradient;
    body._Weight               = weight;
    body._ConstrainPassiveDoFs = incl_passive;
    parallel_for(blocked_range3d<int>(0, ffd->Z(), 0, ffd->Y(), 0, ffd->X()), body);
  }
};


} // namespace LogJacobianConstraintUtils
using namespace LogJacobianConstraintUtils;

// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
LogJacobianConstraint::LogJacobianConstraint(const char *name)
:
  JacobianConstraint(name)
{
}

// -----------------------------------------------------------------------------
LogJacobianConstraint::~LogJacobianConstraint()
{
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
double LogJacobianConstraint::Evaluate()
{
  typedef LogJacobianConstraintEvaluate Body;

  MIRTK_START_TIMING();

  const MultiLevelTransformation *mffd = NULL;
  const FreeFormTransformation   *ffd  = NULL;

  (mffd = MFFD()) || (ffd = FFD());

  double penalty = .0;
  int    ncps    = 0;

  if (mffd) {
    double     *det = _DetJacobian;
    Matrix *adj = _AdjJacobian;
    for (int n = 0; n < mffd->NumberOfLevels(); ++n) {
      if (mffd->LocalTransformationIsActive(n)) {
        ffd = mffd->GetLocalTransformation(n);
        Body::Run(this, ffd, det, adj, penalty, ncps);
        det += ffd->NumberOfCPs();
        adj += ffd->NumberOfCPs();
      }
    }
  } else if (ffd) {
    Body::Run(this, ffd, _DetJacobian, _AdjJacobian, penalty, ncps);
  }

  if (ncps == 0) return .0;
  MIRTK_DEBUG_TIMING(2, "log Jacobian penalty evaluation");
  return penalty / ncps;
}

// -----------------------------------------------------------------------------
void LogJacobianConstraint::EvaluateGradient(double *gradient, double, double weight)
{
  typedef LogJacobianConstraintEvaluateGradient Body;

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
    double *det = _DetJacobian;
    Matrix *adj = _AdjJacobian;
    for (int n = 0; n < mffd->NumberOfLevels(); ++n) {
      if (mffd->LocalTransformationIsActive(n)) {
        ffd = mffd->GetLocalTransformation(n);
        Body::Run(ffd, det, adj, gradient, _Weight / ncps, _ConstrainPassiveDoFs);
        det += ffd->NumberOfCPs();
        adj += ffd->NumberOfCPs();
        gradient += ffd->NumberOfDOFs();
      }
    }
  } else if (ffd) {
    Body::Run(ffd, _DetJacobian, _AdjJacobian, gradient, _Weight / ncps, _ConstrainPassiveDoFs);
  }

  MIRTK_DEBUG_TIMING(2, "log Jacobian penalty gradient computation");
}


} // namespace mirtk
