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
#include "mirtk/MultiLevelTransformation.h"
#include "mirtk/Profiling.h"

#include "mirtk/CommonExport.h" 


namespace mirtk {


// global "debug" flag (cf. mirtk/Options.h)
MIRTK_Common_EXPORT extern int debug;


// =============================================================================
// Auxiliaries
// =============================================================================

namespace {


// -----------------------------------------------------------------------------
/// Body of Update function for parallel execution
struct JacobianConstraintUpdate
{
private:

  const JacobianConstraint     *_This;
  const FreeFormTransformation *_FFD;
  Matrix                       *_AdjJacobian;
  double                       *_DetJacobian;

public:

  void operator()(const blocked_range<int> &re) const
  {
    int    i, j, k, l;
    double x, y, z, t;
    for (int cp = re.begin(); cp != re.end(); ++cp) {
      _FFD->IndexToLattice(cp, i, j, k, l);
      x = i, y = j, z = k, t = _FFD->LatticeToTime(l);
      _FFD->LatticeToWorld(x, y, z);
      _DetJacobian[cp] = _This->Jacobian(_FFD, x, y, z, t, _AdjJacobian[cp]);
    }
  }

  static void Run(const JacobianConstraint     *obj,
                  const FreeFormTransformation *ffd,
                  double                       *det,
                  Matrix                       *adj)
  {
    JacobianConstraintUpdate body;
    body._This        = obj;
    body._FFD         = ffd;
    body._AdjJacobian = adj;
    body._DetJacobian = det;
    blocked_range<int> cps(0, ffd->NumberOfCPs());
    parallel_for(cps, body);
  }
};

// -----------------------------------------------------------------------------
struct JacobianConstraintEvaluate
{
private:

  const JacobianConstraint     *_This;
  const FreeFormTransformation *_FFD;
  const Matrix                 *_AdjJacobian;
  const double                 *_DetJacobian;
  double                        _Penalty;
  int                           _N;

public:

  JacobianConstraintEvaluate() {}

  JacobianConstraintEvaluate(JacobianConstraintEvaluate &lhs, split)
  :
    _This(lhs._This), _FFD(lhs._FFD),
    _AdjJacobian(lhs._AdjJacobian),
    _DetJacobian(lhs._DetJacobian),
    _Penalty(0.), _N(0)
  {}

  void join(JacobianConstraintEvaluate &rhs)
  {
    _Penalty += rhs._Penalty;
    _N       += rhs._N;
  }

  void operator()(const blocked_range<int> &re)
  {
    double logdet;
    for (int cp = re.begin(); cp != re.end(); ++cp) {
      if (_This->ConstrainPassiveDoFs() || _FFD->IsActive(cp)) {
        _Penalty += _This->Penalty(_DetJacobian[cp]);
        ++_N;
      }
    }
  }

  static void Run(const JacobianConstraint     *obj,
                  const FreeFormTransformation *ffd,
                  const double                 *det,
                  const Matrix                 *adj,
                  double                       &penalty,
                  int                          &num)
  {
    JacobianConstraintEvaluate body;
    body._This        = obj;
    body._FFD         = ffd;
    body._AdjJacobian = adj;
    body._DetJacobian = det;
    body._Penalty     = 0.;
    body._N           = 0;
    parallel_reduce(blocked_range<int>(0, ffd->NumberOfCPs()), body);
    penalty += body._Penalty;
    num     += body._N;
  }
};

// -----------------------------------------------------------------------------
struct JacobianConstraintEvaluateGradient
{
  const JacobianConstraint              *_This;
  const BSplineFreeFormTransformation3D *_FFD;
  const double                          *_DetJacobian;
  const Matrix                          *_AdjJacobian;
  double                                *_Gradient;
  double                                 _Weight;

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
      if (_This->ConstrainPassiveDoFs() || _FFD->IsActive(cp)) {
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
          pendrv = _This->DerivativeWrtJacobianDet(_DetJacobian[cp]);

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

        // Note: Not sure why the negative sign...
        pengrad[0] = - pengrad[0];
        pengrad[1] = - pengrad[1];
        pengrad[2] = - pengrad[2];

        n = (i2 - i1 + 1) * (j2 - j1 + 1) * (k2 - k1 + 1) - 1;
        _Gradient[xdof] += _Weight * pengrad[0] / n;
        _Gradient[ydof] += _Weight * pengrad[1] / n;
        _Gradient[zdof] += _Weight * pengrad[2] / n;
      }
    }
  }

  static void Run(const JacobianConstraint     *obj,
                  const FreeFormTransformation *ffd,
                  const double                 *det,
                  const Matrix                 *adj,
                  double                       *gradient,
                  double                        weight)
  {
    JacobianConstraintEvaluateGradient body;
    body._This        = obj;
    body._FFD         = dynamic_cast<const BSplineFreeFormTransformation3D *>(ffd);
    body._DetJacobian = det;
    body._AdjJacobian = adj;
    body._Gradient    = gradient;
    body._Weight      = weight;
    if (body._FFD == nullptr) {
      Throw(ERR_LogicError, "JacobianConstraint::EvaluateGradient", "Only implemented for 3D B-spline (SV) FFD");
    }
    parallel_for(blocked_range3d<int>(0, ffd->Z(), 0, ffd->Y(), 0, ffd->X()), body);
  }
};


} // namespace

// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
JacobianConstraint::JacobianConstraint(const char *name)
:
  TransformationConstraint(name),
  _DetJacobian(nullptr),
  _AdjJacobian(nullptr),
  _NumberOfCPs(0)
{
}

// -----------------------------------------------------------------------------
JacobianConstraint::~JacobianConstraint()
{
  Deallocate(_DetJacobian);
  Deallocate(_AdjJacobian);
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
void JacobianConstraint::Update(bool)
{
  typedef JacobianConstraintUpdate Body;

  const FreeFormTransformation   *ffd  = FFD();
  const MultiLevelTransformation *mffd = MFFD();

  int ncps = 0;
  if      (ffd)  ncps =  ffd->NumberOfCPs();
  else if (mffd) ncps = mffd->NumberOfCPs(true);
  if (ncps == 0) return;

  if (ncps != _NumberOfCPs) {
    Deallocate(_DetJacobian);
    Deallocate(_AdjJacobian);
    _NumberOfCPs = ncps;
    _DetJacobian = Allocate<double>(ncps);
    _AdjJacobian = Allocate<Matrix>(ncps);
  }

  if (mffd) {
    double *det = _DetJacobian;
    Matrix *adj = _AdjJacobian;
    for (int n = 0; n < mffd->NumberOfLevels(); ++n) {
      if (mffd->LocalTransformationIsActive(n)) {
        ffd = mffd->GetLocalTransformation(n);
        Body::Run(this, ffd, det, adj);
        det += ffd->NumberOfCPs();
        adj += ffd->NumberOfCPs();
      }
    }
  } else if (ffd) {
    Body::Run(this, ffd, _DetJacobian, _AdjJacobian);
  }
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
double JacobianConstraint::Evaluate()
{
  typedef JacobianConstraintEvaluate Body;

  MIRTK_START_TIMING();

  const FreeFormTransformation   *ffd  = FFD();
  const MultiLevelTransformation *mffd = MFFD();

  double penalty = .0;
  int    ncps    = 0;

  if (mffd) {
    double *det = _DetJacobian;
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
  MIRTK_DEBUG_TIMING(2, (this->HasName() ? this->Name() : this->DefaultName()) << " evaluation");
  return penalty / ncps;
}

// -----------------------------------------------------------------------------
void JacobianConstraint::EvaluateGradient(double *gradient, double, double weight)
{
  typedef JacobianConstraintEvaluateGradient Body;

  MIRTK_START_TIMING();

  const FreeFormTransformation   *ffd  = FFD();
  const MultiLevelTransformation *mffd = MFFD();

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
        Body::Run(this, ffd, det, adj, gradient, _Weight / ncps);
        det += ffd->NumberOfCPs();
        adj += ffd->NumberOfCPs();
        gradient += ffd->NumberOfDOFs();
      }
    }
  } else if (ffd) {
    Body::Run(this, ffd, _DetJacobian, _AdjJacobian, gradient, _Weight / ncps);
  }

  MIRTK_DEBUG_TIMING(2, (this->HasName() ? this->Name() : this->DefaultName()) << " gradient computation");
}

// =============================================================================
// Debugging
// =============================================================================

// -----------------------------------------------------------------------------
void JacobianConstraint
::WriteJacobian(const char *fname, const FreeFormTransformation *ffd, const double *det) const
{
  GenericImage<double> jacobian(ffd->Attributes(), const_cast<double *>(det));
  jacobian.Write(fname);
}

// -----------------------------------------------------------------------------
void JacobianConstraint::WriteDataSets(const char *p, const char *suffix, bool) const
{
  if (debug < 3) return;

  const int   sz = 1024;
  char        fname[sz];
  string _prefix = Prefix(p);
  const char  *prefix = _prefix.c_str();

  const FreeFormTransformation   *ffd  = FFD();
  const MultiLevelTransformation *mffd = MFFD();

  if (mffd) {
    double *det = _DetJacobian;
    for (int n = 0; n < mffd->NumberOfLevels(); ++n) {
      if (mffd->LocalTransformationIsActive(n)) {
        ffd = mffd->GetLocalTransformation(n);
        if (mffd->NumberOfActiveLevels() == 1) {
          snprintf(fname, sz, "%sjacobian_determinant%s", prefix, suffix);
        } else {
          snprintf(fname, sz, "%sjacobian_determinant_of_ffd_at_level_%d%s", prefix, n+1, suffix);
        }
        WriteJacobian(fname, ffd, det);
        det += ffd->NumberOfCPs();
      }
    }
  } else if (ffd) {
    snprintf(fname, sz, "%sjacobian_determinant%s", prefix, suffix);
    WriteJacobian(fname, ffd, _DetJacobian);
  }
}


} // namespace mirtk
