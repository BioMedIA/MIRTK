/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2017 Imperial College London
 * Copyright 2017 Andreas Schuh
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

#include "mirtk/LinearElasticityConstraint.h"

#include "mirtk/String.h"
#include "mirtk/Memory.h"
#include "mirtk/FreeFormTransformation.h"
#include "mirtk/MultiLevelTransformation.h"
#include "mirtk/BSplineFreeFormTransformation3D.h"


namespace mirtk {


// Register energy term with object factory during static initialization
mirtkAutoRegisterEnergyTermMacro(LinearElasticityConstraint);


// =============================================================================
// Auxiliaries
// =============================================================================

namespace  {


/// Evaluate Jacobian matrices of cubic B-spline FFD at control points
class EvaluateCubicBSplineFFDJacobian
{
  const BSplineFreeFormTransformation3D *_FFD;
  Matrix                                *_Jacobian;
  bool                                   _ConstrainPassiveDoFs;
  bool                                   _RotationInvariant;

public:

  void operator ()(const blocked_range3d<int> &re) const
  {
    Matrix3x3 J;
    FreeFormTransformation::CPStatus status;
    int cp;
    for (int ck = re.pages().begin(); ck != re.pages().end(); ++ck)
    for (int cj = re.rows ().begin(); cj != re.rows ().end(); ++cj)
    for (int ci = re.cols ().begin(); ci != re.cols ().end(); ++ci) {
      cp = _FFD->LatticeToIndex(ci, cj, ck);
      _FFD->GetStatus(cp, status);
      if (_ConstrainPassiveDoFs || (status._x == Active || status._y == Active || status._z == Active)) {
        auto &jac = _Jacobian[cp];
        _FFD->EvaluateJacobian(jac, ci, cj, ck);
        _FFD->JacobianToWorld(jac);
        if (_RotationInvariant) {
          J = jac.To3x3();
          J[0][0] += 1.;
          J[1][1] += 1.;
          J[2][2] += 1.;
          J = J.PolarDecomposition().Inverse() * J;
          J[0][0] -= 1.;
          J[1][1] -= 1.;
          J[2][2] -= 1.;
          jac = J;
        }
      } else {
        _Jacobian[cp].Clear();
      }
    }
  }

  static void Run(const BSplineFreeFormTransformation3D *ffd, Matrix *jac, bool incl_passive_cps, bool rotation_invariant)
  {
    EvaluateCubicBSplineFFDJacobian eval;
    eval._FFD = ffd;
    eval._Jacobian = jac;
    eval._ConstrainPassiveDoFs = incl_passive_cps;
    eval._RotationInvariant = rotation_invariant;
    parallel_for(blocked_range3d<int>(0, ffd->Z(), 0, ffd->Y(), 0, ffd->X()), eval);
  }
};


/// Evaluate linear elastic energy
class EvaluateLinearElasticEnergy
{
  const Matrix *_Jacobian;
  double        _Mu;
  double        _Lambda;
  double        _Energy;
  int           _Count;

  EvaluateLinearElasticEnergy(const Matrix *jac, double mu, double lambda)
  :
    _Jacobian(jac), _Mu(mu), _Lambda(lambda), _Energy(0.), _Count(0)
  {}

  EvaluateLinearElasticEnergy(const EvaluateLinearElasticEnergy &) = default;

public:

  EvaluateLinearElasticEnergy(const EvaluateLinearElasticEnergy &other, split)
  :
    EvaluateLinearElasticEnergy(other)
  {
    _Energy = 0.;
    _Count  = 0;
  }

  void join(const EvaluateLinearElasticEnergy &other)
  {
    _Energy += other._Energy;
    _Count  += other._Count;
  }

  void operator ()(const blocked_range<size_t> &re)
  {
    double a, b;
    for (size_t cp = re.begin(); cp != re.end(); ++cp) {
      const Matrix &jac = _Jacobian[cp];
      const int d = jac.Rows();
      if (d > 0) {
        a = b = 0.;
        for (int j = 0; j < d; ++j) {
          b += jac(j, j);
          for (int k = 0; k < d; ++k) {
            a += pow(.5 * (jac(j, k) + jac(k, j)), 2);
          }
        }
        _Energy += .5 * _Mu * a + _Lambda * b * b;
        ++_Count;
      }
    }
  }

  static double Run(const Matrix *jac, int ncps, double mu, double lambda)
  {
    EvaluateLinearElasticEnergy eval(jac, mu, lambda);
    parallel_reduce(blocked_range<size_t>(0, ncps), eval);
    return (eval._Count > 0 ? eval._Energy / eval._Count : 0.);
  }
};


/// Evaluate gradient of linear elastic energy w.r.t. cubic B-spline FFD coefficients
class AddCubicBSplineFFDGradient
{
  double LookupTable_3D[27][3];

  const BSplineFreeFormTransformation3D *_FFD;
  const Matrix                          *_Jacobian;
  double                                *_Gradient;
  double                                 _Mu;
  double                                 _Lambda;
  bool                                   _ConstrainPassiveDoFs;

  /// Initialize lookup table of cubic B-spline first order derivatives
  void InitializeLookupTable()
  {
    const double *w[2] = {
      BSplineFreeFormTransformation3D::Kernel::LatticeWeights,
      BSplineFreeFormTransformation3D::Kernel::LatticeWeights_I,
    };

    const Matrix &w2l = *_FFD->Attributes()._w2i;

    int n = 0;
    for (int c = 0; c < 3; ++c)
    for (int b = 0; b < 3; ++b)
    for (int a = 0; a < 3; ++a, ++n) {
      double dx = w[1][a] * w[0][b] * w[0][c];
      double dy = w[0][a] * w[1][b] * w[0][c];
      double dz = w[0][a] * w[0][b] * w[1][c];
      LookupTable_3D[n][0] = dx * w2l(0, 0) + dy * w2l(1, 0) + dz * w2l(2, 0);
      LookupTable_3D[n][1] = dx * w2l(0, 1) + dy * w2l(1, 1) + dz * w2l(2, 1);
      LookupTable_3D[n][2] = dx * w2l(0, 2) + dy * w2l(1, 2) + dz * w2l(2, 2);
    }
  }

public:

  void operator ()(const blocked_range3d<int> &re) const
  {
    FreeFormTransformation::CPStatus status;
    Vector3D<double> g1, g2, gb;
    int cp, n, idx, xdof, ydof, zdof;
    double div;

    const Matrix &w2l = *_FFD->Attributes()._w2i;

    const int ncps = _FFD->NumberOfCPs();
    for (int ck = re.pages().begin(); ck != re.pages().end(); ++ck)
    for (int cj = re.rows ().begin(); cj != re.rows ().end(); ++cj)
    for (int ci = re.cols ().begin(); ci != re.cols ().end(); ++ci) {
      cp = _FFD->LatticeToIndex(ci, cj, ck);
      _FFD->GetStatus(cp, status);
      if (_ConstrainPassiveDoFs || (status._x == Active || status._y == Active || status._z == Active)) {
        n = 0, g1 = g2 = 0.;
        for (int k = ck - 1; k <= ck + 1; ++k)
        for (int j = cj - 1; j <= cj + 1; ++j)
        for (int i = ci - 1; i <= ci + 1; ++i, ++n) {
          idx = _FFD->LatticeToIndex(i, j, k);
          if (0 <= idx && idx < ncps) {
            auto &jac = _Jacobian[idx];
            gb._x = LookupTable_3D[n][0];
            gb._y = LookupTable_3D[n][1];
            gb._z = LookupTable_3D[n][2];

            g1._x += gb._x * 2. * jac(0, 0);
            g1._x += gb._y * (jac(0, 1) + jac(1, 0));
            g1._x += gb._z * (jac(0, 2) + jac(2, 0));

            g1._y += gb._x * (jac(0, 1) + jac(1, 0));
            g1._y += gb._y * 2. * jac(1, 1);
            g1._y += gb._z * (jac(2, 1) + jac(1, 2));

            g1._z += gb._x * (jac(0, 2) + jac(2, 0));
            g1._z += gb._y * (jac(2, 1) + jac(1, 2));
            g1._z += gb._z * 2. * jac(2, 2);

            div = jac(0, 0) + jac(1, 1) + jac(2, 2);
            g2._x += gb._x * div;
            g2._y += gb._y * div;
            g2._z += gb._z * div;
          }
        }
        _FFD->IndexToDOFs(cp, xdof, ydof, zdof);
        _Gradient[xdof] -= _Mu * g1._x + 2. * _Lambda * g2._x;
        _Gradient[ydof] -= _Mu * g1._y + 2. * _Lambda * g2._y;
        _Gradient[zdof] -= _Mu * g1._z + 2. * _Lambda * g2._z;
      }
    }
  }

  static void Run(double *gradient,
                  const BSplineFreeFormTransformation3D *ffd, const Matrix *jac,
                  double mu, double lambda, bool incl_passive_cps)
  {
    const int ncps = (incl_passive_cps ? ffd->NumberOfCPs() : ffd->NumberOfActiveCPs());
    AddCubicBSplineFFDGradient eval;
    eval._FFD                  = ffd;
    eval._Jacobian             = jac;
    eval._Mu                   = mu / ncps;
    eval._Lambda               = lambda / ncps;
    eval._ConstrainPassiveDoFs = incl_passive_cps;
    eval._Gradient             = gradient;
    eval.InitializeLookupTable();
    parallel_for(blocked_range3d<int>(0, ffd->Z(), 0, ffd->Y(), 0, ffd->X()), eval);
  }
};


} // anonymous namespace

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
LinearElasticityConstraint::LinearElasticityConstraint(const char *name, double weight)
:
  TransformationConstraint(name, weight),
  _RotationInvariant(false), _Lambda(0.), _Mu(1.)
{
  _ParameterPrefix.push_back("Linear elasticity ");
  _ParameterPrefix.push_back("Linear energy ");
  _ParameterPrefix.push_back("Elasticity ");
}

// =============================================================================
// Parameters
// =============================================================================

// -----------------------------------------------------------------------------
bool LinearElasticityConstraint::SetWithoutPrefix(const char *param, const char *value)
{
  if (strcmp(param, "Rotation invariant") == 0) {
    return FromString(value, _RotationInvariant);
  }
  if (strcmp(param, "Rotation") == 0) {
    if (!FromString(value, _RotationInvariant)) return false;
    _RotationInvariant = !_RotationInvariant;
    return true;
  }
  if (strcmp(param, "Lambda") == 0) {
    return FromString(value, _Lambda);
  }
  if (strcmp(param, "Mu") == 0) {
    return FromString(value, _Mu);
  }
  return TransformationConstraint::SetWithoutPrefix(param, value);
}

// -----------------------------------------------------------------------------
ParameterList LinearElasticityConstraint::Parameter() const
{
  ParameterList params = TransformationConstraint::Parameter();
  InsertWithPrefix(params, "Rotation", !_RotationInvariant);
  InsertWithPrefix(params, "Lambda", _Lambda);
  InsertWithPrefix(params, "Mu", _Mu);
  return params;
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
void LinearElasticityConstraint::Initialize()
{
  // Initialize base class
  TransformationConstraint::Initialize();

  // Allocate memory for Jacobian matrices
  int ncps = 0;

  const MultiLevelTransformation *mffd = nullptr;
  const FreeFormTransformation   *ffd  = nullptr;

  (mffd = MFFD()) || (ffd = FFD());

  if (mffd) {
    for (int l = 0; l < mffd->NumberOfLevels(); ++l) {
      if (!mffd->LocalTransformationIsActive(l)) continue;
      ffd = mffd->GetLocalTransformation(l);
      ncps += ffd->NumberOfCPs();
    }
  } else if (ffd) {
    ncps += ffd->NumberOfCPs();
  }

  _Jacobian.resize(ncps);
}

// -----------------------------------------------------------------------------
void LinearElasticityConstraint::Update(bool gradient)
{
  // Update base class
  TransformationConstraint::Update(gradient);

  // Evaluate Jacobian matrices
  const MultiLevelTransformation *mffd = nullptr;
  const FreeFormTransformation   *ffd  = nullptr;
  const BSplineFreeFormTransformation3D *bffd;

  (mffd = MFFD()) || (ffd = FFD());

  if (mffd) {
    Matrix *jac = _Jacobian.data();
    for (int l = 0; l < mffd->NumberOfLevels(); ++l) {
      if (!mffd->LocalTransformationIsActive(l)) continue;
      ffd = mffd->GetLocalTransformation(l);
      bffd = dynamic_cast<const BSplineFreeFormTransformation3D *>(ffd);
      if (!bffd) {
        Throw(ERR_NotImplemented, __FUNCTION__, "Currently only implemented for cubic BSpline FFD");
      }
      EvaluateCubicBSplineFFDJacobian::Run(bffd, jac, _ConstrainPassiveDoFs, _RotationInvariant);
      jac += ffd->NumberOfCPs();
    }
  } else if (ffd) {
    bffd = dynamic_cast<const BSplineFreeFormTransformation3D *>(ffd);
    if (!bffd) {
      Throw(ERR_NotImplemented, __FUNCTION__, "Currently only implemented for cubic BSpline FFD");
    }
    EvaluateCubicBSplineFFDJacobian::Run(bffd, _Jacobian.data(), _ConstrainPassiveDoFs, _RotationInvariant);
  }
}

// -----------------------------------------------------------------------------
double LinearElasticityConstraint::Evaluate()
{
  double energy = 0.;

  const MultiLevelTransformation *mffd = nullptr;
  const FreeFormTransformation   *ffd  = nullptr;

  (mffd = MFFD()) || (ffd = FFD());

  if (mffd) {
    const Matrix *jac = _Jacobian.data();
    for (int l = 0; l < mffd->NumberOfLevels(); ++l) {
      if (!mffd->LocalTransformationIsActive(l)) continue;
      ffd = mffd->GetLocalTransformation(l);
      energy += EvaluateLinearElasticEnergy::Run(jac, ffd->NumberOfCPs(), _Mu, _Lambda);
      jac += ffd->NumberOfCPs();
    }
  } else if (ffd) {
    energy = EvaluateLinearElasticEnergy::Run(_Jacobian.data(), ffd->NumberOfCPs(), _Mu, _Lambda);
  }

  return energy;
}

// -----------------------------------------------------------------------------
void LinearElasticityConstraint::EvaluateGradient(double *gradient, double, double weight)
{
  const MultiLevelTransformation *mffd = nullptr;
  const FreeFormTransformation   *ffd  = nullptr;
  const BSplineFreeFormTransformation3D *bffd;

  (mffd = MFFD()) || (ffd = FFD());

  double mu     = weight * _Mu;
  double lambda = weight * _Lambda;

  if (mffd) {
    Matrix *jac = _Jacobian.data();
    for (int l = 0; l < mffd->NumberOfLevels(); ++l) {
      if (!mffd->LocalTransformationIsActive(l)) continue;
      ffd = mffd->GetLocalTransformation(l);
      bffd = dynamic_cast<const BSplineFreeFormTransformation3D *>(ffd);
      if (!bffd) {
        Throw(ERR_NotImplemented, __FUNCTION__, "Currently only implemented for cubic BSpline FFD");
      }
      AddCubicBSplineFFDGradient::Run(gradient, bffd, jac, mu, lambda, _ConstrainPassiveDoFs);
      jac += ffd->NumberOfCPs();
      gradient += ffd->NumberOfDOFs();
    }
  } else if (ffd) {
    bffd = dynamic_cast<const BSplineFreeFormTransformation3D *>(ffd);
    if (!bffd) {
      Throw(ERR_NotImplemented, __FUNCTION__, "Currently only implemented for cubic BSpline FFD");
    }
    AddCubicBSplineFFDGradient::Run(gradient, bffd, _Jacobian.data(), mu, lambda, _ConstrainPassiveDoFs);
  }
}

// =============================================================================
// Debugging
// =============================================================================

// -----------------------------------------------------------------------------
void LinearElasticityConstraint::WriteGradient(const char *p, const char *suffix) const
{
  const int   sz = 1024;
  char        fname[sz];
  string _prefix = Prefix(p);
  const char  *prefix = _prefix.c_str();

  const MultiLevelTransformation *mffd = nullptr;
  const FreeFormTransformation   *ffd  = nullptr;
  const BSplineFreeFormTransformation3D *bffd;

  (mffd = MFFD()) || (ffd = FFD());

  if (mffd) {
    const Matrix *jac = _Jacobian.data();
    for (int l = 0; l < mffd->NumberOfLevels(); ++l) {
      if (!mffd->LocalTransformationIsActive(l)) continue;
      ffd = mffd->GetLocalTransformation(l);
      bffd = dynamic_cast<const BSplineFreeFormTransformation3D *>(ffd);
      if (!bffd) {
        Throw(ERR_NotImplemented, __FUNCTION__, "Currently only implemented for cubic BSpline FFD");
      }
      double *gradient = CAllocate<double>(ffd->NumberOfDOFs());
      AddCubicBSplineFFDGradient::Run(gradient, bffd, jac, _Mu, _Lambda, _ConstrainPassiveDoFs);
      if (mffd->NumberOfActiveLevels() == 1) {
        snprintf(fname, sz, "%sgradient%s", prefix, suffix);
      } else {
        snprintf(fname, sz, "%sgradient_of_ffd_at_level_%d%s", prefix, l+1, suffix);
      }
      WriteFFDGradient(fname, ffd, gradient);
      Deallocate(gradient);
      jac += ffd->NumberOfCPs();
    }
  } else if (ffd) {
    bffd = dynamic_cast<const BSplineFreeFormTransformation3D *>(ffd);
    if (!bffd) {
      Throw(ERR_NotImplemented, __FUNCTION__, "Currently only implemented for cubic BSpline FFD");
    }
    snprintf(fname, sz, "%sgradient%s", prefix, suffix);
    double *gradient = CAllocate<double>(ffd->NumberOfDOFs());
    AddCubicBSplineFFDGradient::Run(gradient, bffd, _Jacobian.data(), _Mu, _Lambda, _ConstrainPassiveDoFs);
    WriteFFDGradient(fname, ffd, gradient);
    Deallocate(gradient);
  }
}


} // namespace mirtk
