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
#include "mirtk/Profiling.h"
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


/// Allocate Jacobian matrices in neighborhood of active control points
void AllocateJacobianMatrices(Matrix *jac, const FreeFormTransformation *ffd, bool incl_passive_cps)
{
  const int ncps = ffd->NumberOfCPs();
  if (incl_passive_cps) {
    for (int cp = 0; cp < ncps; ++cp) {
      jac[cp].Initialize(3, 3);
    }
  } else {
    // Mark active control points
    int cp, nc;
    Array<bool> mask(ncps);
    for (cp = 0; cp < ncps; ++cp) {
      mask[cp] = ffd->IsActive(cp);
    }
    // Dilate mask to include control points next to an active control point
    for (int ck = 0; ck < ffd->Z(); ++ck)
    for (int cj = 0; cj < ffd->Y(); ++cj)
    for (int ci = 0; ci < ffd->X(); ++ci) {
      cp = ffd->LatticeToIndex(ci, cj, ck);
      for (int nk = ck - 1; nk <= ck + 1; ++nk)
      for (int nj = cj - 1; nj <= cj + 1; ++nj)
      for (int ni = ci - 1; ni <= ci + 1; ++ni) {
        nc = ffd->LatticeToIndex(ni, nj, nk);
        if (0 <= nc && nc < ncps && mask[nc]) {
          mask[cp] = true;
          nj += 2, nk += 2;
          break;
        }
      }
    }
    // Allocate Jacobian matrices for masked control points
    for (cp = 0; cp < ncps; ++cp) {
      if (mask[cp]) jac[cp].Initialize(3, 3);
      else          jac[cp].Clear();
    }
  }
}


/// Allocate Jacobian matrices in neighborhood of active control points
void AllocateJacobianMatrices(const ImageAttributes *domain, Matrix *jac, const FreeFormTransformation *ffd, bool incl_passive_cps)
{
  const int nvox = domain->NumberOfSpatialPoints();
  if (incl_passive_cps) {
    for (int vox = 0; vox < nvox; ++vox) {
      jac[vox].Initialize(3, 3);
    }
  } else {
    Array<bool> mask(nvox, false);
    int cp, i1, j1, k1, i2, j2, k2;
    for (int ck = 0; ck < ffd->Z(); ++ck)
    for (int cj = 0; cj < ffd->Y(); ++cj)
    for (int ci = 0; ci < ffd->X(); ++ci) {
      cp = ffd->LatticeToIndex(ci, cj, ck);
      if (ffd->IsActive(cp) && ffd->BoundingBox(*domain, cp, i1, j1, k1, i2, j2, k2)) {
        for (int k = k1; k <= k2; ++k)
        for (int j = j1; j <= j2; ++j)
        for (int i = i1; i <= i2; ++i) {
          mask[domain->LatticeToIndex(i, j, k)] = true;
        }
      }
    }
    for (int vox = 0; vox < nvox; ++vox) {
      if (mask[vox]) jac[vox].Initialize(3, 3);
      else           jac[vox].Clear();
    }
  }
}


/// Evaluate Jacobian matrices of cubic B-spline FFD at control points
class EvaluateCubicBSplineFFDJacobianAtCPs
{
  const BSplineFreeFormTransformation3D *_FFD;
  Matrix                                *_Jacobian;

public:

  void operator ()(const blocked_range3d<int> &re) const
  {
    int cp;
    for (int ck = re.pages().begin(); ck != re.pages().end(); ++ck)
    for (int cj = re.rows ().begin(); cj != re.rows ().end(); ++cj)
    for (int ci = re.cols ().begin(); ci != re.cols ().end(); ++ci) {
      cp = _FFD->LatticeToIndex(ci, cj, ck);
      Matrix &jac = _Jacobian[cp];
      if (jac.Rows() > 0) {
        _FFD->EvaluateJacobian(jac, ci, cj, ck);
      }
    }
  }

  static void Run(const BSplineFreeFormTransformation3D *ffd, Matrix *jac)
  {
    EvaluateCubicBSplineFFDJacobianAtCPs eval;
    eval._FFD        = ffd;
    eval._Jacobian   = jac;
    parallel_for(blocked_range3d<int>(0, ffd->Z(), 0, ffd->Y(), 0, ffd->X()), eval);
  }
};


/// Evaluate Jacobian matrices of cubic B-spline FFD at image voxels
class EvaluateCubicBSplineFFDJacobianAtVoxels
{
  const ImageAttributes                 *_Domain;
  const Matrix                          *_CoordMap;
  const BSplineFreeFormTransformation3D *_FFD;
  Matrix                                *_Jacobian;

public:

  void operator ()(const blocked_range3d<int> &re) const
  {
    int vox;
    double u, v, w;
    const Matrix &m = *_CoordMap;
    for (int k = re.pages().begin(); k != re.pages().end(); ++k)
    for (int j = re.rows ().begin(); j != re.rows ().end(); ++j) {
      vox = _Domain->LatticeToIndex(re.cols().begin(), j, k);
      for (int i = re.cols().begin(); i != re.cols().end(); ++i, ++vox) {
        Matrix &jac = _Jacobian[vox];
        if (jac.Rows() > 0) {
          u = m(0, 0) * i + m(0, 1) * j + m(0, 2) * k + m(0, 3);
          v = m(1, 0) * i + m(1, 1) * j + m(1, 2) * k + m(1, 3);
          w = m(2, 0) * i + m(2, 1) * j + m(2, 2) * k + m(2, 3);
          _FFD->EvaluateJacobian(jac, u, v, w);
        }
      }
    }
  }

  static void Run(const ImageAttributes *domain, const BSplineFreeFormTransformation3D *ffd, Matrix *jac)
  {
    Matrix coord_map = ffd->Attributes().GetWorldToLatticeMatrix() * domain->GetLatticeToWorldMatrix();
    EvaluateCubicBSplineFFDJacobianAtVoxels eval;
    eval._Domain     = domain;
    eval._CoordMap   = &coord_map;
    eval._FFD        = ffd;
    eval._Jacobian   = jac;
    parallel_for(blocked_range3d<int>(0, domain->Z(), 0, domain->Y(), 0, domain->X()), eval);
  }
};


/// Apply chain rule to obtain 1st order derivatives w.r.t. world coordinates
class ReorientJacobian
{
  Matrix *_Jacobian;
  const Matrix *_Orient;

  inline void Reorient(double &du, double &dv) const
  {
    const Matrix &m = *_Orient;
    double dx = du * m(0, 0) + dv * m(1, 0);
    double dy = du * m(0, 1) + dv * m(1, 1);
    du = dx, dv = dy;
  }

  inline void Reorient(double &du, double &dv, double &dw) const
  {
    const Matrix &m = *_Orient;
    double dx = du * m(0, 0) + dv * m(1, 0) + dw * m(2, 0);
    double dy = du * m(0, 1) + dv * m(1, 1) + dw * m(2, 1);
    double dz = du * m(0, 2) + dv * m(1, 2) + dw * m(2, 2);
    du = dx, dv = dy, dw = dz;
  }

public:

  void operator ()(const blocked_range<int> &re) const
  {
    for (int i = re.begin(); i != re.end(); ++i) {
      Matrix &jac = _Jacobian[i];
      if (jac.Rows() == 3) {
        Reorient(jac(0, 0), jac(0, 1), jac(0, 2));
        Reorient(jac(1, 0), jac(1, 1), jac(1, 2));
        Reorient(jac(2, 0), jac(2, 1), jac(2, 2));
      } else if (jac.Rows() == 2) {
        Reorient(jac(0, 0), jac(0, 1));
        Reorient(jac(1, 0), jac(1, 1));
      }
    }
  }

  static void Run(int n, Matrix *jac, const Matrix &orient)
  {
    ReorientJacobian eval;
    eval._Jacobian = jac;
    eval._Orient   = &orient;
    parallel_for(blocked_range<int>(0, n), eval);
  }
};


/// Remove rotation components from Jacobian matrices
class RemoveRotation
{
  Matrix *_Jacobian;

public:

  void operator ()(const blocked_range<int> &re) const
  {
    Matrix3x3 J;
    for (auto i = re.begin(); i != re.end(); ++i) {
      Matrix &jac = _Jacobian[i];
      if (jac.Rows() > 0) {
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
    }
  }

  static void Run(Array<Matrix> &jac)
  {
    RemoveRotation eval;
    eval._Jacobian = jac.data();
    parallel_for(blocked_range<int>(0, static_cast<int>(jac.size())), eval);
  }
};


/// Evaluate linear elastic energy at control points
class ApproximateLinearElasticEnergy
{
  const FreeFormTransformation *_FFD;
  bool                          _ConstrainPassiveDoFs;

  const Matrix *_Jacobian;
  double        _Mu;
  double        _Lambda;
  double        _Energy;
  int           _Count;

  ApproximateLinearElasticEnergy() : _Energy(0.), _Count(0) {}
  ApproximateLinearElasticEnergy(const ApproximateLinearElasticEnergy &) = default;

public:

  ApproximateLinearElasticEnergy(const ApproximateLinearElasticEnergy &other, split)
  :
    ApproximateLinearElasticEnergy(other)
  {
    _Energy = 0.;
    _Count  = 0;
  }

  void join(const ApproximateLinearElasticEnergy &other)
  {
    _Energy += other._Energy;
    _Count  += other._Count;
  }

  void operator ()(const blocked_range<int> &re)
  {
    double a, b;
    for (auto cp = re.begin(); cp != re.end(); ++cp) {
      if (_ConstrainPassiveDoFs || _FFD->IsActive(cp)) {
        const Matrix &jac = _Jacobian[cp];
        a = pow(jac(0, 0) + jac(0, 0), 2)
          + pow(jac(0, 1) + jac(1, 0), 2) * 2.
          + pow(jac(0, 2) + jac(2, 0), 2) * 2.
          + pow(jac(1, 1) + jac(1, 1), 2)
          + pow(jac(1, 2) + jac(2, 1), 2) * 2.
          + pow(jac(2, 2) + jac(2, 2), 2);
        b = jac(0, 0) + jac(1, 1) + jac(2, 2);
        _Energy += _Mu * a + _Lambda * b * b;
        ++_Count;
      }
    }
  }

  static double Run(const FreeFormTransformation *ffd, const Matrix *jac,
                    double mu, double lambda, bool incl_passive_cps)
  {
    ApproximateLinearElasticEnergy eval;
    eval._FFD                  = ffd;
    eval._ConstrainPassiveDoFs = incl_passive_cps;
    eval._Jacobian             = jac;
    eval._Mu                   = .5 * mu;
    eval._Lambda               = lambda;
    parallel_reduce(blocked_range<int>(0, ffd->NumberOfCPs()), eval);
    return (eval._Count > 0 ? eval._Energy / eval._Count : 0.);
  }
};


/// Evaluate linear elastic energy at image voxels
class EvaluateLinearElasticEnergy
{
  const Matrix *_Jacobian;
  double        _Mu;
  double        _Lambda;
  double        _Energy;
  int           _Count;

  EvaluateLinearElasticEnergy() : _Energy(0.), _Count(0) {}
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

  void operator ()(const blocked_range<int> &re)
  {
    double a, b;
    for (auto vox = re.begin(); vox != re.end(); ++vox) {
      const Matrix &jac = _Jacobian[vox];
      if (jac.Rows() > 0) {
        a = pow(jac(0, 0) + jac(0, 0), 2)
          + pow(jac(0, 1) + jac(1, 0), 2) * 2.
          + pow(jac(0, 2) + jac(2, 0), 2) * 2.
          + pow(jac(1, 1) + jac(1, 1), 2)
          + pow(jac(1, 2) + jac(2, 1), 2) * 2.
          + pow(jac(2, 2) + jac(2, 2), 2);
        b = jac(0, 0) + jac(1, 1) + jac(2, 2);
        _Energy += _Mu * a + _Lambda * b * b;
        ++_Count;
      }
    }
  }

  static double Run(const ImageAttributes *domain, const Matrix *jac, double mu, double lambda)
  {
    EvaluateLinearElasticEnergy eval;
    eval._Jacobian             = jac;
    eval._Mu                   = .5 * mu;
    eval._Lambda               = lambda;
    parallel_reduce(blocked_range<int>(0, domain->NumberOfSpatialPoints()), eval);
    return (eval._Count > 0 ? eval._Energy / eval._Count : 0.);
  }
};


/// Evaluate gradient of linear elastic energy w.r.t. cubic B-spline FFD coefficients
class AddApproximateCubicBSplineFFDGradient
{
  double LookupTable[27][3];

  const BSplineFreeFormTransformation3D *_FFD;
  const Matrix                          *_Jacobian;
  double                                *_Gradient;
  double                                 _Mu;
  double                                 _Lambda;
  bool                                   _ConstrainRotation;
  bool                                   _ConstrainPassiveDoFs;

  /// Initialize lookup table of cubic B-spline first order derivatives
  void InitializeLookupTable(const Matrix *orient)
  {
    // Note: Order of cubic B-spline pieces is *not* B0, B1, B2, B3!
    //       Therefore, LatticeWeights_I[0] corresponds to an offset
    //       of +1, not -1 as may be expected. Hence the reversed
    //       iteration over the support region.

    const double *w[2] = {
      BSplineFreeFormTransformation3D::Kernel::LatticeWeights,
      BSplineFreeFormTransformation3D::Kernel::LatticeWeights_I,
    };

    int n = 0;
    for (int c = 2; c >= 0; --c)
    for (int b = 2; b >= 0; --b)
    for (int a = 2; a >= 0; --a, ++n) {
      double dx = w[1][a] * w[0][b] * w[0][c];
      double dy = w[0][a] * w[1][b] * w[0][c];
      double dz = w[0][a] * w[0][b] * w[1][c];
      if (orient) {
        const Matrix &m = *orient;
        LookupTable[n][0] = dx * m(0, 0) + dy * m(1, 0) + dz * m(2, 0);
        LookupTable[n][1] = dx * m(0, 1) + dy * m(1, 1) + dz * m(2, 1);
        LookupTable[n][2] = dx * m(0, 2) + dy * m(1, 2) + dz * m(2, 2);
      } else {
        LookupTable[n][0] = dx;
        LookupTable[n][1] = dy;
        LookupTable[n][2] = dz;
      }
    }
  }

public:

  void operator ()(const blocked_range3d<int> &re) const
  {
    double div;
    Vector3D<double> g1, g2, dB;
    int cp, nc, n, xdof, ydof, zdof;

    const int ncps = _FFD->NumberOfCPs();
    for (int ck = re.pages().begin(); ck != re.pages().end(); ++ck)
    for (int cj = re.rows ().begin(); cj != re.rows ().end(); ++cj)
    for (int ci = re.cols ().begin(); ci != re.cols ().end(); ++ci) {
      cp = _FFD->LatticeToIndex(ci, cj, ck);
      if (_ConstrainPassiveDoFs || _FFD->IsActive(cp)) {
        // Sum non-zero gradient contributions within 3x3x3 support region of control point
        n = 0, g1 = g2 = 0.;
        if (_ConstrainRotation) {
          for (int nk = ck - 1; nk <= ck + 1; ++nk)
          for (int nj = cj - 1; nj <= cj + 1; ++nj)
          for (int ni = ci - 1; ni <= ci + 1; ++ni, ++n) {
            nc = _FFD->LatticeToIndex(ni, nj, nk);
            if (0 <= nc && nc < ncps) {
              // First order derivatives at this neighboring control point
              const Matrix &jac = _Jacobian[nc];
              // Get pre-computed derivatives of cubic B-spline kernel centered at current control point (index cp)
              dB._x = LookupTable[n][0];
              dB._y = LookupTable[n][1];
              dB._z = LookupTable[n][2];
              // Apply chain rule to compute derivatives of energy terms w.r.t. coefficients
              g1._x += (jac(0, 0) + jac(0, 0)) * dB._x;
              g1._x += (jac(0, 1) + jac(1, 0)) * dB._y;
              g1._x += (jac(0, 2) + jac(2, 0)) * dB._z;
              g1._y += (jac(1, 0) + jac(0, 1)) * dB._x;
              g1._y += (jac(1, 1) + jac(1, 1)) * dB._y;
              g1._y += (jac(1, 2) + jac(2, 1)) * dB._z;
              g1._z += (jac(2, 0) + jac(0, 2)) * dB._x;
              g1._z += (jac(2, 1) + jac(1, 2)) * dB._y;
              g1._z += (jac(2, 2) + jac(2, 2)) * dB._z;
              div = jac(0, 0) + jac(1, 1) + jac(2, 2);
              g2._x += div * dB._x;
              g2._y += div * dB._y;
              g2._z += div * dB._z;
            }
          }
        } else {
          for (int nk = ck - 1; nk <= ck + 1; ++nk)
          for (int nj = cj - 1; nj <= cj + 1; ++nj)
          for (int ni = ci - 1; ni <= ci + 1; ++ni, ++n) {
            nc = _FFD->LatticeToIndex(ni, nj, nk);
            if (0 <= nc && nc < ncps) {
              // First order derivatives at this neighboring control point
              const Matrix &jac = _Jacobian[nc];
              // Get pre-computed derivatives of cubic B-spline kernel centered at current control point (index cp)
              dB._x = LookupTable[n][0];
              dB._y = LookupTable[n][1];
              dB._z = LookupTable[n][2];
              // Apply chain rule to compute derivatives of energy terms w.r.t. coefficients
              g1._x += 2. * jac(0, 0) * dB._x;
              g1._y += 2. * jac(1, 1) * dB._y;
              g1._z += 2. * jac(2, 2) * dB._z;
              div = jac(0, 0) + jac(1, 1) + jac(2, 2);
              g2._x += div * dB._x;
              g2._y += div * dB._y;
              g2._z += div * dB._z;
            }
          }
        }
        _FFD->IndexToDOFs(cp, xdof, ydof, zdof);
        _Gradient[xdof] += _Mu * g1._x + _Lambda * g2._x;
        _Gradient[ydof] += _Mu * g1._y + _Lambda * g2._y;
        _Gradient[zdof] += _Mu * g1._z + _Lambda * g2._z;
      }
    }
  }

  static void Run(double *gradient,
                  const BSplineFreeFormTransformation3D *ffd, const Matrix *jac,
                  double mu, double lambda, bool incl_rotation, bool incl_passive_cps,
                  const Matrix *orient)
  {
    const int ncps = (incl_passive_cps ? ffd->NumberOfCPs() : ffd->NumberOfActiveCPs());
    AddApproximateCubicBSplineFFDGradient eval;
    eval._FFD                  = ffd;
    eval._Jacobian             = jac;
    eval._Mu                   = 2. * mu / ncps;
    eval._Lambda               = 2. * lambda / ncps;
    eval._ConstrainRotation    = incl_rotation;
    eval._ConstrainPassiveDoFs = incl_passive_cps;
    eval._Gradient             = gradient;
    eval.InitializeLookupTable(orient);
    parallel_for(blocked_range3d<int>(0, ffd->Z(), 0, ffd->Y(), 0, ffd->X()), eval);
  }
};


/// Evaluate gradient of linear elastic energy w.r.t. cubic B-spline FFD coefficients
class AddCubicBSplineFFDGradient
{
  typedef BSplineFreeFormTransformation3D::Kernel Kernel;

  const ImageAttributes                 *_Domain;
  const Matrix                          *_CoordMap;
  const BSplineFreeFormTransformation3D *_FFD;
  const Matrix                          *_Jacobian;
  const Matrix                          *_Orient;
  double                                *_Gradient;
  double                                 _Mu;
  double                                 _Lambda;
  bool                                   _ConstrainRotation;
  bool                                   _ConstrainPassiveDoFs;

public:

  void operator ()(const blocked_range3d<int> &re) const
  {
    Vector3D<double> g1, g2, dB, r;
    int cp, xdof, ydof, zdof, i1, j1, k1, i2, j2, k2, vox;
    double dx, dy, dz, div;

    const Matrix &m = *_CoordMap;
    for (int ck = re.pages().begin(); ck != re.pages().end(); ++ck)
    for (int cj = re.rows ().begin(); cj != re.rows ().end(); ++cj)
    for (int ci = re.cols ().begin(); ci != re.cols ().end(); ++ci) {
      cp = _FFD->LatticeToIndex(ci, cj, ck);
      if ((_ConstrainPassiveDoFs || _FFD->IsActive(cp)) && _FFD->BoundingBox(*_Domain, cp, i1, j1, k1, i2, j2, k2)) {
        // Sum non-zero gradient contributions over support region of control point
        g1 = g2 = 0.;
        if (_ConstrainRotation) {
          for (int k = k1; k <= k2; ++k)
          for (int j = j1; j <= j2; ++j)
          for (int i = i1; i <= i2; ++i) {
            // First order spatial derivatives of displacement field
            vox = _Domain->LatticeToIndex(i, j, k);
            const Matrix &jac = _Jacobian[vox];
            // Derivatives of cubic B-spline kernel centered at this control point
            r._x = m(0, 0) * i + m(0, 1) * j + m(0, 2) * k + m(0, 3) - ci;
            r._y = m(1, 0) * i + m(1, 1) * j + m(1, 2) * k + m(1, 3) - cj;
            r._z = m(2, 0) * i + m(2, 1) * j + m(2, 2) * k + m(2, 3) - ck;
            dB._x = Kernel::B_I(r._x) * Kernel::B(r._y) * Kernel::B(r._z);
            dB._y = Kernel::B(r._x) * Kernel::B_I(r._y) * Kernel::B(r._z);
            dB._z = Kernel::B(r._x) * Kernel::B(r._y) * Kernel::B_I(r._z);
            if (_Orient) {
              const Matrix &m = *_Orient;
              dx = dB._x * m(0, 0) + dB._y * m(1, 0) + dB._z * m(2, 0);
              dy = dB._x * m(0, 1) + dB._y * m(1, 1) + dB._z * m(2, 1);
              dz = dB._x * m(0, 2) + dB._y * m(1, 2) + dB._z * m(2, 2);
              dB._x = dx, dB._y = dy, dB._z = dz;
            }
            // Apply chain rule to compute derivatives of energy terms w.r.t. coefficients
            g1._x += (jac(0, 0) + jac(0, 0)) * dB._x;
            g1._x += (jac(0, 1) + jac(1, 0)) * dB._y;
            g1._x += (jac(0, 2) + jac(2, 0)) * dB._z;
            g1._y += (jac(1, 0) + jac(0, 1)) * dB._x;
            g1._y += (jac(1, 1) + jac(1, 1)) * dB._y;
            g1._y += (jac(1, 2) + jac(2, 1)) * dB._z;
            g1._z += (jac(2, 0) + jac(0, 2)) * dB._x;
            g1._z += (jac(2, 1) + jac(1, 2)) * dB._y;
            g1._z += (jac(2, 2) + jac(2, 2)) * dB._z;
            div = jac(0, 0) + jac(1, 1) + jac(2, 2);
            g2._x += div * dB._x;
            g2._y += div * dB._y;
            g2._z += div * dB._z;
          }
        } else {
          for (int k = k1; k <= k2; ++k)
          for (int j = j1; j <= j2; ++j)
          for (int i = i1; i <= i2; ++i) {
            // First order spatial derivatives of displacement field
            vox = _Domain->LatticeToIndex(i, j, k);
            const Matrix &jac = _Jacobian[vox];
            // Derivatives of cubic B-spline kernel centered at this control point
            r._x = m(0, 0) * i + m(0, 1) * j + m(0, 2) * k + m(0, 3) - ci;
            r._y = m(1, 0) * i + m(1, 1) * j + m(1, 2) * k + m(1, 3) - cj;
            r._z = m(2, 0) * i + m(2, 1) * j + m(2, 2) * k + m(2, 3) - ck;
            dB._x = Kernel::B_I(r._x) * Kernel::B(r._y) * Kernel::B(r._z);
            dB._y = Kernel::B(r._x) * Kernel::B_I(r._y) * Kernel::B(r._z);
            dB._z = Kernel::B(r._x) * Kernel::B(r._y) * Kernel::B_I(r._z);
            if (_Orient) {
              const Matrix &m = *_Orient;
              dx = dB._x * m(0, 0) + dB._y * m(1, 0) + dB._z * m(2, 0);
              dy = dB._x * m(0, 1) + dB._y * m(1, 1) + dB._z * m(2, 1);
              dz = dB._x * m(0, 2) + dB._y * m(1, 2) + dB._z * m(2, 2);
              dB._x = dx, dB._y = dy, dB._z = dz;
            }
            // Apply chain rule to compute derivatives of energy terms w.r.t. coefficients
            g1._x += (jac(0, 0) + jac(0, 0)) * dB._x;
            g1._y += (jac(1, 1) + jac(1, 1)) * dB._y;
            g1._z += (jac(2, 2) + jac(2, 2)) * dB._z;
            div = jac(0, 0) + jac(1, 1) + jac(2, 2);
            g2._x += div * dB._x;
            g2._y += div * dB._y;
            g2._z += div * dB._z;
          }
        }
        _FFD->IndexToDOFs(cp, xdof, ydof, zdof);
        _Gradient[xdof] += _Mu * g1._x + _Lambda * g2._x;
        _Gradient[ydof] += _Mu * g1._y + _Lambda * g2._y;
        _Gradient[zdof] += _Mu * g1._z + _Lambda * g2._z;
      }
    }
  }

  static void Run(double *gradient, const ImageAttributes *domain,
                  const BSplineFreeFormTransformation3D *ffd, const Matrix *jac,
                  double mu, double lambda, bool incl_rotation, bool incl_passive_cps,
                  const Matrix *orient)
  {
    const int npts = domain->NumberOfSpatialPoints();
    Matrix coord_map = ffd->Attributes().GetWorldToLatticeMatrix() * domain->GetLatticeToWorldMatrix();
    AddCubicBSplineFFDGradient eval;
    eval._Domain               = domain;
    eval._CoordMap             = &coord_map;
    eval._FFD                  = ffd;
    eval._Jacobian             = jac;
    eval._Orient               = orient;
    eval._Mu                   = 2. * mu / npts;
    eval._Lambda               = 2. * lambda / npts;
    eval._ConstrainRotation    = incl_rotation;
    eval._ConstrainPassiveDoFs = incl_passive_cps;
    eval._Gradient             = gradient;
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
  _Approximate(true),
  _WithRespectToWorld(true),
  _UseLatticeSpacing(false),
  _ConstrainRotation(true),
  _Lambda(0.), _Mu(1.)
{
  _ParameterPrefix.push_back("Elastic energy ");
  _ParameterPrefix.push_back("Linear energy ");
  _ParameterPrefix.push_back("Linear elasticity ");
  _ParameterPrefix.push_back("Elasticity ");
}

// =============================================================================
// Parameters
// =============================================================================

// -----------------------------------------------------------------------------
bool LinearElasticityConstraint::SetWithoutPrefix(const char *param, const char *value)
{
  if (strcmp(param, "Approximate") == 0) {
    return FromString(value, _Approximate);
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
  if (strcmp(param, "Rotation") == 0 || strcmp(param, "Constrain rotation") == 0) {
    return FromString(value, _ConstrainRotation);
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
  InsertWithPrefix(params, "Approximate", _Approximate);
  InsertWithPrefix(params, "W.r.t. world", _WithRespectToWorld);
  InsertWithPrefix(params, "Use lattice spacing", _UseLatticeSpacing);
  InsertWithPrefix(params, "Rotation", _ConstrainRotation);
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
  if (IsZero(_Weight) || (IsZero(_Mu) && IsZero(_Lambda))) {
    _Jacobian.clear();
    return;
  }

  // Allocate memory for Jacobian matrices
  const MultiLevelTransformation *mffd = nullptr;
  const FreeFormTransformation   *ffd  = nullptr;

  (mffd = MFFD()) || (ffd = FFD());

  if (_Approximate) {
    if (mffd) {
      int n = 0;
      for (int l = 0; l < mffd->NumberOfLevels(); ++l) {
        if (!mffd->LocalTransformationIsActive(l)) continue;
        ffd = mffd->GetLocalTransformation(l);
        n += ffd->NumberOfCPs();
      }
      _Jacobian.resize(n);
      Matrix *jac = _Jacobian.data();
      for (int l = 0; l < mffd->NumberOfLevels(); ++l) {
        if (!mffd->LocalTransformationIsActive(l)) continue;
        ffd = mffd->GetLocalTransformation(l);
        AllocateJacobianMatrices(jac, ffd, _ConstrainPassiveDoFs);
        jac += ffd->NumberOfCPs();
      }
    } else if (ffd) {
      _Jacobian.resize(ffd->NumberOfCPs());
      AllocateJacobianMatrices(_Jacobian.data(), ffd, _ConstrainPassiveDoFs);
    }
  } else {
    if (mffd) {
      int n = 0;
      for (int l = 0; l < mffd->NumberOfLevels(); ++l) {
        if (!mffd->LocalTransformationIsActive(l)) continue;
        n += _Domain.NumberOfSpatialPoints();
      }
      _Jacobian.resize(n);
      Matrix *jac = _Jacobian.data();
      for (int l = 0; l < mffd->NumberOfLevels(); ++l) {
        if (!mffd->LocalTransformationIsActive(l)) continue;
        ffd = mffd->GetLocalTransformation(l);
        AllocateJacobianMatrices(&_Domain, jac, ffd, _ConstrainPassiveDoFs);
        jac += _Domain.NumberOfSpatialPoints();
      }
    } else if (ffd) {
      _Jacobian.resize(_Domain.NumberOfSpatialPoints());
      AllocateJacobianMatrices(&_Domain, _Jacobian.data(), ffd, _ConstrainPassiveDoFs);
    }

  }
}

// -----------------------------------------------------------------------------
void LinearElasticityConstraint::Update(bool gradient)
{
  // Update base class
  TransformationConstraint::Update(gradient);
  if (IsZero(_Weight) || (IsZero(_Mu) && IsZero(_Lambda))) return;
  MIRTK_START_TIMING();

  // Evaluate Jacobian matrices
  int npts;
  const MultiLevelTransformation *mffd = nullptr;
  const FreeFormTransformation   *ffd  = nullptr;

  (mffd = MFFD()) || (ffd = FFD());

  Matrix *jac = _Jacobian.data();
  if (mffd) {
    for (int l = 0; l < mffd->NumberOfLevels(); ++l) {
      if (!mffd->LocalTransformationIsActive(l)) continue;
      ffd = mffd->GetLocalTransformation(l);
      auto bffd = dynamic_cast<const BSplineFreeFormTransformation3D *>(ffd);
      if (!bffd) {
        Throw(ERR_NotImplemented, __FUNCTION__, "Currently only implemented for 3D cubic BSpline (SV)FFD");
      }
      if (_Approximate) {
        npts = bffd->NumberOfCPs();
        EvaluateCubicBSplineFFDJacobianAtCPs::Run(bffd, jac);
      } else {
        npts = _Domain.NumberOfSpatialPoints();
        EvaluateCubicBSplineFFDJacobianAtVoxels::Run(&_Domain, bffd, jac);
      }
      if (_WithRespectToWorld) {
        Matrix orient;
        if (_UseLatticeSpacing) orient = bffd->Attributes().GetWorldToLatticeMatrix();
        else                    orient = bffd->Attributes().GetWorldToLatticeOrientation();
        ReorientJacobian::Run(npts, jac, orient);
      }
      jac += npts;
    }
  } else if (ffd) {
    auto bffd = dynamic_cast<const BSplineFreeFormTransformation3D *>(ffd);
    if (!bffd) {
      Throw(ERR_NotImplemented, __FUNCTION__, "Currently only implemented for 3D cubic BSpline (SV)FFD");
    }
    if (_Approximate) {
      npts = bffd->NumberOfCPs();
      EvaluateCubicBSplineFFDJacobianAtCPs::Run(bffd, jac);
    } else {
      npts = _Domain.NumberOfSpatialPoints();
      EvaluateCubicBSplineFFDJacobianAtVoxels::Run(&_Domain, bffd, jac);
    }
    if (_WithRespectToWorld) {
      Matrix orient;
      if (_UseLatticeSpacing) orient = bffd->Attributes().GetWorldToLatticeMatrix();
      else                    orient = bffd->Attributes().GetWorldToLatticeOrientation();
      ReorientJacobian::Run(npts, jac, orient);
    }
  }

  // Remove rotation components
  if (!_ConstrainRotation) {
    RemoveRotation::Run(_Jacobian);
  }

  MIRTK_DEBUG_TIMING(2, "update of linear energy term");
}

// -----------------------------------------------------------------------------
double LinearElasticityConstraint::Evaluate()
{
  if (IsZero(_Weight) || (IsZero(_Mu) && IsZero(_Lambda))) return 0.;
  MIRTK_START_TIMING();

  double energy = 0.;

  const MultiLevelTransformation *mffd = nullptr;
  const FreeFormTransformation   *ffd  = nullptr;

  (mffd = MFFD()) || (ffd = FFD());

  const Matrix *jac = _Jacobian.data();
  if (mffd) {
    for (int l = 0; l < mffd->NumberOfLevels(); ++l) {
      if (!mffd->LocalTransformationIsActive(l)) continue;
      ffd = mffd->GetLocalTransformation(l);
      if (_Approximate) {
        energy += ApproximateLinearElasticEnergy::Run(ffd, jac, _Mu, _Lambda, _ConstrainPassiveDoFs);
        jac += ffd->NumberOfCPs();
      } else {
        energy += EvaluateLinearElasticEnergy::Run(&_Domain, jac, _Mu, _Lambda);
        jac += _Domain.NumberOfSpatialPoints();
      }
    }
  } else if (ffd) {
    if (_Approximate) {
      energy = ApproximateLinearElasticEnergy::Run(ffd, jac, _Mu, _Lambda, _ConstrainPassiveDoFs);
    } else {
      energy = EvaluateLinearElasticEnergy::Run(&_Domain, jac, _Mu, _Lambda);
    }
  }

  MIRTK_DEBUG_TIMING(2, "evaluation of linear energy");
  return energy;
}

// -----------------------------------------------------------------------------
void LinearElasticityConstraint::EvaluateGradient(double *gradient, double, double weight)
{
  if (IsZero(_Weight) || (IsZero(_Mu) && IsZero(_Lambda))) return;
  MIRTK_START_TIMING();

  const MultiLevelTransformation *mffd = nullptr;
  const FreeFormTransformation   *ffd  = nullptr;

  (mffd = MFFD()) || (ffd = FFD());

  double mu     = weight * _Mu;
  double lambda = weight * _Lambda;

  const Matrix *jac = _Jacobian.data();
  if (mffd) {
    for (int l = 0; l < mffd->NumberOfLevels(); ++l) {
      if (!mffd->LocalTransformationIsActive(l)) continue;
      ffd = mffd->GetLocalTransformation(l);
      auto bffd = dynamic_cast<const BSplineFreeFormTransformation3D *>(ffd);
      if (!bffd) {
        Throw(ERR_NotImplemented, __FUNCTION__, "Currently only implemented for 3D cubic BSpline (SV)FFD");
      }
      UniquePtr<Matrix> orient;
      if (_WithRespectToWorld) {
        if (_UseLatticeSpacing) orient.reset(new Matrix(bffd->Attributes().GetWorldToLatticeMatrix()));
        else                    orient.reset(new Matrix(bffd->Attributes().GetWorldToLatticeOrientation()));
      }
      if (_Approximate) {
        AddApproximateCubicBSplineFFDGradient::Run(gradient, bffd, jac, mu, lambda, _ConstrainRotation, _ConstrainPassiveDoFs, orient.get());
        jac += ffd->NumberOfCPs();
      } else {
        AddCubicBSplineFFDGradient::Run(gradient, &_Domain, bffd, jac, mu, lambda, _ConstrainRotation, _ConstrainPassiveDoFs, orient.get());
        jac += _Domain.NumberOfSpatialPoints();
      }
      gradient += ffd->NumberOfDOFs();
    }
  } else if (ffd) {
    auto bffd = dynamic_cast<const BSplineFreeFormTransformation3D *>(ffd);
    if (!bffd) {
      Throw(ERR_NotImplemented, __FUNCTION__, "Currently only implemented for 3D cubic BSpline (SV)FFD");
    }
    UniquePtr<Matrix> orient;
    if (_WithRespectToWorld) {
      if (_UseLatticeSpacing) orient.reset(new Matrix(bffd->Attributes().GetWorldToLatticeMatrix()));
      else                    orient.reset(new Matrix(bffd->Attributes().GetWorldToLatticeOrientation()));
    }
    if (_Approximate) {
      AddApproximateCubicBSplineFFDGradient::Run(gradient, bffd, jac, mu, lambda, _ConstrainRotation, _ConstrainPassiveDoFs, orient.get());
    } else {
      AddCubicBSplineFFDGradient::Run(gradient, &_Domain, bffd, jac, mu, lambda, _ConstrainRotation, _ConstrainPassiveDoFs, orient.get());
    }
  }

  MIRTK_DEBUG_TIMING(2, "linear energy gradient computation");
}

// =============================================================================
// Debugging
// =============================================================================

// -----------------------------------------------------------------------------
void LinearElasticityConstraint::WriteGradient(const char *p, const char *suffix) const
{
  if (_Jacobian.empty()) return;

  const int   sz = 1024;
  char        fname[sz];
  string _prefix = Prefix(p);
  const char  *prefix = _prefix.c_str();

  const MultiLevelTransformation *mffd = nullptr;
  const FreeFormTransformation   *ffd  = nullptr;
  Array<double> gradient;

  (mffd = MFFD()) || (ffd = FFD());

  const Matrix *jac = _Jacobian.data();
  if (mffd) {
    for (int l = 0; l < mffd->NumberOfLevels(); ++l) {
      if (!mffd->LocalTransformationIsActive(l)) continue;
      ffd = mffd->GetLocalTransformation(l);
      auto bffd = dynamic_cast<const BSplineFreeFormTransformation3D *>(ffd);
      if (bffd) {
        gradient.resize(ffd->NumberOfDOFs());
        memset(gradient.data(), 0, ffd->NumberOfDOFs() * sizeof(double));
        UniquePtr<Matrix> orient;
        if (_WithRespectToWorld) {
          if (_UseLatticeSpacing) orient.reset(new Matrix(bffd->Attributes().GetWorldToLatticeMatrix()));
          else                    orient.reset(new Matrix(bffd->Attributes().GetWorldToLatticeOrientation()));
        }
        if (_Approximate) {
          AddApproximateCubicBSplineFFDGradient::Run(gradient.data(), bffd, jac, _Mu, _Lambda, _ConstrainRotation, _ConstrainPassiveDoFs, orient.get());
          jac += ffd->NumberOfCPs();
        } else {
          AddCubicBSplineFFDGradient::Run(gradient.data(), &_Domain, bffd, jac, _Mu, _Lambda, _ConstrainRotation, _ConstrainPassiveDoFs, orient.get());
          jac += _Domain.NumberOfSpatialPoints();
        }
        if (mffd->NumberOfActiveLevels() == 1) {
          snprintf(fname, sz, "%sgradient%s", prefix, suffix);
        } else {
          snprintf(fname, sz, "%sgradient_of_ffd_at_level_%d%s", prefix, l+1, suffix);
        }
        WriteFFDGradient(fname, ffd, gradient.data());
      }
    }
  } else if (ffd) {
    auto bffd = dynamic_cast<const BSplineFreeFormTransformation3D *>(ffd);
    if (bffd) {
      snprintf(fname, sz, "%sgradient%s", prefix, suffix);
      gradient.resize(ffd->NumberOfDOFs());
      memset(gradient.data(), 0, ffd->NumberOfDOFs() * sizeof(double));
      UniquePtr<Matrix> orient;
      if (_WithRespectToWorld) {
        if (_UseLatticeSpacing) orient.reset(new Matrix(bffd->Attributes().GetWorldToLatticeMatrix()));
        else                    orient.reset(new Matrix(bffd->Attributes().GetWorldToLatticeOrientation()));
      }
      if (_Approximate) {
        AddApproximateCubicBSplineFFDGradient::Run(gradient.data(), bffd, jac, _Mu, _Lambda, _ConstrainRotation, _ConstrainPassiveDoFs, orient.get());
      } else {
        AddCubicBSplineFFDGradient::Run(gradient.data(), &_Domain, bffd, jac, _Mu, _Lambda, _ConstrainRotation, _ConstrainPassiveDoFs, orient.get());
      }
      WriteFFDGradient(fname, ffd, gradient.data());
    }
  }
}


} // namespace mirtk
