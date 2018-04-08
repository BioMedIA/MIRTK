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

#include "mirtk/BSplineFreeFormTransformationSV.h"

#include "mirtk/Math.h"
#include "mirtk/Memory.h"
#include "mirtk/Profiling.h"
#include "mirtk/ScalingAndSquaring.h"
#include "mirtk/DisplacementToVelocityFieldBCH.h"
#include "mirtk/VelocityToDisplacementFieldSS.h"
#include "mirtk/LieBracketImageFilter.h"
#include "mirtk/VoxelFunction.h"
#include "mirtk/NaryVoxelFunction.h"
#include "mirtk/ImageToInterpolationCoefficients.h"
#include "mirtk/HomogeneousTransformation.h"
#include "mirtk/MultiLevelTransformation.h"
#include "mirtk/InterpolateImageFunction.hxx"

#include "FreeFormTransformationIntegration.h"


namespace mirtk {


// =============================================================================
// Integration methods
// =============================================================================

// Tolerance of embedded Runge-Kutta methods with automatic step length control
static const double SVFFD_RKTOL = 1.e-3;

MIRTK_FFDIM2(RKE1,   BSplineFreeFormTransformationSV);
MIRTK_FFDIM2(RKEH12, BSplineFreeFormTransformationSV);
MIRTK_FFDIM2(RKE2,   BSplineFreeFormTransformationSV);
MIRTK_FFDIM2(RKH2,   BSplineFreeFormTransformationSV);
MIRTK_FFDIM2(RKBS23, BSplineFreeFormTransformationSV);
MIRTK_FFDIM2(RK4,    BSplineFreeFormTransformationSV);
MIRTK_FFDIM2(RKF45,  BSplineFreeFormTransformationSV);
MIRTK_FFDIM2(RKCK45, BSplineFreeFormTransformationSV);
MIRTK_FFDIM2(RKDP45, BSplineFreeFormTransformationSV);

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
inline double DefaultMaximumScaledVelocity(double dx, double dy, double dz)
{
  double ds = .0;
  if (dx > .0 && (ds == .0 || dx < ds)) ds = dx;
  if (dy > .0 && (ds == .0 || dy < ds)) ds = dy;
  if (dz > .0 && (ds == .0 || dz < ds)) ds = dz;
  return 0.4 * ds;
}

// -----------------------------------------------------------------------------
BSplineFreeFormTransformationSV
::BSplineFreeFormTransformationSV()
:
  _T                (1.0),
  _TimeUnit         (1.0),
  _NumberOfSteps    (64),
  _MaxScaledVelocity(-1.0),
  _IntegrationMethod(FFDIM_FastSS),
  _UseDenseBCHGrid  (false),
  _LieDerivative    (false),
  _NumberOfBCHTerms (4)
{
  _ExtrapolationMode = Extrapolation_NN;
}

// -----------------------------------------------------------------------------
BSplineFreeFormTransformationSV
::BSplineFreeFormTransformationSV(const ImageAttributes &attr,
                                  double dx, double dy, double dz)
:
  _T                (1.0),
  _TimeUnit         (1.0),
  _NumberOfSteps    (64),
  _MaxScaledVelocity(-1.0),
  _IntegrationMethod(FFDIM_FastSS),
  _UseDenseBCHGrid  (false),
  _LieDerivative    (false),
  _NumberOfBCHTerms (4)
{
  _ExtrapolationMode = Extrapolation_NN;
  Initialize(attr, dx, dy, dz);
}

// -----------------------------------------------------------------------------
BSplineFreeFormTransformationSV
::BSplineFreeFormTransformationSV(const BaseImage &target,
                                  double dx, double dy, double dz)
:
  _T                (1.0),
  _TimeUnit         (1.0),
  _NumberOfSteps    (64),
  _MaxScaledVelocity(-1.0),
  _IntegrationMethod(FFDIM_FastSS),
  _UseDenseBCHGrid  (false),
  _LieDerivative    (false),
  _NumberOfBCHTerms (4)
{
  _ExtrapolationMode = Extrapolation_NN;
  Initialize(target.Attributes(), dx, dy, dz);
}

// -----------------------------------------------------------------------------
BSplineFreeFormTransformationSV
::BSplineFreeFormTransformationSV(const GenericImage<double> &image, bool disp)
:
  _T                (1.0),
  _TimeUnit         (1.0),
  _NumberOfSteps    (64),
  _MaxScaledVelocity(-1.0),
  _IntegrationMethod(FFDIM_FastSS),
  _UseDenseBCHGrid  (false),
  _LieDerivative    (false),
  _NumberOfBCHTerms (4)
{
  Initialize(image, disp);
}

// -----------------------------------------------------------------------------
BSplineFreeFormTransformationSV
::BSplineFreeFormTransformationSV(const BSplineFreeFormTransformationSV &ffd)
:
  BSplineFreeFormTransformation3D(ffd),
  _T                (ffd._T),
  _TimeUnit         (ffd._TimeUnit),
  _NumberOfSteps    (ffd._NumberOfSteps),
  _MaxScaledVelocity(ffd._MaxScaledVelocity),
  _IntegrationMethod(ffd._IntegrationMethod),
  _UseDenseBCHGrid  (ffd._UseDenseBCHGrid),
  _LieDerivative    (ffd._LieDerivative),
  _NumberOfBCHTerms (ffd._NumberOfBCHTerms)
{
}

// -----------------------------------------------------------------------------
BSplineFreeFormTransformationSV
::~BSplineFreeFormTransformationSV()
{
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformationSV::Initialize(const ImageAttributes &attr)
{
  ImageAttributes domain = attr;
  domain._t = 1; // SV FFD can be used to register time series, but not regular FFD
  BSplineFreeFormTransformation3D::Initialize(domain);
  if (attr._t > 0 && !AreEqual(attr._dt, 0.)) {
    _TimeUnit = attr._dt;
    _T = attr._t * attr._dt;
  }
  if (_MaxScaledVelocity < .0) _MaxScaledVelocity = DefaultMaximumScaledVelocity(_dx, _dy, _dz);
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformationSV::Initialize(const ImageAttributes &attr,
                                                 double dx, double dy, double dz,
                                                 const Transformation *dof)
{
  // Initialize free-form deformation (for extended image grid)
  //
  // Ensure that for all target voxels the displacement can be recovered
  // without requiring any extrapolation of the velocity field during
  // computation of the trajectory (integration, i.e., exponentiation)
  this->Initialize(ApproximationDomain(attr, dof), dx, dy, dz);

  // Approximate given transformation
  this->ApproximateAsNew(dof);
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformationSV::Subdivide(bool subdivide_x, bool subdivide_y, bool subdivide_z, bool subdivide_t)
{
  BSplineFreeFormTransformation3D::Subdivide(subdivide_x, subdivide_y, subdivide_z, subdivide_t);
  if (_MaxScaledVelocity > .0) _MaxScaledVelocity /= 2;
}

// =============================================================================
// Auxiliary functions
// =============================================================================

// -----------------------------------------------------------------------------
/// Voxel function used to evaluate Lie bracket at each lattice point using 
/// the Lie derivative definition which is based on the Jacobian matrices
/// of the vector fields. Uses nearest neighbor extrapolation of the velocity field.
class SVFFDEvaluateLieBracket : public VoxelFunction
{
private:

  typedef BSplineFreeFormTransformationSV::Vector Vector;
  typedef BSplineFreeFormTransformationSV::Kernel Kernel;

  const BSplineFreeFormTransformationSV *_FFD;

  double                      _tau; ///< Scaling of left/first vector field
  double                      _eta; ///< Scaling of right/second vector field
  const GenericImage<Vector> &_v;   ///< Left/First vector field
  const GenericImage<Vector> &_w;   ///< Right/Second vector field

public:

  // ---------------------------------------------------------------------------
  /// Constructor
  SVFFDEvaluateLieBracket(const BSplineFreeFormTransformationSV *ffd,
                          const GenericImage<Vector>            &v,
                          const GenericImage<Vector>            &w)
  :
    _FFD(ffd), _tau(1.0), _eta(1.0), _v(v), _w(w)
  {}

  // ---------------------------------------------------------------------------
  /// Constructor
  SVFFDEvaluateLieBracket(const BSplineFreeFormTransformationSV *ffd,
                          double                                 tau,
                          const GenericImage<Vector>            &v,
                          double                                 eta,
                          const GenericImage<Vector>            &w)
  :
    _FFD(ffd), _tau(tau), _eta(eta), _v(v), _w(w)
  {}

  // ---------------------------------------------------------------------------
  /// Evaluate velocity at lattice point
  void Evaluate(Vector &x, double s, const GenericImage<Vector> &v, int i, int j, int k) const
  {
    double B_I, B_J, B_K;
    int    II,  JJ,  KK;
  
    x = .0;
    for (int K = k-1; K <= k+1; K++) {
      B_K = Kernel::LatticeWeights[K - (k-1)];
      if      (K <  0)     KK = 0;
      else if (K >= v.Z()) KK = v.Z()-1;
      else                 KK = K;
      for (int J = j-1; J <= j+1; J++) {
        B_J = Kernel::LatticeWeights[J - (j-1)];
        if      (J <  0)     JJ = 0;
        else if (J >= v.Y()) JJ = v.Y() - 1;
        else                 JJ = J;
        for (int I = i-1; I <= i+1; I++) {
          B_I = Kernel::LatticeWeights[I - (i-1)];
          if      (I <  0)     II = 0;
          else if (I >= v.X()) II = v.X() - 1;
          else                 II = I;
          x += B_I * B_J * B_K * s * v(II, JJ, KK);
        }
      }
    }
  }

  // ---------------------------------------------------------------------------
  /// Evaluate Jacobian of velocity field at lattice point
  void Jacobian(Matrix &jac, double s, const GenericImage<Vector> &v, int i, int j, int k) const
  {
    int    I, J, K, II, JJ, KK;
    double B_I, B_J, B_K, B_I_I, B_J_I, B_K_I;
    Vector dx, dy, dz;

    for (K = k-1; K <= k+1; K++) {
      B_K   = Kernel::LatticeWeights  [K - (k-1)];
      B_K_I = Kernel::LatticeWeights_I[K - (k-1)];
      if      (K <  0)     KK = 0;
      else if (K >= v.Z()) KK = v.Z()-1;
      else                 KK = K;
      for (J = j-1; J <= j+1; J++) {
        B_J   = Kernel::LatticeWeights  [J - (j-1)];
        B_J_I = Kernel::LatticeWeights_I[J - (j-1)];
        if      (J <  0)     JJ = 0;
        else if (J >= v.Y()) JJ = v.Y()-1;
        else                 JJ = J;
        for (I = i-1; I <= i+1; I++) {
          B_I   = Kernel::LatticeWeights  [I - (i-1)];
          B_I_I = Kernel::LatticeWeights_I[I - (i-1)];
          if      (I <  0)     II = 0;
          else if (I >= v.X()) II = v.X()-1;
          else                 II = I;
          const Vector &coeff = v(II, JJ, KK);
          dx += B_I_I * B_J   * B_K   * s * coeff;
          dy += B_I   * B_J_I * B_K   * s * coeff;
          dz += B_I   * B_J   * B_K_I * s * coeff;
        }
      }
    }

    jac.Initialize(3, 3);
    jac(0, 0) = dx._x;
    jac(0, 1) = dy._x;
    jac(0, 2) = dz._x;
    jac(1, 0) = dx._y;
    jac(1, 1) = dy._y;
    jac(1, 2) = dz._y;
    jac(2, 0) = dx._z;
    jac(2, 1) = dy._z;
    jac(2, 2) = dz._z;

    _FFD->JacobianToWorld(jac);
  }

  // ---------------------------------------------------------------------------
  /// Compute product of 3x3 matrix and 3D column vector
  static Vector MatrixProduct(const Matrix &jac, const Vector &vel)
  {
    Vector v;
    v._x = jac(0, 0) * vel._x + jac(0, 1) * vel._y + jac(0, 2) * vel._z;
    v._y = jac(1, 0) * vel._x + jac(1, 1) * vel._y + jac(1, 2) * vel._z;
    v._z = jac(2, 0) * vel._x + jac(2, 1) * vel._y + jac(2, 2) * vel._z;
    return v;
  }

  // ---------------------------------------------------------------------------
  /// Evaluate Lie bracket at given lattice point, u = [v, w]
  void operator ()(int i, int j, int k, int, Vector *u) const
  {
    Matrix jac(3, 3);
    Vector     vel;
    // u = J_w * v
    Jacobian(jac, _eta, _w, i, j, k);
    Evaluate(vel, _tau, _v, i, j, k);
    *u  = MatrixProduct(jac, vel);
    // u = J_w * v - J_v * w
    Jacobian(jac, _tau, _v, i, j, k);
    Evaluate(vel, _eta, _w, i, j, k);
    *u -= MatrixProduct(jac, vel);
  }
};

// -----------------------------------------------------------------------------
template <class TReal>
void EvaluateBCHFormulaDense(int nterms, GenericImage<TReal> &u,
                             const GenericImage<TReal> &v,
                             const GenericImage<TReal> &w,
                             bool minus_v, bool usejac)
{
  MIRTK_START_TIMING();

  GenericImage<TReal> l1, l2, l3, l4;
  if (nterms >= 3) {
    // - [v, w]
    liebracket(&l1, &v, &w, usejac);
    if (nterms >= 4) {
      // - [v, [v, w]]
      liebracket(&l2, &v, &l1, usejac);
      if (nterms >= 5) {
        // - [w, [w, v]] = [[v, w], w]
        liebracket(&l3, &l1, &w, usejac);
        if (nterms >= 6) {
          // - [[v, [v, w]], w]
          liebracket(&l4, &l2, &w, usejac);
          // - [[w, [v, w]], v] == [[v, [v, w]], w]
        }
      }
    }
  }

  // Evaluate BCH formula given all pre-computed terms and their respective weights
  const double weight1[] = {0.0, 1.0, 1.0/2.0, 1.0/12.0, 1.0/12.0, 1.0/48.0, 1.0/48.0};
  const double weight2[] = {1.0, 1.0, 1.0/2.0, 1.0/12.0, 1.0/12.0, 1.0/48.0, 1.0/48.0};

  NaryVoxelFunction::VoxelWiseWeightedSum bch;
  if (minus_v) bch._Weight = weight1;
  else         bch._Weight = weight2;

  switch (nterms) {
    case 1: { if (minus_v) { u = .0; } else { u = v; }                 break; }
    case 2: { ParallelForEachScalar(v, w,                     u, bch); break; }
    case 3: { ParallelForEachScalar(v, w, l1,                 u, bch); break; }
    case 4: { ParallelForEachScalar(v, w, l1, l2,             u, bch); break; }
    case 5: { ParallelForEachScalar(v, w, l1, l2, l3,         u, bch); break; }
    case 6: { ParallelForEachScalar(v, w, l1, l2, l3, l4,     u, bch); break; }
    case 7: { ParallelForEachScalar(v, w, l1, l2, l3, l4, l4, u, bch); break; }
    default:
      cerr << "BSplineFreeFormTransformationSV::EvaluateBCHFormulaDense: Invalid number of terms " << nterms << endl;
      exit(1);
  };
  MIRTK_DEBUG_TIMING(3, "evaluation of BCH formula (dense)");
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformationSV
::EvaluateBCHFormula(int nterms, CPImage &u,
                     double tau, const CPImage &v,
                     double eta, const CPImage &w,
                     bool minus_v) const
{
  MIRTK_START_TIMING();
  GenericImage<Vector> l1, l2, l3, l4;
  const ImageAttributes &lattice = u.Attributes();

  // Calculate required Lie brackets...
  if (_LieDerivative) {
    // ... using Lie derivative
    if (nterms >= 3) {
      // - [v, w]
      l1.Initialize(lattice, 3);
      ParallelForEachVoxel(SVFFDEvaluateLieBracket(this, tau, v, eta, w), lattice, l1);
      ConvertToCubicBSplineCoefficients(l1);
      if (nterms >= 4) {
        // - [v, [v, w]]
        l2.Initialize(lattice, 3);
        ParallelForEachVoxel(SVFFDEvaluateLieBracket(this, tau, v, 1., l1), lattice, l2);
        ConvertToCubicBSplineCoefficients(l2);
        if (nterms >= 5) {
          // - [w, [w, v]] = [[v, w], w]
          l3.Initialize(lattice, 3);
          ParallelForEachVoxel(SVFFDEvaluateLieBracket(this, 1., l1, eta, w), lattice, l3);
          ConvertToCubicBSplineCoefficients(l3);
          if (nterms >= 6) {
            // - [[v, [v, w]], w]
            l4.Initialize(lattice, 3);
            ParallelForEachVoxel(SVFFDEvaluateLieBracket(this, 1., l2, eta, w), lattice, l4);
            ConvertToCubicBSplineCoefficients(l4);
            // - [[w, [v, w]], v] == [[v, [v, w]], w]
          }
        }
      }
    }
  } else {
    // ... using composition of vector fields
    if (nterms >= 3) {
      DifferenceOfCompositionLieBracketImageFilter3D<Vector> lb;
      lb.Interpolation(Interpolation_CubicBSpline);
      lb.Extrapolation(Extrapolation_NN);
      lb.ComputeInterpolationCoefficients(false);
      // - [v, w]
      lb.Input(0, &v), lb.Scaling(0, tau);
      lb.Input(1, &w), lb.Scaling(1, eta);
      lb.Output(&l1);
      lb.Run();
      ConvertToCubicBSplineCoefficients(l1);
      if (nterms >= 4) {
        // - [v, [v, w]]
        lb.Input(0, &v ), lb.Scaling(0, tau);
        lb.Input(1, &l1), lb.Scaling(1, 1.0);
        lb.Output(&l2);
        lb.Run();
        ConvertToCubicBSplineCoefficients(l2);
        if (nterms >= 5) {
          // - [[v, w], w]
          lb.Input(0, &l1), lb.Scaling(0, 1.0);
          lb.Input(1, &w ), lb.Scaling(1, eta);
          lb.Output(&l3);
          lb.Run();
          ConvertToCubicBSplineCoefficients(l3);
          if (nterms >= 6) {
            // - [[v, [v, w]], w]
            lb.Input(0, &l2), lb.Scaling(0, 1.0);
            lb.Input(1, &w ), lb.Scaling(1, eta);
            lb.Output(&l4);
            lb.Run();
            ConvertToCubicBSplineCoefficients(l4);
            // - [[w, [v, w]], v] == [[v, [v, w]], w]
          }
        }
      }
    }
  }

  // Evaluate BCH formula given all pre-computed terms and their respective weights
  const double weight1[] = {0.0, eta, 1.0/2.0, 1.0/12.0, 1.0/12.0, 1.0/48.0, 1.0/48.0};
  const double weight2[] = {tau, eta, 1.0/2.0, 1.0/12.0, 1.0/12.0, 1.0/48.0, 1.0/48.0};

  NaryVoxelFunction::VoxelWiseWeightedSum bch;
  if (minus_v) bch._Weight = weight1;
  else         bch._Weight = weight2;

  switch (nterms) {
    case 1: { if (minus_v) { u = .0; } else { u = v; }                 break; }
    case 2: { ParallelForEachScalar(v, w,                     u, bch); break; }
    case 3: { ParallelForEachScalar(v, w, l1,                 u, bch); break; }
    case 4: { ParallelForEachScalar(v, w, l1, l2,             u, bch); break; }
    case 5: { ParallelForEachScalar(v, w, l1, l2, l3,         u, bch); break; }
    case 6: { ParallelForEachScalar(v, w, l1, l2, l3, l4,     u, bch); break; }
    case 7: { ParallelForEachScalar(v, w, l1, l2, l3, l4, l4, u, bch); break; }
    default:
      cerr << "BSplineFreeFormTransformationSV::EvaluateBCHFormula: Invalid number of terms " << nterms << endl;
      exit(1);
  };
  MIRTK_DEBUG_TIMING(3, "evaluation of BCH formula");
}

// =============================================================================
// Approximation/Interpolation
// =============================================================================

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformationSV
::ApproximateDOFs(const double *x,  const double *y,  const double *z, const double *t,
                  const double *dx, const double *dy, const double *dz, int no)
{
  // FIXME The initial approximation of the displacements on the control point grid
  //       smoothes the displacement field too much and hence introduces quite some
  //       error. Use the overloaded ApproximateAsNew(disp) method when possible.
  //       -as12312

  // Approximate displacements at control points
  double *rx = Allocate<double>(no);
  double *ry = Allocate<double>(no);
  double *rz = Allocate<double>(no);

  memcpy(rx, dx, no * sizeof(double));
  memcpy(ry, dy, no * sizeof(double));
  memcpy(rz, dz, no * sizeof(double));

  BSplineFreeFormTransformation3D::ApproximateDOFs(x, y, z, t, rx, ry, rz, no);

  Deallocate(rx);
  Deallocate(ry);
  Deallocate(rz);

  // Find stationary velocity field which approximates the displacements
  GenericImage<double> disp(this->Attributes(), 3);
  for (int k = 0; k < _z; ++k) {
    for (int j = 0; j < _y; ++j) {
      for (int i = 0; i < _x; ++i) {
        disp(i, j, k, 0) = _CPImage(i, j, k)._x;
        disp(i, j, k, 1) = _CPImage(i, j, k)._y;
        disp(i, j, k, 2) = _CPImage(i, j, k)._z;
      }
    }
  }
  this->ApproximateAsNew(disp);
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformationSV
::ApproximateDOFsGradient(const double *, const double *, const double *, const double *,
                          const double *, const double *, const double *, int,
                          double *, double) const
{
  cerr << this->NameOfClass() << "::ApproximateDOFsGradient: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
ImageAttributes BSplineFreeFormTransformationSV
::ApproximationDomain(const ImageAttributes &attr, const Transformation *dof)
{
  if (!dof) return attr;

  ImageAttributes grid(attr);

  // Ensure that for all target voxels the displacement can be recovered
  // without requiring any extrapolation of the velocity field during
  // computation of the trajectory (integration, i.e., exponentiation)
  const Matrix i2w = grid.GetImageToWorldMatrix();
  const Matrix w2i = grid.GetWorldToImageMatrix();

  double margin_top    = .0;
  double margin_bottom = .0;
  double margin_left   = .0;
  double margin_right  = .0;
  double margin_front  = .0;
  double margin_back   = .0;

  double x, y, z, wx, wy, wz;
  for (int k = 0; k < attr._z; ++k)
  for (int j = 0; j < attr._y; ++j)
  for (int i = 0; i < attr._x; ++i) {
    wx = i2w(0, 0) * i + i2w(0, 1) * j + i2w(0, 2) * k + i2w(0, 3);
    wy = i2w(1, 0) * i + i2w(1, 1) * j + i2w(1, 2) * k + i2w(1, 3);
    wz = i2w(2, 0) * i + i2w(2, 1) * j + i2w(2, 2) * k + i2w(2, 3);
    dof->Transform(wx, wy, wz);
    x = w2i(0, 0) * wx + w2i(0, 1) * wy + w2i(0, 2) * wz + w2i(0, 3);
    y = w2i(1, 0) * wx + w2i(1, 1) * wy + w2i(1, 2) * wz + w2i(1, 3);
    z = w2i(2, 0) * wx + w2i(2, 1) * wy + w2i(2, 2) * wz + w2i(2, 3);
    if (x <  0           && -x              > margin_left)   margin_left   = -x;
    if (y <  0           && -y              > margin_bottom) margin_bottom = -y;
    if (z <  0           && -z              > margin_front)  margin_front  = -z;
    if (x >= grid._x - 1 && x - grid._x - 1 > margin_right)  margin_right  = x - grid._x - 1;
    if (y >= grid._y - 1 && y - grid._y - 1 > margin_top)    margin_top    = y - grid._y - 1;
    if (z >= grid._z - 1 && z - grid._z - 1 > margin_back)   margin_back   = z - grid._z - 1;
  }

  // Disregard small precision errors
  const double eps = 1e-6;
  margin_left   = round(margin_left   / eps) * eps;
  margin_right  = round(margin_right  / eps) * eps;
  margin_bottom = round(margin_bottom / eps) * eps;
  margin_top    = round(margin_top    / eps) * eps;
  margin_front  = round(margin_front  / eps) * eps;
  margin_back   = round(margin_back   / eps) * eps;

  // Account for inter-/extrapolation error on boundary of FFD lattice and
  // therefore make lattice a bit bigger than otherwise needed
  const double margin_safety = 1.5;
  margin_left   = ceil(margin_left   * margin_safety);
  margin_right  = ceil(margin_right  * margin_safety);
  margin_bottom = ceil(margin_bottom * margin_safety);
  margin_top    = ceil(margin_top    * margin_safety);
  margin_front  = ceil(margin_front  * margin_safety);
  margin_back   = ceil(margin_back   * margin_safety);

  // Compute of offsets by which lattice origin must be moved such that
  // the lattice origin is the center of the extended lattice again
  const double ox = (margin_right - margin_left  )  * grid._dx / 2.0;
  const double oy = (margin_top   - margin_bottom) * grid._dy / 2.0;
  const double oz = (margin_back  - margin_front ) * grid._dz / 2.0;

  // Initialize free-form deformation (for extended image grid)
  grid._x       += static_cast<int>(margin_left   + margin_right);
  grid._y       += static_cast<int>(margin_bottom + margin_top  );
  grid._z       += static_cast<int>(margin_front  + margin_back );
  grid._xorigin += grid._xaxis[0] * ox + grid._yaxis[0] * oy + grid._zaxis[0] * oz;
  grid._yorigin += grid._xaxis[1] * ox + grid._yaxis[1] * oy + grid._zaxis[1] * oz;
  grid._zorigin += grid._xaxis[2] * ox + grid._yaxis[2] * oy + grid._zaxis[2] * oz;

  return grid;
}

// -----------------------------------------------------------------------------
double BSplineFreeFormTransformationSV
::ApproximateAsNew(const ImageAttributes &domain, const Transformation *dof,
                   int niter, double max_error)
{
  const HomogeneousTransformation       *lin   = nullptr;
  const BSplineFreeFormTransformationSV *svffd = nullptr;
  const MultiLevelTransformation        *mffd  = nullptr;

  (lin   = dynamic_cast<const HomogeneousTransformation       *>(dof)) ||
  (svffd = dynamic_cast<const BSplineFreeFormTransformationSV *>(dof)) ||
  (mffd  = dynamic_cast<const MultiLevelTransformation        *>(dof));

  if (mffd) {
    if (mffd->NumberOfLevels() == 0) {
      lin = mffd->GetGlobalTransformation();
    } else if (mffd->NumberOfLevels() == 1) {
      if (mffd->GetGlobalTransformation()->IsIdentity()) {
        svffd = dynamic_cast<const BSplineFreeFormTransformationSV *>(mffd->GetLocalTransformation(0));
      }
    }
  }

  // Approximate any other transformation using the base class implementation
  // which simply evaluates the displacement of the transformation at each
  // control point and then calls Interpolate in order to interpolate these
  // control point displacements
  if (!lin && !svffd) {
    return BSplineFreeFormTransformation3D::ApproximateAsNew(domain, dof, niter, max_error);
  }
  // Check if input is the identity transformation
  if ((lin && lin->IsIdentity()) || (svffd && svffd->IsIdentity())) {
    _CPImage = 0.;
    return 0.;
  }

  double error = inf;
  if (niter >= 0) {
    // Compute velocities at control points using log map of affine matrix
    if (lin) {
      EvaluateGlobalSVFFD logA(logm(lin->GetMatrix()), &_CPImage);
      ParallelForEachVoxel(_CPImage.Attributes(), &_CPImage, logA);
    // Evaluate velocities of other SV FFD at control points of this SV FFD
    } else {
      ParallelForEachVoxel(EvaluateBSplineSVFFD(svffd, &_CPImage), _attr, _CPImage);
    }
    // Convert velocities to B-spline coefficients
    ConvertToSplineCoefficients(3, _CPImage);
    // Evaluate approximation error
    error = EvaluateRMSError(domain, dof);
  }

  return error;
}

// -----------------------------------------------------------------------------
//double BSplineFreeFormTransformationSV
//::ApproximateFFDAsNew(GenericImage<double> &disp, double T, int niter)
//{
//  GenericImage<double>  dv;
//  WorldCoordinatesImage wc;
//  disp.ImageToWorld(wc);
//  for (int iter = 0; iter < niter; ++iter) {
//    // Interpolate velocity field
//    const ImageAttributes &grid = d.GetImageAttributes();
//    v->Initialize(grid, d.GetT());
//    switch (d.GetT()) {
//      case 2: { ParallelForEachVoxel(EvaluateBSplineSVFFD2D(this, v), grid, v); break; }
//      case 3: { ParallelForEachVoxel(EvaluateBSplineSVFFD3D(this, v), grid, v); break; }
//      default: {
//        cerr << "BSplineFreeFormTransformationSV::ScalingAndSquaring: Vector field must have 2 or 3 components (_t)" << endl;
//        exit(1);
//      }
//    }
//    // Exponentiate velocity field
//    VelocityToDisplacementFieldSS<double> exp;
//    exp.T(T);
//    exp.NumberOfSteps    (NumberOfStepsForIntervalLength(T));
//    exp.MaxScaledVelocity(_MaxScaledVelocity);
//    exp.Interpolation    (Interpolation_BSpline);
//    exp.Upsample         (false); // better, but too expensive
//    exp.ExternalCache    (cache); // avoid frequent allocation/deallocation
//    exp.SetInput         (0, v);  // velocity field to be exponentiated
//    exp.SetInput         (1, &d); // input displacement field (may be zero)
//    exp.SetOutput        (v);     // result is exp(v) o d
//    exp.Run();
//    d.CopyFrom(*v);
//    // dv = exp(-v) o disp
//    dv = disp;
//    ScalingAndSquaring(dv, T, &wc);
//    // Approximate dv by B-spline coefficients at control points
//    BSplineFreeFormTransformation3D ffd(this->GetFFDAttributes());
//    ffd.Approximate(wc.GetPointerToVoxels(0, 0, 0, 0),
//                    wc.GetPointerToVoxels(0, 0, 0, 1),
//                    wc.GetPointerToVoxels(0, 0, 0, 2),
//                    dv.GetPointerToVoxels(0, 0, 0, 0),
//                    dv.GetPointerToVoxels(0, 0, 0, 1),
//                    dv.GetPointerToVoxels(0, 0, 0, 2));
//    dv.Initialize(_x, _y, _z, 1, 3);
//    for (int k = 0; k < _z; ++k) {
//      for (int j = 0; j < _y; ++j) {
//        for (int i = 0; i < _x; ++i) {
//          dv(i, j, k, 0) = ffd._data[k][j][i]._x;
//          dv(i, j, k, 1) = ffd._data[k][j][i]._x;
//          dv(i, j, k, 2) = ffd._data[k][j][i]._x;
//        }
//      }
//    }
//    // Calculate required Lie brackets
//    if (_NumberOfBCHTerms >= 3) {
//      DifferenceOfCompositionLieBracketImageFilter3D<double> lb;
//      lb.Interpolation(Interpolation_CubicBSpline);
//      lb.Extrapolation(Extrapolation_NN);
//      lb.ComputeInterpolationCoefficients(false);
//      // - [v, dv]
//      lb.SetInput (0, v);
//      lb.SetInput (1, d);
//      lb.SetOutput(l1);
//      lb.Run();
//      ConvertToCubicBSplineCoefficients(*l1);
//      if (_NumberOfBCHTerms >= 4) {
//        // - [v, [v, dv]]
//        lb.SetInput (0, v);
//        lb.SetInput (1, l1);
//        lb.SetOutput(l2);
//        lb.Run();
//        ConvertToCubicBSplineCoefficients(*l2);
//        if (_NumberOfBCHTerms >= 5) {
//          // - [dv, [v, dv]]
//          lb.SetInput (0, d);
//          lb.SetInput (1, l1);
//          lb.SetOutput(l3);
//          lb.Run();
//          ConvertToCubicBSplineCoefficients(*l3);
//          if (_NumberOfBCHTerms >= 6) {
//            // - [dv, [v, [v, dv]]]
//            lb.SetInput (0, d);
//            lb.SetInput (1, l2);
//            lb.SetOutput(l4);
//            lb.Run();
//            ConvertToCubicBSplineCoefficients(*l4);
//          }
//        }
//      }
//    }
//    // Update velocity coefficients using BCH formula
//    NaryVoxelFunction::EvaluateBCHFormula bch;
//    switch (_NumberOfBCHTerms) {
//      case 2: { /* nothing else to do */                              break; }
//      case 3: { ParallelForEachScalar(v, dv, l1,             v, bch); break; }
//      case 4: { ParallelForEachScalar(v, dv, l1, l2,         v, bch); break; }
//      case 5: { ParallelForEachScalar(v, dv, l1, l2, l3,     v, bch); break; }
//      case 6: { ParallelForEachScalar(v, dv, l1, l2, l3, l4, v, bch); break; }
//    };
//  }
//}

// -----------------------------------------------------------------------------
double BSplineFreeFormTransformationSV
::ApproximateAsNew(GenericImage<double> &disp, int niter, double max_error)
{
  return this->ApproximateAsNew(disp, false, 3, niter * 8);
}

// -----------------------------------------------------------------------------
double BSplineFreeFormTransformationSV
::ApproximateAsNew(GenericImage<double> &disp, bool smooth, int nterms, int niter)
{
  // TODO: Refactor/review implementation again after update of velocities
  //       from gradient is now implemented and working.

  ImageAttributes grid = this->Attributes();
  grid._t  = 3;
  grid._dt = disp.GetTSize(); // ignore difference in _dt

  // Sample displacement field at control points using linear interpolation
  GenericImage<double> *d = NULL;

  if (disp.Attributes() == grid) {
    d = &disp;
  } else {    
    double x, y, z, vec[3] = {.0, .0, .0};

    UniquePtr<InterpolateImageFunction> f;
    f.reset(InterpolateImageFunction::New(Interpolation_Linear, Extrapolation_NN, &disp));
    f->Input(&disp);
    f->Initialize();

    d = new GenericImage<double>(grid);

    for (int k = 0; k < grid._z; ++k)
    for (int j = 0; j < grid._y; ++j)
    for (int i = 0; i < grid._x; ++i) {
      x = i, y = j, z = k;
      d->ImageToWorld(x, y, z);
      disp.WorldToImage(x, y, z);
      f->Evaluate(vec, x, y, z);
      d->Put(i, j, k, 0, vec[0]);
      d->Put(i, j, k, 1, vec[1]);
      d->Put(i, j, k, 2, vec[2]);
    }
  }

  // Compute stationary velocity field at control points
  GenericImage<double>                   v;
  DisplacementToVelocityFieldBCH<double> dtov;

  dtov.Input (d);
  dtov.Output(&v);

  dtov.UpperIntegrationLimit(UpperIntegrationLimit(0, 1));
  dtov.NumberOfIterations(niter);
  dtov.NumberOfTerms     (nterms);
  dtov.NumberOfSteps     (NumberOfStepsForIntervalLength(dtov.UpperIntegrationLimit()));
  dtov.SmoothVelocities  (smooth);

  dtov.Run();

  // Free temporary displacement field
  if (d != &disp) {
    delete d;
    d = NULL;
  }

  // Interpolate velocities by B-spline function
  BSplineFreeFormTransformation3D::Interpolate(v.Data(0, 0, 0, 0),
                                               v.Data(0, 0, 0, 1),
                                               v.Data(0, 0, 0, 2));

  // Evaluate RMS of approximation error
  double error = .0;

  v = .0;
  this->Displacement(v, 0, 1);

  for (int k = 0; k < disp.Z(); ++k)
  for (int j = 0; j < disp.Y(); ++j)
  for (int i = 0; i < disp.X(); ++i) {
    disp(i, j, k, 0) -= v(i, j, k, 0);
    disp(i, j, k, 1) -= v(i, j, k, 1);
    disp(i, j, k, 2) -= v(i, j, k, 2);
    error += sqrt(disp(i, j ,k, 0) * disp(i, j, k, 0) +
                  disp(i, j ,k, 1) * disp(i, j, k, 1) +
                  disp(i, j ,k, 2) * disp(i, j, k, 2));
  }
  error /= static_cast<double>(disp.X() * disp.Y() * disp.Z());

  return error;
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformationSV
::Interpolate(const double *, const double *, const double *)
{
  cerr << this->NameOfClass() << "::Interpolate: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformationSV
::InterpolateVelocities(const double *vx, const double *vy, const double *vz)
{
  BSplineFreeFormTransformation3D::Interpolate(vx, vy, vz);
}

// -----------------------------------------------------------------------------
double BSplineFreeFormTransformationSV
::ApproximateVelocitiesAsNew(GenericImage<double> &v)
{
  double error = 0.;

  auto * const vx = v.Data(0, 0, 0, 0);
  auto * const vy = v.Data(0, 0, 0, 1);
  auto * const vz = v.Data(0, 0, 0, 2);

  if (_attr.EqualInSpace(v.Attributes())) {

    BSplineFreeFormTransformation3D::Interpolate(vx, vy, vz);

  } else {

    // Compute world coordinates of lattice points
    const int no = v.NumberOfSpatialVoxels();
    double *x = Allocate<double>(no);
    double *y = Allocate<double>(no);
    double *z = Allocate<double>(no);
    double *t = Allocate<double>(no);
    v.Attributes().LatticeToWorld(x, y, z, t);

    // Approximate velocities
    BSplineFreeFormTransformation3D::ApproximateDOFs(x, y, z, t, vx, vy, vz, no);

    // Evaluate error of approximation
    GenericImage<double> u(v.Attributes(), 3);
    this->Velocity(u);
    for (int k = 0; k < v.Z(); ++k)
    for (int j = 0; j < v.Y(); ++j)
    for (int i = 0; i < v.X(); ++i) {
      const int idx = v.VoxelToIndex(i, j, k);
      vx[idx] -= u(i, j, k, 0);
      vy[idx] -= u(i, j, k, 1);
      vz[idx] -= u(i, j, k, 2);
      error += pow(vx[idx], 2) + pow(vy[idx], 2) + pow(vz[idx], 2);
    }
    error = sqrt(error / no);

    // Free memory
    Deallocate(x);
    Deallocate(y);
    Deallocate(z);
    Deallocate(t);
  }

  return error;
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformationSV::CombineWith(const Transformation *dof)
{
  // Convert transformation into SV FFD
  const BSplineFreeFormTransformationSV *svffd = NULL;
  svffd = dynamic_cast<const BSplineFreeFormTransformationSV *>(dof);
  if (!svffd) {
    BSplineFreeFormTransformationSV *tmp;
    tmp = new BSplineFreeFormTransformationSV();
    tmp->Initialize(this->Attributes());
    tmp->ApproximateAsNew(dof);
    svffd = tmp;
  }
  // Compute coefficients of composite SV FFD using BCH formula
  this->CombineWith(svffd);
  // Delete temporary SV FFD
  if (svffd != dof) delete svffd;
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformationSV::CombineWith(const BSplineFreeFormTransformationSV *svffd)
{
  EvaluateBCHFormula(_NumberOfBCHTerms, _CPImage, _T, _CPImage, svffd->T(), svffd->_CPImage);
  if (_T < 0.) this->Invert();
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformationSV::ScaleVelocities(double e)
{
  _CPImage *= e;
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformationSV::Invert()
{
  this->ScaleVelocities(-1.);
}

// =============================================================================
// Parameters (non-DoFs)
// =============================================================================

// -----------------------------------------------------------------------------
bool BSplineFreeFormTransformationSV::Set(const char *name, const char *value)
{
  if (strcmp(name, "Cross-sectional time interval") == 0 ||
      strcmp(name, "Cross sectional time interval") == 0) {
    return FromString(value, _T);
  } else if (strcmp(name, "Time unit of integration interval") == 0) {
    return FromString(value, _TimeUnit);
  } else if (strcmp(name,    "No. of integration steps") == 0 ||
             strcmp(name, "Number of integration steps") == 0) {
    return FromString(value, _NumberOfSteps) && _NumberOfSteps > 0;
  } else if (strcmp(name,    "No. of squaring steps")    == 0 ||
             strcmp(name, "Number of squaring steps") == 0) {
    if (!FromString(value, _NumberOfSteps) || _NumberOfSteps <= 0) return false;
    _NumberOfSteps = static_cast<int>(pow(2, _NumberOfSteps));
    if (_IntegrationMethod != FFDIM_SS && _IntegrationMethod != FFDIM_FastSS) {
      _IntegrationMethod = FFDIM_FastSS;
    }
  } else if (strcmp(name, "Maximum scaled velocity") == 0) {
    return FromString(value, _MaxScaledVelocity);
  } else if (strcmp(name, "Use Lie derivative") == 0) {
    return FromString(value, _LieDerivative);
  } else if (strcmp(name, "Use dense BCH lattice") == 0) {
    return FromString(value, _UseDenseBCHGrid);
  } else if (strcmp(name,    "No. of BCH terms")    == 0 ||
             strcmp(name, "Number of BCH terms") == 0) {
    return FromString(value, _NumberOfBCHTerms) && _NumberOfBCHTerms <= 6;
  } else if (strcmp(name, "Integration method") == 0) {
    return FromString(value, _IntegrationMethod) && _IntegrationMethod != FFDIM_Unknown;
  // deprecated parameters
  } else if (strcmp(name, "Use scaling and squaring") == 0) {
    bool useSS = false;
    if (!FromString(value, useSS)) return false;
    if (useSS) {
      if (_IntegrationMethod != FFDIM_SS && _IntegrationMethod != FFDIM_FastSS) {
        _IntegrationMethod = FFDIM_FastSS;
      }
    } else {
      if (_IntegrationMethod == FFDIM_SS || _IntegrationMethod == FFDIM_FastSS) {
        _IntegrationMethod = FFDIM_RKE1;
      }
    }
    return true;
  } else if (strcmp(name, "Fast scaling and squaring") == 0) {
    bool fastSS;
    if (!FromString(value, fastSS)) return false;
    if (_IntegrationMethod == FFDIM_SS && fastSS) _IntegrationMethod = FFDIM_FastSS;
    return true;
  }
  return BSplineFreeFormTransformation3D::Set(name, value);
}

// -----------------------------------------------------------------------------
ParameterList BSplineFreeFormTransformationSV::Parameter() const
{
  ParameterList params = BSplineFreeFormTransformation3D::Parameter();
  Insert(params, "Integration method",                _IntegrationMethod);
  Insert(params, "Cross-sectional time interval",     _T);
  Insert(params, "Time unit of integration interval", _TimeUnit);
  Insert(params, "No. of integration steps",          _NumberOfSteps);
  Insert(params, "Maximum scaled velocity",           _MaxScaledVelocity);
  Insert(params, "Use Lie derivative",                _LieDerivative);
  Insert(params, "Use dense BCH lattice",             _UseDenseBCHGrid);
  Insert(params, "No. of BCH terms",                  _NumberOfBCHTerms);
  return params;
}

// =============================================================================
// Point transformation
// =============================================================================

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformationSV::IntegrateVelocities(double &x, double &y, double &z, double T) const
{
  const double dt = StepLengthForIntervalLength(T);
  if (dt) {
    if      (_IntegrationMethod == FFDIM_FastSS ||
             _IntegrationMethod == FFDIM_SS     ||
             _IntegrationMethod == FFDIM_RKE1)   RKE1  ::Transform(this, x, y, z, .0, T, dt);
    else if (_IntegrationMethod == FFDIM_RKE2)   RKE2  ::Transform(this, x, y, z, .0, T, dt);
    else if (_IntegrationMethod == FFDIM_RKH2)   RKH2  ::Transform(this, x, y, z, .0, T, dt);
    else if (_IntegrationMethod == FFDIM_RK4)    RK4   ::Transform(this, x, y, z, .0, T, dt);
    else if (_IntegrationMethod == FFDIM_RKEH12) RKEH12::Transform(this, x, y, z, .0, T, 0.5 * dt, 2.0 * dt, SVFFD_RKTOL);
    else if (_IntegrationMethod == FFDIM_RKBS23) RKBS23::Transform(this, x, y, z, .0, T, 0.5 * dt, 2.0 * dt, SVFFD_RKTOL);
    else if (_IntegrationMethod == FFDIM_RKF45)  RKF45 ::Transform(this, x, y, z, .0, T, 0.5 * dt, 2.0 * dt, SVFFD_RKTOL);
    else if (_IntegrationMethod == FFDIM_RKCK45) RKCK45::Transform(this, x, y, z, .0, T, 0.5 * dt, 2.0 * dt, SVFFD_RKTOL);
    else if (_IntegrationMethod == FFDIM_RKDP45) RKDP45::Transform(this, x, y, z, .0, T, 0.5 * dt, 2.0 * dt, SVFFD_RKTOL);
    else {
      cerr << "BSplineFreeFormTransformationSV::IntegrateVelocities: Unknown integration method: " << _IntegrationMethod << endl;
      exit(1);
    }
  }
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformationSV
::LocalTransform(double &x, double &y, double &z, double t, double t0) const
{
  IntegrateVelocities(x, y, z, + UpperIntegrationLimit(t, t0));
}

// -----------------------------------------------------------------------------
bool BSplineFreeFormTransformationSV
::LocalInverse(double &x, double &y, double &z, double t, double t0) const
{
  IntegrateVelocities(x, y, z, - UpperIntegrationLimit(t, t0));
  return true;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void BSplineFreeFormTransformationSV
::ScalingAndSquaring(GenericImage<VoxelType> *d,
                     double T, const WorldCoordsImage *wc) const
{
  ScalingAndSquaring<VoxelType>(d->Attributes(), d, nullptr, nullptr, nullptr, T, wc);
}

template void BSplineFreeFormTransformationSV
::ScalingAndSquaring(GenericImage<float> *d,
                     double T, const WorldCoordsImage *) const;

template void BSplineFreeFormTransformationSV
::ScalingAndSquaring(GenericImage<double> *d,
                     double T, const WorldCoordsImage *) const;

// -----------------------------------------------------------------------------
template <class VoxelType>
void BSplineFreeFormTransformationSV
::ScalingAndSquaring(const ImageAttributes   &a,
                     GenericImage<VoxelType> *d,
                     GenericImage<VoxelType> *dx,
                     GenericImage<VoxelType> *dj,
                     GenericImage<VoxelType> *lj,
                     double T, const WorldCoordsImage *) const
{
  // Attributes of output images
  ImageAttributes attr(a);
  if (!attr) {
    if      (d ) attr = d ->Attributes();
    else if (dx) attr = dx->Attributes();
    else if (dj) attr = dj->Attributes();
    else if (lj) attr = lj->Attributes();
  }
  attr._t = 1, attr._dt = .0;
  if (!attr) return;
  // Copy input displacement field
  UniquePtr<GenericImage<VoxelType> > din(d ? new GenericImage<VoxelType>(*d) : nullptr);
  // TODO: The runtime of the ScalingAndSquaring filter has been greatly improved
  //       to almost match the old VelocityToDisplacementFieldSS implementation.
  //       However, the latter is still about 15% faster...
  if (d && !dx && !dj && !lj) {
    GenericImage<VoxelType> v;
    if (_IntegrationMethod == FFDIM_FastSS) {
      v.Initialize(this->Attributes(), 3);
      ParallelForEachVoxel(EvaluateBSplineSVFFD3D(this, &v), this->Attributes(), v);
    } else {
      v.Initialize(attr, 3);
      ParallelForEachVoxel(EvaluateBSplineSVFFD3D(this, &v), attr, v);
    }
    // Exponentiate velocity field
    VelocityToDisplacementFieldSS<VoxelType> exp;
    exp.UpperIntegrationLimit(T);
    exp.NumberOfSteps(NumberOfStepsForIntervalLength(T));
    exp.MaxScaledVelocity(static_cast<VoxelType>(_MaxScaledVelocity));
    exp.Interpolation(Interpolation_Linear);
    exp.Upsample(false);     // better, but too expensive
    exp.Input(0, &v);        // velocity field to be exponentiated
    exp.Input(1, din.get()); // input displacement field (may be nullptr)
    exp.Output(d);           // result is exp(v) o d
    exp.Run();
  } else {
    // Copy B-spline coefficients of velocity field
    GenericImage<VoxelType> v(this->Attributes(), 3);
    VoxelType *vx = v.Data(0, 0, 0, 0);
    VoxelType *vy = v.Data(0, 0, 0, 1);
    VoxelType *vz = v.Data(0, 0, 0, 2);
    const Vector *vp = _CPImage.Data();
    for (int idx = 0; idx < _CPImage.NumberOfVoxels(); ++idx, ++vx, ++vy, ++vz, ++vp) {
      *vx = static_cast<VoxelType>(vp->_x);
      *vy = static_cast<VoxelType>(vp->_y);
      *vz = static_cast<VoxelType>(vp->_z);
    }
    // Exponentiate velocity field
    mirtk::ScalingAndSquaring<VoxelType> exp;
    exp.UpperIntegrationLimit(T);
    exp.NumberOfSteps(NumberOfStepsForIntervalLength(T));
    exp.MaxScaledVelocity(_MaxScaledVelocity);
    exp.InterimAttributes(_IntegrationMethod == FFDIM_FastSS ? this->Attributes() : attr);
    exp.OutputAttributes(attr);
    exp.Upsample(false);                         // better, but too expensive
    exp.ComputeInterpolationCoefficients(false); // v contains B-spline coefficients
    exp.InputVelocity(&v);                       // velocity field to be exponentiated
    exp.InputDisplacement(din.get());            // input displacement field (may be zero)
    exp.OutputDisplacement(d);                   // i.e., d = exp(v) o din
    exp.OutputJacobian(dx);                      // i.e., Jacobian
    exp.OutputDetJacobian(dj);                   // i.e., det(Jacobian)
    exp.OutputLogJacobian(lj);                   // i.e., log(det(Jacobian)
    exp.Run();
  }
}

template void BSplineFreeFormTransformationSV
::ScalingAndSquaring(const ImageAttributes &,
                     GenericImage<float> *,
                     GenericImage<float> *,
                     GenericImage<float> *,
                     GenericImage<float> *,
                     double, const WorldCoordsImage *) const;

template void BSplineFreeFormTransformationSV
::ScalingAndSquaring(const ImageAttributes &,
                     GenericImage<double> *,
                     GenericImage<double> *,
                     GenericImage<double> *,
                     GenericImage<double> *,
                     double, const WorldCoordsImage *) const;

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformationSV
::Displacement(GenericImage<float> &d, double t, double t0, const WorldCoordsImage *wc) const
{
  double T;
  if ((T = + UpperIntegrationLimit(t, t0))) {
    MIRTK_START_TIMING();
    // Use scaling and squaring method to efficiently compute displacements when possible
    if ((_IntegrationMethod == FFDIM_SS || _IntegrationMethod == FFDIM_FastSS) &&
        ((_z <= 1 && d.Z() <= 1) || (_z > 1 && d.Z() > 1))) {
      ScalingAndSquaring(&d, T, wc);
    // Evaluate transformation at each voxel separately using numerical integration
    } else {
      Transformation::Displacement(d, t, t0, wc);
    }
    MIRTK_DEBUG_TIMING(3, "computation of exp(" << T << "*v)");
  }
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformationSV
::Displacement(GenericImage<double> &d, double t, double t0, const WorldCoordsImage *wc) const
{
  double T;
  if ((T = + UpperIntegrationLimit(t, t0))) {
    MIRTK_START_TIMING();
    // Use scaling and squaring method to efficiently compute displacements when possible
    if ((_IntegrationMethod == FFDIM_SS || _IntegrationMethod == FFDIM_FastSS) &&
        ((_z <= 1 && d.Z() <= 1) || (_z > 1 && d.Z() > 1))) {
      ScalingAndSquaring(&d, T, wc);
    // Evaluate transformation at each voxel separately using numerical integration
    } else {
      Transformation::Displacement(d, t, t0, wc);
    }
    MIRTK_DEBUG_TIMING(3, "computation of exp(" << T << "*v)");
  }
}

// -----------------------------------------------------------------------------
int BSplineFreeFormTransformationSV
::InverseDisplacement(GenericImage<float> &d, double t, double t0, const WorldCoordsImage *wc) const
{
  double T;
  if ((T = - UpperIntegrationLimit(t, t0))) {
    MIRTK_START_TIMING();
    // Use scaling and squaring method to efficiently compute displacements when possible
    if ((_IntegrationMethod == FFDIM_SS || _IntegrationMethod == FFDIM_FastSS) &&
        ((_z <= 1 && d.Z() <= 1) || (_z > 1 && d.Z() > 1))) {
      ScalingAndSquaring(&d, T, wc);
    // Evaluate transformation at each voxel separately using numerical integration
    } else {
      Transformation::InverseDisplacement(d, t, t0, wc);
    }
    MIRTK_DEBUG_TIMING(3, "computation of exp(" << T << "*v)");
  }
  return 0;
}

// -----------------------------------------------------------------------------
int BSplineFreeFormTransformationSV
::InverseDisplacement(GenericImage<double> &d, double t, double t0, const WorldCoordsImage *wc) const
{
  double T;
  if ((T = - UpperIntegrationLimit(t, t0))) {
    MIRTK_START_TIMING();
    // Use scaling and squaring method to efficiently compute displacements when possible
    if ((_IntegrationMethod == FFDIM_SS || _IntegrationMethod == FFDIM_FastSS) &&
        ((_z <= 1 && d.Z() <= 1) || (_z > 1 && d.Z() > 1))) {
      ScalingAndSquaring(&d, T, wc);
    // Evaluate transformation at each voxel separately uing numerical integration
    } else {
      Transformation::InverseDisplacement(d, t, t0, wc);
    }
    MIRTK_DEBUG_TIMING(3, "computation of exp(" << T << "*v)");
  }
  return 0;
}

// =============================================================================
// Derivatives
// =============================================================================

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformationSV::LocalJacobian(Matrix &jac, double x, double y, double z, double t, double t0) const
{
  jac.Initialize(3, 3);
  jac.Ident();
  double dt, T;
  if ((dt = StepLengthForIntervalLength(T = UpperIntegrationLimit(t, t0)))) {
    if      (_IntegrationMethod == FFDIM_SS     ||
             _IntegrationMethod == FFDIM_FastSS ||
             _IntegrationMethod == FFDIM_RKE1)   RKE1  ::Jacobian(this, jac, x, y, z, .0, T, dt);
    else if (_IntegrationMethod == FFDIM_RKE2)   RKE2  ::Jacobian(this, jac, x, y, z, .0, T, dt);
    else if (_IntegrationMethod == FFDIM_RKH2)   RKH2  ::Jacobian(this, jac, x, y, z, .0, T, dt);
    else if (_IntegrationMethod == FFDIM_RK4)    RK4   ::Jacobian(this, jac, x, y, z, .0, T, dt);
    else if (_IntegrationMethod == FFDIM_RKEH12) RKEH12::Jacobian(this, jac, x, y, z, .0, T, 0.5 * dt, 2.0 * dt, SVFFD_RKTOL);
    else if (_IntegrationMethod == FFDIM_RKBS23) RKBS23::Jacobian(this, jac, x, y, z, .0, T, 0.5 * dt, 2.0 * dt, SVFFD_RKTOL);
    else if (_IntegrationMethod == FFDIM_RKF45)  RKF45 ::Jacobian(this, jac, x, y, z, .0, T, 0.5 * dt, 2.0 * dt, SVFFD_RKTOL);
    else if (_IntegrationMethod == FFDIM_RKCK45) RKCK45::Jacobian(this, jac, x, y, z, .0, T, 0.5 * dt, 2.0 * dt, SVFFD_RKTOL);
    else if (_IntegrationMethod == FFDIM_RKDP45) RKDP45::Jacobian(this, jac, x, y, z, .0, T, 0.5 * dt, 2.0 * dt, SVFFD_RKTOL);
    else {
      cerr << "BSplineFreeFormTransformationSV::Jacobian: Unknown integration method: " << _IntegrationMethod << endl;
      exit(1);
    }
  }
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformationSV::LocalHessian(Matrix [3], double, double, double, double, double) const
{
  cerr << this->NameOfClass() << "::LocalHessian: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformationSV::JacobianDOFs(Matrix &jac, int cp, double x, double y, double z, double t, double t0) const
{
  jac.Initialize(3, 3);
  jac.Ident();
  double dt, T;
  if ((dt = StepLengthForIntervalLength(T = UpperIntegrationLimit(t, t0)))) {
    int ci, cj, ck, cl;
    this->IndexToLattice(cp, ci, cj, ck, cl);
    if      (_IntegrationMethod == FFDIM_SS     ||
             _IntegrationMethod == FFDIM_FastSS ||
             _IntegrationMethod == FFDIM_RKE1)   RKE1  ::JacobianDOFs(this, jac, ci, cj, ck, cl, x, y, z, .0, T, dt);
    else if (_IntegrationMethod == FFDIM_RKE2)   RKE2  ::JacobianDOFs(this, jac, ci, cj, ck, cl, x, y, z, .0, T, dt);
    else if (_IntegrationMethod == FFDIM_RKH2)   RKH2  ::JacobianDOFs(this, jac, ci, cj, ck, cl, x, y, z, .0, T, dt);
    else if (_IntegrationMethod == FFDIM_RK4)    RK4   ::JacobianDOFs(this, jac, ci, cj, ck, cl, x, y, z, .0, T, dt);
    else if (_IntegrationMethod == FFDIM_RKEH12) RKEH12::JacobianDOFs(this, jac, ci, cj, ck, cl, x, y, z, .0, T, 0.5 * dt, 2.0 * dt, SVFFD_RKTOL);
    else if (_IntegrationMethod == FFDIM_RKBS23) RKBS23::JacobianDOFs(this, jac, ci, cj, ck, cl, x, y, z, .0, T, 0.5 * dt, 2.0 * dt, SVFFD_RKTOL);
    else if (_IntegrationMethod == FFDIM_RKF45)  RKF45 ::JacobianDOFs(this, jac, ci, cj, ck, cl, x, y, z, .0, T, 0.5 * dt, 2.0 * dt, SVFFD_RKTOL);
    else if (_IntegrationMethod == FFDIM_RKCK45) RKCK45::JacobianDOFs(this, jac, ci, cj, ck, cl, x, y, z, .0, T, 0.5 * dt, 2.0 * dt, SVFFD_RKTOL);
    else if (_IntegrationMethod == FFDIM_RKDP45) RKDP45::JacobianDOFs(this, jac, ci, cj, ck, cl, x, y, z, .0, T, 0.5 * dt, 2.0 * dt, SVFFD_RKTOL);
    else {
      cerr << "BSplineFreeFormTransformationSV::JacobianDOFs: Unknown integration method: " << _IntegrationMethod << endl;
      exit(1);
    }
  }
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformationSV::JacobianDOFs(double jac[3], int dof, double x, double y, double z, double t, double t0) const
{
  Matrix dTdp(3, 3);
  this->JacobianDOFs(dTdp, dof / 3, x, y, z, t, t0);
  const int c = dof % 3;
  jac[0] = dTdp(0, c);
  jac[1] = dTdp(1, c);
  jac[2] = dTdp(2, c);
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformationSV
::EvaluateJacobianDOFs(TransformationJacobian &jac, double x, double y) const
{
  int i = ifloor(x);
  int j = ifloor(y);

  const int A = Kernel::VariableToIndex(x - i);
  const int B = Kernel::VariableToIndex(y - j);

  --i, --j;

  double wxy, wy;
  int    ci, cj, xdof, ydof;

  for (int b = 0; b <= 3; ++b) {
    cj = j + b;
    if (cj < 0 || cj >= _y) continue;
    wy = Kernel::LookupTable[B][b];
    for (int a = 0; a <= 3; ++a) {
      ci = i + a;
      if (ci < 0 || ci >= _x) continue;
      wxy = Kernel::LookupTable[A][a] * wy;
      IndexToDOFs(LatticeToIndex(ci, cj), xdof, ydof);
      jac(xdof)._x = jac(ydof)._y = wxy;
    }
  }
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformationSV
::EvaluateJacobianDOFs(TransformationJacobian &jac, double x, double y, double z) const
{
  int i = ifloor(x);
  int j = ifloor(y);
  int k = ifloor(z);

  const int A = Kernel::VariableToIndex(x - i);
  const int B = Kernel::VariableToIndex(y - j);
  const int C = Kernel::VariableToIndex(z - k);

  --i, --j, --k;

  double wxyz, wyz, wz;
  int    ci, cj, ck, xdof, ydof, zdof;

  for (int c = 0; c <= 3; ++c) {
    ck = k + c;
    if (ck < 0 || ck >= _z) continue;
    wz = Kernel::LookupTable[C][c];
    for (int b = 0; b <= 3; ++b) {
      cj = j + b;
      if (cj < 0 || cj >= _y) continue;
      wyz = Kernel::LookupTable[B][b] * wz;
      for (int a = 0; a <= 3; ++a) {
        ci = i + a;
        if (ci < 0 || ci >= _x) continue;
        wxyz = Kernel::LookupTable[A][a] * wyz;
        IndexToDOFs(LatticeToIndex(ci, cj, ck), xdof, ydof, zdof);
        jac(xdof)._x = jac(ydof)._y = jac(zdof)._z = wxyz;
      }
    }
  }
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformationSV::JacobianDOFs(TransformationJacobian &jac, double x, double y, double z, double t, double t0) const
{
  jac.Clear();
  double dt, T;
  if ((dt = StepLengthForIntervalLength(T = UpperIntegrationLimit(t, t0)))) {
    if      (_IntegrationMethod == FFDIM_SS     ||
             _IntegrationMethod == FFDIM_FastSS ||
             _IntegrationMethod == FFDIM_RKE1)   RKE1  ::JacobianDOFs(this, jac, x, y, z, .0, T, dt);
    else if (_IntegrationMethod == FFDIM_RKE2)   RKE2  ::JacobianDOFs(this, jac, x, y, z, .0, T, dt);
    else if (_IntegrationMethod == FFDIM_RKH2)   RKH2  ::JacobianDOFs(this, jac, x, y, z, .0, T, dt);
    else if (_IntegrationMethod == FFDIM_RK4)    RK4   ::JacobianDOFs(this, jac, x, y, z, .0, T, dt);
    else if (_IntegrationMethod == FFDIM_RKEH12) RKEH12::JacobianDOFs(this, jac, x, y, z, .0, T, 0.5 * dt, 2.0 * dt, SVFFD_RKTOL);
    else if (_IntegrationMethod == FFDIM_RKBS23) RKBS23::JacobianDOFs(this, jac, x, y, z, .0, T, 0.5 * dt, 2.0 * dt, SVFFD_RKTOL);
    else if (_IntegrationMethod == FFDIM_RKF45)  RKF45 ::JacobianDOFs(this, jac, x, y, z, .0, T, 0.5 * dt, 2.0 * dt, SVFFD_RKTOL);
    else if (_IntegrationMethod == FFDIM_RKCK45) RKCK45::JacobianDOFs(this, jac, x, y, z, .0, T, 0.5 * dt, 2.0 * dt, SVFFD_RKTOL);
    else if (_IntegrationMethod == FFDIM_RKDP45) RKDP45::JacobianDOFs(this, jac, x, y, z, .0, T, 0.5 * dt, 2.0 * dt, SVFFD_RKTOL);
    else {
      cerr << "BSplineFreeFormTransformationSV::JacobianDOFs: Unknown integration method: " << _IntegrationMethod << endl;
      exit(1);
    }
  }
}

// -----------------------------------------------------------------------------
inline void BSplineFreeFormTransformationSV
::JacobianDetDerivative(double [3], const Matrix &, int, double, double, double, double, double, bool, bool) const
{
  Throw(ERR_NotImplemented, __FUNCTION__, "Not implemented");
}

namespace BSplineFreeFormTransformationSVUtils {

// -----------------------------------------------------------------------------
struct MultiplyDerivatives : public VoxelFunction
{
  const int x, y, z, xx, xy, xz, yx, yy, yz, zx, zy, zz; // offsets

  MultiplyDerivatives(int n)
  :
    x(0), y(x+n), z(y+n),
    xx(   0), xy(xx+n), xz(xy+n),
    yx(xz+n), yy(yx+n), yz(yy+n),
    zx(yz+n), zy(zx+n), zz(zy+n)
  {}

  void operator ()(const GenericImage<double> &, int, const double *in, const double *d, double *out) const
  {
    // Attention: out can be equal d, therefore use temporary variables
    double gx = in[x] * d[xx] + in[x] * d[xy] + in[x] * d[xz];
    double gy = in[y] * d[yx] + in[y] * d[yy] + in[y] * d[yz];
    double gz = in[z] * d[zx] + in[z] * d[zy] + in[z] * d[zz];
    out[x] = gx, out[y] = gy, out[z] = gz;
  }
};


} // namespace BSplineFreeFormTransformationSVUtils
using namespace BSplineFreeFormTransformationSVUtils;

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformationSV
::ParametricGradient(const GenericImage<double> *in, double *out,
                     const WorldCoordsImage *i2w, const WorldCoordsImage *wc,
                     double t, double t0, double w) const
{
  // Upper integration limit for given interval
  const double T = UpperIntegrationLimit(t, t0);
  if (IsZero(T)) return;

  // ---------------------------------------------------------------------------
  // BCH based velocity update computation
  if (_NumberOfBCHTerms > 1) {

    MIRTK_START_TIMING();
    if (_UseDenseBCHGrid && _NumberOfBCHTerms > 2) {
      GenericImage<double> v(in->Attributes(), 3);
      this->Velocity(v);
      EvaluateBCHFormulaDense<double>(_NumberOfBCHTerms, v, v, *in, true, _LieDerivative);
      BSplineFreeFormTransformation3D::ParametricGradient(&v, out, i2w, wc, t0, w / T);
   } else {
      // Compute spline coefficients of update velocity field
      CPImage u(this->Attributes());
      DOFValue * const grd = reinterpret_cast<DOFValue *>(u.Data());
      BSplineFreeFormTransformation3D::ParametricGradient(in, grd, i2w, wc, t0, 1.);
      // Approximate velocity spline coefficients of composite transformation
      // using Baker-Campbell-Hausdorff (BCH) formula and subtract current coefficients
      EvaluateBCHFormula(_NumberOfBCHTerms, u, T, _CPImage, 1., u, true);
      // Adjust weight as update field is computed for tau * v, i.e.,
      //   exp(tau * v_{i+1}) = exp(tau v_i) o exp(\delta u)
      //   ==> v_{i+1} = log(exp(tau * v_{i+1})) / tau
      w /= T;
      // Add weighted gradient to total energy gradient
      for (int dof = 0; dof < this->NumberOfDOFs(); ++dof) {
        if (_Status[dof] == Active) out[dof] += w * grd[dof];
      }
    }
    MIRTK_DEBUG_TIMING(2, "parametric gradient computation (BCH, dense=" << (_UseDenseBCHGrid ? "yes" : "no") << ")");

  // ---------------------------------------------------------------------------
  // Scaling and squaring based gradient computation
  } else if (_IntegrationMethod == FFDIM_SS || _IntegrationMethod == FFDIM_FastSS) {

    MIRTK_START_TIMING();

    // Copy B-spline coefficients of velocity field
    GenericImage<double> v(this->Attributes(), 3);
    double *vx = v.Data(0, 0, 0, 0);
    double *vy = v.Data(0, 0, 0, 1);
    double *vz = v.Data(0, 0, 0, 2);
    const Vector *vp = _CPImage.Data();
    for (int idx = 0; idx < _CPImage.NumberOfVoxels(); ++idx, ++vx, ++vy, ++vz, ++vp) {
      *vx = vp->_x;
      *vy = vp->_y;
      *vz = vp->_z;
    }

    // Exponentiate non-parametric gradient using inverse mapping
    // (cf. Modat et al., "Parametric non-rigid registration using a stationary velocity field", MMBIA, 145–150, 2012)
    GenericImage<double> dv;
    mirtk::ScalingAndSquaring<double> exp;
    exp.UpperIntegrationLimit(T);
    exp.ComputeInverse(true);
    exp.NumberOfSteps(NumberOfStepsForIntervalLength(T));
    exp.MaxScaledVelocity(_MaxScaledVelocity);
    exp.InterimAttributes(_IntegrationMethod == FFDIM_FastSS ? this->Attributes() : in->Attributes());
    exp.OutputAttributes(exp.InterimAttributes());
    exp.InputVelocity(&v);
    exp.ComputeInterpolationCoefficients(false);
    exp.InputGradient(in);
    exp.OutputGradient(&dv);
    exp.Run();

    // Free copy of velocity coefficients
    v.Clear();

    // Multiply resulting vectors by derivative of v w.r.t. the DoFs
    BSplineFreeFormTransformation3D::ParametricGradient(&dv, out, i2w, wc, t0, w);

    MIRTK_DEBUG_TIMING(2, "parametric gradient computation" << (_IntegrationMethod == FFDIM_FastSS ? " (FastSS)" : " (SS)"));

  // ---------------------------------------------------------------------------
  // Runge-Kutta integration based gradient computation similar to TD FFD
  // transformation parameterized by non-stationary velocity field
  } else {

    // Note: T = in->GetTOrigin() - t0
    FreeFormTransformation::ParametricGradient(in, out, i2w, wc, in->GetTOrigin() - T, w);

  }
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformationSV
::ParametricGradient(const PointSet &pos, const Vector3D<double> *in,
                     double *out, double t, double t0, double w) const
{
  // Runge-Kutta integration based gradient computation similar to TD FFD
  // transformation parameterized by non-stationary velocity field
  FreeFormTransformation::ParametricGradient(pos, in, out, t, t0, w);
}

// =============================================================================
// I/O
// =============================================================================

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformationSV::Print(ostream &os, Indent indent) const
{
  os << indent << "B-spline SV FFD:" << endl;
  indent++;
  // Print FFD attributes
  FreeFormTransformation3D::Print(os, indent);
  // Change output stream settings
  const streamsize    w = os.width    (0);
  const streamsize    p = os.precision(2);
  const ios::fmtflags f = os.flags();
  cout.unsetf(ios::floatfield);
  // Print SV FFD parameters
  os << indent << "Integration method:                " << setw(6) << ToString(_IntegrationMethod) << endl;
  os << indent << "Cross-sectional time interval:     " << setw(6) << _T << endl;
  os << indent << "Time unit of integration interval: " << setw(6) << _TimeUnit << endl;
  os << indent << "Maximum scaled velocity:           " << setw(6) << _MaxScaledVelocity << endl;
  os << indent << "No. of integration steps per unit: " << setw(6) << _NumberOfSteps << endl;
  os << indent << "No. of cross-sectional steps:      " << setw(6) << NumberOfStepsForIntervalLength(_T) << endl;
  os << indent << "No. of BCH terms:                  " << setw(6) << _NumberOfBCHTerms << endl;
  os << indent << "Use Lie derivative:                " << setw(6) << ToString(_LieDerivative) << endl;
  // Restore output stream settings
  os.width    (w);
  os.precision(p);
  os.flags    (f);
}

// -----------------------------------------------------------------------------
bool BSplineFreeFormTransformationSV::CanRead(TransformationType format) const
{
  switch (format) {
    case TRANSFORMATION_BSPLINE_FFD_SV_v1:
    case TRANSFORMATION_BSPLINE_FFD_SV_v2:
    case TRANSFORMATION_BSPLINE_FFD_SV_v3:
    case TRANSFORMATION_BSPLINE_FFD_SV_v4:
    case TRANSFORMATION_BSPLINE_FFD_SV_v5:
    case TRANSFORMATION_BSPLINE_FFD_SV_v6:
    case TRANSFORMATION_BSPLINE_FFD_SV_v7:
    case TRANSFORMATION_BSPLINE_FFD_SV_v8:
      return true;
    default:
      return false;
  }
}

// -----------------------------------------------------------------------------
Cifstream &BSplineFreeFormTransformationSV::ReadDOFs(Cifstream &from, TransformationType format)
{
  // Read FFD data
  switch (format) {
    case TRANSFORMATION_BSPLINE_FFD_SV_v1:
    case TRANSFORMATION_BSPLINE_FFD_SV_v2:
    case TRANSFORMATION_BSPLINE_FFD_SV_v3:
    case TRANSFORMATION_BSPLINE_FFD_SV_v4:
    case TRANSFORMATION_BSPLINE_FFD_SV_v5:
    case TRANSFORMATION_BSPLINE_FFD_SV_v6:
      BSplineFreeFormTransformation3D::ReadDOFs(from, TRANSFORMATION_BSPLINE_FFD_3D_v2);
    case TRANSFORMATION_BSPLINE_FFD_SV_v7:
      BSplineFreeFormTransformation3D::ReadDOFs(from, TRANSFORMATION_BSPLINE_FFD_3D_v3);
    default:
      BSplineFreeFormTransformation3D::ReadDOFs(from, TRANSFORMATION_BSPLINE_FFD_3D);
  }

  // Read number of integration steps
  from.ReadAsInt(&_NumberOfSteps, 1);

  if (format == TRANSFORMATION_BSPLINE_FFD_SV_v1) return from;

  // Read upper integration limit
  from.ReadAsDouble(&_T, 1);
  // Read number of BCH terms to use for update
  from.ReadAsInt(&_NumberOfBCHTerms, 1);

  if (format == TRANSFORMATION_BSPLINE_FFD_SV_v2) return from;

  // Read time unit of integration interval
  from.ReadAsDouble(&_TimeUnit, 1);

  if (format == TRANSFORMATION_BSPLINE_FFD_SV_v3) return from;

  if (format <= TRANSFORMATION_BSPLINE_FFD_SV_v6) {

    // Whether to use scaling and squaring
    char useSS;
    from.ReadAsChar(&useSS, 1);

    // Maximum scaled velocity
    from.ReadAsDouble(&_MaxScaledVelocity, 1);

    if (format == TRANSFORMATION_BSPLINE_FFD_SV_v4) return from;

    // Whether to use fast scaling and squaring
    char fastSS;
    from.ReadAsChar(&fastSS, 1);

    // Set integration method
    _IntegrationMethod = (useSS ? (fastSS ? FFDIM_FastSS : FFDIM_SS) : FFDIM_RKE1);

  } else {

    // Integration method
    unsigned int integration_method;
    from.ReadAsUInt(&integration_method, 1);
    _IntegrationMethod = static_cast<FFDIntegrationMethod>(integration_method);

    // Maximum scaled velocity
    from.ReadAsDouble(&_MaxScaledVelocity, 1);

  }

  return from;
}

// -----------------------------------------------------------------------------
Cofstream &BSplineFreeFormTransformationSV::WriteDOFs(Cofstream &to) const
{
  // Write FFD data
  BSplineFreeFormTransformation3D::WriteDOFs(to);

  // Write number of integration steps
  to.WriteAsInt(&_NumberOfSteps, 1);
  // Write upper integration limit
  to.WriteAsDouble(&_T, 1);
  // Write number of BCH terms to use for update
  to.WriteAsInt(&_NumberOfBCHTerms, 1);
  // Write time unit of integration interval
  to.WriteAsDouble(&_TimeUnit, 1);
  // Integration method
  const unsigned int integration_method = _IntegrationMethod;
  to.WriteAsUInt(&integration_method, 1);
  // Maximum scaled velocity
  to.WriteAsDouble(&_MaxScaledVelocity, 1);

  return to;
}


} // namespace mirtk
