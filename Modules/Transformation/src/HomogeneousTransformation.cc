/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2015 Imperial College London
 * Copyright 2008-2013 Daniel Rueckert, Julia Schnabel
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

#include "mirtk/HomogeneousTransformation.h"

#include "mirtk/Math.h"
#include "mirtk/Memory.h"
#include "mirtk/PointSamples.h"
#include "mirtk/AdaptiveLineSearch.h"
#include "mirtk/ConjugateGradientDescent.h"
#include "mirtk/TransformationApproximationError.h"
#include "mirtk/Profiling.h"


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
// Note: DoFs are parameters (e.g., rotation) from which 4x4 matrix is constructed
HomogeneousTransformation::HomogeneousTransformation(int ndofs)
:
  Transformation(ndofs),
  _matrix (4, 4),
  _inverse(4, 4)
{
  _matrix .Ident();
  _inverse.Ident();
}

// -----------------------------------------------------------------------------
// Note: DoFs are parameters (e.g., rotation) from which 4x4 matrix is constructed
HomogeneousTransformation::HomogeneousTransformation(const HomogeneousTransformation &t, int ndofs)
:
  Transformation(t, ndofs),
  _matrix (t._matrix),
  _inverse(t._inverse)
{
}

// -----------------------------------------------------------------------------
// Note: DoFs are the matrix elements
HomogeneousTransformation::HomogeneousTransformation()
:
  Transformation(16),
  _inverse(4, 4)
{
  // Use memory allocated by Transformation when possible
  if (sizeof(DOFValue) == sizeof(Matrix::ElementType)) {
    _matrix.Initialize(4, 4, _Param);
  } else {
    _matrix.Initialize(4, 4);
  }
  _matrix .Ident();
  _inverse.Ident();
  // Last row should always be "0 0 0 1"
  _Status[12] = Passive;
  _Status[13] = Passive;
  _Status[14] = Passive;
  _Status[15] = Passive;
}

// -----------------------------------------------------------------------------
// Note: DoFs are the matrix elements
HomogeneousTransformation::HomogeneousTransformation(const Matrix &matrix)
:
  Transformation(16),
  _matrix (4, 4, _Param), // uses memory allocated by Transformation
  _inverse(4, 4)
{
  _matrix  = matrix;
  _inverse = matrix.Inverse();
  // Last row should always be "0 0 0 1"
  _Status[12] = Passive;
  _Status[13] = Passive;
  _Status[14] = Passive;
  _Status[15] = Passive;
}

// -----------------------------------------------------------------------------
// Note: DoFs are either matrix elements if other transformation is an instance
//       of HomogeneousTransformation itself and not a subclass or
//       parameters (e.g., rotation) from which 4x4 matrix is constructed
HomogeneousTransformation::HomogeneousTransformation(const HomogeneousTransformation &t)
:
  Transformation(t)
{
  if (t.GetMatrix().GetPointerToElements() == t._Param) {
    _matrix.Initialize(4, 4, _Param);
  } else {
    _matrix = t._matrix;
  }
  _inverse = t._inverse;
}

// -----------------------------------------------------------------------------
HomogeneousTransformation::~HomogeneousTransformation()
{
}

// =============================================================================
// Approximation
// =============================================================================

// -----------------------------------------------------------------------------
double HomogeneousTransformation
::Approximate(const ImageAttributes &domain, double *dx, double *dy, double *dz,
              int niter, double max_error)
{
  // Check input arguments
  const int no = domain.NumberOfPoints();
  if (no <= 0) return .0;
  double error = numeric_limits<double>::infinity();
  if (niter < 1) return error;

  // Compute world coordinates of lattice points
  double *x = Allocate<double>(no);
  double *y = Allocate<double>(no);
  double *z = Allocate<double>(no);
  double *t = Allocate<double>(no);
  domain.LatticeToWorld(x, y, z, t);

  // Allocate memory for transformed world coordinates
  double *wx = Allocate<double>(no);
  double *wy = Allocate<double>(no);
  double *wz = Allocate<double>(no);

  // Copy original displacements
  double *tx = Allocate<double>(no);
  double *ty = Allocate<double>(no);
  double *tz = Allocate<double>(no);
  memcpy(tx, dx, no * sizeof(double));
  memcpy(ty, dy, no * sizeof(double));
  memcpy(tz, dz, no * sizeof(double));

  // Evaluate approximation error and residual displacements
  if (this->RequiresCachingOfDisplacements()) {
    error = EvaluateRMSError(domain, dx, dy, dz);
  } else {
    error = EvaluateRMSError(x, y, z, t, dx, dy, dz, no);
  }

  // Repeat approximation n times or until error drops below a threshold
  for (int iter = 0; iter < niter && error > max_error; ++iter) {

    // Transform world coordinates
    memcpy(wx, x, no * sizeof(double));
    memcpy(wy, y, no * sizeof(double));
    memcpy(wz, z, no * sizeof(double));
    this->Transform(no, wx, wy, wz, t);

    // Approximate residual displacements by new parameters
    const Matrix matrix = _matrix;
    this->ApproximateDOFs(wx, wy, wz, t, dx, dy, dz, no);

    // Compose with previous transformation
    this->PutMatrix(_matrix * matrix);

    // Evaluate approximation error and residual displacements
    memcpy(dx, tx, no * sizeof(double));
    memcpy(dy, ty, no * sizeof(double));
    memcpy(dz, tz, no * sizeof(double));
    if (this->RequiresCachingOfDisplacements()) {
      error = EvaluateRMSError(domain, dx, dy, dz);
    } else {
      error = EvaluateRMSError(x, y, z, t, dx, dy, dz, no);
    }
  }

  // Free memory
  Deallocate(tx);
  Deallocate(ty);
  Deallocate(tz);
  Deallocate(wx);
  Deallocate(wy);
  Deallocate(wz);
  Deallocate(x);
  Deallocate(y);
  Deallocate(z);
  Deallocate(t);

  return error;
}

// -----------------------------------------------------------------------------
double HomogeneousTransformation
::Approximate(const double *x,  const double *y,  const double *z,
              double       *dx, double       *dy, double       *dz, int no,
              int niter, double max_error)
{
  // Check input arguments
  if (no <= 0) return .0;
  double error = numeric_limits<double>::infinity();
  if (niter < 1) return error;

  // Compute (fixed) time coordinates
  const double t0 = .0;
  double *t = CAllocate<double>(no, &t0);

  // Allocate memory for transformed world coordinates
  double *wx = Allocate<double>(no);
  double *wy = Allocate<double>(no);
  double *wz = Allocate<double>(no);

  // Copy original displacements
  double *tx = Allocate<double>(no);
  double *ty = Allocate<double>(no);
  double *tz = Allocate<double>(no);
  memcpy(tx, dx, no * sizeof(double));
  memcpy(ty, dy, no * sizeof(double));
  memcpy(tz, dz, no * sizeof(double));

  // Evaluate approximation error and residual displacements
  error = EvaluateRMSError(x, y, z, t0, dx, dy, dz, no);

  // Repeat approximation n times or until error drops below a threshold
  for (int iter = 0; iter < niter && error > max_error; ++iter) {

    // Transform world coordinates
    memcpy(wx, x, no * sizeof(double));
    memcpy(wy, y, no * sizeof(double));
    memcpy(wz, z, no * sizeof(double));
    this->Transform(no, wx, wy, wz, t0);

    // Approximate residual displacements by new parameters
    const Matrix matrix = _matrix;
    this->ApproximateDOFs(wx, wy, wz, t, dx, dy, dz, no);

    // Compose with previous transformation
    this->PutMatrix(_matrix * matrix);

    // Evaluate error of approximation and residual displacements
    memcpy(dx, tx, no * sizeof(double));
    memcpy(dy, ty, no * sizeof(double));
    memcpy(dz, tz, no * sizeof(double));
    error = EvaluateRMSError(x, y, z, t0, dx, dy, dz, no);
  }

  // Free memory
  Deallocate(tx);
  Deallocate(ty);
  Deallocate(tz);
  Deallocate(wx);
  Deallocate(wy);
  Deallocate(wz);
  Deallocate(t);

  return error;
}

// -----------------------------------------------------------------------------
double HomogeneousTransformation
::Approximate(const double *x,  const double *y,  const double *z,  const double *t,
              double       *dx, double       *dy, double       *dz, int no,
              int niter, double max_error)
{
  // Check input arguments
  if (no <= 0) return .0;
  double error = numeric_limits<double>::infinity();
  if (niter < 1) return error;

  // Allocate memory for transformed world coordinates
  double *wx = Allocate<double>(no);
  double *wy = Allocate<double>(no);
  double *wz = Allocate<double>(no);

  // Copy original displacements
  double *tx = Allocate<double>(no);
  double *ty = Allocate<double>(no);
  double *tz = Allocate<double>(no);
  memcpy(tx, dx, no * sizeof(double));
  memcpy(ty, dy, no * sizeof(double));
  memcpy(tz, dz, no * sizeof(double));

  // Evaluate error of approximation and residual displacements
  error = EvaluateRMSError(x, y, z, t, dx, dy, dz, no);

  // Repeat approximation n times or until error drops below a threshold
  for (int iter = 0; iter < niter; ++iter) {

    // Transform world coordinates
    memcpy(wx, x, no * sizeof(double));
    memcpy(wy, y, no * sizeof(double));
    memcpy(wz, z, no * sizeof(double));
    this->Transform(no, wx, wy, wz, t);

    // Approximate residual displacements by new parameters
    const Matrix matrix = _matrix;
    this->ApproximateDOFs(wx, wy, wz, t, dx, dy, dz, no);

    // Compose with previous transformation
    this->PutMatrix(_matrix * matrix);

    // Evaluate error of approximation and residual displacements
    memcpy(dx, tx, no * sizeof(double));
    memcpy(dy, ty, no * sizeof(double));
    memcpy(dz, tz, no * sizeof(double));
    error = EvaluateRMSError(x, y, z, t, dx, dy, dz, no);
  }

  // Free memory
  Deallocate(tx);
  Deallocate(ty);
  Deallocate(tz);
  Deallocate(wx);
  Deallocate(wy);
  Deallocate(wz);

  return error;
}

// -----------------------------------------------------------------------------
double HomogeneousTransformation::Approximate(const Matrix &matrix)
{
  return this->ApproximateAsNew(matrix);
}

// -----------------------------------------------------------------------------
double HomogeneousTransformation::ApproximateAsNew(const Matrix &matrix)
{
  if (this->NumberOfActiveDOFs() == 12) {
    this->PutMatrix(matrix);
    return .0;
  }

  const int    nsamples = 100;
  const double t0       =  .0;

  PointSamples samples(nsamples, /*fixed seed=*/ 0);
  samples.SampleSphere(.0, 100.0);

  double *x  =  Allocate<double>(nsamples);
  double *y  =  Allocate<double>(nsamples);
  double *z  =  Allocate<double>(nsamples);
  double *t  = CAllocate<double>(nsamples, &t0);
  double *dx =  Allocate<double>(nsamples);
  double *dy =  Allocate<double>(nsamples);
  double *dz =  Allocate<double>(nsamples);

  for (int i = 0; i < nsamples; ++i) {
    x [i] = samples(i)._x;
    y [i] = samples(i)._y;
    z [i] = samples(i)._z;
    dx[i] = matrix(0, 0) * x[i] + matrix(0, 1) * y[i] + matrix(0, 2) * z[i] + matrix(0, 3) - x[i];
    dy[i] = matrix(1, 0) * x[i] + matrix(1, 1) * y[i] + matrix(1, 2) * z[i] + matrix(1, 3) - y[i];
    dz[i] = matrix(2, 0) * x[i] + matrix(2, 1) * y[i] + matrix(2, 2) * z[i] + matrix(2, 3) - z[i];
  }

  double rms = this->ApproximateAsNew(x, y, z, t, dx, dy, dz, nsamples);

  Deallocate(x);
  Deallocate(y);
  Deallocate(z);
  Deallocate(t);
  Deallocate(dx);
  Deallocate(dy);
  Deallocate(dz);

  return rms;
}

// -----------------------------------------------------------------------------
void HomogeneousTransformation
::ApproximateDOFs(const double *x,  const double *y,  const double *z, const double *t,
                  const double *dx, const double *dy, const double *dz, int no)
{
  MIRTK_START_TIMING();

  // Note: Do not reset transformation parameters as this function is usually
  //       called by subclass implementations after these were initialized!
  Matrix pre(4, 4), post(4, 4);

  // Mean squared error function
  TransformationApproximationError error(this, x, y, z, t, dx, dy, dz, no);

  // Center point sets
  if (_Status[TX] == Active && _Status[TY] == Active && _Status[TZ] == Active) {
    error.CenterPoints();
  }

  // Adjust transformation
  pre.Ident();
  pre(0, 3)  = + error.TargetCenter()._x;
  pre(1, 3)  = + error.TargetCenter()._y;
  pre(2, 3)  = + error.TargetCenter()._z;

  post.Ident();
  post(0, 3) = - error.SourceCenter()._x;
  post(1, 3) = - error.SourceCenter()._y;
  post(2, 3) = - error.SourceCenter()._z;

  this->PutMatrix(post * this->GetMatrix() * pre);

  // Optimization method
  AdaptiveLineSearch linesearch;
  linesearch.MaxRejectedStreak(0);
  linesearch.StrictStepLengthRange(false);
  linesearch.ReusePreviousStepLength(true);
  linesearch.NumberOfIterations(20);
  linesearch.MinStepLength(1e-6);
  linesearch.MaxStepLength(10.0);

  ConjugateGradientDescent optimizer;
  optimizer.Function(&error);
  optimizer.LineSearch(&linesearch);
  optimizer.NumberOfSteps(100);
  optimizer.Epsilon(1e-6);
  optimizer.Delta(1e-12);

  // Find transformation parameters which minimize the approximation error
  optimizer.Run();
  optimizer.ConjugateGradientOff();
  optimizer.Run();

  // Include centering transformations in final transformation
  pre.Ident();
  pre(0, 3)  = - error.TargetCenter()._x;
  pre(1, 3)  = - error.TargetCenter()._y;
  pre(2, 3)  = - error.TargetCenter()._z;

  post.Ident();
  post(0, 3) = + error.SourceCenter()._x;
  post(1, 3) = + error.SourceCenter()._y;
  post(2, 3) = + error.SourceCenter()._z;

  this->PutMatrix(post * this->GetMatrix() * pre);
  MIRTK_DEBUG_TIMING(5, "HomogeneousTransformation::ApproximateDOFs");
}

// =============================================================================
// Transformation parameters
// =============================================================================

// -----------------------------------------------------------------------------
void HomogeneousTransformation::UpdateMatrix()
{
  if (_matrix.GetPointerToElements() != _Param) {
    if (NumberOfDOFs() != 16) {
      cerr << "HomogeneousTransformation::UpdateMatrix: Override in subclass!" << endl;
      exit(1);
    }
    memcpy(_matrix.GetPointerToElements(), _Param, 16 * sizeof(double));
  }
  _inverse = _matrix.Inverse();
}

// -----------------------------------------------------------------------------
void HomogeneousTransformation::UpdateDOFs()
{
  if (_matrix.GetPointerToElements() != _Param) {
    if (NumberOfDOFs() != 16) {
      cerr << "HomogeneousTransformation::UpdateDOFs: Override in subclass!" << endl;
      exit(1);
    }
    memcpy(_Param, _matrix.GetPointerToElements(), 16 * sizeof(double));
  }
}

// ---------------------------------------------------------------------------
bool HomogeneousTransformation::CopyFrom(const Transformation *other)
{
  const HomogeneousTransformation *linear;
  if ((linear = dynamic_cast<const HomogeneousTransformation *>(other))) {
    this->Reset();
    const int ndofs = min(_NumberOfDOFs, other->NumberOfDOFs());
    for (int dof = 0; dof < ndofs; ++dof) {
      if (_Status[dof] == Active) _Param[dof] = other->Get(dof);
    }
    this->Update(MATRIX);
    return true;
  } else {
    return false;
  }
}

// -----------------------------------------------------------------------------
void HomogeneousTransformation::Reset()
{
  _matrix.Ident();
  _inverse.Ident();
  this->Update(DOFS);
}

// -----------------------------------------------------------------------------
void HomogeneousTransformation::Invert()
{
  Matrix matrix = _matrix;
  _matrix  = _inverse;
  _inverse = matrix;
  this->Update(DOFS);
}

// =============================================================================
// Properties
// =============================================================================

// -----------------------------------------------------------------------------
bool HomogeneousTransformation::IsIdentity() const
{
  for (int i = 0; i < _matrix.Rows(); ++i)
  for (int j = 0; j < _matrix.Cols(); ++j) {
    if (_matrix(i, j) != static_cast<double>(i == j)) return false;
  }
  return true;
}

// =============================================================================
// I/O
// =============================================================================

// -----------------------------------------------------------------------------
void HomogeneousTransformation::Print(ostream &os, Indent indent) const
{
  _matrix.Print(os, indent);
}

// ---------------------------------------------------------------------------
Cifstream &HomogeneousTransformation::ReadDOFs(Cifstream &from, TransformationType format)
{
  Transformation::ReadDOFs(from, format);
  this->Update(MATRIX);
  return from;
}


} // namespace mirtk
