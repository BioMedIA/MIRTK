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

#include "mirtk/TransformationApproximationError.h"

#include "mirtk/Math.h"
#include "mirtk/Memory.h"
#include "mirtk/PointSet.h"
#include "mirtk/Vector3D.h"
#include "mirtk/Parallel.h"
#include "mirtk/Transformation.h"
#include "mirtk/OrderedSet.h"


namespace mirtk {


// =============================================================================
// Auxiliary functors
// =============================================================================

namespace TransformationApproximationErrorUtils {


// -----------------------------------------------------------------------------
/// Evaluate mean squared error of current approximation
struct Evaluate
{
  PointSet *_Target;
  PointSet *_Source;
  double    _Sum;

  Evaluate() : _Sum(.0) {}

  Evaluate(const Evaluate &other, split)
  :
    _Target(other._Target),
    _Source(other._Source),
    _Sum(.0)
  {}

  void join(const Evaluate &other)
  {
    _Sum += other._Sum;
  }

  void operator ()(const blocked_range<int> &re)
  {
    double dx, dy, dz;
    for (int n = re.begin(); n != re.end(); ++n) {
      const Point &p1 = _Target->GetPoint(n);
      const Point &p2 = _Source->GetPoint(n);
      dx = p1._x - p2._x;
      dy = p1._y - p2._y;
      dz = p1._z - p2._z;
      _Sum += dx * dx + dy * dy + dz * dz;
    }
  }
};

// -----------------------------------------------------------------------------
/// Evaluate non-parametric gradient of mean squared error of current approximation
struct EvaluateGradient
{
  PointSet         *_Target;
  PointSet         *_Source;
  Vector3D<double> *_Gradient;
  double            _Norm;

  void operator ()(const blocked_range<int> &re) const
  {
    for (int n = re.begin(); n != re.end(); ++n) {
      const Point &p1 = _Target->GetPoint(n);
      const Point &p2 = _Source->GetPoint(n);
      _Gradient[n]._x = _Norm * (p1._x - p2._x);
      _Gradient[n]._y = _Norm * (p1._y - p2._y);
      _Gradient[n]._z = _Norm * (p1._z - p2._z);
    }
  }
};


} // namespace TransformationApproximationErrorUtils

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
TransformationApproximationError
::TransformationApproximationError(class Transformation *transformation,
    const double *x,  const double *y,  const double *z, const double *t,
    const double *dx, const double *dy, const double *dz, int no)
:
  _Transformation(transformation),
  _NumberOfTimePoints(0),
  _NumberOfPoints(no),
  _Gradient(NULL)
{
  // Determine distinct time points
  OrderedSet<double> times;
  for (int n = 0; n < no; ++n) times.insert(t[n]);
  _NumberOfTimePoints = static_cast<int>(times.size());

  // Initialize input point sets
  _TargetTime.resize(_NumberOfTimePoints);
  _Target    .resize(_NumberOfTimePoints);
  _Current   .resize(_NumberOfTimePoints);
  _Source    .resize(_NumberOfTimePoints);

  int max_no = 0, cur_no;
  OrderedSet<double>::const_iterator time = times.begin();
  for (int l = 0; l < _NumberOfTimePoints; ++l, ++time) {
    _TargetTime[l] = *time;
    cur_no = 0;
    for (int n = 0; n < no; ++n) {
      if (fequal(t[n], _TargetTime[l])) ++cur_no;
    }
    max_no = max(max_no, cur_no);
    _Target[l].Reserve(cur_no);
    _Source[l].Reserve(cur_no);
    for (int n = 0; n < no; ++n) {
      if (fequal(t[n], _TargetTime[l])) {
        _Target[l].Add(Point(x[n],         y[n],         z[n]));
        _Source[l].Add(Point(x[n] + dx[n], y[n] + dy[n], z[n] + dz[n]));
      }
    }
  }

  // Initialize transformed target point set (cf. Update)
  _Current.resize(_NumberOfTimePoints);
  if (_NumberOfTimePoints > 0) _Transformation->Changed(true);

  // Allocate memory for non-parametric gradient
  if (max_no > 0) Allocate(_Gradient, max_no);
}

// -----------------------------------------------------------------------------
TransformationApproximationError::~TransformationApproximationError()
{
  Deallocate(_Gradient);
}

// -----------------------------------------------------------------------------
void TransformationApproximationError::CenterPoints()
{
  _TargetCenter = _SourceCenter = Point(.0, .0, .0);
  for (int l = 0; l < _NumberOfTimePoints; ++l) {
    for (int n = 0; n < _Source[l].Size(); ++n) {
      _TargetCenter += _Target[l](n) / _NumberOfPoints;
      _SourceCenter += _Source[l](n) / _NumberOfPoints;
    }
  }
  for (int l = 0; l < _NumberOfTimePoints; ++l) {
    for (int n = 0; n < _Source[l].Size(); ++n) {
      _Target[l](n) -= _TargetCenter;
      _Source[l](n) -= _SourceCenter;
    }
  }
}

// =============================================================================
// Parameters (DoFs)
// =============================================================================

// -----------------------------------------------------------------------------
int TransformationApproximationError::NumberOfDOFs() const
{
  return _Transformation->NumberOfDOFs();
}

// -----------------------------------------------------------------------------
void TransformationApproximationError::Put(const double *x)
{
  _Transformation->Put(x);
}

// -----------------------------------------------------------------------------
double TransformationApproximationError::Get(int i) const
{
  return _Transformation->Get(i);
}

// -----------------------------------------------------------------------------
void TransformationApproximationError::Get(double *x) const
{
  _Transformation->Get(x);
}

// -----------------------------------------------------------------------------
double TransformationApproximationError::Step(double *dx)
{
  return _Transformation->Update(dx);
}

// -----------------------------------------------------------------------------
void TransformationApproximationError::Update(bool)
{
  if (_Transformation->Changed()) {
    for (int l = 0; l < _NumberOfTimePoints; ++l) {
      _Current[l] = _Target[l];
      _Transformation->Transform(_Current[l], _TargetTime[l]);
    }
    _Transformation->Changed(false);
  }
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
double TransformationApproximationError::Value()
{
  double sum = .0;
  int    num =  0;
  for (int l = 0; l < _NumberOfTimePoints; ++l) {
    TransformationApproximationErrorUtils::Evaluate eval;
    eval._Target = &_Current[l];
    eval._Source = &_Source [l];
    parallel_reduce(blocked_range<int>(0, _Source[l].Size()), eval);
    sum += eval._Sum;
    num += _Source[l].Size();
  }
  return (num > 0 ? sum / num : .0);
}

// -----------------------------------------------------------------------------
void TransformationApproximationError::Gradient(double *dx, double, bool *)
{
  memset(dx, 0, this->NumberOfDOFs() * sizeof(double));
  TransformationApproximationErrorUtils::EvaluateGradient eval;
  eval._Gradient = _Gradient;
  eval._Norm     = 2.0 / _NumberOfPoints;
  for (int l = 0; l < _NumberOfTimePoints; ++l) {
    eval._Target = &_Current[l];
    eval._Source = &_Source [l];
    parallel_for(blocked_range<int>(0, _Target[l].Size()), eval);
    _Transformation->ParametricGradient(_Target[l], _Gradient, dx, _TargetTime[l]);
  }
}

// -----------------------------------------------------------------------------
double TransformationApproximationError::GradientNorm(const double *dx) const
{
  return _Transformation->DOFGradientNorm(dx);
}


} // namespace mirtk
