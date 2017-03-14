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

#include "mirtk/LocalOptimizer.h"
#include "mirtk/StoppingCriterion.h"

#include "mirtk/Math.h"
#include "mirtk/Array.h"
#include "mirtk/Memory.h"


namespace mirtk {


// =============================================================================
// Factory method
// =============================================================================

// Define singleton factory class for instantiation of local optimizers
mirtkDefineObjectFactory(OptimizationMethod, LocalOptimizer);

// -----------------------------------------------------------------------------
LocalOptimizer::FactoryType &LocalOptimizer::Factory()
{
  return LocalOptimizerFactory::Instance();
}

// -----------------------------------------------------------------------------
LocalOptimizer *LocalOptimizer::New(enum OptimizationMethod type, ObjectiveFunction *f)
{
  LocalOptimizer *optimizer = Factory().New(type);
  if (!optimizer) {
    cerr << "LocalOptimizer::New: Unknown optimizer or optimizer not available: "
         << ToString(type) << " (code " << type << ")" << endl;
    exit(1);
  }
  optimizer->Function(f);
  return optimizer;
}

// =============================================================================
// Line fitting
// =============================================================================

// -----------------------------------------------------------------------------
/// Compute slope of least squares fit of line to last n objective function values
///
/// \see https://www.che.udel.edu/pdf/FittingData.pdf
/// \see https://en.wikipedia.org/wiki/1_%2B_2_%2B_3_%2B_4_%2B_%E2%8B%AF
/// \see https://proofwiki.org/wiki/Sum_of_Sequence_of_Squares
double SlopeOfLeastSquaresFit(const Deque<double> &values)
{
  // sum_x1 divided by n as a slight modified to reduce no. of operations,
  // i.e., the other terms are divided by n as well by dropping one factor n
  const int n = static_cast<int>(values.size());
  if (n <  2) return NaN;
  if (n == 2) return values.back() - values.front();
  double sum_x1 = double(n + 1) / 2.;                     // sum of x / n
  double sum_x2 = n * double(n + 1) * (2. * n + 1) / 6.;  // sum of x^2
  double sum_y1 = 0.;
  double sum_xy = 0.;
  double x = 1.;
  for (auto y : values) {
    sum_y1 += y;
    sum_xy += x * y;
    x += 1.;
  }
  return (sum_xy - sum_x1 * sum_y1) / (sum_x2 - n * sum_x1 * sum_x1);
}

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
LocalOptimizer::LocalOptimizer(ObjectiveFunction *f)
:
  _Function(f),
  _NumberOfSteps(100),
  _Epsilon(1e-4),
  _Delta(1e-12),
  _NumberOfLastValues(2),
  _Converged(false)
{
}

// -----------------------------------------------------------------------------
void LocalOptimizer::CopyAttributes(const LocalOptimizer &other)
{
  _Function           = other._Function;
  _NumberOfSteps      = other._NumberOfSteps;
  _Epsilon            = other._Epsilon;
  _Delta              = other._Delta;
  _LastValues         = other._LastValues;
  _NumberOfLastValues = other._NumberOfLastValues;
  _Converged          = other._Converged;

  ClearStoppingCriteria();
  for (int i = 0; i < other.NumberOfStoppingCriteria(); ++i) {
    this->AddStoppingCriterion(other.StoppingCriterion(i)->New());
  }
}

// -----------------------------------------------------------------------------
LocalOptimizer::LocalOptimizer(const LocalOptimizer &other)
:
  Observable(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
LocalOptimizer &LocalOptimizer::operator =(const LocalOptimizer &other)
{
  if (this != &other) {
    Observable::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
LocalOptimizer::~LocalOptimizer()
{
  ClearStoppingCriteria();
}

// -----------------------------------------------------------------------------
void LocalOptimizer::Initialize()
{
  // Check objective function
  if (!_Function) {
    cerr << this->NameOfClass() << "::Initialize: Objective function not set" << endl;
    exit(1);
  }
  if (_Function->NumberOfDOFs() <= 0) {
    cerr << this->NameOfClass() << "::Initialize: Objective function has no free parameters (DoFs)" << endl;
    exit(1);
  }

  // Initialize container for last objective function values
  _LastValues.clear();
  _LastValuesSlope = NaN;
  if (_NumberOfLastValues < 2) {
    // Need at least the previous and current function value for _Epsilon stopping criterion
    _NumberOfLastValues = 2;
  }

  // Initialize stopping criteria
  for (size_t n = 0; n < _StoppingCriteria.size(); ++n) {
    class StoppingCriterion *criterion = _StoppingCriteria[n];
    criterion->Optimizer(this);
    criterion->Function(_Function);
    criterion->Initialize();
  }
  _Converged = false;
}

// =============================================================================
// Stopping criteria
// =============================================================================

// -----------------------------------------------------------------------------
int LocalOptimizer::NumberOfStoppingCriteria() const
{
  return static_cast<int>(_StoppingCriteria.size());
}

// -----------------------------------------------------------------------------
void LocalOptimizer::AddStoppingCriterion(class StoppingCriterion *criterion)
{
  Array<class StoppingCriterion *>::iterator it;
  for (it = _StoppingCriteria.begin(); it != _StoppingCriteria.end(); ++it) {
    if (*it == criterion) return;
  }
  _StoppingCriteria.push_back(criterion);
}

// -----------------------------------------------------------------------------
void LocalOptimizer::RemoveStoppingCriterion(class StoppingCriterion *criterion)
{
  Array<class StoppingCriterion *>::iterator it;
  for (it = _StoppingCriteria.begin(); it != _StoppingCriteria.end(); ++it) {
    if (*it == criterion) {
      _StoppingCriteria.erase(it);
      break;
    }
  }
}

// -----------------------------------------------------------------------------
StoppingCriterion *LocalOptimizer::StoppingCriterion(int n)
{
  return _StoppingCriteria[n];
}

// -----------------------------------------------------------------------------
const StoppingCriterion *LocalOptimizer::StoppingCriterion(int n) const
{
  return _StoppingCriteria[n];
}

// -----------------------------------------------------------------------------
void LocalOptimizer::ClearStoppingCriteria()
{
  for (size_t n = 0; n < _StoppingCriteria.size(); ++n) {
    Delete(_StoppingCriteria[n]);
  }
  _StoppingCriteria.clear();
}

// -----------------------------------------------------------------------------
bool LocalOptimizer::Converged(int iter, double value, const double *dx)
{
  // Detect numerical problem of objective function implementation
  if (IsNaN(value)) {
    cerr << this->NameOfClass() << ": NaN objective function value!" << endl;
    exit(1);
  }

  // Update history of last n function values
  double prev = NaN;
  if (!_LastValues.empty()) {
    prev = _LastValues.back();
  }
  _LastValues.push_back(value);
  if (_LastValues.size() > static_cast<size_t>(_NumberOfLastValues)) {
    _LastValues.pop_front();
  }

  // Objective function value changed from < inf to inf although the value may either
  // only be always inf in case of a deformable surface model, or always < inf otherwise.
  if (!IsNaN(prev) && !IsInf(prev) && IsInf(value)) {
    return true;
  }

  // Check convergence towards local extremum of objective function
  const double epsilon = (_Epsilon < 0. ? abs(_Epsilon * value) : _Epsilon);
  _LastValuesSlope = SlopeOfLeastSquaresFit(_LastValues);
  if (abs(_LastValuesSlope) < epsilon) return true;

  // Check minimum change requirement
  if (_Function->GradientNorm(dx) <= _Delta) return true;

  // Test other stopping criteria
  const int ncriteria = static_cast<int>(_StoppingCriteria.size());
  for (int n = 0; n < ncriteria; ++n) {
    if (StoppingCriterion(n)->Fulfilled(iter, value, dx)) return true;
  }
  return false;
}

// =============================================================================
// Parameters
// =============================================================================

// -----------------------------------------------------------------------------
bool LocalOptimizer::Set(const char *name, const char *value)
{
  if (strcmp(name, "Maximum no. of iterations")        == 0 ||
      strcmp(name, "Maximum number of iterations")     == 0 ||
      strcmp(name, "Maximum no. of gradient steps")    == 0 ||
      strcmp(name, "Maximum number of gradient steps") == 0 ||
      strcmp(name,    "No. of iterations")             == 0 ||
      strcmp(name, "Number of iterations")             == 0 ||
      strcmp(name,    "No. of gradient steps")         == 0 ||
      strcmp(name, "Number of gradient steps")         == 0) {
    return FromString(value, _NumberOfSteps);
  }
  if (strcmp(name, "Epsilon") == 0) {
    return FromString(value, _Epsilon);
  }
  if (strcmp(name, "Delta") == 0) {
    return FromString(value, _Delta) && _Delta >= .0;
  }
  if (strcmp(name, "No. of last function values") == 0 ||
      strcmp(name, "No. of past function values") == 0) {
    return FromString(value, _NumberOfLastValues);
  }
  return false;
}

// -----------------------------------------------------------------------------
ParameterList LocalOptimizer::Parameter() const
{
  ParameterList params;
  Insert(params, "Maximum no. of iterations", _NumberOfSteps);
  Insert(params, "Epsilon", _Epsilon);
  Insert(params, "Delta", _Delta);
  Insert(params, "No. of last function values", _NumberOfLastValues);
  return params;
}


} // namespace mirtk
