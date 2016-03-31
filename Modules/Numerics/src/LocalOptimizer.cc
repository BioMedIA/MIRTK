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

#include <mirtkLocalOptimizer.h>
#include <mirtkStoppingCriterion.h>

#include <mirtkMath.h>
#include <mirtkArray.h>
#include <mirtkMemory.h>
#include <mirtkObjectFactory.h>


namespace mirtk {


// =============================================================================
// Factory method
// =============================================================================

// -----------------------------------------------------------------------------
LocalOptimizer *LocalOptimizer::New(enum OptimizationMethod type, ObjectiveFunction *f)
{
  typedef ObjectFactory<enum OptimizationMethod, LocalOptimizer> Factory;
  LocalOptimizer *optimizer = Factory::Instance().New(type);
  if (!optimizer) {
    cerr << "LocalOptimizer::New: Unknown optimizer or optimizer not available: "
         << ToString(type) << " (code " << type << ")" << endl;
    exit(1);
  }
  optimizer->Function(f);
  return optimizer;
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
  _Delta(1e-12)
{
}

// -----------------------------------------------------------------------------
void LocalOptimizer::CopyAttributes(const LocalOptimizer &other)
{
  _Function      = other._Function;
  _NumberOfSteps = other._NumberOfSteps;
  _Epsilon       = other._Epsilon;
  _Delta         = other._Delta;

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

  // Initialize stopping criteria
  for (size_t n = 0; n < _StoppingCriteria.size(); ++n) {
    class StoppingCriterion *criterion = _StoppingCriteria[n];
    criterion->Optimizer(this);
    criterion->Function(_Function);
    criterion->Initialize();
  }
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
bool LocalOptimizer::Converged(int iter, double prev, double value, const double *dx)
{
  // Detect numerical problem of objective function implementation
  if (IsNaN(value)) {
    cerr << this->NameOfClass() << ": NaN objective function value!" << endl;
    exit(1);
  }

  // Check convergence towards local extremum of objective function
  if (!IsInf(prev)) {
    if (IsInf(value) || abs(prev - value) <= _Epsilon) return true;
  }

  // Check minimum change requirement
  if (_Function->GradientNorm(dx) <= _Delta) return true;

  // Test other stopping criteria
  const int ncriteria = static_cast<int>(_StoppingCriteria.size());
  for (int n = 0; n < ncriteria; ++n) {
    if (StoppingCriterion(n)->Fulfilled(iter, prev, value, dx)) return true;
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
  return false;
}

// -----------------------------------------------------------------------------
ParameterList LocalOptimizer::Parameter() const
{
  ParameterList params;
  Insert(params, "Maximum no. of iterations", _NumberOfSteps);
  Insert(params, "Epsilon", _Epsilon);
  Insert(params, "Delta", _Delta);
  return params;
}


} // namespace mirtk
