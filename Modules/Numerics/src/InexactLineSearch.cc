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

#include "mirtk/InexactLineSearch.h"

#include "mirtk/Memory.h"
#include "mirtk/String.h"


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
InexactLineSearch::InexactLineSearch(ObjectiveFunction *f)
:
  LineSearch(f),
  _MaxRejectedStreak      (-1),
  _ReusePreviousStepLength(true),
  _StrictStepLengthRange  (1),
  _AllowSignChange        (NULL),
  _CurrentDoFValues       (NULL),
  _ScaledDirection        (NULL)
{
  if (f) {
    Allocate(_CurrentDoFValues, f->NumberOfDOFs());
    Allocate(_ScaledDirection,  f->NumberOfDOFs());
  }
}

// -----------------------------------------------------------------------------
void InexactLineSearch::CopyAttributes(const InexactLineSearch &other)
{
  Deallocate(_CurrentDoFValues);
  Deallocate(_ScaledDirection);

  _MaxRejectedStreak       = other._MaxRejectedStreak;
  _ReusePreviousStepLength = other._ReusePreviousStepLength;
  _StrictStepLengthRange   = other._StrictStepLengthRange;
  _AllowSignChange         = other._AllowSignChange;

  if (Function()) {
    Allocate(_CurrentDoFValues, Function()->NumberOfDOFs());
    Allocate(_ScaledDirection,  Function()->NumberOfDOFs());
  }
}

// -----------------------------------------------------------------------------
InexactLineSearch::InexactLineSearch(const InexactLineSearch &other)
:
  LineSearch(other),
  _CurrentDoFValues(NULL),
  _ScaledDirection (NULL)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
InexactLineSearch &InexactLineSearch::operator =(const InexactLineSearch &other)
{
  if (this != &other) {
    LineSearch::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
InexactLineSearch::~InexactLineSearch()
{
  Deallocate(_CurrentDoFValues);
  Deallocate(_ScaledDirection);
}

// -----------------------------------------------------------------------------
void InexactLineSearch::Function(ObjectiveFunction *f)
{
  if (Function() != f) {
    Deallocate(_CurrentDoFValues);
    Deallocate(_ScaledDirection);
    LineSearch::Function(f);
    if (f) {
      Allocate(_CurrentDoFValues, f->NumberOfDOFs());
      Allocate(_ScaledDirection,  f->NumberOfDOFs());
    }
  }
}

// =============================================================================
// Parameters
// =============================================================================

// -----------------------------------------------------------------------------
bool InexactLineSearch::Set(const char *name, const char *value)
{
  // Maximum number of consectutive rejections
  if (strcmp(name, "Maximum streak of rejected steps") == 0) {
    return FromString(value, _MaxRejectedStreak);
  // Whether to start new search using step length of previous search
  }
  if (strcmp(name, "Reuse previous step length") == 0) {
    return FromString(value, _ReusePreviousStepLength);
  // Whether [min, max] step length range is strict
  }
  if (strcmp(name, "Strict step length range")             == 0 ||
      strcmp(name, "Strict incremental step length range") == 0) {
    bool limit_increments;
    if (!FromString(value, limit_increments)) return false;
    if (limit_increments) {
      _StrictStepLengthRange |= 1;
    } else {
      _StrictStepLengthRange &= ~1;
    }
    return true;
  }
  if (strcmp(name, "Strict total step length range")       == 0 ||
      strcmp(name, "Strict accumulated step length range") == 0) {
    bool limit_step;
    if (!FromString(value, limit_step)) return false;
    if (limit_step) {
      _StrictStepLengthRange |= 2;
    } else {
      _StrictStepLengthRange &= ~2;
    }
    return true;
  }
  return LineSearch::Set(name, value);
}

// -----------------------------------------------------------------------------
ParameterList InexactLineSearch::Parameter() const
{
  ParameterList params = LineSearch::Parameter();
  Insert(params, "Maximum streak of rejected steps", _MaxRejectedStreak);
  Insert(params, "Reuse previous step length",       _ReusePreviousStepLength);
  Insert(params, "Strict incremental step length range", (_StrictStepLengthRange & 1) != 0);
  Insert(params, "Strict total step length range", (_StrictStepLengthRange & 2) != 0);
  return params;
}

// =============================================================================
// Optimization
// =============================================================================

// -----------------------------------------------------------------------------
double InexactLineSearch::Advance(double alpha)
{
  if (_StepLengthUnit == .0 || alpha == .0) return .0;
  // Backup current function parameter values
  Function()->Get(_CurrentDoFValues);
  // Compute gradient for given step length
  alpha /= _StepLengthUnit;
  if (_Revert) alpha *= -1.0;
  const int ndofs = Function()->NumberOfDOFs();
  for (int dof = 0; dof < ndofs; ++dof) {
    _ScaledDirection[dof] = alpha * _Direction[dof];
  }
  // Set scaled gradient to the negative of the current parameter value if sign
  // changes are not allowed for this parameter s.t. updated value is zero
  //
  // Note: This is used for the optimization of the registration cost
  //       function with L1-norm sparsity constraint on the multi-level
  //       free-form deformation parameters (Wenzhe et al.'s Sparse FFD).
  if (_AllowSignChange) {
    double next_value;
    for (int dof = 0; dof < ndofs; ++dof) {
      if (_AllowSignChange[dof]) continue;
      next_value = _CurrentDoFValues[dof] + _ScaledDirection[dof];
      if ((_CurrentDoFValues[dof] * next_value) <= .0) {
        _ScaledDirection[dof] = - _CurrentDoFValues[dof];
      }
    }
  }
  // Update all parameters at once to only trigger a single modified event
  return Function()->Step(_ScaledDirection);
}

// -----------------------------------------------------------------------------
void InexactLineSearch::Retreat(double alpha)
{
  if (_StepLengthUnit == .0 || alpha == .0) return;
  Function()->Put(_CurrentDoFValues);
}

// -----------------------------------------------------------------------------
double InexactLineSearch::Value(double alpha, double *delta)
{
  const double max_delta = Advance(alpha);
  if (delta) *delta = max_delta;
  Function()->Update(false);
  const double value = Function()->Value();
  Retreat(alpha);
  return value;
}


} // namespace mirtk
