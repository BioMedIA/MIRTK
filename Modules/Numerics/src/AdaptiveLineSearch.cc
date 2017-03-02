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

#include "mirtk/AdaptiveLineSearch.h"

#include <algorithm>


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
AdaptiveLineSearch::AdaptiveLineSearch(ObjectiveFunction *f)
:
  InexactLineSearch(f),
  _StepLengthRise(1.1),
  _StepLengthDrop(0.5)
{
}

// -----------------------------------------------------------------------------
void AdaptiveLineSearch::CopyAttributes(const AdaptiveLineSearch &other)
{
  _StepLengthRise = other._StepLengthRise;
  _StepLengthDrop = other._StepLengthDrop;
}

// -----------------------------------------------------------------------------
AdaptiveLineSearch::AdaptiveLineSearch(const AdaptiveLineSearch &other)
:
  InexactLineSearch(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
AdaptiveLineSearch &AdaptiveLineSearch::operator =(const AdaptiveLineSearch &other)
{
  if (this != &other) {
    InexactLineSearch::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
AdaptiveLineSearch::~AdaptiveLineSearch()
{
}

// =============================================================================
// Parameters
// =============================================================================

// -----------------------------------------------------------------------------
bool AdaptiveLineSearch::Set(const char *name, const char *value)
{
  if (strcmp(name, "Step length rise") == 0) {
    return FromString(value, _StepLengthRise) && _StepLengthRise > .0;
  }
  if (strcmp(name, "Step length drop") == 0) {
    return FromString(value, _StepLengthDrop) && _StepLengthDrop > .0;
  }
  return InexactLineSearch::Set(name, value);
}

// -----------------------------------------------------------------------------
ParameterList AdaptiveLineSearch::Parameter() const
{
  ParameterList params = InexactLineSearch::Parameter();
  Insert(params, "Step length rise", _StepLengthRise);
  Insert(params, "Step length drop", _StepLengthDrop);
  return params;
}

// =============================================================================
// Optimization
// =============================================================================

// -----------------------------------------------------------------------------
void AdaptiveLineSearch::Initialize()
{
  // Initialize base class
  InexactLineSearch::Initialize();

  // Set initial step length
  if (_StrictStepLengthRange >= 2) {
    _StepLength = min(_MinStepLength + .5 * (_MaxStepLength - _MinStepLength), .5 * _MaxStepLength);
  } else {
    _StepLength = _MaxStepLength;
  }
}

// -----------------------------------------------------------------------------
double AdaptiveLineSearch::Run()
{
  const double min_length   = _MinStepLength;
  const double max_length   = _MaxStepLength;
  const bool   strict_range = (_StrictStepLengthRange != 0 || fequal(_MinStepLength, _MaxStepLength));
  const bool   strict_total = (_StrictStepLengthRange >= 2);
  const double rise         = (strict_total ? .5 : _StepLengthRise);
  const double drop         = (strict_total ? .5 : _StepLengthDrop);
  // Note: Binary search-like steps in case of strict total step length range

  // Number of consecutively rejected steps after at least one accepted step
  int rejected_streak = -1;

  // Get current value of objective function
  LineSearch::Run();

  // Define "aliases" for event data for prettier code below
  LineSearchStep step;
  step._Info        = (strict_total ? nullptr : (strict_range ? "incremental" : "unbound"));
  step._Direction   = _Direction;
  step._Unit        = _StepLengthUnit;
  step._MinLength   = min_length;
  step._MaxLength   = max_length;
  step._TotalLength = .0;
  step._TotalDelta  = .0;
  double &value     = step._Value;
  double &current   = step._Current;
  double &alpha     = step._Length;
  double &total     = step._TotalLength;
  double &delta     = step._Delta;
  double &Delta     = step._TotalDelta;

  // Start line search
  total   = .0;  // accumulated length of accepted steps
  value   = _CurrentValue;
  current = _CurrentValue;

  if (strict_total) {
    alpha = min(min_length + .5 * (max_length - min_length), .5 * max_length);
  } else {
    alpha = max_length;
  }

  // Continue search using total step length of previous search
  if (_ReusePreviousStepLength && _StepLength > .0) alpha = _StepLength;

  // Limit incremental step length strictly to [min_length, max_length]
  if (strict_range) {
    if      (alpha > max_length) alpha = max_length;
    else if (alpha < min_length) alpha = min_length;
  }

  // Notify observers about start of line search
  Broadcast(LineSearchStartEvent, &step);

  // Walk along search direction until no further improvement
  Iteration iteration(0, _NumberOfIterations);
  while (iteration.Next()) {

    // Notify observers about start of iteration
    Broadcast(LineSearchIterationStartEvent, &iteration);

    // Take step along search direction
    delta = Advance(alpha);

    // Check minimum maximum change of DoFs convergence criterium
    if (delta <= _Delta) {

      Retreat(alpha);
      break;

    } else {

      // Re-evaluate objective function
      Function()->Update(false);
      value = Function()->Value();

      // Check minimum change of objective function value convergence criterium
      if (IsImprovement(current, value)) {

        // Notify observers about progress
        Broadcast(AcceptedStepEvent, &step);

        // New current value
        current = value;

        // Accumulate accepted steps
        total += alpha;
        Delta += delta;

        // Increase step length and continue search from new point
        // If the step length range is not strict, keep possibly
        // greater previous step length as it was accepted again
        if (strict_range || alpha < max_length) {
          alpha = min(alpha * rise, max_length);
          if (strict_total && total + alpha > max_length) {
            alpha = max_length - total;
            if (alpha < 1e-12) break;
          }
        }

        // Reset counter of consecutive rejections
        rejected_streak = 0;

      } else {

        // Notify observers about progress
        Broadcast(RejectedStepEvent, &step);

        // Reject step and start over at previous point
        Retreat(alpha);

        // Stop search if maximum number of consecutive rejections exceeded
        if (rejected_streak != -1) {
          ++rejected_streak;
          if (_MaxRejectedStreak >= 0 && rejected_streak > _MaxRejectedStreak) break;
        }

        // Decrease step length and continue if possible
        if (alpha == min_length) break;
        // If previous step length was reused even though it was outside
        // the considered [min, max] step length range, start over with
        // step length limited by the original range as it would have been
        // if previous length had not been reused.
        if (alpha > max_length) alpha = max_length;
        else alpha = max(alpha * drop, min_length);
      }
    }

    // Notify observers about end of iteration
    Broadcast(LineSearchIterationEndEvent, &iteration);

  }

  // Notify observers about end of line search
  Broadcast(LineSearchEndEvent, &step);

  // Re-use final step length upon next line search
  _StepLength = total;

  // Return new value of objective function
  return current;
}


} // namespace mirtk
