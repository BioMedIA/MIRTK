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

#include "mirtk/MaxStepLineSearch.h"

#include "mirtk/Event.h"


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
MaxStepLineSearch::MaxStepLineSearch(ObjectiveFunction *f)
:
  InexactLineSearch(f)
{
}

// -----------------------------------------------------------------------------
MaxStepLineSearch::MaxStepLineSearch(const MaxStepLineSearch &other)
:
  InexactLineSearch(other)
{
}

// -----------------------------------------------------------------------------
MaxStepLineSearch &MaxStepLineSearch::operator =(const MaxStepLineSearch &other)
{
  InexactLineSearch::operator =(other);
  return *this;
}

// -----------------------------------------------------------------------------
MaxStepLineSearch::~MaxStepLineSearch()
{
}

// =============================================================================
// Optimization
// =============================================================================

// -----------------------------------------------------------------------------
double MaxStepLineSearch::Run()
{
  // Get current value of objective function
  LineSearch::Run();

  // Define "aliases" for event data for prettier code below
  LineSearchStep step;
  step._Direction   = _Direction;
  step._Unit        = _StepLengthUnit;
  step._Length      = _MaxStepLength;
  step._MinLength   = _MaxStepLength;
  step._MaxLength   = _MaxStepLength;
  step._TotalLength = .0;
  step._Value       = _CurrentValue;
  step._Current     = _CurrentValue;
  step._Delta       = .0;

  double &value     = step._Value;
  double &current   = step._Current;
  double &alpha     = step._Length;
  double &total     = step._TotalLength;
  double &delta     = step._Delta;

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

        // Accumulate accepted steps
        total += alpha;
        current = value;

      } else {

        // Notify observers about progress
        Broadcast(RejectedStepEvent, &step);

        // Reject step and start over at previous point
        Retreat(alpha);

        break;
      }
    }

    // Notify observers about end of iteration
    Broadcast(LineSearchIterationEndEvent, &iteration);
  }

  // Notify observers about end of line search
  Broadcast(LineSearchEndEvent, &step);

  // Return new value of objective function
  return current;
}


} // namespace mirtk
