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

#include "mirtk/BrentLineSearch.h"

#include "mirtk/Math.h"
#include "mirtk/Memory.h"
#include "mirtk/String.h"
#include "mirtk/Profiling.h"

#include <algorithm>


namespace mirtk {


// =============================================================================
// Constants
// =============================================================================

const double BrentLineSearch::_MinificationRatio     = 0.3819660;
const double BrentLineSearch::_MagnificationRatio    = 1.6180340;
const double BrentLineSearch::_MaxMagnificationRatio = 100.0;

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
BrentLineSearch::BrentLineSearch(ObjectiveFunction *f)
:
  InexactLineSearch(f),
  _Tolerance(.01)
{
}

// -----------------------------------------------------------------------------
void BrentLineSearch::CopyAttributes(const BrentLineSearch &other)
{
  _Tolerance = other._Tolerance;
}

// -----------------------------------------------------------------------------
BrentLineSearch::BrentLineSearch(const BrentLineSearch &other)
:
  InexactLineSearch(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
BrentLineSearch &BrentLineSearch::operator =(const BrentLineSearch &other)
{
  if (this != &other) {
    InexactLineSearch::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
BrentLineSearch::~BrentLineSearch()
{
}

// =============================================================================
// Parameters
// =============================================================================

// -----------------------------------------------------------------------------
bool BrentLineSearch::Set(const char *name, const char *value)
{
  if (strcmp(name, "Brent's line search tolerance") == 0) {
    return FromString(value, _Tolerance) && _Tolerance >= .0;
  }
  return InexactLineSearch::Set(name, value);
}

// -----------------------------------------------------------------------------
ParameterList BrentLineSearch::Parameter() const
{
  ParameterList params = InexactLineSearch::Parameter();
  Insert(params, "Brent's line search tolerance", _Tolerance);
  return params;
}

// =============================================================================
// Optimization
// =============================================================================

// -----------------------------------------------------------------------------
double BrentLineSearch::Run()
{
  double m, v, w, fv, fw, r, q, p, d, tol1, tol2, etemp, e = .0, x = .0;
  Iteration iteration(0, _NumberOfIterations);

  // Number of consecutively rejected steps after at least one accepted step
  int rejected_streak = -1;

  // Get current value of objective function
  LineSearch::Run();

  // Notify observers about start of line search
  LineSearchStep step;
  step._Info        = "Brent's method,";
  step._Direction   = _Direction;
  step._Unit        = _StepLengthUnit;
  step._MinLength   = _MinStepLength;
  step._MaxLength   = _MaxStepLength;
  step._Length      = _StepLength;
  step._Value       = _CurrentValue;
  step._Current     = _CurrentValue;
  step._TotalLength = .0;

  double &a  = step._MinLength; // current minimum step length
  double &b  = step._MaxLength; // current maximum step length
  double &u  = step._Length;    // trial step length
  double &fx = step._Current;   // current best value
  double &fu = step._Value;     // function value for trial step length

  // Bracket extremum for step length x in interval [a, b]
  if (!_StrictStepLengthRange) {
    fu = BracketExtremum(a, u, b, &step._Delta);
    step._Info = "Brent's method and bracketed";
  }

  // Notify observers about start of line search
  Broadcast(LineSearchStartEvent, &step);

  // Allocate memory for parameters corresponding to last accepted step
  double *params = Allocate<double>(Function()->NumberOfDOFs());
  Function()->Get(params);

  // Initial iteration
  if (_ReusePreviousStepLength || !_StrictStepLengthRange) {
    if (iteration.Next()) {

      // Notify observers about start of iteration
      Broadcast(LineSearchIterationStartEvent, &iteration);

      // Try previous step length if not done so already during bracketing
      if (_StrictStepLengthRange && _ReusePreviousStepLength) {
        u = _StepLength;
        // Use Advance/Retreat instead of Value such that the objective
        // function need not be updated to report the individual term values
        // upon the following AcceptedStepEvent or RejectedStepEvent
        step._Delta = Advance(u);
        Function()->Update(false);
        fu = Function()->Value();
      }

      // Notify observers about progress
      if (IsImprovement(fx, fu)) {
        Broadcast(AcceptedStepEvent, &step);
        Function()->Get(params);
        step._TotalLength = x = u, fx = fu;
        rejected_streak = 0;
      } else {
        Broadcast(RejectedStepEvent, &step);
      }

      // Reset objective function parameters
      Retreat(u);

      // Notify observers about end of iteration
      Broadcast(LineSearchIterationEndEvent, &iteration);
    }
  }

  // Iterations of Brent's method
  u  = v  = w  = x;
  fu = fv = fw = fx;
  while (iteration.Next()) {

    // Current tolerance values
    tol1 = _Tolerance * abs(x) + 1e-10;
    tol2 = 2.0 * tol1;

    // Get midpoint of current step length interval
    m = .5 * (a + b);

    // Test if we are done
    if (x >= a && (x - m) <= (tol2 - .5 * (b - a))) {
      break;
    }

    // Notify observers about start of iteration
    Broadcast(LineSearchIterationStartEvent, &iteration);

    // Construct a trial parabolic fit
    if (abs(e) > tol1) {
      r = (x - w) * (fx - fv);
      q = (x - v) * (fx - fw);
      p = (x - v) * q - (x - w) * r;
      q = 2.0 * (q - r);
      if (q > .0) p = -p;
      q = abs(q);
      etemp = e;
      e = d;
      if (abs(p) >= abs(.5 * q * etemp) || p <= q * (a - x) || p >= q * (b - x)) {
        e = (x >= m ? a : b) - x;
        d = _MinificationRatio * e;
      } else {
        d = p / q;
        u = x + d;
        if (u - a < tol2 || b - u < tol2) {
          d = copysign(tol1, m - x);
        }
      }
    } else {
      e = (x >= m ? a : b) - x;
      d = _MinificationRatio * e;
    }

    if (abs(d) < tol1) d = copysign(tol1, d);
    u = x + d;

    // Use Advance/Retreat instead of Value such that the objective
    // function need not be updated to report the individual term values
    // upon the following AcceptedStepEvent or RejectedStepEvent
    step._Delta = Advance(u);
    Function()->Update(false);
    fu = Function()->Value();

		if (IsImprovement(fx, fu)) {

      // Notify observers about progress
      Broadcast(AcceptedStepEvent, &step);

      // Store function parameters
      Function()->Get(params);

      // Update bracketing triplet and current optimal step
      if (u >= x) a = x;
      else        b = x;
      v = w, fv = fw;
      w = x, fw = fx;
      x = u, fx = fu;

      step._TotalLength = x;

      // Reset counter of consecutive rejections
      rejected_streak = 0;

    } else {

      // Notify observers about progress
      Broadcast(RejectedStepEvent, &step);

      // Stop search if maximum number of consecutive rejections exceeded
      if (rejected_streak != -1) {
        ++rejected_streak;
        if (_MaxRejectedStreak >= 0 && rejected_streak > _MaxRejectedStreak) break;
      }

      // Update bracketing triplet
      if (u < x) a = u;
      else       b = u;

      if (IsImprovement(fw, fu) || w == x) {
        v = w, fv = fw;
        w = u, fw = fu;
      } else if (IsImprovement(fv, fu) || v == x || v == w) {
        v = u, fv = fu;
      }

    }

    // Reset objective function parameters
    Retreat(u);

    // Notify observers about end of iteration
    Broadcast(LineSearchIterationEndEvent, &iteration);

  }

  // Restore parameters corresponding to optimal step
  Function()->Put(params);
  Deallocate(params);

  // Notify observers about end of line search
  Broadcast(LineSearchEndEvent, &step);

  // Re-use optimal step length upon next line search
  _StepLength = x;

  // Return new value of objective function
  return fx;
}

// -----------------------------------------------------------------------------
double BrentLineSearch::BracketExtremum(double &a, double &b, double &c, double *delta)
{
  double r, q, d, u, ulim, fu, da, db, dc, du;

  MIRTK_START_TIMING();

  // Initial guess
  a = _MinStepLength;
  if (_ReusePreviousStepLength) {
    b = _StepLength;
    c = b + _MagnificationRatio * (b - a);
  } else {
    c = _MaxStepLength;
    b = (c - _MagnificationRatio * a) / (1.0 + _MagnificationRatio);
    b = max(a, b);
  }

  double fa = Value(a, &da);
  double fb = Value(b, &db);
  double fc = Value(c, &dc);

  // Iterate until an extremum at b is bracketed in [a, c]
  while (IsImprovement(fb, fc)) {
    // Parabolic extrapolation
    r = (b - a) * (fb - fc);
    q = (b - c) * (fb - fa);
    d = q - r;
    if (abs(d) < 1e-20) d = copysign(1e-20, d);
    u = b - ((b - c) * q - (b - a) * r) / (2.0 * d);
    // When b < u < c, try it
    if ((b - u) * (u - c) > .0) {
      fu = Value(u, &du);
      // Got a minimum between b and c
      if (IsImprovement(fc, fu)) {
        a = b, fa = fb, da = db;
        b = u, fb = fu, db = du;
        break;
      // Got a minimum between a and u
      } else if (!IsImprovement(fb, fu)) {
        c = u, fc = fu, dc = du;
        break;
      }
      // Parabolic fit was of no use, use default magnification
      u  = c + _MagnificationRatio * (c - b);
      fu = Value(u, &du);
    } else {
      ulim = b + _MaxMagnificationRatio * (c - b);
      // When c < u < ulim
      if ((c - u) * (u - ulim) > .0) {
        fu = Value(u, &du);
        if (IsImprovement(fc, fu)) {
          b = c, fb = fc, db = dc;
          c = u, fc = fu, dc = du;
          u = c + _MagnificationRatio * (c - b);
          fu = Value(u, &du);
        }
      // Limit u to maximum allowed value
      } else if ((u - ulim) * (ulim - c) >= .0) {
        u = ulim, fu = Value(u, &du);
      // Reject u, use default magnification
      } else {
        u  = c + _MagnificationRatio * (c - b);
        fu = Value(u, &du);
      }
    }
    // Eliminate oldest point and continue
    a = b, fa = fb, da = db;
    b = c, fb = fc, db = dc;
    c = u, fc = fu, dc = du;
  }

  MIRTK_DEBUG_TIMING(4, "bracketing of extremum");

  if (delta) *delta = db;
  return fb;
}


} // namespace mirtk
