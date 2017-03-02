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

#include "mirtk/LineSearch.h"

#include "mirtk/String.h"

#include "mirtk/MaxStepLineSearch.h"
#include "mirtk/AdaptiveLineSearch.h"
#include "mirtk/BrentLineSearch.h"

#include <limits>


namespace mirtk {


// =============================================================================
// Factory method
// =============================================================================

// -----------------------------------------------------------------------------
LineSearch *LineSearch::New(LineSearchStrategy &strategy, ObjectiveFunction *f)
{
  switch (strategy) {
    case LS_None:     return new MaxStepLineSearch (f);
    case LS_Adaptive: return new AdaptiveLineSearch(f);
    case LS_Brent:    return new BrentLineSearch(f);
    default:
      cerr << "LineSearch::New: Unknown line search strategy: " << strategy << endl;
      exit(1);
  }
  return NULL;
}

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
LineSearch::LineSearch(ObjectiveFunction *f)
:
  LocalOptimizer(f),
  _Direction         (NULL),
  _Revert            (false),
  _CurrentValue      (numeric_limits<double>::quiet_NaN()),
  _NumberOfIterations(20),
  _MinStepLength     ( 0.1),
  _MaxStepLength     (10.0),
  _StepLengthUnit    ( 1.0),
  _StepLength        ( 0.0)
{
}

// -----------------------------------------------------------------------------
void LineSearch::CopyAttributes(const LineSearch &other)
{
  _Direction          = other._Direction;
  _Revert             = other._Revert;
  _CurrentValue       = other._CurrentValue;
  _NumberOfIterations = other._NumberOfIterations;
  _MinStepLength      = other._MinStepLength;
  _MaxStepLength      = other._MaxStepLength;
  _StepLengthUnit     = other._StepLengthUnit;
  _StepLength         = other._StepLength;
}

// -----------------------------------------------------------------------------
LineSearch::LineSearch(const LineSearch &other)
:
  LocalOptimizer(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
LineSearch &LineSearch::operator =(const LineSearch &other)
{
  if (this != &other) {
    LocalOptimizer::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
LineSearch::~LineSearch()
{
}

// =============================================================================
// Parameters
// =============================================================================

// -----------------------------------------------------------------------------
bool LineSearch::Set(const char *name, const char *value)
{
  // Number of line search iterations
  if (strcmp(name, "Maximum no. of line search iterations")    == 0 ||
      strcmp(name, "Maximum number of line search iterations") == 0 ||
      strcmp(name, "Maximum no. of line iterations")           == 0 ||
      strcmp(name, "Maximum number of line iterations")        == 0 ||
      strcmp(name,    "No. of line search iterations")         == 0 ||
      strcmp(name, "Number of line search iterations")         == 0 ||
      strcmp(name,    "No. of line iterations")                == 0 ||
      strcmp(name, "Number of line iterations")                == 0) {
    return FromString(value, _NumberOfIterations) && _NumberOfIterations > 0;
  // Minimum length of step
  } else if (strcmp(name, "Minimum length of steps") == 0) {
    return FromString(value, _MinStepLength) && _MinStepLength > .0;
  // Maximum length of step
  } else if (strcmp(name, "Maximum length of steps") == 0) {
    return FromString(value, _MaxStepLength) && _MaxStepLength > .0;
  // Length of steps
  } else if (strcmp(name, "Length of steps") == 0) {
    if (FromString(value, _MaxStepLength) && _MaxStepLength > .0) {
      _MinStepLength = _MaxStepLength;
    }
  // Range of step lengths (e.g., "[0.1 1]" or "0.1 1")
  } else if (strcmp(name, "Range of steps")                  == 0 ||
             strcmp(name, "Range of step lengths")           == 0 ||
             strcmp(name, "Min-/Maximum length of steps")    == 0 ||
             strcmp(name, "Min-/maximum length of steps")    == 0 ||
             strcmp(name, "Minimum/Maximum length of steps") == 0 ||
             strcmp(name, "Minimum/maximum length of steps") == 0) {
    istringstream ss(value);
    char c1, c2;
    string str;
    ss >> c1;
    if (c1 != '[') ss.putback(c1);
    ss >> str;
    if (!FromString(str.c_str(), _MinStepLength) || _MinStepLength < .0) return false;
    ss >> str;
    if (!FromString(str.c_str(), _MaxStepLength) || _MaxStepLength < _MinStepLength) return false;
    ss >> c2;
    if (c1 == '[' && c2 != ']') return false;
    return true;
  }
  return false;
}

// -----------------------------------------------------------------------------
ParameterList LineSearch::Parameter() const
{
  ParameterList params;
  Insert(params, "Maximum no. of line iterations", ToString(_NumberOfIterations));
  if (_MinStepLength == _MaxStepLength) {
    Insert(params, "Length of steps", ToString(_MaxStepLength));
  } else {
    Insert(params, "Minimum length of steps", ToString(_MinStepLength));
    Insert(params, "Maximum length of steps", ToString(_MaxStepLength));
  }
  return params;
}

// =============================================================================
// Optimization
// =============================================================================

// -----------------------------------------------------------------------------
void LineSearch::Initialize()
{
  // Initialize base class
  LocalOptimizer::Initialize();

  // Check parameters
  if (_MinStepLength <= 0.) {
    Throw(ERR_InvalidArgument, __FUNCTION__, "Minimum length of steps must be positive, but is ", _MinStepLength);
  }
  if (_MaxStepLength < _MinStepLength) {
    Throw(ERR_InvalidArgument, __FUNCTION__, "Maximum length of steps must be greater or equal minimum length, but is ", _MaxStepLength);
  }

  // Set initial step length
  _StepLength = _MaxStepLength;
}


} // namespace mirtk
