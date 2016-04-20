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

#include "mirtk/StoppingCriterion.h"
#include "mirtk/LocalOptimizer.h"


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
StoppingCriterion::StoppingCriterion(const ObjectiveFunction *f)
:
  _Optimizer(NULL),
  _Function(f)
{
}

// -----------------------------------------------------------------------------
void StoppingCriterion::CopyAttributes(const StoppingCriterion &other)
{
  _Optimizer = other._Optimizer;
  _Function  = other._Function;
}

// -----------------------------------------------------------------------------
StoppingCriterion::StoppingCriterion(const LocalOptimizer *optimizer)
:
  _Optimizer(optimizer),
  _Function(optimizer->Function())
{
}

// -----------------------------------------------------------------------------
StoppingCriterion::StoppingCriterion(const StoppingCriterion &other)
:
  Object(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
StoppingCriterion &StoppingCriterion::operator =(const StoppingCriterion &other)
{
  if (this != &other) {
    Object::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
StoppingCriterion::~StoppingCriterion()
{
}

// -----------------------------------------------------------------------------
void StoppingCriterion::Initialize()
{
  if (_Optimizer) {
    if (_Function == NULL) _Function = _Optimizer->Function();
  }
  if (!_Function) {
    cerr << "StoppingCriterion::Initialize: Objective function not set" << endl;
    exit(1);
  }
  if (_Optimizer && _Function != _Optimizer->Function()) {
    cerr << "StoppingCriterion::Initialize: Objective function of optimizer and stopping criterion differ" << endl;
    exit(1);
  }
}

// =============================================================================
// Logging
// =============================================================================

// -----------------------------------------------------------------------------
void StoppingCriterion::Print(ostream &) const
{
}


} // namespace mirtk
