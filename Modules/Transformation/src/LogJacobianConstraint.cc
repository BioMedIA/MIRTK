/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2017 Imperial College London
 * Copyright 2013-2017 Andreas Schuh
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

#include "mirtk/LogJacobianConstraint.h"


namespace mirtk {


// Register energy term with object factory during static initialization
mirtkAutoRegisterEnergyTermMacro(LogJacobianConstraint);


// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
LogJacobianConstraint::LogJacobianConstraint(const char *name, bool constrain_spline)
:
  JacobianConstraint(name, constrain_spline),
  _Epsilon(.005)
{
  _ParameterPrefix.push_back("Jacobian penalty ");
}

// -----------------------------------------------------------------------------
LogJacobianConstraint::~LogJacobianConstraint()
{
}

// =============================================================================
// Parameters
// =============================================================================

// -----------------------------------------------------------------------------
bool LogJacobianConstraint::SetWithPrefix(const char *param, const char *value)
{
  if (strcmp(param, "Jacobian determinant epsilon") == 0) {
    return FromString(value, _Epsilon);
  }
  return JacobianConstraint::SetWithPrefix(param, value);
}

// -----------------------------------------------------------------------------
bool LogJacobianConstraint::SetWithoutPrefix(const char *param, const char *value)
{
  if (strcmp(param, "Epsilon") == 0) {
    return FromString(value, _Epsilon);
  }
  return JacobianConstraint::SetWithoutPrefix(param, value);
}

// -----------------------------------------------------------------------------
ParameterList LogJacobianConstraint::Parameter() const
{
  ParameterList params = JacobianConstraint::Parameter();
  InsertWithPrefix(params, "Epsilon", _Epsilon);
  return params;
}


} // namespace mirtk
