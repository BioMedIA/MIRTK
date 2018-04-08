/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2017 Imperial College London
 * Copyright 2017 Andreas Schuh
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

#include "mirtk/NegJacobianConstraint.h"


namespace mirtk {


// Register energy term with object factory during static initialization
mirtkAutoRegisterEnergyTermMacro(NegJacobianConstraint);


// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
NegJacobianConstraint::NegJacobianConstraint(const char *name, bool constrain_spline)
:
  JacobianConstraint(name, constrain_spline),
  _Epsilon(.01), _Gamma(.5), _Power(2)
{
  _SubDomain = SD_Shifted;
}

// -----------------------------------------------------------------------------
NegJacobianConstraint::~NegJacobianConstraint()
{
}

// =============================================================================
// Parameters
// =============================================================================

// -----------------------------------------------------------------------------
bool NegJacobianConstraint::SetWithPrefix(const char *param, const char *value)
{
  if (strcmp(param, "Jacobian determinant epsilon") == 0) {
    return FromString(value, _Epsilon);
  }
  return JacobianConstraint::SetWithPrefix(param, value);
}

// -----------------------------------------------------------------------------
bool NegJacobianConstraint::SetWithoutPrefix(const char *param, const char *value)
{
  if (strcmp(param, "Gamma") == 0 || strcmp(param, "Threshold") == 0) {
    return FromString(value, _Gamma);
  }
  if (strcmp(param, "Epsilon") == 0) {
    return FromString(value, _Epsilon);
  }
  if (strcmp(param, "Power") == 0) {
    return FromString(value, _Power);
  }
  return JacobianConstraint::SetWithoutPrefix(param, value);
}

// -----------------------------------------------------------------------------
ParameterList NegJacobianConstraint::Parameter() const
{
  ParameterList params = JacobianConstraint::Parameter();
  InsertWithPrefix(params, "Gamma",   _Gamma);
  InsertWithPrefix(params, "Epsilon", _Epsilon);
  InsertWithPrefix(params, "Power",   _Power);
  return params;
}


} // namespace mirtk
