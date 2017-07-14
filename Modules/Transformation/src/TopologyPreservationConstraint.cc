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

#include "mirtk/TopologyPreservationConstraint.h"


namespace mirtk {


// Register energy term with object factory during static initialization
mirtkAutoRegisterEnergyTermMacro(TopologyPreservationConstraint);


// -----------------------------------------------------------------------------
TopologyPreservationConstraint::TopologyPreservationConstraint(const char *name)
:
  JacobianConstraint(name),
  _Gamma(.3)
{
  _ConstrainPassiveDoFs = true;
}


// =============================================================================
// Parameters
// =============================================================================

// -----------------------------------------------------------------------------
bool TopologyPreservationConstraint::SetWithoutPrefix(const char *param, const char *value)
{
  if (strcmp(param, "Threshold") == 0 ||
      strcmp(param, "Gamma")     == 0) {
    return FromString(value, _Gamma);
  }
  return JacobianConstraint::SetWithoutPrefix(param, value);
}

// -----------------------------------------------------------------------------
ParameterList TopologyPreservationConstraint::Parameter() const
{
  ParameterList params = JacobianConstraint::Parameter();
  InsertWithPrefix(params, "Threshold", _Gamma);
  return params;
}


} // namespace mirtk
