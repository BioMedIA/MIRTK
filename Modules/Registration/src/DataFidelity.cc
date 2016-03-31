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

#include "mirtk/DataFidelity.h"


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
DataFidelity::DataFidelity(const char *name, double weight)
:
  EnergyTerm(name, weight)
{
  _DivideByInitialValue = true;
}

// -----------------------------------------------------------------------------
DataFidelity::DataFidelity(const DataFidelity &other)
:
  EnergyTerm(other)
{
}

// -----------------------------------------------------------------------------
DataFidelity &DataFidelity::operator =(const DataFidelity &other)
{
  EnergyTerm::operator =(other);
  return *this;
}

// -----------------------------------------------------------------------------
DataFidelity::~DataFidelity()
{
}

// =============================================================================
// Parameters
// =============================================================================

// -----------------------------------------------------------------------------
bool DataFidelity::SetWithPrefix(const char *param, const char *value)
{
  if (strcmp(param, "Divide data term by initial value")           == 0 ||
      strcmp(param, "Divide data terms by initial value")          == 0 ||
      strcmp(param, "Divide data fidelity term by initial value")  == 0 ||
      strcmp(param, "Divide data fidelity terms by initial value") == 0) {
    return FromString(value, _DivideByInitialValue);
  }
  return EnergyTerm::SetWithPrefix(param, value);
}


} // namespace mirtk
