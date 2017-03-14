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

#include "mirtk/EnergyThreshold.h"

#include "mirtk/Math.h"


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
EnergyThreshold::EnergyThreshold(const ObjectiveFunction *f)
:
  StoppingCriterion(f),
  _Threshold(.0)
{
}

// -----------------------------------------------------------------------------
void EnergyThreshold::CopyAttributes(const EnergyThreshold &other)
{
  _Threshold = other._Threshold;
}

// -----------------------------------------------------------------------------
EnergyThreshold::EnergyThreshold(const EnergyThreshold &other)
:
  StoppingCriterion(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
EnergyThreshold &EnergyThreshold::operator =(const EnergyThreshold &other)
{
  if (this != &other) {
    StoppingCriterion::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
StoppingCriterion *EnergyThreshold::New() const
{
  return new EnergyThreshold(*this);
}

// -----------------------------------------------------------------------------
EnergyThreshold::~EnergyThreshold()
{
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
bool EnergyThreshold::Fulfilled(int, double value, const double *)
{
  return value <= _Threshold;
}

// =============================================================================
// Logging
// =============================================================================

// -----------------------------------------------------------------------------
void EnergyThreshold::Print(ostream &os) const
{
  const ios::fmtflags fmt = os.flags();
  os << "target value = ";
  if (_Threshold != .0 && (abs(_Threshold) < 1.e-5 || abs(_Threshold) >= 1.e5)) {
    os << scientific << setprecision(5) << setw(8) << _Threshold; // e-0x
  } else os << fixed << setprecision(5) << setw(8) << _Threshold;
  os.flags(fmt);
}


} // namespace mirtk
