/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2016 Imperial College London
 * Copyright 2016 Andreas Schuh
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

#include "mirtk/PeakSignalToNoiseRatio.h"

#include "mirtk/Math.h"
#include "mirtk/ObjectFactory.h"


namespace mirtk {


// Register energy term with object factory during static initialization
mirtkAutoRegisterEnergyTermMacro(PeakSignalToNoiseRatio);

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
PeakSignalToNoiseRatio::PeakSignalToNoiseRatio(const char *name)
:
  SumOfSquaredIntensityDifferences(name)
{
}

// -----------------------------------------------------------------------------
PeakSignalToNoiseRatio::PeakSignalToNoiseRatio(const PeakSignalToNoiseRatio &other)
:
  SumOfSquaredIntensityDifferences(other)
{
}

// -----------------------------------------------------------------------------
PeakSignalToNoiseRatio &PeakSignalToNoiseRatio::operator =(const PeakSignalToNoiseRatio &other)
{
  SumOfSquaredIntensityDifferences::operator =(other);
  return *this;
}

// -----------------------------------------------------------------------------
PeakSignalToNoiseRatio::~PeakSignalToNoiseRatio()
{
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
double PeakSignalToNoiseRatio::Evaluate()
{
  return - (20.0 * log10(_MaxTargetIntensity) -
            10.0 * log10(_SumSqDiff / _NumberOfForegroundVoxels));
}

// -----------------------------------------------------------------------------
double PeakSignalToNoiseRatio::RawValue(double value) const
{
  return - ImageSimilarity::RawValue(value);
}


} // namespace mirtk
