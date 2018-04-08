/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2016-2017 Imperial College London
 * Copyright 2016-2017 Andreas Schuh
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

#include "mirtk/LabelConsistency.h"

#include "mirtk/ObjectFactory.h"


namespace mirtk {


// Register energy term with object factory during static initialization
mirtkAutoRegisterEnergyTermMacro(LabelConsistency);


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
LabelConsistency::LabelConsistency(const char *name)
:
  HistogramImageSimilarity(name)
{
}

// -----------------------------------------------------------------------------
LabelConsistency::LabelConsistency(const LabelConsistency &other)
:
  HistogramImageSimilarity(other)
{
}

// -----------------------------------------------------------------------------
LabelConsistency::~LabelConsistency()
{
}

// =============================================================================
// Initialization
// =============================================================================

// -----------------------------------------------------------------------------
void LabelConsistency::Initialize()
{
  // Number of target and source bins must be equal
  if (_Samples && !_SamplesOwner) {
    if (_Samples->NumberOfBinsX() != _Samples->NumberOfBinsY()) {
      cerr << this->NameOfType() << "::Initialize: External joint histogram must"
          " have equal number of bins in both dimensions" << endl;
      exit(1);
    }
  } else {
    if (_NumberOfTargetBins <= 0) {
      double tmin, tmax;
      Target()->InputImage()->GetMinMaxAsDouble(&tmin, &tmax);
      _NumberOfTargetBins = DefaultNumberOfBins(Target()->InputImage(), tmin, tmax);
    }
    if (_NumberOfSourceBins <= 0) {
      double smin, smax;
      Source()->InputImage()->GetMinMaxAsDouble(&smin, &smax);
      _NumberOfSourceBins = DefaultNumberOfBins(Source()->InputImage(), smin, smax);
    }
    _NumberOfTargetBins = _NumberOfSourceBins = (_NumberOfTargetBins + _NumberOfSourceBins) / 2;
  }

  // Initialize base class
  HistogramImageSimilarity::Initialize();
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
double LabelConsistency::Evaluate()
{
  return _Histogram.LabelConsistency();
}


} // namespace mirtk
