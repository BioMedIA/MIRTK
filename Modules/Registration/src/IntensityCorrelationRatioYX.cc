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

#include "mirtk/IntensityCorrelationRatioYX.h"

#include "mirtk/ObjectFactory.h"


namespace mirtk {


// Register energy term with object factory during static initialization
mirtkAutoRegisterEnergyTermMacro(IntensityCorrelationRatioYX);


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
IntensityCorrelationRatioYX
::IntensityCorrelationRatioYX(const char *name)
:
  HistogramImageSimilarity(name)
{
}

// -----------------------------------------------------------------------------
IntensityCorrelationRatioYX
::IntensityCorrelationRatioYX(const IntensityCorrelationRatioYX &other)
:
  HistogramImageSimilarity(other)
{
}

// -----------------------------------------------------------------------------
IntensityCorrelationRatioYX::~IntensityCorrelationRatioYX()
{
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
double IntensityCorrelationRatioYX::Evaluate()
{
  return _Histogram.CorrelationRatioYX();
}


} // namespace mirtk
