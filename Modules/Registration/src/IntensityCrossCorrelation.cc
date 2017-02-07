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

#include "mirtk/IntensityCrossCorrelation.h"

#include "mirtk/ObjectFactory.h"


namespace mirtk {


// Register energy term with object factory during static initialization
mirtkAutoRegisterEnergyTermMacro(IntensityCrossCorrelation);


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
IntensityCrossCorrelation
::IntensityCrossCorrelation(const char *name)
:
  NormalizedIntensityCrossCorrelation(name)
{
}

// -----------------------------------------------------------------------------
IntensityCrossCorrelation
::IntensityCrossCorrelation(const IntensityCrossCorrelation &other)
:
  NormalizedIntensityCrossCorrelation(other)
{
}

// -----------------------------------------------------------------------------
IntensityCrossCorrelation::~IntensityCrossCorrelation()
{
}


} // namespace mirtk
