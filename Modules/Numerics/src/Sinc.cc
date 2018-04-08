/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2015 Imperial College London
 * Copyright 2013-2018 Andreas Schuh
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

#include "mirtk/Sinc.h"
#include "mirtk/Math.h"


namespace mirtk {


// -----------------------------------------------------------------------------
#ifndef __clang__
template <class TReal>
MIRTK_Numerics_EXPORT TReal *Sinc<TReal>::LookupTable = nullptr;
#endif // !defined(__clang__)

// -----------------------------------------------------------------------------
template <class TReal>
void Sinc<TReal>::Initialize()
{
  if (!LookupTable) {
    // Allocate lookup table
    const int N = iround(LookupTableSize * (Radius + 0.5));
    LookupTable = new TReal[N];
    // Value at zero distance
    LookupTable[0] = TReal(1);
    // Fill remaining fields (symmetric, so we only use positive distances)
    TReal alpha;
    for (int i = 1; i < N; ++i) {
      alpha          = TReal((pi * i) / LookupTableSize);
      LookupTable[i] = TReal(0.5) * (TReal(1) + cos(alpha/Radius)) * sin(alpha) / alpha;
    }
  }
}

// -----------------------------------------------------------------------------
template class Sinc<float>;
template class Sinc<double>;


} // namespace mirtk
