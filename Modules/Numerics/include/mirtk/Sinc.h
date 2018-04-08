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

#ifndef MIRTK_Sinc_H
#define MIRTK_Sinc_H

#include "mirtk/Config.h"
#include "mirtk/Math.h"

#include "mirtk/NumericsExport.h"


namespace mirtk {


/**
 * Sinc function
 */
template <class TReal>
class Sinc
{
public:

  /// Floating point precision type
  typedef TReal Real;
 
  /// Radius of Sinc kernel
  /// (for use as Kernel::Radius, where Kernel is a typedef of Sinc<T>)
  MIRTKCU_API static const int Radius = 6;

  /// Size of Sinc kernel
  MIRTKCU_API static const int KernelSize = 2 * Radius + 1;

  /// Size of Sinc function lookup table
  /// Note that the actual array size is: LookupTableSize * (Radius + 0.5)
  MIRTKCU_API static const int LookupTableSize = 1000000;

  /// Lookup table of Sinc function values
  MIRTK_Numerics_EXPORT MIRTKCU_API static Real *LookupTable;

  /// Initialize lookup table of Sinc function values
  MIRTKCU_API static void Initialize();

  /// Lookup Sinc function value
  MIRTKCU_API static Real Lookup(TReal);

};

// -----------------------------------------------------------------------------
// Fix Clang warning regarding missing definition of static template member
#ifdef __clang__
template <class TReal> TReal *Sinc<TReal>::LookupTable = nullptr;
#endif // defined(__clang__)

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
template <class TReal>
inline TReal Sinc<TReal>::Lookup(TReal x)
{
  return LookupTable[iround(abs(x) * LookupTableSize)];
}


} // namespace mirtk

#endif // MIRTK_Sinc_H
