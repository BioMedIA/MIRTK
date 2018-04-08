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

#include "mirtk/BSpline.h"

namespace mirtk {


#ifndef __clang__

// Lookup table of B-spline function values
template <class TReal> MIRTK_Numerics_EXPORT TReal BSpline<TReal>::WeightLookupTable[BSpline<TReal>::LookupTableSize];

// Lookup table of B-spline basis function values
template <class TReal> MIRTK_Numerics_EXPORT TReal BSpline<TReal>::LookupTable   [BSpline<TReal>::LookupTableSize][4];
template <class TReal> MIRTK_Numerics_EXPORT TReal BSpline<TReal>::LookupTable_I [BSpline<TReal>::LookupTableSize][4];
template <class TReal> MIRTK_Numerics_EXPORT TReal BSpline<TReal>::LookupTable_II[BSpline<TReal>::LookupTableSize][4];

// Wether lookup tables of B-spline kernel were initialized
template <class TReal> MIRTK_Numerics_EXPORT bool BSpline<TReal>::_initialized = false;

#endif // !defined(__clang__)


// Explicit template instantiations
template class BSpline<float>;
template class BSpline<double>;


} // namespace mirtk
