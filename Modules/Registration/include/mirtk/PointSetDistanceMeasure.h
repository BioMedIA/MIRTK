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

#ifndef MIRTK_PointSetDistanceMeasure_H
#define MIRTK_PointSetDistanceMeasure_H

#include "mirtk/EnergyMeasure.h"


namespace mirtk {


// -----------------------------------------------------------------------------
/// Enumeration of available point set distance measures
///
/// \note This enumeration constains only a subset of all EnergyMeasure
///       enumeration values, whereby the integer value of corresponding
///       enumeration values is equal.
///
/// \sa EnergyMeasure
enum PointSetDistanceMeasure
{
  PDM_Unknown                = PDM_Begin,
  PDM_FRE                    = EM_FRE,
  PDM_CorrespondenceDistance = EM_CorrespondenceDistance,
  PDM_CurrentsDistance       = EM_CurrentsDistance,
  PDM_VarifoldDistance       = EM_VarifoldDistance,
};

// -----------------------------------------------------------------------------
template <>
inline string ToString(const PointSetDistanceMeasure &pdm, int w, char c, bool left)
{
  EnergyMeasure em = static_cast<EnergyMeasure>(pdm);
  if (em <= PDM_Begin || em >= PDM_End) return ToString("Unknown", w, c, left);
  return ToString(em, w, c, left);
}

// -----------------------------------------------------------------------------
template <>
inline bool FromString(const char *str, PointSetDistanceMeasure &pdm)
{
  EnergyMeasure em = EM_Unknown;
  if (FromString(str, em) && PDM_Begin < em && em < PDM_End) {
    pdm = static_cast<PointSetDistanceMeasure>(em);
    return true;
  } else {
    pdm = PDM_Unknown;
    return false;
  }
}


} // namespace mirtk

#endif // MIRTK_PointSetDistanceMeasure_H
