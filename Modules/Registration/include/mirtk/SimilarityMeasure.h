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

#ifndef MIRTK_SimilarityMeasure_H
#define MIRTK_SimilarityMeasure_H

#include "mirtk/EnergyMeasure.h"


namespace mirtk {


// -----------------------------------------------------------------------------
/// Enumeration of available similarity measures
///
/// \note This enumeration constains only a subset of all EnergyMeasure
///       enumeration values, whereby the integer value of corresponding
///       enumeration values is equal.
///
/// \sa EnergyMeasure
enum SimilarityMeasure
{
  SIM_Unknown = SIM_Begin,  ///< Unknown/invalid image (dis-)similarity measure
  SIM_JE      = EM_JE,      ///< Joint entropy
  SIM_CC      = EM_CC,      ///< Cross-correlation
  SIM_MI      = EM_MI,      ///< Mutual information
  SIM_NMI     = EM_NMI,     ///< Normalized mutual information
  SIM_SSD     = EM_SSD,     ///< Sum of squared differences
  SIM_CR_XY   = EM_CR_XY,   ///< Correlation ratio
  SIM_CR_YX   = EM_CR_YX,   ///< Correlation ratio
  SIM_LC      = EM_LC,      ///< Label consistency
  SIM_K       = EM_K,       ///< Kappa statistic
  SIM_ML      = EM_ML,      ///< Maximum likelihood
  SIM_NGF_COS = EM_NGF_COS, ///< Cosine of normalzed gradient field
  SIM_NCC     = EM_NCC,     ///< Normalized cross-correlation
  SIM_LNCC    = EM_LNCC,    ///< Local normalized cross-correlation
  SIM_CoVar   = EM_CoVar,   ///< Covariance
  SIM_PSNR    = EM_PSNR     ///< Peak signal-to-noise ratio
};

// -----------------------------------------------------------------------------
template <>
inline string ToString(const SimilarityMeasure &sim, int w, char c, bool left)
{
  EnergyMeasure em = static_cast<EnergyMeasure>(sim);
  if (em <= SIM_Begin || em >= SIM_End) return ToString("Unknown", w, c, left);
  return ToString(em, w, c, left);
}

// -----------------------------------------------------------------------------
inline string ToPrettyString(const SimilarityMeasure &sim, int w = 0, char c = ' ', bool left = true)
{
  EnergyMeasure em = static_cast<EnergyMeasure>(sim);
  if (em <= SIM_Begin || em >= SIM_End) return ToString("Unknown", w, c, left);
  return ToPrettyString(em, w, c, left);
}

// -----------------------------------------------------------------------------
template <>
inline bool FromString(const char *str, SimilarityMeasure &sim)
{
  EnergyMeasure em = EM_Unknown;
  if (FromString(str, em) && SIM_Begin < em && em < SIM_End) {
    sim = static_cast<SimilarityMeasure>(em);
    return true;
  } else {
    sim = SIM_Unknown;
    return false;
  }
}


} // namespace mirtk

#endif // MIRTK_SimilarityMeasure_H
