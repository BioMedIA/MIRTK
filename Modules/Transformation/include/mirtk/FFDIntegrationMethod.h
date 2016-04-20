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

#ifndef MIRTK_FreeFormTransformationIntegrationMethod_H
#define MIRTK_FreeFormTransformationIntegrationMethod_H

#include "mirtk/String.h"


namespace mirtk {


// ---------------------------------------------------------------------------
/// Enumeration of implemented numerical integration methods
enum FFDIntegrationMethod
{
  FFDIM_Unknown = 0,                                                   // keep at begin
  FFDIM_RKE1,   FFDIM_RKE2,   FFDIM_RKH2,  FFDIM_RK4,                  // explicit RK
  FFDIM_RKEH12, FFDIM_RKBS23, FFDIM_RKF45, FFDIM_RKDP45, FFDIM_RKCK45, // embedded RK
  FFDIM_SS, FFDIM_FastSS,                                              // scaling and squaring
  FFDIM_Last                                                           // keep at end
};

typedef FFDIntegrationMethod FFDIM;

// ---------------------------------------------------------------------------
/// Convert FFD integration method enumeration value to string
template <>
inline string ToString(const FFDIM &m, int w, char c, bool left)
{
  const char *str;
  switch (m) {
    case FFDIM_SS:     str = "SS";      break;
    case FFDIM_FastSS: str = "FastSS";  break;
    case FFDIM_RKE1:   str = "RKE1";    break;
    case FFDIM_RKE2:   str = "RKE2";    break;
    case FFDIM_RKH2:   str = "RKH2";    break;
    case FFDIM_RK4:    str = "RK4";     break;
    case FFDIM_RKEH12: str = "RKEH12";  break;
    case FFDIM_RKBS23: str = "RKBS23";  break;
    case FFDIM_RKF45:  str = "RKF45";   break;
    case FFDIM_RKDP45: str = "RKDP45";  break;
    case FFDIM_RKCK45: str = "RKCK45";  break;
    default:           str = "Unknown"; break;
  }
  return ToString(str, w, c, left);
}

// ---------------------------------------------------------------------------
/// Convert FFD integration method string to enumeration value
template <>
inline bool FromString(const char *str, FFDIM &value)
{
  // Aliases of methods
  const string lstr = ToLower(str);
  if      (lstr == "scalingandsquaring")        value = FFDIM_SS;
  else if (lstr == "scaling and squaring")      value = FFDIM_SS;
  else if (lstr == "fastscalingandsquaring")    value = FFDIM_FastSS;
  else if (lstr == "fast scaling and squaring") value = FFDIM_FastSS;
  else if (lstr == "euler")                     value = FFDIM_RKE1;
  else if (lstr == "forwardeuler")              value = FFDIM_RKE1;
  else if (lstr == "forward euler")             value = FFDIM_RKE1;
  else if (lstr == "fwdeuler")                  value = FFDIM_RKE1;
  else if (lstr == "modifiedeuler")             value = FFDIM_RKE2;
  else if (lstr == "modified euler")            value = FFDIM_RKE2;
  else if (lstr == "modeuler")                  value = FFDIM_RKE2;
  else if (lstr == "heun")                      value = FFDIM_RKH2;
  else if (lstr == "improvedeuler")             value = FFDIM_RKH2;
  else if (lstr == "improved euler")            value = FFDIM_RKH2;
  else if (lstr == "impeuler")                  value = FFDIM_RKH2;
  else if (lstr == "rk1")                       value = FFDIM_RKE1;
  else if (lstr == "rk2")                       value = FFDIM_RKH2;
  else if (lstr == "rk12")                      value = FFDIM_RKEH12;
  else if (lstr == "rk23")                      value = FFDIM_RKBS23;
  else if (lstr == "rk45")                      value = FFDIM_RKCK45;
  else                                          value = FFDIM_Unknown;
  // Default names
  if (value == FFDIM_Unknown) {
    value = static_cast<FFDIM>(FFDIM_Last - 1);
    while (value != FFDIM_Unknown) {
      if (iequal(ToString(value).c_str(), str)) break;
      value = static_cast<FFDIM>(value - 1);
    }
  }
  // Return whether conversion was successful
  return (value != FFDIM_Unknown);
}


} // namespace mirtk

#endif // MIRTK_FreeFormTransformationIntegrationMethod_H
