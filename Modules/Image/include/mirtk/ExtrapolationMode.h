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

#ifndef MIRTK_ExtrapolationMode_H
#define MIRTK_ExtrapolationMode_H

#include "mirtk/String.h"


namespace mirtk {


// ----------------------------------------------------------------------------
/// Image extrapolation modes
enum ExtrapolationMode {
  Extrapolation_Default, // i.e., not specified
  Extrapolation_None,    // i.e., partial interpolation or default value
  Extrapolation_Const,
  Extrapolation_NN,
  Extrapolation_Repeat,  // i.e., periodic
  Extrapolation_Mirror,
  Extrapolation_ConstWithPeriodicTime,
  // Add new enumeration values above
  Extrapolation_Last
};

// ----------------------------------------------------------------------------
template <>
inline string ToString(const ExtrapolationMode &m, int w, char c, bool left)
{
  const char *str;
  switch(m) {
    case Extrapolation_Default:               str = "Default"; break;
    case Extrapolation_None:                  str = "None"; break;
    case Extrapolation_Const:                 str = "Const"; break;
    case Extrapolation_NN:                    str = "NN"; break;
    case Extrapolation_Repeat:                str = "Repeat"; break;
    case Extrapolation_Mirror:                str = "Mirror"; break;
    case Extrapolation_ConstWithPeriodicTime: str = "ConstWithPeriodicTime"; break;
    default:                                  str = "Unknown"; break;
  }
  return ToString(str, w, c, left);
}

// ----------------------------------------------------------------------------
/// Get corresponding extrapolation with periodic time
inline ExtrapolationMode ExtrapolationWithPeriodicTime(ExtrapolationMode m)
{
  switch(m) {
    case Extrapolation_ConstWithPeriodicTime:
    case Extrapolation_Const:  return Extrapolation_ConstWithPeriodicTime;
    case Extrapolation_Repeat: return Extrapolation_Repeat;
    default:                   return Extrapolation_None;
  }
}

// ----------------------------------------------------------------------------
/// Get corresponding extrapolation without periodic time
inline ExtrapolationMode ExtrapolationWithoutPeriodicTime(ExtrapolationMode m)
{
  switch(m) {
    case Extrapolation_Const:
    case Extrapolation_ConstWithPeriodicTime: return Extrapolation_Const;
    // Note: Extrapolation_Repeat remains periodic in all dimensions
    default: return m;
  }
}

// ----------------------------------------------------------------------------
template <>
inline bool FromString(const char *str, ExtrapolationMode &m)
{
  string lstr = ToLower(str);
  if      (lstr == "default"                 ) m = Extrapolation_Default;
  else if (lstr == "none"                    ) m = Extrapolation_None;
  else if (lstr == "const"                   ) m = Extrapolation_Const;
  else if (lstr == "nn"                      ) m = Extrapolation_NN;
  else if (lstr == "repeat" || lstr == "tile") m = Extrapolation_Repeat;
  else if (lstr == "mirror"                  ) m = Extrapolation_Mirror;
  else if (lstr == "constwithperiodictime"   ) m = Extrapolation_ConstWithPeriodicTime;
  else return false;
  return true;
}


} // namespace mirtk

#endif // MIRTK_ExtrapolationMode_H
