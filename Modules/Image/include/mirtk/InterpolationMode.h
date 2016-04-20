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

#ifndef MIRTK_InterpolationMode_H
#define MIRTK_InterpolationMode_H

#include "mirtk/String.h"


namespace mirtk {


// ----------------------------------------------------------------------------
/// Image interpolation modes
enum InterpolationMode {
  Interpolation_NN,
  Interpolation_Linear,
  Interpolation_FastLinear,
  Interpolation_BSpline,
  Interpolation_CubicBSpline,
  Interpolation_FastCubicBSpline,
  Interpolation_CSpline,
  Interpolation_SBased,
  Interpolation_Sinc,
  Interpolation_Gaussian,
  Interpolation_NNWithPadding,
  Interpolation_LinearWithPadding,
  Interpolation_FastLinearWithPadding,
  Interpolation_BSplineWithPadding,
  Interpolation_CubicBSplineWithPadding,
  Interpolation_FastCubicBSplineWithPadding,
  Interpolation_CSplineWithPadding,
  Interpolation_SBasedWithPadding,
  Interpolation_SincWithPadding,
  Interpolation_GaussianWithPadding,
  // Add new enumeration values above
  Interpolation_Last
};

// ----------------------------------------------------------------------------
template <>
inline string ToString(const InterpolationMode &m, int w, char c, bool left)
{
  const char *str;
  switch(m) {
    case Interpolation_NN:                          str = "NN"; break;
    case Interpolation_Linear:                      str = "Linear"; break;
    case Interpolation_FastLinear:                  str = "Fast linear"; break;
    case Interpolation_BSpline:                     str = "BSpline"; break;
    case Interpolation_CSpline:                     str = "CSpline"; break;
    case Interpolation_CubicBSpline:                str = "Cubic BSpline"; break;
    case Interpolation_FastCubicBSpline:            str = "Fast cubic BSpline"; break;
    case Interpolation_SBased:                      str = "SBased"; break;
    case Interpolation_Sinc:                        str = "Sinc"; break;
    case Interpolation_Gaussian:                    str = "Gaussian"; break;
    case Interpolation_NNWithPadding:               str = "NN with padding"; break;
    case Interpolation_LinearWithPadding:           str = "Linear with padding"; break;
    case Interpolation_FastLinearWithPadding:       str = "Fast linear with padding"; break;
    case Interpolation_BSplineWithPadding:          str = "BSpline with padding"; break;
    case Interpolation_CubicBSplineWithPadding:     str = "Cubic BSpline with padding"; break;
    case Interpolation_FastCubicBSplineWithPadding: str = "Fast cubic BSpline with padding"; break;
    case Interpolation_CSplineWithPadding:          str = "CSpline with padding"; break;
    case Interpolation_SBasedWithPadding:           str = "SBased with padding"; break;
    case Interpolation_SincWithPadding:             str = "Sinc with padding"; break;
    case Interpolation_GaussianWithPadding:         str = "Gaussian with padding"; break;
    default:                                        str = "Unknown"; break;
  }
  return ToString(str, w, c, left);
}

// ----------------------------------------------------------------------------
template <>
inline bool FromString(const char *str, InterpolationMode &m)
{
  const string lstr = ToLower(str);
  if      (lstr == "nn") m = Interpolation_NN;
  else if (lstr == "linear") m = Interpolation_Linear;
  else if (lstr == "fast linear") m = Interpolation_FastLinear;
  else if (lstr == "bspline") m = Interpolation_BSpline;
  else if (lstr == "cspline") m = Interpolation_CSpline;
  else if (lstr == "cubic bspline") m = Interpolation_CubicBSpline;
  else if (lstr == "fast cubic bspline") m = Interpolation_FastCubicBSpline;
  else if (lstr == "sbased") m = Interpolation_SBased;
  else if (lstr == "sinc") m = Interpolation_Sinc;
  else if (lstr == "gaussian") m = Interpolation_Gaussian;
  else if (lstr == "nn with padding") m = Interpolation_NNWithPadding;
  else if (lstr == "nn (with padding)") m = Interpolation_NNWithPadding;
  else if (lstr == "linear with padding") m = Interpolation_LinearWithPadding;
  else if (lstr == "linear (with padding)") m = Interpolation_LinearWithPadding;
  else if (lstr == "fast linear with padding") m = Interpolation_FastLinearWithPadding;
  else if (lstr == "fast linear (with padding)") m = Interpolation_FastLinearWithPadding;
  else if (lstr == "bspline with padding") m = Interpolation_BSplineWithPadding;
  else if (lstr == "bspline (with padding)") m = Interpolation_BSplineWithPadding;
  else if (lstr == "cubic bspline with padding") m = Interpolation_CubicBSplineWithPadding;
  else if (lstr == "cubic bspline (with padding)") m = Interpolation_CubicBSplineWithPadding;
  else if (lstr == "fast cubic bspline with padding") m = Interpolation_FastCubicBSplineWithPadding;
  else if (lstr == "fast cubic bspline (with padding)") m = Interpolation_FastCubicBSplineWithPadding;
  else if (lstr == "cspline with padding") m = Interpolation_CSplineWithPadding;
  else if (lstr == "cspline (with padding)") m = Interpolation_CSplineWithPadding;
  else if (lstr == "sbased with padding") m = Interpolation_SBasedWithPadding;
  else if (lstr == "sbased (with padding)") m = Interpolation_SBasedWithPadding;
  else if (lstr == "sinc with padding") m = Interpolation_SincWithPadding;
  else if (lstr == "sinc (with padding)") m = Interpolation_SincWithPadding;
  else if (lstr == "gaussian with padding") m = Interpolation_GaussianWithPadding;
  else if (lstr == "gaussian (with padding)") m = Interpolation_GaussianWithPadding;
  else return false;
  return true;
}

// ----------------------------------------------------------------------------
/// Get corresponding interpolation with padding
inline InterpolationMode InterpolationWithPadding(InterpolationMode m)
{
  switch(m) {
    case Interpolation_NN:               return Interpolation_NNWithPadding;
    case Interpolation_Linear:           return Interpolation_LinearWithPadding;
    case Interpolation_FastLinear:       return Interpolation_FastLinearWithPadding;
    case Interpolation_BSpline:          return Interpolation_BSplineWithPadding;
    case Interpolation_CubicBSpline:     return Interpolation_CubicBSplineWithPadding;
    case Interpolation_FastCubicBSpline: return Interpolation_FastCubicBSplineWithPadding;
    case Interpolation_CSpline:          return Interpolation_CSplineWithPadding;
    case Interpolation_SBased:           return Interpolation_SBasedWithPadding;
    case Interpolation_Sinc:             return Interpolation_SincWithPadding;
    case Interpolation_Gaussian:         return Interpolation_GaussianWithPadding;
    default:                             return m;
  }
}

// ----------------------------------------------------------------------------
/// Get corresponding interpolation without padding
inline InterpolationMode InterpolationWithoutPadding(InterpolationMode m)
{
  switch(m) {
    case Interpolation_NNWithPadding:               return Interpolation_NN;
    case Interpolation_LinearWithPadding:           return Interpolation_Linear;
    case Interpolation_FastLinearWithPadding:       return Interpolation_FastLinear;
    case Interpolation_BSplineWithPadding:          return Interpolation_BSpline;
    case Interpolation_CubicBSplineWithPadding:     return Interpolation_CubicBSpline;
    case Interpolation_FastCubicBSplineWithPadding: return Interpolation_FastCubicBSpline;
    case Interpolation_CSplineWithPadding:          return Interpolation_CSpline;
    case Interpolation_SBasedWithPadding:           return Interpolation_SBased;
    case Interpolation_SincWithPadding:             return Interpolation_Sinc;
    case Interpolation_GaussianWithPadding:         return Interpolation_Gaussian;
    default:                                        return m;
  }
}


} // namespace mirtk

#endif // MIRTK_InterpolationMode_H
