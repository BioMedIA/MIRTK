/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2015 Imperial College London
 * Copyright 2013-2015 Stefan Pszczolkowski Parraguez, Andreas Schuh
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

#include "mirtk/RadialErrorFunction.h"

#include "mirtk/String.h"

#include "mirtk/DistanceErrorFunction.h"
#include "mirtk/SquaredErrorFunction.h"
#include "mirtk/GaussianErrorFunction.h"
#include "mirtk/CharbonnierErrorFunction.h"
#include "mirtk/PeronaMalikErrorFunction.h"


namespace mirtk {


// ----------------------------------------------------------------------------
RadialErrorFunction::RadialErrorFunction()
{
}

// ----------------------------------------------------------------------------
RadialErrorFunction::~RadialErrorFunction()
{
}

// ----------------------------------------------------------------------------
RadialErrorFunction *RadialErrorFunction::New(TypeId type)
{
  switch (type) {
    case Distance:    return new DistanceErrorFunction();
    case Squared:     return new SquaredErrorFunction();
    case Gaussian:    return new GaussianErrorFunction();
    case Charbonnier: return new CharbonnierErrorFunction();
    case PeronaMalik: return new PeronaMalikErrorFunction();
    default:
      cerr << "RadialErrorFunction::New: Unknown type = " << type << endl;
      exit(1);
  }
}

// ----------------------------------------------------------------------------
RadialErrorFunction *RadialErrorFunction::New(const char *type_name)
{
  TypeId type = Unknown;
  if (!FromString(type_name, type)) {
    cerr << "RadialErrorFunction::New: Unknown type = " << type_name << endl;
    exit(1);
  }
  return New(type);
}

// -----------------------------------------------------------------------------
template <> bool FromString(const char *str, RadialErrorFunction::TypeId &type)
{
  if (strcmp(str, "Distance") == 0) {
    type = RadialErrorFunction::Distance;
  } else if (strcmp(str, "Square")           == 0 ||
             strcmp(str, "Squared")          == 0 ||
             strcmp(str, "Squared distance") == 0) {
    type = RadialErrorFunction::Squared;
  } else if (strcmp(str, "Gaussian") == 0) {
    type = RadialErrorFunction::Gaussian;
  } else if (strcmp(str, "Charbonnier") == 0) {
    type = RadialErrorFunction::Charbonnier;
  } else if (strcmp(str, "PeronaMalik")      == 0 ||
             strcmp(str, "Perona Malik")     == 0 ||
             strcmp(str, "Perona-Malik")     == 0 ||
             strcmp(str, "Perona and Malik") == 0) {
    type = RadialErrorFunction::PeronaMalik;
  } else {
    type = RadialErrorFunction::Unknown;
  }
  return (type != RadialErrorFunction::Unknown);
}

// -----------------------------------------------------------------------------
template <> string ToString(const RadialErrorFunction::TypeId &type, int w, char c, bool left)
{
  string str;
  switch (type) {
    case RadialErrorFunction::Distance:    str = "Distance";         break;
    case RadialErrorFunction::Squared:     str = "Squared distance"; break;
    case RadialErrorFunction::Gaussian:    str = "Gaussian";         break;
    case RadialErrorFunction::Charbonnier: str = "Charbonnier";      break;
    case RadialErrorFunction::PeronaMalik: str = "Perona-Malik";     break;
    default:                               str = "Unknown";          break;
  }
  return ToString(str, w, c, left);
}


} // namespace mirtk
