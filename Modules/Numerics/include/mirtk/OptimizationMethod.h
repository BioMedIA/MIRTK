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

#ifndef MIRTK_OptimizationMethod_H
#define MIRTK_OptimizationMethod_H

#include "mirtk/String.h"

#include <functional>


namespace mirtk {


// -----------------------------------------------------------------------------
/// Enumeration of available optimization methods
enum OptimizationMethod
{
  OM_DownhillDescent,
  OM_GradientDescent,
  OM_GradientDescentConstrained,
  OM_SteepestGradientDescent,
  OM_ConjugateGradientDescent,
  OM_ClosedForm,
  OM_LBFGS,
  OM_LineSearch,
  OM_EulerMethod,                ///< Explicit Euler method for deformable surface models
  OM_EulerMethodWithDamping,     ///< Explicit Euler method with momentum for deformable surface models
  OM_EulerMethodWithMomentum     ///< Explicit Euler method with momentum for deformable surface models
};

// -----------------------------------------------------------------------------
template <>
inline string ToString(const OptimizationMethod &m, int w, char c, bool left)
{
  const char *str;
  switch (m) {
    case OM_DownhillDescent:            str = "DownhillDescent"; break;
    case OM_GradientDescent:            str = "GradientDescent"; break;
    case OM_GradientDescentConstrained: str = "GradientDescentConstrained"; break;
    case OM_SteepestGradientDescent:    str = "SteepestGradientDescent"; break;
    case OM_ConjugateGradientDescent:   str = "ConjugateGradientDescent"; break;
    case OM_ClosedForm:                 str = "ClosedForm"; break;
    case OM_LBFGS:                      str = "LBFGS"; break;
    case OM_LineSearch:                 str = "LineSearch"; break;
    case OM_EulerMethod:                str = "EulerMethod"; break;
    case OM_EulerMethodWithDamping:     str = "EulerMethodWithDamping"; break;
    case OM_EulerMethodWithMomentum:    str = "EulerMethodWithMomentum"; break;
    default:                            str = "Unknown"; break;
  }
  return ToString(str, w, c, left);
}

// -----------------------------------------------------------------------------
template <>
inline bool FromString(const char *str, OptimizationMethod &m)
{
  if (strcmp(str, "DownhillDescent") == 0 ||
      strcmp(str, "Downhill")        == 0) {
    m = OM_DownhillDescent;
  } else if (strcmp(str, "GradientDescent") == 0 ||
           strcmp(str, "GD")              == 0) {
    m = OM_GradientDescent;
  } else if (strcmp(str, "ConstrainedGradientDescent") == 0 ||
             strcmp(str, "GradientDescentConstrained") == 0) {
    m = OM_GradientDescentConstrained;
  } else if (strcmp(str, "SteepestGradientDescent") == 0 ||
             strcmp(str, "SteepestGradient")        == 0 ||
             strcmp(str, "SGD")                     == 0) {
    m = OM_SteepestGradientDescent;
  } else if (strcmp(str, "ConjugateGradientDescent") == 0 ||
             strcmp(str, "ConjugateGradient")        == 0 ||
             strcmp(str, "CGD")                      == 0) {
    m = OM_ConjugateGradientDescent;
  } else if (strcmp(str, "ClosedForm") == 0) {
    m = OM_ClosedForm;
  } else if (strcmp(str, "LBFGS") == 0) {
    m = OM_LBFGS;
  } else if (strcmp(str, "LineSearch") == 0) {
    m = OM_LineSearch;
  } else if (strcmp(str, "EulerMethod") == 0) {
    m = OM_EulerMethod;
  } else if (strcmp(str, "EulerMethodWithDamping") == 0) {
    m = OM_EulerMethodWithDamping;
  } else if (strcmp(str, "EulerMethodWithMomentum") == 0) {
    m = OM_EulerMethodWithMomentum;
  } else {
    return false;
  }
  return true;
}

} // namespace mirtk


namespace std {

/// Compute hash value from OptimizationMethod enumeration value
template<>
struct hash<mirtk::OptimizationMethod> {
    size_t operator()(const mirtk::OptimizationMethod &enum_value) const {
        return std::hash<int>()(enum_value);
    }
};


} // namespace std


#endif // MIRTK_OptimizationMethod_H
