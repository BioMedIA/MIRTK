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

#ifndef MIRTK_CharbonnierErrorFunction_H
#define MIRTK_CharbonnierErrorFunction_H

#include "mirtk/RadialErrorFunction.h"

#include "mirtk/Math.h"


namespace mirtk {


/**
 * Charbonnier fiducial registration error function
 */
class CharbonnierErrorFunction : public RadialErrorFunction
{
  mirtkObjectMacro(CharbonnierErrorFunction);

  /// Squared fiducial registration error threshold
  mirtkPublicAttributeMacro(double, SquaredThreshold);
 
public:

  /// Constructor
  CharbonnierErrorFunction(double threshold = 1.0)
  :
    _SquaredThreshold(threshold * threshold)
  {}

  /// Copy constructor
  CharbonnierErrorFunction(const CharbonnierErrorFunction &other)
  :
    _SquaredThreshold(other._SquaredThreshold)
  {}

  /// Copy construct a new instance
  virtual RadialErrorFunction *NewInstance() const
  {
    return new CharbonnierErrorFunction(*this);
  }

  /// Type enumeration value
  virtual TypeId Type() const
  {
    return Charbonnier;
  }

  /// Set parameter value from string
  bool Set(const char *name, const char *value)
  {
    if (strcmp(name, "Threshold") == 0) {
      double threshold = .0;
      if (!FromString(value, threshold) || threshold <= .0) return false;
      _SquaredThreshold = threshold * threshold;
      return true;
    } else if (strcmp(name, "Squared threshold") == 0) {
      return FromString(value, _SquaredThreshold) && _SquaredThreshold <= .0;
    }
    return false;
  }

  // Import other overloads
  using RadialErrorFunction::Parameter;

  /// Get parameter key/value as string map
  ParameterList Parameter() const
  {
    ParameterList params;
    Insert(params, "Threshold", ToString(sqrt(_SquaredThreshold)));
    return params;
  }

  /// Evaluate radial registration error
  virtual double Value(double d) const
  {
    return 2.0 * _SquaredThreshold * (sqrt(1.0 + d / _SquaredThreshold) - 1.0);
  }

  /// Evaluate derivative of radial registration error
  virtual double Derivative(double d) const
  {
    return 1.0 / sqrt(1.0 + d / _SquaredThreshold);
  }

};


} // namespace mirtk

#endif // MIRTK_CharbonnierErrorFunction_H
