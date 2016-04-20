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

#ifndef MIRTK_PeronaMalikErrorFunction_H
#define MIRTK_PeronaMalikErrorFunction_H

#include "mirtk/RadialErrorFunction.h"
#include "mirtk/Math.h"


namespace mirtk {


/**
 * Perona and Malik fiducial registration error function
 */
class PeronaMalikErrorFunction : public RadialErrorFunction
{
  mirtkObjectMacro(PeronaMalikErrorFunction);

  /// Squared fiducial registration error threshold
  mirtkPublicAttributeMacro(double, SquaredThreshold);

public:

  /// Constructor
  PeronaMalikErrorFunction(double threshold = 1.0)
  :
    _SquaredThreshold(threshold * threshold)
  {}

  /// Copy constructor
  PeronaMalikErrorFunction(const PeronaMalikErrorFunction &other)
  :
    _SquaredThreshold(other._SquaredThreshold)
  {}

  /// Destructor
  ~PeronaMalikErrorFunction() {}

  /// Copy construct a new instance
  virtual RadialErrorFunction *NewInstance() const
  {
    return new PeronaMalikErrorFunction(*this);
  }

  /// Type enumeration value
  virtual TypeId Type() const
  {
    return PeronaMalik;
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
    return _SquaredThreshold * log(1.0 + d / _SquaredThreshold);
  }

  /// Evaluate derivative of radial registration error
  virtual double Derivative(double d) const
  {
    return 1.0 / (1.0 + d / _SquaredThreshold);
  }

};


} // namespace mirtk

#endif // MIRTK_PeronaMalikErrorFunction_H
