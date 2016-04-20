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

#ifndef MIRTK_DistanceErrorFunction_H
#define MIRTK_DistanceErrorFunction_H

#include "mirtk/RadialErrorFunction.h"
#include "mirtk/Math.h"


namespace mirtk {


/**
 * Euclidean distance registration error function
 */
class DistanceErrorFunction : public RadialErrorFunction
{
  mirtkObjectMacro(DistanceErrorFunction);

public:

  /// Constructor
  DistanceErrorFunction() {}

  /// Copy constructor
  DistanceErrorFunction(const DistanceErrorFunction &) {}

  /// Destructor
  ~DistanceErrorFunction() {}

  /// Copy construct a new instance
  virtual RadialErrorFunction *NewInstance() const
  {
    return new DistanceErrorFunction(*this);
  }

  /// Type enumeration value
  virtual TypeId Type() const
  {
    return Distance;
  }

  /// Evaluate radial registration error
  virtual double Value(double d) const
  {
    return sqrt(d);
  }

  /// Evaluate derivative of radial registration error
  virtual double Derivative(double d) const
  {
    return (d == .0 ? .0 : 0.5 / sqrt(d));
  }

};


} // namespace mirtk

#endif // MIRTK_DistanceErrorFunction_H
