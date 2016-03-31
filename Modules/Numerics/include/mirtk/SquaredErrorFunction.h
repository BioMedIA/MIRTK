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

#ifndef MIRTK_SquaredErrorFunction_H
#define MIRTK_SquaredErrorFunction_H

#include "mirtk/RadialErrorFunction.h"


namespace mirtk {


/**
 * Squared Euclidean distance registration error function
 */
class SquaredErrorFunction : public RadialErrorFunction
{
  mirtkObjectMacro(SquaredErrorFunction);

public:

  /// Constructor
  SquaredErrorFunction() {}

  /// Copy constructor
  SquaredErrorFunction(const SquaredErrorFunction &) {}

  /// Destructor
  ~SquaredErrorFunction() {}

  /// Copy construct a new instance
  virtual RadialErrorFunction *NewInstance() const
  {
    return new SquaredErrorFunction(*this);
  }

  /// Type enumeration value
  virtual TypeId Type() const
  {
    return Squared;
  }

  /// Evaluate radial registration error
  virtual double Value(double d) const
  {
    return d;
  }

  /// Evaluate derivative of radial registration error
  virtual double Derivative(double d) const
  {
    return 1.0;
  }

};


} // namespace mirtk

#endif // MIRTK_SquaredErrorFunction_H
