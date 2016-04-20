/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2015 Imperial College London
 * Copyright 2008-2015 Daniel Rueckert, Julia Schnabel
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

#ifndef MIRTK_ScalarFunction_H
#define MIRTK_ScalarFunction_H

#include "mirtk/Object.h"


namespace mirtk {


/**
 * Base class of scalar functions.
 */
class ScalarFunction : public Object
{
  mirtkAbstractMacro(ScalarFunction);

public:

  /// Virtual destructor
  virtual ~ScalarFunction() {}

  /// Evaluation function (pure virtual)
  virtual double Evaluate(double x, double y, double z) = 0;

};


} // namespace mirtk

#endif // MIRTK_ScalarFunction_H
