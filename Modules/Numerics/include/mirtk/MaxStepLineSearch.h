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

#ifndef MIRTK_MaxStepLineSearch_H
#define MIRTK_MaxStepLineSearch_H

#include "mirtk/InexactLineSearch.h"


namespace mirtk {


/**
 * Dummy line search which always takes the maximum step
 *
 * This line search implements the LS_None line search strategy.
 */
class MaxStepLineSearch : public InexactLineSearch
{
  mirtkLineSearchMacro(MaxStepLineSearch, LS_None);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
public:

  /// Constructor
  MaxStepLineSearch(ObjectiveFunction * = NULL);

  /// Copy constructor
  MaxStepLineSearch(const MaxStepLineSearch &);

  /// Assignment operator
  MaxStepLineSearch &operator =(const MaxStepLineSearch &);

  /// Destructor
  virtual ~MaxStepLineSearch();

  // ---------------------------------------------------------------------------
  // Optimization

  /// Make optimal step along search direction
  virtual double Run();

};


} // namespace mirtk

#endif // MIRTK_MaxStepLineSearch_H
