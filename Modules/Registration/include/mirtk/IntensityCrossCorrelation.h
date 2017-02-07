/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2017 Imperial College London
 * Copyright 2013-2017 Andreas Schuh
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

#ifndef MIRTK_IntensityCrossCorrelation_H
#define MIRTK_IntensityCrossCorrelation_H

#include "mirtk/NormalizedIntensityCrossCorrelation.h"


namespace mirtk {


/**
 * Normalized cross correlation image similarity measure
 *
 * \deprecated Use NormalizedIntensityCrossCorrelation with (default) window size zero instead.
 */
class IntensityCrossCorrelation : public NormalizedIntensityCrossCorrelation
{
  mirtkEnergyTermMacro(IntensityCrossCorrelation, EM_CC);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
public:

  /// Constructor
  IntensityCrossCorrelation(const char * = "");

  /// Copy constructor
  IntensityCrossCorrelation(const IntensityCrossCorrelation &);

  /// Assignment operator
  IntensityCrossCorrelation &operator =(const IntensityCrossCorrelation &);

  /// Destructor
  virtual ~IntensityCrossCorrelation();

};


} // namespace mirtk

#endif // MIRTK_IntensityCrossCorrelation_H
