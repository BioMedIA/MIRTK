/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2016 Imperial College London
 * Copyright 2016 Andreas Schuh
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

#ifndef MIRTK_IntensityCorrelationRatioYX_H
#define MIRTK_IntensityCorrelationRatioYX_H

#include "mirtk/HistogramImageSimilarity.h"


namespace mirtk {


/**
 * Correlation ratio of target (X) and source (Y) image intensities
 */
class IntensityCorrelationRatioYX : public HistogramImageSimilarity
{
  mirtkEnergyTermMacro(IntensityCorrelationRatioYX, EM_CR_YX);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
public:

  /// Constructor
  IntensityCorrelationRatioYX(const char * = "");

  /// Copy constructor
  IntensityCorrelationRatioYX(const IntensityCorrelationRatioYX &);

  /// Assignment operator
  IntensityCorrelationRatioYX &operator =(const IntensityCorrelationRatioYX &);

  /// Destructor
  ~IntensityCorrelationRatioYX();

  // ---------------------------------------------------------------------------
  // Evaluation
protected:

  /// Evaluate similarity of images
  virtual double Evaluate();

};


} // namespace mirtk

#endif // MIRTK_IntensityCorrelationRatioYX_H
