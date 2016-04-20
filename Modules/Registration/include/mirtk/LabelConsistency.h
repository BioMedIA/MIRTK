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

#ifndef MIRTK_LabelConsistency_H
#define MIRTK_LabelConsistency_H

#include "mirtk/HistogramImageSimilarity.h"


namespace mirtk {


/**
 * Label consistency image similarity measure
 */
class LabelConsistency : public HistogramImageSimilarity
{
  mirtkEnergyTermMacro(LabelConsistency, EM_LC);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
public:

  /// Constructor
  LabelConsistency(const char * = "");

  /// Copy constructor
  LabelConsistency(const LabelConsistency &);

  /// Assignment operator
  LabelConsistency &operator =(const LabelConsistency &);

  /// Destructor
  ~LabelConsistency();

  // ---------------------------------------------------------------------------
  // Initialization

  /// Initialize similarity measure once input and parameters have been set
  virtual void Initialize();

  // ---------------------------------------------------------------------------
  // Evaluation
protected:

  /// Evaluate similarity of images
  virtual double Evaluate();

};


} // namespace mirtk

#endif // MIRTK_LabelConsistency_H
