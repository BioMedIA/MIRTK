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

#ifndef MIRTK_PeakSignalToNoiseRatio_H
#define MIRTK_PeakSignalToNoiseRatio_H

#include "mirtk/SumOfSquaredIntensityDifferences.h"


namespace mirtk {


/**
 * Peak signal to noise ratio of two registered images
 */
class PeakSignalToNoiseRatio : public SumOfSquaredIntensityDifferences
{
  mirtkEnergyTermMacro(PeakSignalToNoiseRatio, EM_PSNR);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
public:

  /// Constructor
  PeakSignalToNoiseRatio(const char * = "");

  /// Copy constructor
  PeakSignalToNoiseRatio(const PeakSignalToNoiseRatio &);

  /// Assignment operator
  PeakSignalToNoiseRatio &operator =(const PeakSignalToNoiseRatio &);

  /// Destructor
  ~PeakSignalToNoiseRatio();

  // ---------------------------------------------------------------------------
  // Evaluation
protected:

  /// Evaluate similarity of images
  virtual double Evaluate();

  // ---------------------------------------------------------------------------
  // Debugging
public:

  /// Return unweighted and unnormalized raw energy term value
  /// \remarks Use for progress reporting only.
  virtual double RawValue(double) const;

};


} // namespace mirtk

#endif // MIRTK_PeakSignalToNoiseRatio_H 
