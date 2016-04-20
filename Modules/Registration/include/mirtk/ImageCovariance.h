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

#ifndef MIRTK_ImageCovariance_H
#define MIRTK_ImageCovariance_H

#include "mirtk/HistogramImageSimilarity.h"


namespace mirtk {


/**
 * Joint entropy image similarity measure
 */
class ImageCovariance : public HistogramImageSimilarity
{
  mirtkEnergyTermMacro(ImageCovariance, EM_CoVar);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
public:

  /// Constructor
  ImageCovariance(const char * = "");

  /// Copy constructor
  ImageCovariance(const ImageCovariance &);

  /// Assignment operator
  ImageCovariance &operator =(const ImageCovariance &);

  /// Destructor
  ~ImageCovariance();

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

#endif // MIRTK_ImageCovariance_H
