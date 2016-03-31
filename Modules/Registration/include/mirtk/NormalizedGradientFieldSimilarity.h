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

#ifndef MIRTK_NormalizedGradientFieldSimilarity_H
#define MIRTK_NormalizedGradientFieldSimilarity_H

#include "mirtk/GradientFieldSimilarity.h"


namespace mirtk {


/**
 * Base class for normalized gradient image similarity measures
 *
 * Subclasses of this intensity-based image similarity measure compute similarity
 * from the normalized gradient field of the images.
 */
class NormalizedGradientFieldSimilarity : public GradientFieldSimilarity
{
  mirtkAbstractMacro(NormalizedGradientFieldSimilarity);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Noise parameter for target image
  mirtkPublicAttributeMacro(double, TargetNoise);

  /// Noise parameter for source image
  mirtkPublicAttributeMacro(double, SourceNoise);

  /// Normalization constant for target image
  mirtkPublicAttributeMacro(double, TargetNormalization);

  /// Normalization constant for source image
  mirtkPublicAttributeMacro(double, SourceNormalization);

  /// Compute normalization constant for given image and noise value
  double NormalizationFactor(RegisteredImage *, double = 1.0) const;

  // ---------------------------------------------------------------------------
  // Construction/Destruction

protected:

  /// Constructor
  NormalizedGradientFieldSimilarity(const char * = "", double = 1.0);

  /// Copy constructor
  NormalizedGradientFieldSimilarity(const NormalizedGradientFieldSimilarity &);

  /// Destructor
  virtual ~NormalizedGradientFieldSimilarity();

  // ---------------------------------------------------------------------------
  // Parameters

protected:

  /// Set parameter value from string
  virtual bool SetWithPrefix(const char *, const char *);

public:

  // Import other overloads
  using GradientFieldSimilarity::Parameter;

  /// Get parameter name/value map
  virtual ParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Evaluation

  /// Update moving input image(s) and internal state of similarity measure
  virtual void Update(bool = true);

};


} // namespace mirtk

#endif // MIRTK_NormalizedGradientFieldSimilarity_H
