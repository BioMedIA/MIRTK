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

#ifndef MIRTK_CosineOfNormalizedGradientField_H
#define MIRTK_CosineOfNormalizedGradientField_H

#include "mirtk/NormalizedGradientFieldSimilarity.h"


namespace mirtk {


/**
 * Cosine of normalized gradient field similarity measure
 */
class CosineOfNormalizedGradientField : public NormalizedGradientFieldSimilarity
{
  mirtkEnergyTermMacro(CosineOfNormalizedGradientField, EM_NGF_COS);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Power (exponent) of cosine function
  mirtkPublicAttributeMacro(int, Power);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
public:

  /// Constructor
  CosineOfNormalizedGradientField(const char * = "");

  /// Copy constructor
  CosineOfNormalizedGradientField(const CosineOfNormalizedGradientField &);

  /// Assignment operator
  CosineOfNormalizedGradientField &operator =(const CosineOfNormalizedGradientField &);

  /// Destructor
  ~CosineOfNormalizedGradientField();

  // ---------------------------------------------------------------------------
  // Parameters
protected:

  /// Set parameter value from string
  virtual bool SetWithoutPrefix(const char *, const char *);

public:

  // Import other overloads
  using ImageSimilarity::Parameter;

  /// Get parameter key/value as string map
  virtual ParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Evaluation
protected:

  /// Evaluate similarity of images
  virtual double Evaluate();

  /// Evaluate non-parametric similarity gradient w.r.t the given image
  virtual bool NonParametricGradient(const RegisteredImage *, GradientImageType *);

  // ---------------------------------------------------------------------------
  // Debugging
public:

  /// Return unweighted and unnormalized raw energy term value
  /// \remarks Use for progress reporting only.
  virtual double RawValue(double) const;

  /// Write input of data fidelity term
  virtual void WriteDataSets(const char *, const char *, bool = true) const;

};


} // namespace mirtk

#endif // MIRTK_CosineOfNormalizedGradientField_H
