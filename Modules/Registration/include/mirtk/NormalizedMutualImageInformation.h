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

#ifndef MIRTK_NormalizedMutualImageInformation_H
#define MIRTK_NormalizedMutualImageInformation_H

#include "mirtk/HistogramImageSimilarity.h"

#include "mirtk/Histogram1D.h"


namespace mirtk {


/**
 * Normalized mutual information image similarity measure
 */
class NormalizedMutualImageInformation : public HistogramImageSimilarity
{
  mirtkEnergyTermMacro(NormalizedMutualImageInformation, EM_NMI);

  // ---------------------------------------------------------------------------
  // Types

  /// Type of marginal histograms
  typedef Histogram1D<BinType> MarginalHistogram;

  // ---------------------------------------------------------------------------
  // Attributes

  /// Marginal histogram w.r.t. target image intensities
  mirtkAttributeMacro(MarginalHistogram, TargetHistogram);

  /// Marginal histogram w.r.t. source image intensities
  mirtkAttributeMacro(MarginalHistogram, SourceHistogram);

  /// Marginal entropy w.r.t. target image intensities
  mirtkAttributeMacro(double, TargetEntropy);

  /// Marginal entropy w.r.t. source image intensities
  mirtkAttributeMacro(double, SourceEntropy);

  /// Joint entropy of target and source image intensities
  mirtkAttributeMacro(double, JointEntropy);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const NormalizedMutualImageInformation &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
public:

  /// Constructor
  NormalizedMutualImageInformation(const char * = "");

  /// Copy constructor
  NormalizedMutualImageInformation(const NormalizedMutualImageInformation &);

  /// Assignment operator
  NormalizedMutualImageInformation &operator =(const NormalizedMutualImageInformation &);

  /// Destructor
  ~NormalizedMutualImageInformation();

  // ---------------------------------------------------------------------------
  // Evaluation

  /// Update moving image and internal state of similarity measure
  virtual void Update(bool = true);

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

};


} // namespace mirtk

#endif // MIRTK_NormalizedMutualImageInformation_H
