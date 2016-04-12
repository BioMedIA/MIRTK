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

#ifndef MIRTK_SumOfSquaredIntensityDifferences_H
#define MIRTK_SumOfSquaredIntensityDifferences_H

#include "mirtk/ImageSimilarity.h"


namespace mirtk {


/**
 * Sum of squared differences image similarity measure
 */
class SumOfSquaredIntensityDifferences : public ImageSimilarity
{
  mirtkEnergyTermMacro(SumOfSquaredIntensityDifferences, EM_SSD);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Minimum input target image intensity value
  mirtkReadOnlyAttributeMacro(double, MinTargetIntensity);

  /// Maximum input target image intensity value
  mirtkReadOnlyAttributeMacro(double, MaxTargetIntensity);

  /// Minimum input source image intensity value
  mirtkReadOnlyAttributeMacro(double, MinSourceIntensity);

  /// Maximum input source image intensity value
  mirtkReadOnlyAttributeMacro(double, MaxSourceIntensity);

  /// Maximum squared intensity difference used for normalization
  mirtkReadOnlyAttributeMacro(double, MaxSqDiff);

  /// Sum of squared intensity difference value
  mirtkReadOnlyAttributeMacro(double, SumSqDiff);

  /// Number of foreground voxels for which similarity is evaluated
  mirtkReadOnlyAttributeMacro(int, NumberOfForegroundVoxels);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const SumOfSquaredIntensityDifferences &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
public:

  /// Constructor
  SumOfSquaredIntensityDifferences(const char * = "");

  /// Copy constructor
  SumOfSquaredIntensityDifferences(const SumOfSquaredIntensityDifferences &);

  /// Assignment operator
  SumOfSquaredIntensityDifferences &operator =(const SumOfSquaredIntensityDifferences &);

  /// Destructor
  ~SumOfSquaredIntensityDifferences();

  // ---------------------------------------------------------------------------
  // Initialization

  /// Initialize similarity measure
  virtual void Initialize();

  /// Update moving input image(s) and internal state of similarity measure
  virtual void Update(bool = true);

  // ---------------------------------------------------------------------------
  // Evaluation
protected:

  /// Exclude region from similarity evaluation
  ///
  /// Called by ApproximateGradient \b before the registered image region of
  /// the transformed image is updated.
  virtual void Exclude(const blocked_range3d<int> &);

  /// Include region in similarity evaluation
  ///
  /// Called by ApproximateGradient \b after the registered image region of
  /// the transformed image is updated.
  virtual void Include(const blocked_range3d<int> &);

  /// Evaluate similarity of images
  virtual double Evaluate();

  /// Evaluate non-parametric similarity gradient w.r.t the given image
  virtual bool NonParametricGradient(const RegisteredImage *, GradientImageType *);

};


} // namespace mirtk

#endif // MIRTK_SumOfSquaredIntensityDifferences_H 
