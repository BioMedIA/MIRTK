/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2017 Imperial College London
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

#ifndef MIRTK_HistogramMatching_H
#define MIRTK_HistogramMatching_H

#include "mirtk/ImageToImage.h"


namespace mirtk {


/**
 * Match histogram of image to match distribution of reference image
 *
 *   Nyul, L.G., Udupa, J.K., Xuan Zhang, "New variants of a method of MRI scale standardization",
 *   IEEE TMI, vol.19, no.2, pp.143-150, Feb. 2000 (http://dx.doi.org/10.1109/42.836373)
 *
 * \note Source code adapted from an implementation by Vladimir Fonov in EZminc (https://github.com/vfonov/EZminc)
 * \note Only the foreground intensities of the input images are considered.
 */
template <class TVoxel>
class HistogramMatching : public ImageToImage<TVoxel>
{
  mirtkInPlaceImageFilterMacro(HistogramMatching, TVoxel);

  /// Reference image
  mirtkPublicAggregateMacro(const ImageType, Reference);

  /// Number of histogram bins
  mirtkPublicAttributeMacro(int, NumberOfBins);

  /// Number of linear histogram rescaling intervals
  mirtkPublicAttributeMacro(int, NumberOfSteps);

  /// Cut off value
  mirtkPublicAttributeMacro(double, CutOff);

public:

  /// Constructor
  HistogramMatching();

  /// Run filter on entire image
  virtual void Run();

protected:

  /// Initialize the filter. This function must be called by any derived
  /// filter class to perform some initialize tasks.
  virtual void Initialize();
};


} // namespace mirtk

#endif // MIRTK_HistogramMatching_H
