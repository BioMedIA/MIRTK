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

#ifndef MIRTK_ImageSurfaceStatistics_H
#define MIRTK_ImageSurfaceStatistics_H

#include "mirtk/SurfaceFilter.h"

#include "mirtk/Array.h"
#include "mirtk/Memory.h"
#include "mirtk/DataStatistics.h"
#include "mirtk/InterpolateImageFunction.h"


namespace mirtk {


/**
 * Extract image patches at surface mesh points and compute patch statistics
 *
 * This filter samples the values of a given image at a local patch centered at
 * the points of a surface mesh and (optionally) computes common image statistics
 * of these local image samples. The rectangular image patches can be either all
 * aligned with the global world coordinate axes, the global image coordinate axes,
 * or the local tangent space coordinate axes at each surface point.
 */
class ImageSurfaceStatistics : public SurfaceFilter
{
  mirtkObjectMacro(ImageSurfaceStatistics);

  // ---------------------------------------------------------------------------
  // Types
public:

  enum Space
  {
    ImageSpace,
    WorldSpace,
    TangentSpace
  };

  // ---------------------------------------------------------------------------
  // Attributes
private:

  /// Continuous scalar image
  mirtkPublicAttributeMacro(SharedPtr<const InterpolateImageFunction>, Image);

  /// Name of output point data array
  mirtkPublicAttributeMacro(string, ArrayName);

  /// Coordinate system within which to extract image patches
  mirtkPublicAttributeMacro(Space, PatchSpace);

  /// Size of each image patch centered at the mesh points
  mirtkPublicAttributeMacro(int3, PatchSize);

  /// Spacing of patch samples in world units
  mirtkPublicAttributeMacro(double3, PatchSpacing);

  /// Whether to include patch samples in output
  mirtkPublicAttributeMacro(bool, PatchSamples);

  /// Whether to subtract mean of patch samples from sample values
  mirtkPublicAttributeMacro(bool, DemeanSamples);

  /// Whether to divide sample values by standard deviation of patch samples
  mirtkPublicAttributeMacro(bool, WhitenSamples);

  /// Functor objects used to compute desired patch statistics
  mirtkPublicAttributeMacro(Array<SharedPtr<data::Statistic> >, Statistics);

  /// Copy attributes of this filter from another instance
  void CopyAttributes(const ImageSurfaceStatistics &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Default constructor
  ImageSurfaceStatistics();

  /// Copy constructor
  ImageSurfaceStatistics(const ImageSurfaceStatistics &);

  /// Assignment operator
  ImageSurfaceStatistics &operator =(const ImageSurfaceStatistics &);

  /// Destructor
  virtual ~ImageSurfaceStatistics();

  // ---------------------------------------------------------------------------
  // Execution

protected:

  /// Execute filter
  virtual void Execute();

private:

  /// Helper used for implementation of Execute function
  template <class Body> void ExecuteInParallel(Body &);

};


} // namespace mirtk

#endif // MIRTK_ImageSurfaceStatistics_H
