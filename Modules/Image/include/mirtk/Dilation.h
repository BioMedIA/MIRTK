/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2015 Imperial College London
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

#ifndef MIRTK_Dilation_H
#define MIRTK_Dilation_H

#include "mirtk/ImageToImage.h"

#include "mirtk/Assert.h"
#include "mirtk/NeighborhoodOffsets.h"


namespace mirtk {


/**
 * morphological erosion of images.
 */
template <class TVoxel>
class Dilation : public ImageToImage<TVoxel>
{
  mirtkImageFilterMacro(Dilation, TVoxel);

  /// What connectivity to assume when running the filter.
  mirtkPublicAttributeMacro(ConnectivityType, Connectivity);

  /// List of voxel offsets of the neighborhood
  mirtkAttributeMacro(NeighborhoodOffsets, Offsets);

public:

  /// Constructor
  Dilation();

  /// Destructor
  virtual ~Dilation();

  /// Run erosion
  virtual void Run();

protected:

  /// Initialize the filter
  virtual void Initialize();

};


// -----------------------------------------------------------------------------
template <class TVoxel>
void Dilate(BaseImage *image, int iterations, ConnectivityType connectivity)
{
  GenericImage<TVoxel> * const im = dynamic_cast<GenericImage<TVoxel> *>(image);
  mirtkAssert(im != nullptr, "template function called with correct type");
  Dilation<TVoxel> dilation;
  dilation.Connectivity(connectivity);
  dilation.Input (im);
  dilation.Output(im);
  for (int i = 0; i < iterations; ++i) dilation.Run();
}


} // namespace mirtk

#endif // MIRTK_Dilation_H
