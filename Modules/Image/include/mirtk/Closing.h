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

#ifndef MIRTK_Closing_H
#define MIRTK_Closing_H

#include "mirtk/ImageToImage.h"

#include "mirtk/Assert.h"
#include "mirtk/NeighborhoodOffsets.h"


namespace mirtk {


/**
 * morphological closing of images.
 */
template <class TVoxel>
class Closing : public ImageToImage<TVoxel>
{
  mirtkInPlaceImageFilterMacro(Closing, TVoxel);

  /// What connectivity to assume when running the filter.
  mirtkPublicAttributeMacro(ConnectivityType, Connectivity);

  /// Number of dilation/erosion iterations
  mirtkPublicAttributeMacro(int, NumberOfIterations);

public:

  /// Constructor
  Closing();

  /// Destructor
  virtual ~Closing();

  /// Run erosion
  virtual void Run();

protected:

  /// Initialize filter
  virtual void Initialize();

};


// -----------------------------------------------------------------------------
template <class TVoxel>
void Close(BaseImage *image, int iterations, ConnectivityType connectivity)
{
  GenericImage<TVoxel> * const im = dynamic_cast<GenericImage<TVoxel> *>(image);
  mirtkAssert(im != nullptr, "template function called with correct type");
  Closing<TVoxel> closing;
  closing.Connectivity(connectivity);
  closing.NumberOfIterations(iterations);
  closing.Input (im);
  closing.Output(im);
  closing.Run();
}


} // namespace mirtk

#endif // MIRTK_Closing_H
