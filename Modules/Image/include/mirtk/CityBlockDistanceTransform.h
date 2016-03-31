/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2011-2016 Imperial College London
 * Copyright 2011      Paul Aljabar
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

#ifndef MIRTK_CityBlockDistanceTransform_H
#define MIRTK_CityBlockDistanceTransform_H

#include "mirtk/ImageToImage.h"

#include "mirtk/NeighborhoodOffsets.h"


namespace mirtk {


/**
 * City block image distance transform
 *
 * Find the City Block (Manhattan, L1) distance for all object voxels in an
 * image from the boundary. Object voxels have a value greater than zero and
 * background voxels are the rest. The distance map is initialised to zero.
 * In each iteration, the distance map is incremented by 1 for all object
 * voxels. The border voxels are removed and the the next iteration starts.
 * If any dimension is a singleton, a 2D version is applied.
 */
template <class TVoxel>
class CityBlockDistanceTransform : public ImageToImage<TVoxel>
{
  mirtkInPlaceImageFilterMacro(CityBlockDistanceTransform, TVoxel);

  // ---------------------------------------------------------------------------
  // Types

  /// In the 2D case, a flip may be necessary so that the singleton dimension
  /// is the z-direction. This makes processing easier.
  enum FlipType { FlipNone, FlipXY,  FlipXZ,  FlipYZ };

  // ---------------------------------------------------------------------------
  // Attriubutes

  /// Storage for the object voxels that can be updated.
  GreyImage _data;

  /// List of voxel offsets of the 6 neighbourhood of a voxel. Enables checking to
  /// see if there is a face-neighbour labelled background during each iteration.
  NeighborhoodOffsets _offsets;

  FlipType _flipType;

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Default constructor
  CityBlockDistanceTransform();

  /// Destructor
  virtual ~CityBlockDistanceTransform() {};

  // ---------------------------------------------------------------------------
  // Execution

  /// Run distance transform
  virtual void Run();

protected:

  /// Run distance transform in 2D
  virtual void Run2D();

  /// Run distance transform in 3D
  virtual void Run3D();

  /// Initialize the filter
  virtual void Initialize();

  /// Initialize the filter
  virtual void Initialize2D();

  /// Initialize the filter
  virtual void Initialize3D();

  /// Finalize the filter
  virtual void Finalize();

};


} // namespace mirtk

#endif // MIRTK_CityBlockDistanceTransform_H
