/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2015 Imperial College London
 * Copyright 2008-2015 Daniel Rueckert, Julia Schnabel
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

#ifndef MIRTK_ScalarFunctionToImage_H
#define MIRTK_ScalarFunctionToImage_H

#include "mirtk/Object.h"
#include "mirtk/GenericImage.h"
#include "mirtk/ScalarFunction.h"


namespace mirtk {


/**
 * Scalar function to image filter.
 *
 * This class uses a scalar function to produce an image as output. The
 * filter loops through each voxel of the output image and calculates its
 * intensity as the value of the scalar function as a function of spatial
 * location.
 */
template <class VoxelType>
class ScalarFunctionToImage : public Object
{
  mirtkObjectMacro(ScalarFunctionToImage);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Input for the filter
  mirtkPublicAggregateMacro(ScalarFunction, Input);

  /// Output for the filter
  mirtkPublicAggregateMacro(GenericImage<VoxelType>, Output);

  /// Whether to print debug information
  mirtkPublicAttributeMacro(bool, DebugFlag);

  /// Flag to use world or image coordinates for scalar function evaluation
  mirtkPublicAttributeMacro(bool, UseWorldCoordinates);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Constructor (using world coordinates by default)
  ScalarFunctionToImage(bool = true);

  /// Deconstuctor
  virtual ~ScalarFunctionToImage();

  /// Run the filter on entire image
  virtual void Run();

};


} // namespace mirtk

#endif // MIRTK_ScalarFunctionToImage_H
