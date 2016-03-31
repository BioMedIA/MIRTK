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

#ifndef MIRTK_ConstExtrapolateImageFunctionWithPeriodicTime_H
#define MIRTK_ConstExtrapolateImageFunctionWithPeriodicTime_H

#include "mirtk/BaseImage.h"
#include "mirtk/ConstExtrapolateImageFunction.h"
#include "mirtk/RepeatExtrapolateImageFunction.h"


namespace mirtk {


/**
 * Extrapolation of generic image by padding it with a constant value
 * in the spatial domain, but a periodic extension in the temporal dimension
 */
template <class TImage>
class GenericConstExtrapolateImageFunctionWithPeriodicTime
: public GenericConstExtrapolateImageFunction<TImage>
{
  mirtkExtrapolatorMacro(
    GenericConstExtrapolateImageFunctionWithPeriodicTime,
    Extrapolation_ConstWithPeriodicTime
  );

public:

  typedef TImage                        ImageType; ///< Input image type
  typedef typename ImageType::VoxelType VoxelType; ///< Input voxel type
  typedef typename ImageType::RealType  RealType;  ///< Compatible floating-point type

  /// Constructor
  GenericConstExtrapolateImageFunctionWithPeriodicTime(double padding_value = .0)
  :
    GenericConstExtrapolateImageFunction<TImage>(padding_value)
  {}

  /// Destructor
  virtual ~GenericConstExtrapolateImageFunctionWithPeriodicTime() {}

  // Import overloaded non-virtual member functions from base class
  using GenericConstExtrapolateImageFunction<TImage>::Get;

  /// Get image value at an arbitrary discrete image location
  virtual VoxelType Get(int i, int j, int k = 0, int l = 0) const
  {
    RepeatExtrapolateImageFunction::Apply(l, this->T() - 1);
    return GenericConstExtrapolateImageFunction<TImage>::Get(i, j, k, l);
  }

};


/**
 * Extrapolation of any image by padding it with a constant value
 * in the spatial domain, but a periodic extension in the temporal dimension
 */
class ConstExtrapolateImageFunctionWithPeriodicTime
: public GenericConstExtrapolateImageFunctionWithPeriodicTime<BaseImage>
{
  mirtkObjectMacro(ConstExtrapolateImageFunctionWithPeriodicTime);

public:

  /// Constructor
  ConstExtrapolateImageFunctionWithPeriodicTime() {}

  /// Destructor
  virtual ~ConstExtrapolateImageFunctionWithPeriodicTime() {}

};


} // namespace mirtk

#endif // MIRTK_ConstExtrapolateImageFunctionWithPeriodicTime_H
