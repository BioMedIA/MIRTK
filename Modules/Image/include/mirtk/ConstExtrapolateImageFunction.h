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

#ifndef MIRTK_ConstExtrapolateImageFunction_H
#define MIRTK_ConstExtrapolateImageFunction_H

#include "mirtk/ExtrapolateImageFunction.h"

#include "mirtk/VoxelCast.h"
#include "mirtk/BaseImage.h"


namespace mirtk {


/**
 * Constant extrapolation (i.e., padding) of generic image
 *
 * This discrete image extrapolation function assumes a constant image value
 * for image grid points outside the discrete domain on which the discrete image
 * function is defined.
 */
template <class TImage>
class GenericConstExtrapolateImageFunction
: public GenericExtrapolateImageFunction<TImage>
{
  mirtkGenericExtrapolatorMacro(
    GenericConstExtrapolateImageFunction,
    Extrapolation_Const
  );

public:

  /// Default constructor
  GenericConstExtrapolateImageFunction(double padding_value = .0)
  {
    this->_DefaultValue = padding_value;
  }

  /// Destructor
  virtual ~GenericConstExtrapolateImageFunction() {}

  /// Get image value at an arbitrary discrete image location
  virtual VoxelType Get(int i, int j, int k = 0, int l = 0) const
  {
    if (0 <= i && i < this->X() &&
        0 <= j && j < this->Y() &&
        0 <= k && k < this->Z() &&
        0 <= l && l < this->T()) {
      return this->Input()->Get(i, j, k, l);
    } else {
      return voxel_cast<VoxelType>(this->_DefaultValue);
    }
  }

};


/**
 * Constant extrapolation (i.e., padding) of any image
 */
class ConstExtrapolateImageFunction
: public GenericConstExtrapolateImageFunction<BaseImage>
{
  mirtkObjectMacro(ConstExtrapolateImageFunction);

public:

  /// Constructor
  ConstExtrapolateImageFunction(double padding_value = .0)
  :
    GenericConstExtrapolateImageFunction<BaseImage>(padding_value)
  {}

  /// Destructor
  virtual ~ConstExtrapolateImageFunction() {}

};


} // namespace mirtk

#endif // MIRTK_ConstExtrapolateImageFunction_H
