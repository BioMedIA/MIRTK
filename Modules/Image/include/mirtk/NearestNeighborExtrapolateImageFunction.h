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

#ifndef MIRTK_NearestNeighorExtrapolateImageFunction_H
#define MIRTK_NearestNeighorExtrapolateImageFunction_H

#include "mirtk/ExtrapolateImageFunction.h"
#include "mirtk/BaseImage.h"


namespace mirtk {


/**
 * Nearest neighbor extrapolation of generic image
 *
 * The nearest neighor extrapolation of a discrete image corresponds to a
 * Neumann boundary condition with the derivative value outside the image domain
 * set to zero.
 */
template <class TImage>
class GenericNearestNeighborExtrapolateImageFunction
: public IndexExtrapolateImageFunction<TImage>
{
  mirtkExtrapolatorMacro(
    GenericNearestNeighborExtrapolateImageFunction,
    Extrapolation_NN
  );

public:

  /// Constructor
  GenericNearestNeighborExtrapolateImageFunction() {}

  /// Destructor
  virtual ~GenericNearestNeighborExtrapolateImageFunction() {}

  /// Transform index such that it is inside the range [0, max]
  virtual void TransformIndex(int &index, int max) const
  {
    if      (index < 0  ) index = 0;
    else if (index > max) index = max;
  }

};


/**
 * Nearest neighbor extrapolation of any image
 *
 * The nearest neighor extrapolation of a discrete image corresponds to a
 * Neumann boundary condition with the derivative value outside the image domain
 * set to zero.
 */
class NearestNeighborExtrapolateImageFunction
: public GenericNearestNeighborExtrapolateImageFunction<BaseImage>
{
  mirtkObjectMacro(NearestNeighborExtrapolateImageFunction);

public:

  /// Constructor
  NearestNeighborExtrapolateImageFunction() {}

  /// Destructor
  virtual ~NearestNeighborExtrapolateImageFunction() {}

};


} // namespace mirtk

#endif // MIRTK_NearestNeighorExtrapolateImageFunction_H
