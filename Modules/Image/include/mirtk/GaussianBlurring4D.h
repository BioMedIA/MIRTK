/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2015 Imperial College London
 * Copyright 2008-2013 Daniel Rueckert, Julia Schnabel
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

#ifndef MIRTK_GaussianBlurring4D_H
#define MIRTK_GaussianBlurring4D_H

#include "mirtk/GaussianBlurring.h"


namespace mirtk {


/**
 * Class for Gaussian blurring of image sequences
 *
 * This class defines and implements the Gaussian blurring of image sequences.
 * The blurring is implemented by four successive 1D convolutions with a 1D
 * Gaussian kernel.
 */
template <class TVoxel>
class GaussianBlurring4D : public GaussianBlurring<TVoxel>
{
  mirtkInPlaceImageFilterMacro(GaussianBlurring4D, TVoxel);

public:

  /// Constructor
  GaussianBlurring4D(double = 1.0);

  /// Constructor
  GaussianBlurring4D(double, double, double, double);

  /// Destructor
  ~GaussianBlurring4D();

};


} // namespace mirtk

#endif // MIRTK_GaussianBlurring4D_H
