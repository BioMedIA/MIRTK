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

#include "mirtk/GaussianBlurringWithPadding2D.h"


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
GaussianBlurringWithPadding2D<VoxelType>
::GaussianBlurringWithPadding2D(double sigma, VoxelType padding_value)
:
  GaussianBlurringWithPadding<VoxelType>(sigma, sigma, .0, .0, padding_value)
{
}

// -----------------------------------------------------------------------------
template <class VoxelType>
GaussianBlurringWithPadding2D<VoxelType>
::GaussianBlurringWithPadding2D(double xsigma, double ysigma, VoxelType padding_value)
:
  GaussianBlurringWithPadding<VoxelType>(xsigma, ysigma, .0, .0, padding_value)
{
}

// =============================================================================
// Explicit template instantiations
// =============================================================================

template class GaussianBlurringWithPadding2D<unsigned char>;
template class GaussianBlurringWithPadding2D<short>;
template class GaussianBlurringWithPadding2D<unsigned short>;
template class GaussianBlurringWithPadding2D<float>;
template class GaussianBlurringWithPadding2D<double>;


} // namespace mirtk
