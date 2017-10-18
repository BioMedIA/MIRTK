/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2017 Imperial College London
 * Copyright 2008-2013 Daniel Rueckert, Julia Schnabel
 * Copyright 2013-2017 Andreas Schuh
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

#include "mirtk/GaussianBlurringWithPadding.h"


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType> GaussianBlurringWithPadding<VoxelType>
::GaussianBlurringWithPadding(double sigma, VoxelType padding_value)
:
  GaussianBlurringWithPadding<VoxelType>(sigma, sigma, sigma, sigma, padding_value)
{
}

// -----------------------------------------------------------------------------
template <class VoxelType> GaussianBlurringWithPadding<VoxelType>
::GaussianBlurringWithPadding(double xsigma, double ysigma, VoxelType padding_value)
:
  GaussianBlurringWithPadding<VoxelType>(xsigma, ysigma, 0., 0., padding_value)
{
}

// -----------------------------------------------------------------------------
template <class VoxelType> GaussianBlurringWithPadding<VoxelType>
::GaussianBlurringWithPadding(double xsigma, double ysigma, double zsigma, VoxelType padding_value)
:
  GaussianBlurringWithPadding<VoxelType>(xsigma, ysigma, zsigma, 0., padding_value)
{
}

// -----------------------------------------------------------------------------
template <class VoxelType> GaussianBlurringWithPadding<VoxelType>
::GaussianBlurringWithPadding(double xsigma, double ysigma, double zsigma, double tsigma, VoxelType padding_value)
:
  GaussianBlurring<VoxelType>(xsigma, ysigma, zsigma, tsigma)
{
  this->_UseBackgroundMask  = false;
  this->_UseBackgroundValue = false;
  this->_UsePaddingValue    = true;
  this->_PaddingValue       = padding_value;
}

// =============================================================================
// Explicit template instantiations
// =============================================================================

template class GaussianBlurringWithPadding<unsigned char>;
template class GaussianBlurringWithPadding<short>;
template class GaussianBlurringWithPadding<unsigned short>;
template class GaussianBlurringWithPadding<float>;
template class GaussianBlurringWithPadding<double>;


} // namespace mirtk
