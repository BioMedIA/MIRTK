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

#ifndef MIRTK_LinearInterpolateImageFunction2D_HXX
#define MIRTK_LinearInterpolateImageFunction2D_HXX

#include "mirtk/LinearInterpolateImageFunction2D.h"
#include "mirtk/LinearInterpolateImageFunction.hxx"


namespace mirtk {


// -----------------------------------------------------------------------------
template <class TImage>
GenericLinearInterpolateImageFunction2D<TImage>
::GenericLinearInterpolateImageFunction2D()
{
  this->NumberOfDimensions(2);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericLinearInterpolateImageFunction2D<TImage>::VoxelType
GenericLinearInterpolateImageFunction2D<TImage>
::Get(double x, double y, double z, double t) const
{
  return this->Get2D(x, y, z, t);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericLinearInterpolateImageFunction2D<TImage>::VoxelType
GenericLinearInterpolateImageFunction2D<TImage>
::GetWithPadding(double x, double y, double z, double t) const
{
  return this->GetWithPadding2D(x, y, z, t);
}

// -----------------------------------------------------------------------------
template <class TImage> template <class TOtherImage>
inline typename TOtherImage::VoxelType
GenericLinearInterpolateImageFunction2D<TImage>
::Get(const TOtherImage *input, double x, double y, double z, double t) const
{
  return this->Get2D(input, x, y, z, t);
}

// -----------------------------------------------------------------------------
template <class TImage> template <class TOtherImage>
inline typename TOtherImage::VoxelType
GenericLinearInterpolateImageFunction2D<TImage>
::GetWithPadding(const TOtherImage *input, double x, double y, double z, double t) const
{
  return this->GetWithPadding2D(input, x, y, z, t);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericLinearInterpolateImageFunction2D<TImage>::VoxelType
GenericLinearInterpolateImageFunction2D<TImage>
::GetInside(double x, double y, double z, double t) const
{
  // Use faster pixel access than Get2D(Input(), x, y, z, t), which requires
  // a 4D array lookup with three indirections instead of a direct 1D lookup
  return this->GetInside2D(x, y, z, t);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericLinearInterpolateImageFunction2D<TImage>::VoxelType
GenericLinearInterpolateImageFunction2D<TImage>
::GetOutside(double x, double y, double z, double t) const
{
  if (this->Extrapolator()) {
    return Get(this->Extrapolator(), x, y, z, t);
  } else {
    return Get(x, y, z, t);
  }
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericLinearInterpolateImageFunction2D<TImage>::VoxelType
GenericLinearInterpolateImageFunction2D<TImage>
::GetWithPaddingInside(double x, double y, double z, double t) const
{
  return this->GetWithPadding(this->Input(), x, y, z, t);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericLinearInterpolateImageFunction2D<TImage>::VoxelType
GenericLinearInterpolateImageFunction2D<TImage>
::GetWithPaddingOutside(double x, double y, double z, double t) const
{
  if (this->Extrapolator()) {
    return GetWithPadding(this->Extrapolator(), x, y, z, t);
  } else {
    return GetWithPadding(x, y, z, t);
  }
}


} // namespace mirtk

#endif // MIRTK_LinearInterpolateImageFunction2D_HXX
