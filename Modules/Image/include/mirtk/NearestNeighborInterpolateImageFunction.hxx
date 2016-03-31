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

#ifndef MIRTK_NearestNeighorInterpolateImageFunction_HXX
#define MIRTK_NearestNeighorInterpolateImageFunction_HXX

#include "mirtk/NearestNeighborInterpolateImageFunction.h"

#include "mirtk/InterpolateImageFunction.hxx"
#include "mirtk/Math.h"


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
template <class TImage>
GenericNearestNeighborInterpolateImageFunction<TImage>
::GenericNearestNeighborInterpolateImageFunction()
{
}

// -----------------------------------------------------------------------------
template <class TImage>
GenericNearestNeighborInterpolateImageFunction<TImage>
::~GenericNearestNeighborInterpolateImageFunction()
{
}

// =============================================================================
// Domain checks
// =============================================================================

// -----------------------------------------------------------------------------
template <class TImage>
void GenericNearestNeighborInterpolateImageFunction<TImage>
::BoundingInterval(double x, int &i, int &I) const
{
  i = I = iround(x);
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericNearestNeighborInterpolateImageFunction<TImage>::VoxelType
GenericNearestNeighborInterpolateImageFunction<TImage>
::Get(double x, double y, double z, double t) const
{
  const int i = iround(x);
  const int j = iround(y);
  const int k = iround(z);
  const int l = iround(t);

  if (this->Input()->IsInside(i, j, k, l)) {
    return this->Input()->Get(i, j, k, l);
  } else {
    return voxel_cast<VoxelType>(this->DefaultValue());
  }
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericNearestNeighborInterpolateImageFunction<TImage>::VoxelType
GenericNearestNeighborInterpolateImageFunction<TImage>
::GetWithPadding(double x, double y, double z, double t) const
{
  const int i = iround(x);
  const int j = iround(y);
  const int k = iround(z);
  const int l = iround(t);

  if (this->Input()->IsInsideForeground(i, j, k, l)) {
    return this->Input()->Get(i, j, k, l);
  } else {
    return voxel_cast<VoxelType>(this->DefaultValue());
  }
}

// -----------------------------------------------------------------------------
template <class TImage> template <class TOtherImage>
inline typename TOtherImage::VoxelType
GenericNearestNeighborInterpolateImageFunction<TImage>
::Get(const TOtherImage *input, double x, double y, double z, double t) const
{
  return input->Get(iround(x), iround(y), iround(z), iround(t));
}

// -----------------------------------------------------------------------------
template <class TImage> template <class TOtherImage>
inline typename TOtherImage::VoxelType
GenericNearestNeighborInterpolateImageFunction<TImage>
::GetWithPadding(const TOtherImage *input, double x, double y, double z, double t) const
{
  const int i = iround(x);
  const int j = iround(y);
  const int k = iround(z);
  const int l = iround(t);

  if (input->IsForeground(i, j, k, l)) {
    return input->Get(i, j, k, l);
  } else {
    return voxel_cast<VoxelType>(this->DefaultValue());
  }
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericNearestNeighborInterpolateImageFunction<TImage>::VoxelType
GenericNearestNeighborInterpolateImageFunction<TImage>
::GetInside(double x, double y, double z, double t) const
{
  return Get(this->Input(), x, y, z, t);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericNearestNeighborInterpolateImageFunction<TImage>::VoxelType
GenericNearestNeighborInterpolateImageFunction<TImage>
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
inline typename GenericNearestNeighborInterpolateImageFunction<TImage>::VoxelType
GenericNearestNeighborInterpolateImageFunction<TImage>
::GetWithPaddingInside(double x, double y, double z, double t) const
{
  return GetWithPadding(this->Input(), x, y, z, t);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericNearestNeighborInterpolateImageFunction<TImage>::VoxelType
GenericNearestNeighborInterpolateImageFunction<TImage>
::GetWithPaddingOutside(double x, double y, double z, double t) const
{
  if (this->Extrapolator()) {
    return GetWithPadding(this->Extrapolator(), x, y, z, t);
  } else {
    return GetWithPadding(x, y, z, t);
  }
}


} // namespace mirtk

#endif // MIRTK_NearestNeighorInterpolateImageFunction_HXX
