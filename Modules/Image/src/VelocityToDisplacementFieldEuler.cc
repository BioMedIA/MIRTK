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

#include "mirtk/VelocityToDisplacementFieldEuler.h"

#include "mirtk/Parallel.h"
#include "mirtk/NaryVoxelFunction.h"


namespace mirtk {


// ===========================================================================
// Construction/Destruction
// ===========================================================================

// ---------------------------------------------------------------------------
template <class VoxelType>
VelocityToDisplacementFieldEuler<VoxelType>::VelocityToDisplacementFieldEuler()
:
  VelocityToDisplacementField<VoxelType>(),
  _VelocityInterpolator(NULL)
{
}

// ---------------------------------------------------------------------------
template <class VoxelType>
VelocityToDisplacementFieldEuler<VoxelType>::~VelocityToDisplacementFieldEuler()
{
  delete _VelocityInterpolator;
}

// ===========================================================================
// Filter implementation
// ===========================================================================

// ---------------------------------------------------------------------------
template <class VoxelType>
void VelocityToDisplacementFieldEuler<VoxelType>::Initialize()
{
  // Initialize base class
  VelocityToDisplacementField<VoxelType>::Initialize();

  // Initialize interpolator
  if (_VelocityInterpolator) delete _VelocityInterpolator;
  _VelocityInterpolator = InterpolateImageFunction::New(this->Interpolation(),
                                                        this->Extrapolation(),
                                                        this->Input());
  _VelocityInterpolator->Input(this->Input());
  _VelocityInterpolator->Initialize();
}

// ---------------------------------------------------------------------------
template <class VoxelType>
void VelocityToDisplacementFieldEuler<VoxelType>::Run()
{
  this->Initialize();

  using NaryVoxelFunction::ExpVelocityFieldEuler2D;
  using NaryVoxelFunction::ExpVelocityFieldEuler3D;
  const ImageAttributes &grid = this->Output()->Attributes();

  if (this->Output()->Z() > 1) {
    ExpVelocityFieldEuler3D<> exp(_VelocityInterpolator, this->NumberOfSteps(), this->UpperIntegrationLimit());
    if (this->Input(1)) ParallelForEachVoxel(grid, this->Input(1), this->Output(), exp);
    else                ParallelForEachVoxel(grid,                 this->Output(), exp);
  } else {
    ExpVelocityFieldEuler2D<> exp(_VelocityInterpolator, this->NumberOfSteps(), this->UpperIntegrationLimit());
    if (this->Input(1)) ParallelForEachVoxel(grid, this->Input(1), this->Output(), exp);
    else                ParallelForEachVoxel(grid,                 this->Output(), exp);
  }

  this->Finalize();
}

// ===========================================================================
// Explicit template instantiations
// ===========================================================================

template class VelocityToDisplacementFieldEuler<float>;
template class VelocityToDisplacementFieldEuler<double>;


} // namespace mirtk
