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

#include "mirtk/VelocityToDisplacementField.h"


namespace mirtk {


// ===========================================================================
// Construction/Destruction
// ===========================================================================

// ---------------------------------------------------------------------------
template <class VoxelType>
VelocityToDisplacementField<VoxelType>::VelocityToDisplacementField()
:
  _Interpolation                   (Interpolation_Linear),
  _Extrapolation                   (Extrapolation_NN),
  _ComputeInterpolationCoefficients(true),
  _ComputeInverse                  (false),
  _NumberOfSteps                   (64), // i.e., 6 squaring steps in case of SS method
  _UpperIntegrationLimit           (1.0),
  _InputDisplacementField          (NULL)
{
}

// ---------------------------------------------------------------------------
template <class VoxelType>
VelocityToDisplacementField<VoxelType>::~VelocityToDisplacementField()
{
}

// ===========================================================================
// Filter implementation
// ===========================================================================

// --------------------------------------------------------------------------
template <class VoxelType>
void VelocityToDisplacementField<VoxelType>::Input(int i, const ImageType *image)
{
  if      (i == 0) Baseclass::Input(image);
  else if (i == 1) _InputDisplacementField = image;
  else {
    cerr << this->NameOfClass() << "::Input: Input index out of range: " << i << endl;
    exit(1);
  }
}

// --------------------------------------------------------------------------
template <class VoxelType>
const typename VelocityToDisplacementField<VoxelType>::ImageType *
VelocityToDisplacementField<VoxelType>::Input(int i)
{
  if      (i == 0) return this->Input();
  else if (i == 1) return _InputDisplacementField;
  else {
    cerr << this->NameOfClass() << "::Input: Input index out of range: " << i << endl;
    exit(1);
  }
}

// ---------------------------------------------------------------------------
template <class VoxelType>
void VelocityToDisplacementField<VoxelType>::Initialize()
{
  // Initialize base class
  Baseclass::Initialize();

  // Check input
  if (this->Input()->Z() < 1) {
    // not allowed as otherwise a for loop over z from 0 to _z-1 would never
    // be executed, not even once as it should be in case of 2 dimensions
    cerr << this->NameOfClass() << "::Initialize: Size of z dimension must be 1 in case of 2D vector field" << endl;
    exit(1);
  }

  if ((this->Input()->Z() <= 1 && this->Input()->T() != 2 && this->Input()->T() != 3) ||
      (this->Input()->Z() >  1 && this->Input()->T() != 3)) {
    cerr << this->NameOfClass() << "::Initialize: Input must be a 2D or 3D vector field" << endl;
    exit(1);
  }

  // Check parameters
  if (_NumberOfSteps < 1) {
    cerr << this->NameOfClass() << "::Initialize: Number of integration steps must be positive" << endl;
    exit(1);
  }
  if (_UpperIntegrationLimit == 0) {
    cerr << this->NameOfClass() << "::Initialize: Upper integration limit (T) must not be zero" << endl;
    exit(1);
  }
  if (_ComputeInverse) _UpperIntegrationLimit = -_UpperIntegrationLimit; // Inverse <=> backward integration

  // Output is vector field
  this->Output()->PutTSize(.0);
}

// ---------------------------------------------------------------------------
template <class VoxelType>
void VelocityToDisplacementField<VoxelType>::Finalize()
{
  // Reset upper integration limit attribute
  if (_ComputeInverse) _UpperIntegrationLimit = -_UpperIntegrationLimit;
  // Finalize base class
  Baseclass::Finalize();
}

// ===========================================================================
// Explicit template instantiations
// ===========================================================================

template class VelocityToDisplacementField<float>;
template class VelocityToDisplacementField<double>;


} // namespace mirtk
