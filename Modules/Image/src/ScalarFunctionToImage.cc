/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2015 Imperial College London
 * Copyright 2008-2015 Daniel Rueckert, Julia Schnabel
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

#include "mirtk/ScalarFunctionToImage.h"

#include "mirtk/Math.h"


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
ScalarFunctionToImage<VoxelType>::ScalarFunctionToImage(bool use_world_coords)
:
  _Input(NULL),
  _Output(NULL),
  _UseWorldCoordinates(use_world_coords)
{
}

// -----------------------------------------------------------------------------
template <class VoxelType>
ScalarFunctionToImage<VoxelType>::~ScalarFunctionToImage()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
void ScalarFunctionToImage<VoxelType>::Run()
{
  const double eps = static_cast<double>(numeric_limits<float>::min());

  double x, y, z;
  for (int l = 0; l < _Output->T(); ++l)
  for (int k = 0; k < _Output->Z(); ++k)
  for (int j = 0; j < _Output->Y(); ++j)
  for (int i = 0; i < _Output->X(); ++i) {
    x = i, y = j, z = k;
    if (_UseWorldCoordinates) _Output->ImageToWorld(x, y, z);
    _Output->PutAsDouble(i, j, k, l, _Input->Evaluate(x, y, z));
    if (abs(_Output->GetAsDouble(i, j, k, l)) < eps) {
      _Output->Put(i, j, k, l, 0);
    }
  }
}

// =============================================================================
// Explicit template instantiations
// =============================================================================

template class ScalarFunctionToImage<BytePixel>;
template class ScalarFunctionToImage<GreyPixel>;
template class ScalarFunctionToImage<RealPixel>;


} // namespace mirtk
