/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2016 Imperial College London
 * Copyright 2016 Andreas Schuh
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

#include "mirtk/Closing.h"
#include "mirtk/Dilation.h"
#include "mirtk/Erosion.h"


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
Closing<VoxelType>::Closing()
:
  _Connectivity(ConnectivityType::CONNECTIVITY_26),
  _NumberOfIterations(1)
{
}

// -----------------------------------------------------------------------------
template <class VoxelType>
Closing<VoxelType>::~Closing()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
void Closing<VoxelType>::Initialize()
{
  ImageToImage<VoxelType>::Initialize(false);
  *this->_Output = *this->_Input;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void Closing<VoxelType>::Run()
{
  this->Initialize();

  ImageType * const output = this->_Output;

  const double bg_val = output->GetBackgroundValueAsDouble();
  const bool   bg_set = output->HasBackgroundValue();

  int i1, j1, k1, i2, j2, k2;
  output->PutBackgroundValueAsDouble(.0);
  output->BoundingBox(i1, j1, k1, i2, j2, k2);

  if (bg_set) {
    output->PutBackgroundValueAsDouble(bg_val);
  } else {
    output->ClearBackgroundValue();
  }

  i1 -= _NumberOfIterations;
  j1 -= _NumberOfIterations;
  k1 -= _NumberOfIterations;
  i2 += _NumberOfIterations;
  j2 += _NumberOfIterations;
  k2 += _NumberOfIterations;

  if (i1 < 0 || i2 >= output->X() ||
      j1 < 0 || j2 >= output->Y() ||
      k1 < 0 || k2 >= output->Z()) {

    GenericImage<VoxelType> padded(i2 - i1 + 1, j2 - j1 + 1, k2 - k1 + 1);

    i1 += _NumberOfIterations;
    j1 += _NumberOfIterations;
    k1 += _NumberOfIterations;
    i2 -= _NumberOfIterations;
    j2 -= _NumberOfIterations;
    k2 -= _NumberOfIterations;

    for (int k = k1; k <= k2; ++k)
    for (int j = j1; j <= j2; ++j)
    for (int i = i1; i <= i2; ++i) {
      padded(i - i1 + _NumberOfIterations,
             j - j1 + _NumberOfIterations,
             k - k1 + _NumberOfIterations)
          = static_cast<VoxelType>(output->GetAsDouble(i, j, k));
    }

    Dilate<VoxelType>(&padded, _NumberOfIterations, _Connectivity);
    Erode <VoxelType>(&padded, _NumberOfIterations, _Connectivity);

    for (int k = k1; k <= k2; ++k)
    for (int j = j1; j <= j2; ++j)
    for (int i = i1; i <= i2; ++i) {
      output->PutAsDouble(i, j, k, static_cast<double>(padded(i - i1 + _NumberOfIterations,
                                                              j - j1 + _NumberOfIterations,
                                                              k - k1 + _NumberOfIterations)));
    }

  } else {

    Dilate<VoxelType>(output, _NumberOfIterations, _Connectivity);
    Erode <VoxelType>(output, _NumberOfIterations, _Connectivity);

  }

  this->Finalize();
}

// =============================================================================
// Explicit template instantiations
// =============================================================================

template class Closing<BytePixel>;
template class Closing<GreyPixel>;
template class Closing<RealPixel>;


} // namespace mirtk
