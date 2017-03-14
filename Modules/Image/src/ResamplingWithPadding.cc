/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2017 Imperial College London
 * Copyright 2008-2015 Daniel Rueckert, Julia Schnabel
 * Copyright 2017      Andreas Schuh
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

#include "mirtk/ResamplingWithPadding.h"

#include "mirtk/Math.h"
#include "mirtk/Parallel.h"
#include "mirtk/Profiling.h"
#include "mirtk/InterpolateImageFunction.h"


namespace mirtk {


// =============================================================================
// Auxiliaries
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
class MultiThreadedResamplingWithPadding
{
  /// Input image interpolator
  const InterpolateImageFunction *_Input;

  /// Output image
  BaseImage *_Output;

  /// Time frame to resample
  double _t;

public:

  MultiThreadedResamplingWithPadding(const InterpolateImageFunction *input, BaseImage *output, int t)
  :
    _Input(input),
    _Output(output),
    _t(t)
  {}

  void operator ()(const blocked_range<int> &r) const
  {
    double x, y, z;
    for (int k = r.begin(); k != r.end(); ++k) {
      for (int j = 0; j < _Output->Y(); ++j)
      for (int i = 0; i < _Output->X(); ++i) {
        x = i, y = j, z = k;
        _Output->ImageToWorld(x, y, z);
        _Input ->WorldToImage(x, y, z);
        _Output->PutAsDouble(i, j, k, _Input->EvaluateWithPadding(x, y, z, _t));
      }
    }
  }
};

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
ResamplingWithPadding<VoxelType>
::ResamplingWithPadding(double dx, double dy, double dz, VoxelType padding)
:
  Resampling<VoxelType>(dx, dy, dz),
  _PaddingValue(padding)
{
}

// -----------------------------------------------------------------------------
template <class VoxelType>
ResamplingWithPadding<VoxelType>
::ResamplingWithPadding(int x, int y, int z, VoxelType padding)
:
  Resampling<VoxelType>(x, y, z),
  _PaddingValue(padding)
{
}

// -----------------------------------------------------------------------------
template <class VoxelType>
ResamplingWithPadding<VoxelType>
::ResamplingWithPadding(int x, int y, int z, double dx, double dy, double dz, VoxelType padding)
:
  Resampling<VoxelType>(x, y, z, dx, dy, dz),
  _PaddingValue(padding)
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
void ResamplingWithPadding<VoxelType>::Initialize()
{
  // Initialize base class
  Resampling<VoxelType>::Initialize();

  // Initialize output image
  this->InitializeOutput();
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void ResamplingWithPadding<VoxelType>::Run()
{
  MIRTK_START_TIMING();

  this->Initialize();

  ImageType input(this->_Input->Attributes(), const_cast<VoxelType *>(this->_Input->Data()));
  input.PutBackgroundValueAsDouble(static_cast<double>(_PaddingValue), true);
  this->_Interpolator->Input(&input);

  for (int l = 0; l < this->_Output->T(); ++l) {
    MultiThreadedResamplingWithPadding<VoxelType> body(this->_Interpolator, this->_Output, l);
    parallel_for(blocked_range<int>(0, this->_Output->Z(), 1), body);
  }

  this->_Interpolator->Input(this->_Input);

  this->Finalize();

  MIRTK_DEBUG_TIMING(5, this->NameOfClass());
}

// =============================================================================
// Explicit template instantiations
// =============================================================================

template class ResamplingWithPadding<char>;
template class ResamplingWithPadding<unsigned char>;
template class ResamplingWithPadding<short>;
template class ResamplingWithPadding<unsigned short>;
template class ResamplingWithPadding<int>;
template class ResamplingWithPadding<float>;
template class ResamplingWithPadding<double>;


} // namespace mirtk
