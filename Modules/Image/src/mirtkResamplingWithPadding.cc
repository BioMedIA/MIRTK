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

#include <mirtkResamplingWithPadding.h>

#include <mirtkMath.h>
#include <mirtkParallel.h>
#include <mirtkProfiling.h>


namespace mirtk {


// =============================================================================
// Auxiliaries
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
class MultiThreadedResamplingWithPadding
{
  /// Pointer to image transformation class
  ResamplingWithPadding<VoxelType> *_Filter;

  /// Background padding value
  double _PaddingValue;

  /// Time frame to transform
  int _t;

public:

  MultiThreadedResamplingWithPadding(ResamplingWithPadding<VoxelType> *filter, int t)
  :
    _Filter(filter),
    _PaddingValue(filter->PaddingValue()),
    _t(t)
  {}

  void operator ()(const blocked_range<int> &r) const
  {
    int    u, v, w, pad;
    double val, sum, w1, w2, w3, w4, w5, w6, w7, w8, dx, dy, dz, x, y, z;

    const GenericImage<VoxelType> *input  = _Filter->Input();
    GenericImage<VoxelType>       *output = _Filter->Output();

    // TODO: Use InterpolateImageFunction::GetWithPadding instead

    int l = _t;
    for (int k = r.begin(); k != r.end();    ++k)
    for (int j = 0;         j < output->Y(); ++j)
    for (int i = 0;         i < output->X(); ++i) {
      x = i, y = j, z = k;
      output->ImageToWorld(x, y, z);
      input ->WorldToImage(x, y, z);

      // Calculate integer fraction of points
      u = ifloor(x);
      v = ifloor(y);
      w = ifloor(z);

      // Calculate floating point fraction of points
      dx = x - u;
      dy = y - v;
      dz = z - w;

      // Calculate weights for trilinear interpolation
      w1 = (1 - dx) * (1 - dy) * (1 - dz);
      w2 = (1 - dx) * (1 - dy) * dz;
      w3 = (1 - dx) * dy * (1 - dz);
      w4 = (1 - dx) * dy * dz;
      w5 = dx * (1 - dy) * (1 - dz);
      w6 = dx * (1 - dy) * dz;
      w7 = dx * dy * (1 - dz);
      w8 = dx * dy * dz;

      // Calculate trilinear interpolation, ignoring padded values
      val = 0;
      pad = 8;
      sum = 0;
      if ((u >= 0) && (u < input->X()) &&
          (v >= 0) && (v < input->Y()) &&
          (w >= 0) && (w < input->Z())) {
        if (input->Get(u, v, w, l) != _PaddingValue) {
          pad--;
          val += input->Get(u, v, w, l) * w1;
          sum += w1;
        }
      } else {
        pad--;
      }
      if ((u   >= 0) && (u   < input->X()) &&
          (v   >= 0) && (v   < input->Y()) &&
          (w+1 >= 0) && (w+1 < input->Z())) {
        if (input->Get(u, v, w+1, l) != _PaddingValue) {
          pad--;
          val += input->Get(u, v, w+1, l) * w2;
          sum += w2;
        }
      } else {
        pad--;
      }
      if ((u   >= 0) && (u   < input->X()) &&
          (v+1 >= 0) && (v+1 < input->Y()) &&
          (w   >= 0) && (w   < input->Z())) {
        if (input->Get(u, v+1, w, l) != _PaddingValue) {
          pad--;
          val += input->Get(u, v+1, w, l) * w3;
          sum += w3;
        }
      } else {
        pad--;
      }
      if ((u   >= 0) && (u   < input->X()) &&
          (v+1 >= 0) && (v+1 < input->Y()) &&
          (w+1 >= 0) && (w+1 < input->Z())) {
        if (input->Get(u, v+1, w+1, l) != _PaddingValue) {
          pad--;
          val += input->Get(u, v+1, w+1, l) * w4;
          sum += w4;
        }
      } else {
        pad--;
      }
      if ((u+1 >= 0) && (u+1 < input->X()) &&
          (v   >= 0) && (v   < input->Y()) &&
          (w   >= 0) && (w   < input->Z())) {
        if (input->Get(u+1, v, w, l) != _PaddingValue) {
          pad--;
          val += input->Get(u+1, v, w, l) * w5;
          sum += w5;
        }
      } else {
        pad--;
      }
      if ((u+1 >= 0) && (u+1 < input->X()) &&
          (v   >= 0) && (v   < input->Y()) &&
          (w+1 >= 0) && (w+1 < input->Z())) {
        if (input->Get(u+1, v, w+1, l) != _PaddingValue) {
          pad--;
          val += input->Get(u+1, v, w+1, l) * w6;
          sum += w6;
        }
      } else {
        pad--;
      }
      if ((u+1 >= 0) && (u+1 < input->X()) &&
          (v+1 >= 0) && (v+1 < input->Y()) &&
          (w   >= 0) && (w   < input->Z())) {
        if (input->Get(u+1, v+1, w, l) != _PaddingValue) {
          pad--;
          val += input->Get(u+1, v+1, w, l) * w7;
          sum += w7;
        }
      } else {
        pad--;
      }
      if ((u+1 >= 0) && (u+1 < input->X()) &&
          (v+1 >= 0) && (v+1 < input->Y()) &&
          (w+1 >= 0) && (w+1 < input->Z())) {
        if (input->Get(u+1, v+1, w+1, l) != _PaddingValue) {
          pad--;
          val += input->Get(u+1, v+1, w+1, l) * w8;
          sum += w8;
        }
      } else {
        pad--;
      }
      if (pad < 4) {
        if (sum > 0) {
          output->PutAsDouble(i, j, k, l, val / sum);
        } else {
          output->Put(i, j, k, l, _PaddingValue);
        }
      } else {
        output->Put(i, j, k, l, _PaddingValue);
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
  // Not Resampling::Initialize which requires _Interpolator to be set!
  ImageToImage<VoxelType>::Initialize();

  // Initialize output image
  this->InitializeOutput();
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void ResamplingWithPadding<VoxelType>::Run()
{
  MIRTK_START_TIMING();

  this->Initialize();

  for (int l = 0; l < this->_Output->T(); ++l) {
    MultiThreadedResamplingWithPadding<VoxelType> body(this, l);
    parallel_for(blocked_range<int>(0, this->_Output->Z(), 1), body);
  }

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
