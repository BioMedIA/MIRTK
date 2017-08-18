/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2017 Imperial College London
 * Copyright 2017 Andreas Schuh
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

#include "mirtk/SeparableConvolution.h"

#include "mirtk/ConvolutionFunction.h"
#include "mirtk/Profiling.h"


namespace mirtk {


// -----------------------------------------------------------------------------
template <class TVoxel, class TKernel>
SeparableConvolution<TVoxel, TKernel>::SeparableConvolution(const KernelType *k)
:
  SeparableConvolution<TVoxel, TKernel>(k, k, k, k)
{
}

// -----------------------------------------------------------------------------
template <class TVoxel, class TKernel>
SeparableConvolution<TVoxel, TKernel>
::SeparableConvolution(const KernelType *kx, const KernelType *ky, const KernelType *kz, const KernelType *kt)
:
  _KernelX(kx),
  _KernelY(ky),
  _KernelZ(kz),
  _KernelT(kt),
  _Normalize(true),
  _UseBackgroundMask(false),
  _UseBackgroundValue(false),
  _UsePaddingValue(false),
  _PaddingValue(NaN),
  _Components(0)
{
}

// -----------------------------------------------------------------------------
template <class TVoxel, class TKernel>
SeparableConvolution<TVoxel, TKernel>::~SeparableConvolution()
{
}

// -----------------------------------------------------------------------------
template <class TVoxel, class TKernel>
void SeparableConvolution<TVoxel, TKernel>::Kernel(const KernelType *k)
{
  _KernelX = _KernelY = _KernelZ = _KernelT = k;
}

// -----------------------------------------------------------------------------
template <class TVoxel, class TKernel>
bool SeparableConvolution<TVoxel, TKernel>::CheckKernel(const KernelType *k) const
{
  if (k == nullptr) return false;
  if (k->NumberOfVoxels() == 0) return false;
  return (k->Y() == 1 || k->Z() == 1 || k->T() == 1);
}

// -----------------------------------------------------------------------------
template <class TVoxel, class TKernel>
void SeparableConvolution<TVoxel, TKernel>::Initialize()
{
  // Initialize base class
  ImageToImage<TVoxel>::Initialize(false);

  // Allocate output image
  if (this->Input() != this->Output()) {
    this->Output()->Initialize(this->Input()->Attributes(), _Components);
  }

  // Check arguments
  if (_KernelX && !CheckKernel(_KernelX)) {
    Throw(ERR_InvalidArgument, __FUNCTION__, "Convolution kernel along x axis is invalid");
  }
  if (_KernelY && !CheckKernel(_KernelY)) {
    Throw(ERR_InvalidArgument, __FUNCTION__, "Convolution kernel along y axis is invalid");
  }
  if (_KernelZ && !CheckKernel(_KernelZ)) {
    Throw(ERR_InvalidArgument, __FUNCTION__, "Convolution kernel along z axis is invalid");
  }
  if (_KernelT && !CheckKernel(_KernelT)) {
    Throw(ERR_InvalidArgument, __FUNCTION__, "Convolution kernel along t axis is invalid");
  }
}

// -----------------------------------------------------------------------------
template <class TVoxel, class TKernel>
void SeparableConvolution<TVoxel, TKernel>::Run()
{
  MIRTK_START_TIMING();

  // Do the initial set up
  this->Initialize();

  ImageType *input  = const_cast<ImageType *>(this->Input());
  ImageType *output = this->Output();

  ImageAttributes attr = input->Attributes();

  // Number of vector components
  int N = 1;
  if (_Components > 0) {
    N = _Components;
    attr._dt = 0.;
  } else if (attr._t == 1) {
    if (attr._z > 1 && AreEqual(attr._dz, 0.)) {
      N = attr._z;
    }
  } else if (attr._t > 1 && AreEqual(attr._dt, 0.)) {
    N = attr._t;
  }

  int    padmode = 0;
  double padding = _PaddingValue;
  if (_UseBackgroundMask && input->HasMask()) {
    padmode = 1;
  } else if (_UseBackgroundValue && input->HasBackgroundValue()) {
    padding = input->GetBackgroundValueAsDouble();
    padmode = 2;
  } else if (_UsePaddingValue) {
    padmode = 2;
  }

  // Convolve along x axis
  if (_KernelX && input->X() > 1) {
    if (output == this->Input()) {
      output = new ImageType(attr);
    }
    for (int n = 0; n < N; ++n) {
      switch (padmode) {
        case 1: {
          typedef ConvolutionFunction::ConvolveForegroundInX<TKernel> ConvFunc;
          ConvFunc conv(input, _KernelX->Data(), _KernelX->X(), _Normalize, n);
          ParallelForEachVoxel(attr, input, output, conv);
        } break;
        case 2: {
          typedef ConvolutionFunction::ConvolveTruncatedForegroundInX<TKernel> ConvFunc;
          ConvFunc conv(input, padding, _KernelX->Data(), _KernelX->X(), _Normalize, n);
          ParallelForEachVoxel(attr, input, output, conv);
        } break;
        default: {
          typedef ConvolutionFunction::ConvolveInX<TKernel> ConvFunc;
          ConvFunc conv(input, _KernelX->Data(), _KernelX->X(), _Normalize, n);
          ParallelForEachVoxel(attr, input, output, conv);
        } break;
      }
    }
    swap(input, output);
  }

  // Convolve along y axis
  if (_KernelY && input->Y() > 1) {
    if (output == this->Input()) {
      output = new ImageType(attr);
    }
    for (int n = 0; n < N; ++n) {
      switch (padmode) {
        case 1: {
          typedef ConvolutionFunction::ConvolveForegroundInY<TKernel> ConvFunc;
          ConvFunc conv(input, _KernelY->Data(), _KernelY->X(), _Normalize, n);
          ParallelForEachVoxel(attr, input, output, conv);
        } break;
        case 2: {
          typedef ConvolutionFunction::ConvolveTruncatedForegroundInY<TKernel> ConvFunc;
          ConvFunc conv(input, padding, _KernelY->Data(), _KernelY->X(), _Normalize, n);
          ParallelForEachVoxel(attr, input, output, conv);
        } break;
        default: {
          typedef ConvolutionFunction::ConvolveInY<TKernel> ConvFunc;
          ConvFunc conv(input, _KernelY->Data(), _KernelY->X(), _Normalize, n);
          ParallelForEachVoxel(attr, input, output, conv);
        } break;
      }
    }
    swap(input, output);
  }

  // Convolve along z axis
  if (_KernelZ && input->Z() > 1 && !AreEqual(input->ZSize(), 0.)) {
    if (output == this->Input()) {
      output = new ImageType(attr);
    }
    for (int n = 0; n < N; ++n) {
      switch (padmode) {
        case 1: {
          typedef ConvolutionFunction::ConvolveForegroundInZ<TKernel> ConvFunc;
          ConvFunc conv(input, _KernelZ->Data(), _KernelZ->X(), _Normalize, n);
          ParallelForEachVoxel(attr, input, output, conv);
        } break;
        case 2: {
          typedef ConvolutionFunction::ConvolveTruncatedForegroundInZ<TKernel> ConvFunc;
          ConvFunc conv(input, padding, _KernelZ->Data(), _KernelZ->X(), _Normalize, n);
          ParallelForEachVoxel(attr, input, output, conv);
        } break;
        default: {
          typedef ConvolutionFunction::ConvolveInZ<TKernel> ConvFunc;
          ConvFunc conv(input, _KernelZ->Data(), _KernelZ->X(), _Normalize, n);
          ParallelForEachVoxel(attr, input, output, conv);
        } break;
      }
    }
    swap(input, output);
  }

  // Convolve along t axis
  if (_KernelT && input->T() > 1 && !AreEqual(input->TSize(), 0.)) {
    if (output == this->Input()) {
      output = new ImageType(attr);
    }
    switch (padmode) {
      case 1: {
        typedef ConvolutionFunction::ConvolveForegroundInT<TKernel> ConvFunc;
        ConvFunc conv(input, _KernelT->Data(), _KernelT->X(), _Normalize);
        ParallelForEachVoxel(attr, input, output, conv);
      } break;
      case 2: {
        typedef ConvolutionFunction::ConvolveTruncatedForegroundInT<TKernel> ConvFunc;
        ConvFunc conv(input, padding, _KernelT->Data(), _KernelT->X(), _Normalize);
        ParallelForEachVoxel(attr, input, output, conv);
      } break;
      default: {
        typedef ConvolutionFunction::ConvolveInT<TKernel> ConvFunc;
        ConvFunc conv(input, _KernelT->Data(), _KernelT->X(), _Normalize);
        ParallelForEachVoxel(attr, input, output, conv);
      } break;
    }
    swap(input, output);
  }

  // Copy result if last output (i.e., after swap input pointer) was not filter
  // output image and delete possibly additionally allocated image
  if (input != output) {
    if (input == this->Output()) {
      if (output != this->Input()) {
        Delete(output);
      }
    } else {
      this->Output()->CopyFrom(input->Data());
      if (input != this->Input()) {
        Delete(input);
      }
    }
  }

  // Do the final cleaning up
  this->Finalize();

  MIRTK_DEBUG_TIMING(5, this->NameOfClass());
}

// -----------------------------------------------------------------------------
template <class TVoxel, class TKernel>
void SeparableConvolution<TVoxel, TKernel>::Finalize()
{
  // Finalize base class
  ImageToImage<TVoxel>::Finalize();

  // Copy background information
  if (this->Input() != this->Output()) {
    if (this->Input()->HasBackgroundValue()) {
      this->Output()->PutBackgroundValueAsDouble(this->Input()->GetBackgroundValueAsDouble());
    }
    if (this->Input()->HasMask()) {
      this->Output()->PutMask(new BinaryImage(*this->Input()->GetMask()), true);
    }
  }
}

// -----------------------------------------------------------------------------
template <class TVoxel, class TKernel>
void SeparableConvolution<TVoxel, TKernel>::RunX()
{
  const KernelType * const ky = _KernelY;
  const KernelType * const kz = _KernelZ;
  const KernelType * const kt = _KernelT;
  _KernelY = _KernelZ = _KernelT = nullptr;

  this->Run();

  _KernelY = ky;
  _KernelZ = kz;
  _KernelT = kt;
}

// -----------------------------------------------------------------------------
template <class TVoxel, class TKernel>
void SeparableConvolution<TVoxel, TKernel>::RunY()
{
  const KernelType * const kx = _KernelX;
  const KernelType * const kz = _KernelZ;
  const KernelType * const kt = _KernelT;
  _KernelX = _KernelZ = _KernelT = nullptr;

  this->Run();

  _KernelX = kx;
  _KernelZ = kz;
  _KernelT = kt;
}

// -----------------------------------------------------------------------------
template <class TVoxel, class TKernel>
void SeparableConvolution<TVoxel, TKernel>::RunZ()
{
  const KernelType * const kx = _KernelX;
  const KernelType * const ky = _KernelY;
  const KernelType * const kt = _KernelT;
  _KernelX = _KernelY = _KernelT = nullptr;

  this->Run();

  _KernelX = kx;
  _KernelY = ky;
  _KernelT = kt;
}

// -----------------------------------------------------------------------------
template <class TVoxel, class TKernel>
void SeparableConvolution<TVoxel, TKernel>::RunT()
{
  const KernelType * const kx = _KernelX;
  const KernelType * const ky = _KernelY;
  const KernelType * const kz = _KernelZ;
  _KernelX = _KernelY = _KernelZ = nullptr;

  this->Run();

  _KernelX = kx;
  _KernelY = ky;
  _KernelZ = kz;
}

// =============================================================================
// Explicit template instantiations
// =============================================================================

template class SeparableConvolution<unsigned char, float>;
template class SeparableConvolution<short, float>;
template class SeparableConvolution<unsigned short, float>;
template class SeparableConvolution<float, float>;
template class SeparableConvolution<double, float>;

template class SeparableConvolution<unsigned char, double>;
template class SeparableConvolution<short, double>;
template class SeparableConvolution<unsigned short, double>;
template class SeparableConvolution<float, double>;
template class SeparableConvolution<double, double>;


} // namespace mirtk
