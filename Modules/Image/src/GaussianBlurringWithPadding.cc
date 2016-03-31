/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2015 Imperial College London
 * Copyright 2008-2013 Daniel Rueckert, Julia Schnabel
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

#include <mirtkGaussianBlurringWithPadding.h>

#include <mirtkMemory.h>
#include <mirtkGenericImage.h>
#include <mirtkProfiling.h>
#include <mirtkParallel.h>
#include <mirtkConvolutionFunction.h>


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType> GaussianBlurringWithPadding<VoxelType>
::GaussianBlurringWithPadding(double sigma, VoxelType padding_value)
:
  GaussianBlurring<VoxelType>(sigma),
  _PaddingValue(padding_value)
{
}

// -----------------------------------------------------------------------------
template <class VoxelType> GaussianBlurringWithPadding<VoxelType>
::GaussianBlurringWithPadding(double xsigma, double ysigma, VoxelType padding_value)
:
  GaussianBlurring<VoxelType>(xsigma, ysigma, .0, .0),
  _PaddingValue(padding_value)
{
}

// -----------------------------------------------------------------------------
template <class VoxelType> GaussianBlurringWithPadding<VoxelType>
::GaussianBlurringWithPadding(double xsigma, double ysigma, double zsigma, VoxelType padding_value)
:
  GaussianBlurring<VoxelType>(xsigma, ysigma, zsigma, .0),
  _PaddingValue(padding_value)
{
}

// -----------------------------------------------------------------------------
template <class VoxelType> GaussianBlurringWithPadding<VoxelType>
::GaussianBlurringWithPadding(double xsigma, double ysigma, double zsigma, double tsigma, VoxelType padding_value)
:
  GaussianBlurring<VoxelType>(xsigma, ysigma, zsigma, tsigma),
  _PaddingValue(padding_value)
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
void GaussianBlurringWithPadding<VoxelType>::Run()
{
  MIRTK_START_TIMING();

  // Do the initial set up
  this->Initialize();

  GenericImage<VoxelType> *input  = const_cast<GenericImage<VoxelType> *>(this->Input());
  GenericImage<VoxelType> *output = this->Output();

  const ImageAttributes &attr = input->Attributes();
  const int N = ((attr._dt == .0) ? attr._t : 1);

  // Blur along x axis
  if (this->_SigmaX != .0 && input->X() > 1) {
    if (output == this->Input()) output = new GenericImage<VoxelType>(attr);
    this->InitializeKernel(this->_SigmaX / input->GetXSize());
    typedef ConvolutionFunction::ConvolveTruncatedForegroundInX<RealPixel> Conv;
    for (int n = 0; n < N; ++n) {
      Conv conv(input, _PaddingValue, this->_Kernel->Data(), this->_Kernel->X(), true, n);
      ParallelForEachVoxel(attr, input, output, conv);
    }
    swap(input, output);
  }

  // Blur along y axis
  if (this->_SigmaY != .0 && input->Y() > 1) {
    if (output == this->Input()) output = new GenericImage<VoxelType>(attr);
    this->InitializeKernel(this->_SigmaY / input->GetYSize());
    typedef ConvolutionFunction::ConvolveTruncatedForegroundInY<RealPixel> Conv;
    for (int n = 0; n < N; ++n) {
      Conv conv(input, _PaddingValue, this->_Kernel->Data(), this->_Kernel->X(), true, n);
      ParallelForEachVoxel(attr, input, output, conv);
    }
    swap(input, output);
  }

  // Blur along z axis
  if (this->_SigmaZ != .0 && input->Z() > 1) {
    if (output == this->Input()) output = new GenericImage<VoxelType>(attr);
    this->InitializeKernel(this->_SigmaZ / input->GetZSize());
    typedef ConvolutionFunction::ConvolveTruncatedForegroundInZ<RealPixel> Conv;
    for (int n = 0; n < N; ++n) {
      Conv conv(input, _PaddingValue, this->_Kernel->Data(), this->_Kernel->X(), true, n);
      ParallelForEachVoxel(attr, input, output, conv);
    }
    swap(input, output);
  }

  // Blur along t axis
  if (this->_SigmaT != .0 && input->T() > 1) {
    if (output == this->Input()) output = new GenericImage<VoxelType>(attr);
    this->InitializeKernel(this->_SigmaT / input->GetTSize());
    typedef ConvolutionFunction::ConvolveTruncatedForegroundInT<RealPixel> Conv;
    Conv conv(input, _PaddingValue, this->_Kernel->Data(), this->_Kernel->X(), true);
    ParallelForEachVoxel(attr, input, output, conv);
    swap(input, output);
  }

  // Copy result if last output (i.e., after swap input pointer) was not filter
  // output image and delete possibly additionally allocated image
  if (input != output) {
    if (input == this->Output()) {
      if (output != this->Input()) Delete(output);
    } else {
      this->Output()->CopyFrom(input->Data());
      if (input != this->Input ()) Delete(input);
    }
  }

  // Do the final cleaning up
  this->Finalize();

  MIRTK_DEBUG_TIMING(5, this->NameOfClass());
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
