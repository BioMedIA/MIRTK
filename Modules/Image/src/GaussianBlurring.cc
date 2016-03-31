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

#include <mirtkGaussianBlurring.h>

#include <mirtkMath.h>
#include <mirtkMemory.h>
#include <mirtkGenericImage.h>
#include <mirtkConvolutionFunction.h>
#include <mirtkScalarFunctionToImage.h>
#include <mirtkScalarGaussian.h>
#include <mirtkProfiling.h>
#include <mirtkDeallocate.h>


namespace mirtk {


// -----------------------------------------------------------------------------
template <class VoxelType>
GaussianBlurring<VoxelType>::GaussianBlurring(double sigma)
:
  _SigmaX(sigma),
  _SigmaY(sigma),
  _SigmaZ(sigma),
  _SigmaT(.0),
  _Kernel(NULL)
{
}

// -----------------------------------------------------------------------------
template <class VoxelType>
GaussianBlurring<VoxelType>::GaussianBlurring(double xsigma, double ysigma, double zsigma, double tsigma)
:
  _SigmaX(xsigma),
  _SigmaY(ysigma),
  _SigmaZ(zsigma),
  _SigmaT(tsigma),
  _Kernel(NULL)
{
}

// -----------------------------------------------------------------------------
template <class VoxelType>
GaussianBlurring<VoxelType>::~GaussianBlurring()
{
  Delete(_Kernel);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void GaussianBlurring<VoxelType>::SetSigma(double sigma)
{
  _SigmaX = _SigmaY = _SigmaZ = sigma;
  _SigmaT = .0;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void GaussianBlurring<VoxelType>::SetSigma(double xsigma, double ysigma, double zsigma, double tsigma)
{
  _SigmaX = xsigma;
  _SigmaY = ysigma;
  _SigmaZ = zsigma;
  _SigmaT = tsigma;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void GaussianBlurring<VoxelType>::Initialize()
{
  // Initialize base class
  ImageToImage<VoxelType>::Initialize();

  // Instantiate convolution kernel of proper image type
  _Kernel = new GenericImage<RealPixel>;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
int GaussianBlurring<VoxelType>::KernelSize(double sigma)
{
  return 2 * ifloor(3.0 * sigma) + 1;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void GaussianBlurring<VoxelType>::InitializeKernel(double sigma)
{
  // Create filter kernel for 1D Gaussian function
  _Kernel->Initialize(KernelSize(sigma), 1, 1);

  // Create scalar function which corresponds to a 1D Gaussian function
  ScalarGaussian gaussian(sigma,            .0, .0,  // stddev in x, y, z
                          _Kernel->X() / 2, .0, .0); // center in x, y, z

  // Convert from scalar function to filter kernel
  RealPixel *p = _Kernel->Data();
  for (int i = 0; i < _Kernel->X(); ++i, ++p) {
    (*p) = gaussian.Evaluate(i);
  }
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void GaussianBlurring<VoxelType>::Run()
{
  MIRTK_START_TIMING();

  // Do the initial set up
  this->Initialize();

  GenericImage<VoxelType> *input  = const_cast<GenericImage<VoxelType> *>(this->Input());
  GenericImage<VoxelType> *output = this->Output();

  const ImageAttributes &attr = input->Attributes();
  const int N = ((attr._dt == .0) ? attr._t : 1);

  // Blur along x axis
  if (_SigmaX != .0 && input->X() > 1) {
    if (output == this->Input()) output = new GenericImage<VoxelType>(attr);
    this->InitializeKernel(_SigmaX / input->GetXSize());
    for (int n = 0; n < N; ++n) {
      using ConvolutionFunction::ConvolveInX;
      ConvolveInX<RealPixel> conv(input, _Kernel->Data(), _Kernel->X(), true, n);
      ParallelForEachVoxel(attr, input, output, conv);
    }
    swap(input, output);
  }

  // Blur along y axis
  if (_SigmaY != .0 && input->Y() > 1) {
    if (output == this->Input()) output = new GenericImage<VoxelType>(attr);
    this->InitializeKernel(_SigmaY / input->GetYSize());
    for (int n = 0; n < N; ++n) {
      using ConvolutionFunction::ConvolveInY;
      ConvolveInY<RealPixel> conv(input, _Kernel->Data(), _Kernel->X(), true, n);
      ParallelForEachVoxel(attr, input, output, conv);
    }
    swap(input, output);
  }

  // Blur along z axis
  if (_SigmaZ != .0 && input->Z() > 1) {
    if (output == this->Input()) output = new GenericImage<VoxelType>(attr);
    this->InitializeKernel(_SigmaZ / input->GetZSize());
    for (int n = 0; n < N; ++n) {
      using ConvolutionFunction::ConvolveInZ;
      ConvolveInZ<RealPixel> conv(input, _Kernel->Data(), _Kernel->X(), true, n);
      ParallelForEachVoxel(attr, input, output, conv);
    }
    swap(input, output);
  }

  // Blur along t axis
  if (_SigmaT != .0 && input->T() > 1) {
    if (output == this->Input()) output = new GenericImage<VoxelType>(attr);
    this->InitializeKernel(_SigmaT / input->GetTSize());
    using ConvolutionFunction::ConvolveInT;
    ConvolveInT<RealPixel> conv(input, _Kernel->Data(), _Kernel->X(), true);
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

// -----------------------------------------------------------------------------
template <class VoxelType>
void GaussianBlurring<VoxelType>::Finalize()
{
  // Finalize base class
  ImageToImage<VoxelType>::Finalize();

  // Free convolution kernel
  Delete(_Kernel);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void GaussianBlurring<VoxelType>::RunX()
{
  const double ysigma = _SigmaY;
  const double zsigma = _SigmaZ;
  const double tsigma = _SigmaT;

  _SigmaY = _SigmaZ = _SigmaT = .0;
  this->Run();

  _SigmaY = ysigma;
  _SigmaZ = zsigma;
  _SigmaT = tsigma;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void GaussianBlurring<VoxelType>::RunY()
{
  const double xsigma = _SigmaX;
  const double zsigma = _SigmaZ;
  const double tsigma = _SigmaT;

  _SigmaX = _SigmaZ = _SigmaT = .0;
  this->Run();

  _SigmaX = xsigma;
  _SigmaZ = zsigma;
  _SigmaT = tsigma;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void GaussianBlurring<VoxelType>::RunZ()
{
  const double xsigma = _SigmaX;
  const double ysigma = _SigmaY;
  const double tsigma = _SigmaT;

  _SigmaX = _SigmaY = _SigmaT = .0;
  this->Run();

  _SigmaX = xsigma;
  _SigmaY = ysigma;
  _SigmaT = tsigma;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void GaussianBlurring<VoxelType>::RunT()
{
  const double xsigma = _SigmaX;
  const double ysigma = _SigmaY;
  const double zsigma = _SigmaZ;
  const double tsigma = _SigmaT;

  _SigmaX = _SigmaY = _SigmaZ = .0;
  if (_SigmaT == .0) _SigmaT = _SigmaX;
  this->Run();

  _SigmaX = xsigma;
  _SigmaY = ysigma;
  _SigmaZ = zsigma;
  _SigmaT = tsigma;
}

// =============================================================================
// Explicit template instantiations
// =============================================================================

template class GaussianBlurring<unsigned char>;
template class GaussianBlurring<short>;
template class GaussianBlurring<unsigned short>;
template class GaussianBlurring<float>;
template class GaussianBlurring<double>;


} // namespace mirtk
