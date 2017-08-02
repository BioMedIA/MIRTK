/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2017 Imperial College London
 * Copyright 2008-2013 Daniel Rueckert, Julia Schnabel
 * Copyright 2013-2017 Andreas Schuh
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

#include "mirtk/GaussianBlurring.h"

#include "mirtk/Math.h"
#include "mirtk/ScalarFunctionToImage.h"
#include "mirtk/ScalarGaussian.h"


namespace mirtk {


// -----------------------------------------------------------------------------
template <class VoxelType>
GaussianBlurring<VoxelType>::GaussianBlurring(double sigma)
:
  GaussianBlurring<VoxelType>::GaussianBlurring(sigma, sigma, sigma, 0.)
{
}

// -----------------------------------------------------------------------------
template <class VoxelType>
GaussianBlurring<VoxelType>::GaussianBlurring(double xsigma, double ysigma, double zsigma, double tsigma)
:
  _SigmaX(xsigma),
  _SigmaY(ysigma),
  _SigmaZ(zsigma),
  _SigmaT(tsigma)
{
}

// -----------------------------------------------------------------------------
template <class VoxelType>
GaussianBlurring<VoxelType>::~GaussianBlurring()
{
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void GaussianBlurring<VoxelType>::SetSigma(double sigma)
{
  _SigmaX = _SigmaY = _SigmaZ = sigma;
  _SigmaT = 0.;
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
int GaussianBlurring<VoxelType>::KernelSize(double sigma)
{
  return 2 * ifloor(3.0 * sigma) + 1;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
UniquePtr<typename GaussianBlurring<VoxelType>::KernelType>
GaussianBlurring<VoxelType>::InitializeKernel(double sigma)
{
  // Create filter kernel for 1D Gaussian function
  UniquePtr<KernelType> kernel(new KernelType(KernelSize(sigma), 1, 1));

  // Create scalar function which corresponds to a 1D Gaussian function
  ScalarGaussian gaussian(sigma,           0., 0.,  // stddev in x, y, z
                          kernel->X() / 2, 0., 0.); // center in x, y, z

  // Convert from scalar function to filter kernel
  for (int i = 0; i < kernel->X(); ++i) {
    kernel->PutAsDouble(i, gaussian.Evaluate(i));
  }

  return kernel;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void GaussianBlurring<VoxelType>::Initialize()
{
  // Initialize base class
  SeparableConvolution<VoxelType>::Initialize();

  // Instantiate convolution kernels
  const ImageType *input  = this->Input();
  if (!AreEqual(_SigmaX, 0.) && input->X() > 1) {
    _GaussianKernel[0] = this->InitializeKernel(_SigmaX / input->XSize());
  } else {
    _GaussianKernel[0].reset();
  }
  if (!AreEqual(_SigmaY, 0.) && input->Y() > 1) {
    _GaussianKernel[1] = this->InitializeKernel(_SigmaY / input->YSize());
  } else {
    _GaussianKernel[1].reset();
  }
  if (!AreEqual(_SigmaZ, 0.) && input->Z() > 1 && !AreEqual(input->ZSize(), 0.)) {
    _GaussianKernel[2] = this->InitializeKernel(_SigmaZ / input->ZSize());
  } else {
    _GaussianKernel[2].reset();
  }
  if (!AreEqual(_SigmaT, 0.) && input->T() > 1 && !AreEqual(input->TSize(), 0.)) {
    _GaussianKernel[3] = this->InitializeKernel(_SigmaT / input->TSize());
  } else {
    _GaussianKernel[3].reset();
  }

  // Set kernels
  this->_KernelX = _GaussianKernel[0].get();
  this->_KernelY = _GaussianKernel[1].get();
  this->_KernelZ = _GaussianKernel[2].get();
  this->_KernelT = _GaussianKernel[3].get();
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
