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

#include "mirtk/CubicBSplineConvolution.h"

#include "mirtk/BSpline.h"


namespace mirtk {


// -----------------------------------------------------------------------------
template <class VoxelType>
CubicBSplineConvolution<VoxelType>::CubicBSplineConvolution(double r)
:
  CubicBSplineConvolution<VoxelType>::CubicBSplineConvolution(r, r, r, 0.)
{
}

// -----------------------------------------------------------------------------
template <class VoxelType>
CubicBSplineConvolution<VoxelType>::CubicBSplineConvolution(double rx, double ry, double rz, double rt)
:
  _RadiusX(rx),
  _RadiusY(ry),
  _RadiusZ(rz),
  _RadiusT(rt)
{
}

// -----------------------------------------------------------------------------
template <class VoxelType>
CubicBSplineConvolution<VoxelType>::~CubicBSplineConvolution()
{
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void CubicBSplineConvolution<VoxelType>::Radius(double r)
{
  _RadiusX = _RadiusY = _RadiusZ = r;
  _RadiusT = 0.;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void CubicBSplineConvolution<VoxelType>::Radius(double rx, double ry, double rz, double rt)
{
  _RadiusX = rx;
  _RadiusY = ry;
  _RadiusZ = rz;
  _RadiusT = rt;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
int CubicBSplineConvolution<VoxelType>::KernelSize(double rv)
{
  int r = ifloor(2. * rv);
  if (r % 2 == 0) ++r;
  return r;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
UniquePtr<typename CubicBSplineConvolution<VoxelType>::KernelType>
CubicBSplineConvolution<VoxelType>::InitializeKernel(double rv)
{
  UniquePtr<KernelType> kernel(new KernelType(KernelSize(rv), 1, 1));
  double dt = 4. / kernel->X();
  double t  = - (kernel->X() / 2) * dt;
  for (int i = 0; i < kernel->X(); ++i, t += dt) {
    kernel->Put(i, BSpline<typename KernelType::VoxelType>::B(t));
  }
  return kernel;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void CubicBSplineConvolution<VoxelType>::Initialize()
{
  // Initialize base class
  SeparableConvolution<VoxelType>::Initialize();

  // Instantiate convolution kernels
  const ImageType *input  = this->Input();
  if (!AreEqual(_RadiusX, 0.) && input->X() > 1) {
    _BSplineKernel[0] = this->InitializeKernel(_RadiusX < 0. ? -_RadiusX : _RadiusX / input->XSize());
  } else {
    _BSplineKernel[0].reset();
  }
  if (!AreEqual(_RadiusY, 0.) && input->Y() > 1) {
    _BSplineKernel[1] = this->InitializeKernel(_RadiusY < 0. ? -_RadiusY : _RadiusY / input->YSize());
  } else {
    _BSplineKernel[1].reset();
  }
  if (!AreEqual(_RadiusZ, 0.) && input->Z() > 1 && !AreEqual(input->ZSize(), 0.)) {
    _BSplineKernel[2] = this->InitializeKernel(_RadiusZ < 0. ? -_RadiusZ : _RadiusZ / input->ZSize());
  } else {
    _BSplineKernel[2].reset();
  }
  if (!AreEqual(_RadiusT, 0.) && input->T() > 1 && !AreEqual(input->TSize(), 0.)) {
    _BSplineKernel[3] = this->InitializeKernel(_RadiusT < 0. ? -_RadiusT : _RadiusT / input->TSize());
  } else {
    _BSplineKernel[3].reset();
  }

  // Set kernels
  this->_KernelX = _BSplineKernel[0].get();
  this->_KernelY = _BSplineKernel[1].get();
  this->_KernelZ = _BSplineKernel[2].get();
  this->_KernelT = _BSplineKernel[3].get();
}

// =============================================================================
// Explicit template instantiations
// =============================================================================

template class CubicBSplineConvolution<unsigned char>;
template class CubicBSplineConvolution<short>;
template class CubicBSplineConvolution<unsigned short>;
template class CubicBSplineConvolution<float>;
template class CubicBSplineConvolution<double>;


} // namespace mirtk
