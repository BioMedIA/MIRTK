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

#include "mirtk/Downsampling.h"
#include "mirtk/GenericImage.h"
#include "mirtk/UnaryVoxelFunction.h"
#include "mirtk/ConvolutionFunction.h"
#include "mirtk/ScalarFunction.h"
#include "mirtk/ScalarFunctionToImage.h"
#include "mirtk/Profiling.h"


namespace mirtk {


// ---------------------------------------------------------------------------
template <class VoxelType>
Downsampling<VoxelType>::Downsampling(int m)
:
  _DownsampleFactorX(m),
  _DownsampleFactorY(m),
  _DownsampleFactorZ(m),
  _Kernel           (NULL),
  _KernelSize       (0),
  _NormalizeKernel  (true)
{
}

// ---------------------------------------------------------------------------
template <class VoxelType>
Downsampling<VoxelType>::Downsampling(int mx, int my, int mz)
:
  _DownsampleFactorX(mx),
  _DownsampleFactorY(my),
  _DownsampleFactorZ(mz),
  _Kernel           (NULL),
  _KernelSize       (0),
  _NormalizeKernel  (true)
{
}

// ---------------------------------------------------------------------------
template <class VoxelType>
void Downsampling<VoxelType>::Kernel(ScalarFunction *kernel, int r, bool normalize)
{
  _Kernel          = kernel;
  _KernelSize      = 2 * r + 1;
  _NormalizeKernel = normalize;
}

// ---------------------------------------------------------------------------
template <class VoxelType>
void Downsampling<VoxelType>::Initialize()
{
  // Initialize base class, but not output image except of creating a
  // temporary image instance if set output equals input image instance
  ImageToImage<VoxelType>::Initialize(false);
}

// ---------------------------------------------------------------------------
template <class VoxelType>
void Downsampling<VoxelType>::DiscreteKernel(GenericImage<RealPixel> &kernel, double dx)
{
  // Initialize discrete kernel
  kernel.Initialize  (_KernelSize, 1,   1  );
  kernel.PutPixelSize(         dx, 1.0, 1.0);
  // Sample continuous kernel function
  ScalarFunctionToImage<RealPixel> generator;
  generator.Input (_Kernel);
  generator.Output(&kernel);
  generator.Run();
  // Normalize discrete kernel weights
  if (_NormalizeKernel) {
    RealPixel *v, sum = .0;
    v = kernel.Data();
    for (int i = 0; i < _KernelSize; ++i) sum += (*v++);
    v = kernel.Data();
    for (int i = 0; i < _KernelSize; ++i) (*v++) /= sum;
  }
}

// ---------------------------------------------------------------------------
template <class VoxelType>
void Downsampling<VoxelType>::Run()
{
  using namespace UnaryVoxelFunction;
  using namespace ConvolutionFunction;

  MIRTK_START_TIMING();

  // Do the initial set up
  this->Initialize();

  GenericImage<RealPixel>        kernel;
  const GenericImage<VoxelType> &input  = *this->Input();
  GenericImage<VoxelType>       &output = *this->Output();
  GenericImage<VoxelType>        image;

  if (input.HasBackgroundValue()) {
    image .PutBackgroundValueAsDouble(input.GetBackgroundValueAsDouble());
    output.PutBackgroundValueAsDouble(input.GetBackgroundValueAsDouble());
  }

  // Initial output image attributes
  ImageAttributes attr = input.Attributes();

  // Downsample x dimension
  attr._x  /= _DownsampleFactorX;
  attr._dx *= _DownsampleFactorX;
  output.Initialize(attr);

  if (_Kernel) {
    DiscreteKernel(kernel, input.GetXSize());
    DownsampleConvolvedMirroredForegroundInX<VoxelType, RealPixel> downsample(&input, _DownsampleFactorX, &kernel);
    ParallelForEachVoxel(attr, output, downsample);
  } else {
    DownsampleX<VoxelType> downsample(&input, _DownsampleFactorX);
    ParallelForEachVoxel(attr, output, downsample);
  }

  // Downsample y dimension
  attr._y  /= _DownsampleFactorY;
  attr._dy *= _DownsampleFactorY;
  image.Initialize(attr);

  if (_Kernel) {
    DiscreteKernel(kernel, input.GetYSize());
    DownsampleConvolvedMirroredForegroundInY<VoxelType, RealPixel> downsample(&output, _DownsampleFactorY, &kernel);
    ParallelForEachVoxel(attr, image, downsample);
  } else {
    DownsampleY<VoxelType> downsample(&output, _DownsampleFactorY);
    ParallelForEachVoxel(attr, image, downsample);
  }

  // Downsample z dimension
  attr._z  /= _DownsampleFactorZ;
  attr._dz *= _DownsampleFactorZ;
  output.Initialize(attr);

  if (_Kernel) {
    DiscreteKernel(kernel, input.GetZSize());
    DownsampleConvolvedMirroredForegroundInZ<VoxelType, RealPixel> downsample(&image, _DownsampleFactorZ, &kernel);
    ParallelForEachVoxel(attr, output, downsample);
  } else {
    DownsampleZ<VoxelType> downsample(&image, _DownsampleFactorZ);
    ParallelForEachVoxel(attr, output, downsample);
  }

  // Adjust origin if necessary
  //
  // Note: The first voxel of the output is offset such that the downsampled image
  //       region is centered around the origin of the input image as best as
  //       possible already. Only if the left and right voxel margins are not
  //       identical, the origin of the output must be adjusted by half an input voxel.
  double x = (input.X() - 1 - (input.X() % _DownsampleFactorX) % 2) / 2.0;
  double y = (input.Y() - 1 - (input.Y() % _DownsampleFactorY) % 2) / 2.0;
  double z = (input.Z() - 1 - (input.Z() % _DownsampleFactorZ) % 2) / 2.0;
  input .ImageToWorld(x, y, z);
  output.PutOrigin   (x, y, z);

  // Do the final cleaning up
  this->Finalize();

  MIRTK_DEBUG_TIMING(3, "Downsampling::Run");
}

// ---------------------------------------------------------------------------
// Explicit template instantiations
template class Downsampling<char>;
template class Downsampling<unsigned char>;
template class Downsampling<short>;
template class Downsampling<unsigned short>;
template class Downsampling<int>;
template class Downsampling<float>;
template class Downsampling<double>;


} // namespace mirtk
