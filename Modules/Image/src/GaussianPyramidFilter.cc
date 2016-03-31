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

#include "mirtk/GaussianPyramidFilter.h"

#include "mirtk/GenericImage.h"
#include "mirtk/UnaryVoxelFunction.h"
#include "mirtk/ConvolutionFunction.h"
#include "mirtk/Profiling.h"

#include <iostream>


namespace mirtk {

using namespace UnaryVoxelFunction;
using namespace ConvolutionFunction;


// =============================================================================
// Construction/Destruction
// =============================================================================

// ---------------------------------------------------------------------------
template <class VoxelType>
GaussianPyramidFilter<VoxelType>::GaussianPyramidFilter(int j)
:
  _InputLevel (0),
  _OutputLevel(j)
{
}

// ---------------------------------------------------------------------------
template <class VoxelType>
GaussianPyramidFilter<VoxelType>::GaussianPyramidFilter(int i, int j)
:
  _InputLevel (i),
  _OutputLevel(j)
{
}

// =============================================================================
// Execution
// =============================================================================

// ---------------------------------------------------------------------------
template <class VoxelType>
void GaussianPyramidFilter<VoxelType>::Initialize()
{
  // Initialize base class, but not output image except of creating a
  // temporary image instance if set output equals input image instance
  ImageToImage<VoxelType>::Initialize(false);
}

// ---------------------------------------------------------------------------
template <class VoxelType>
void GaussianPyramidFilter<VoxelType>::Downsample()
{
  const GenericImage<VoxelType> *input  = this->Input();
  GenericImage<VoxelType>       *temp   = this->Output();
  GenericImage<VoxelType>       *output = new GenericImage<VoxelType>();

  if (input->HasBackgroundValue()) {
    temp  ->PutBackgroundValueAsDouble(input->GetBackgroundValueAsDouble());
    output->PutBackgroundValueAsDouble(input->GetBackgroundValueAsDouble());
  }

  // Type of downsampling voxel function
  typedef DownsampleConvolvedExtendedForegroundInX<VoxelType, RealPixel> DownsampleInX;
  typedef DownsampleConvolvedExtendedForegroundInY<VoxelType, RealPixel> DownsampleInY;
  typedef DownsampleConvolvedExtendedForegroundInZ<VoxelType, RealPixel> DownsampleInZ;

  // FIXME: The GaussianPyramidFilter shifts the image origin incorrectly!
  cerr << "WARNING: The GaussianPyramidFilter shifts the image origin" << endl;

  // Smoothing kernel used for downsampling
  const int    size         = 5;
  RealPixel    kernel[size] = {1, 4, 6, 4, 1};
  const double norm         = 1.0 / 16.0;

  for (int i = _InputLevel; i < _OutputLevel; ++i) {
    // Set output of previous iteration as input
    if (i != _InputLevel) {
      input  = output;
      output = new GenericImage<VoxelType>();
    }

    // Input image attributes
    ImageAttributes attr = input->Attributes();
    const int n = attr._dt ? 1 : attr._t;
    attr._dt = 1.0; // s.t. ParallelForEachVoxel downsamples each component

    // Downsample x dimension
    if (attr._x > 1) {
      DownsampleInX downsampleX(input, 2, kernel, size, norm);
      attr._x  /= 2;
      attr._dx *= 2;
      attr._xorigin += attr._xaxis[0] * downsampleX._Offset * attr._dx;
      attr._yorigin += attr._yaxis[0] * downsampleX._Offset * attr._dx;
      attr._zorigin += attr._zaxis[0] * downsampleX._Offset * attr._dx;
      output->Initialize(attr, n);
      ParallelForEachVoxel(attr, output, downsampleX);
    }

    // Downsample y dimension
    if (attr._y > 1) {
      DownsampleInY downsampleY(output, 2, kernel, size, norm);
      attr._y  /= 2;
      attr._dy *= 2;
      attr._xorigin += attr._xaxis[1] * downsampleY._Offset * attr._dy;
      attr._yorigin += attr._yaxis[1] * downsampleY._Offset * attr._dy;
      attr._zorigin += attr._zaxis[1] * downsampleY._Offset * attr._dy;
      temp->Initialize(attr, n);
      ParallelForEachVoxel(attr, temp, downsampleY);
    }

    // Downsample z dimension
    if (attr._z > 1) {
      DownsampleInZ downsampleZ(temp, 2, kernel, size, norm);
      attr._z  /= 2;
      attr._dz *= 2;
      attr._xorigin += attr._xaxis[2] * downsampleZ._Offset * attr._dz;
      attr._yorigin += attr._yaxis[2] * downsampleZ._Offset * attr._dz;
      attr._zorigin += attr._zaxis[2] * downsampleZ._Offset * attr._dz;
      output->Initialize(attr, n);
      ParallelForEachVoxel(attr, output, downsampleZ);
    }

    // Delete intermediate input
    if (input != this->Input()) {
      delete input;
      input = NULL;
    }
  }

  // Copy final output to actual output
  *temp = *output;
  delete output;
}

// ---------------------------------------------------------------------------
template <class VoxelType>
void GaussianPyramidFilter<VoxelType>::Upsample()
{
  const GenericImage<VoxelType> &input  = *this->Input();
  GenericImage<VoxelType>       &output = *this->Output();
  GenericImage<VoxelType>        temp;

  // Input image attributes
  ImageAttributes attr = input.Attributes();

  attr._x *= 2; attr._dx /= 2.0;
  attr._y *= 2; attr._dy /= 2.0;
  attr._z *= 2; attr._dz /= 2.0;
  temp  .Initialize(attr);
  output.Initialize(attr);

  cerr << "GaussianPyramidFilter::Upsample: Not implemented" << endl;
  exit(1);
}

// ---------------------------------------------------------------------------
template <class VoxelType>
void GaussianPyramidFilter<VoxelType>::Run()
{
  MIRTK_START_TIMING();

  // Do the initial set up
  this->Initialize();

  // Down-/upsample image
  if (_InputLevel == _OutputLevel) {
    this->Output()->Initialize(this->Input()->Attributes());
    this->Output()->CopyFrom  (this->Input()->Data());
  } else if (_InputLevel < _OutputLevel) {
    this->Downsample();
  } else {
    this->Upsample();
  }

  // Do the final cleaning up
  this->Finalize();

  MIRTK_DEBUG_TIMING(5, "GaussianPyramidFilter::Run");
}

// =============================================================================
// Explicit template instantiations
// =============================================================================

template class GaussianPyramidFilter<char>;
template class GaussianPyramidFilter<unsigned char>;
template class GaussianPyramidFilter<short>;
template class GaussianPyramidFilter<unsigned short>;
template class GaussianPyramidFilter<int>;
template class GaussianPyramidFilter<unsigned int>;
template class GaussianPyramidFilter<float>;
template class GaussianPyramidFilter<double>;


} // namespace mirtk
