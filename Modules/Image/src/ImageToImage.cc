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

#include "mirtk/ImageToImage.h"

#include "mirtk/Memory.h"
#include "mirtk/Parallel.h"
#include "mirtk/Profiling.h"


namespace mirtk {

// =============================================================================
// Multi-threaded execution
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
class MultiThreadedImageToImage
{
  /// Time frame to transform
  int _t;

  /// Pointer to image transformation class
  ImageToImage<VoxelType> *_Filter;

public:

  MultiThreadedImageToImage(ImageToImage<VoxelType> *filter, int t) {
    _t = t;
    _Filter = filter;
  }

  void operator()(const blocked_range<int> &r) const
  {
    for (int k = r.begin(); k != r.end(); ++k)
    for (int j = 0; j < _Filter->Input()->Y(); ++j)
    for (int i = 0; i < _Filter->Input()->X(); ++i) {
      _Filter->Output()->PutAsDouble(i, j, k, _t, _Filter->Run(i, j, k, _t));
    }
  }
};

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
ImageToImage<VoxelType>::ImageToImage()
:
  _Input(NULL),
  _Output(NULL),
  _Buffer(NULL)
{
}

// -----------------------------------------------------------------------------
template <class VoxelType>
ImageToImage<VoxelType>::~ImageToImage()
{
}

// =============================================================================
// Execution
// =============================================================================

// ---------------------------------------------------------------------------
template <class VoxelType>
bool ImageToImage<VoxelType>::RequiresBuffering() const
{
  return true;
}

// ---------------------------------------------------------------------------
template <class VoxelType>
void ImageToImage<VoxelType>::Initialize(bool initialize_output)
{
  // Check inputs and outputs
  if (_Input == NULL) {
    cerr << this->NameOfClass() << "::Initialize: Filter has no input" << endl;
    exit(1);
  }

  if (_Output == NULL) {
    cerr << this->NameOfClass() << "::Initialize: Filter has no output" << endl;
    exit(1);
  }

  if (_Input->IsEmpty()) {
    cerr << this->NameOfClass() << "::Initialize: Input is empty" << endl;
    exit(1);
  }

  // Allocate buffer output if needed
  if (this->RequiresBuffering()) {
    if (_Input == _Output) {
      _Buffer = _Output;
      _Output = new GenericImage<VoxelType>;
    } else {
      _Buffer = NULL;
    }
  }

  // Make sure that output has the correct dimensions
  if (initialize_output && _Input != _Output) {
    _Output->Initialize(_Input->Attributes());
  }
}

// ---------------------------------------------------------------------------
template <class VoxelType>
void ImageToImage<VoxelType>::Initialize()
{
  Initialize(true);
}

// ---------------------------------------------------------------------------
template <class VoxelType>
void ImageToImage<VoxelType>::Finalize()
{
  if (_Buffer) {
    *_Buffer = *_Output;
    swap(_Buffer, _Output);
    Delete(_Buffer);
  }
}

// ---------------------------------------------------------------------------
template <class VoxelType>
double ImageToImage<VoxelType>::Run(int, int, int, int)
{
  cerr << "Filter " << this->NameOfClass() << " has no Run(int, int, int) ";
  cerr << "member function. Using ImageToImage::Run." << endl;
  return 0;
}

// ---------------------------------------------------------------------------
template <class VoxelType>
void ImageToImage<VoxelType>::Run()
{
  MIRTK_START_TIMING();

  this->Initialize();

  for (int t = 0; t < _Input->T(); ++t) {
    MultiThreadedImageToImage<VoxelType> eval(this, t);
    parallel_for(blocked_range<int>(0, _Output->Z(), 1), eval);
  }

  this->Finalize();

  MIRTK_DEBUG_TIMING(5, this->NameOfClass());
}

// =============================================================================
// Explicit template instantiations
// =============================================================================

template class ImageToImage<char>;
template class ImageToImage<unsigned char>;
template class ImageToImage<short>;
template class ImageToImage<unsigned short>;
template class ImageToImage<int>;
template class ImageToImage<unsigned int>;
template class ImageToImage<float>;
template class ImageToImage<double>;
template class ImageToImage<Float3>;
template class ImageToImage<Double3>;
template class ImageToImage<Float9>;
template class ImageToImage<Double9>;


} // namespace mirtk
