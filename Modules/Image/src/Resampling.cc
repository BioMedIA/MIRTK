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

#include "mirtk/Resampling.h"

#include "mirtk/Parallel.h"
#include "mirtk/Profiling.h"
#include "mirtk/InterpolateImageFunction.h"
#include "mirtk/NearestNeighborInterpolateImageFunction.h"


namespace mirtk {


// ---------------------------------------------------------------------------
template <class VoxelType>
Resampling<VoxelType>::Resampling(double dx, double dy, double dz)
{
  _X     = 0;
  _Y     = 0;
  _Z     = 0;
  _XSize = dx;
  _YSize = dy;
  _ZSize = dz;
  _Interpolator = NULL;
}

// ---------------------------------------------------------------------------
template <class VoxelType>
Resampling<VoxelType>::Resampling(int x, int y, int z)
{
  _X     = x;
  _Y     = y;
  _Z     = z;
  _XSize = 0;
  _YSize = 0;
  _ZSize = 0;
  _Interpolator = NULL;
}

// ---------------------------------------------------------------------------
template <class VoxelType>
Resampling<VoxelType>::Resampling(int x, int y, int z, double dx, double dy, double dz)
{
  _X     = x;
  _Y     = y;
  _Z     = z;
  _XSize = dx;
  _YSize = dy;
  _ZSize = dz;
  _Interpolator = NULL;
}

// ---------------------------------------------------------------------------
template <class VoxelType>
void Resampling<VoxelType>::Initialize()
{
  // Initialize base class
  ImageToImage<VoxelType>::Initialize(false);

  // Set up interpolator
  if (_Interpolator == nullptr) {
    Throw(ERR_InvalidArgument, __FUNCTION__, "Interpolator must be set");
  }
  _Interpolator->Input(this->_Input);
  _Interpolator->Initialize();

  // Initialize output image
  this->InitializeOutput();
}

// ---------------------------------------------------------------------------
template <class VoxelType>
void Resampling<VoxelType>::InitializeOutput()
{
  ImageAttributes attr = this->_Input->Attributes();

  if (_X > 0 && _XSize > 0) {
    attr._x  = _X;
    attr._dx = _XSize;
  } else if (_X > 0) {
    attr._x  = _X;
    attr._dx = this->_Input->X() * this->_Input->XSize() / static_cast<double>(_X);
  } else if (_XSize > 0 && this->_Input->XSize() > 0. && this->_Input->X() > 1) {
    attr._x  = iceil(this->_Input->X() * this->_Input->XSize() / _XSize);
    attr._dx = _XSize;
  }

  if (_Y > 0 && _YSize > 0) {
    attr._y  = _Y;
    attr._dy = _YSize;
  } else if (_Y > 0) {
    attr._y  = _Y;
    attr._dy = this->_Input->Y() * this->_Input->YSize() / static_cast<double>(_Y);
  } else if (_YSize > 0 && this->_Input->YSize() > 0. && this->_Input->Y() > 1) {
    attr._y  = iceil(this->_Input->Y() * this->_Input->YSize() / _YSize);
    attr._dy = _YSize;
  }
  
  if (_Z > 0 && _ZSize > 0) {
    attr._z  = _Z;
    attr._dz = _ZSize;
  } else if (_Z > 0) {
    attr._z  = _Z;
    attr._dz = this->_Input->Z() * this->_Input->ZSize() / static_cast<double>(_Z);
  } else if (_ZSize > 0 && this->_Input->ZSize() > 0. && this->_Input->Z() > 1) {
    attr._z  = iceil(this->_Input->Z() * this->_Input->ZSize() / _ZSize);
    attr._dz = _ZSize;
  }

  this->_Output->Initialize(attr);
  if (this->_Input->HasBackgroundValue()) {
    this->_Output->PutBackgroundValueAsDouble(this->_Input->GetBackgroundValueAsDouble());
  }
  if (this->_Input->HasMask()) {
    const BinaryImage *other = this->_Input->GetMask();
    UniquePtr<BinaryImage> mask(new BinaryImage(attr, other->T()));
    GenericNearestNeighborExtrapolateImageFunction<BinaryImage> nn;
    nn.Input(other);
    nn.Initialize();
    double x, y, z;
    for (int l = 0; l < mask->T(); ++l)
    for (int k = 0; k < mask->Z(); ++k)
    for (int j = 0; j < mask->Y(); ++j)
    for (int i = 0; i < mask->X(); ++i) {
      x = i, y = j, z = k;
      mask->ImageToWorld(x, y, z);
      other->WorldToImage(x, y, z);
      mask->PutAsDouble(i, j, k, l, nn.Evaluate(x, y, z, l));
    }
    this->_Output->PutMask(mask.release(), true);
  }
}

// ---------------------------------------------------------------------------
// Serial/Multi-threaded body of Resampling::Run method
class ResamplingRun
{
public:
  // -------------------------------------------------------------------------
  // Data members
  const BaseImage          *_Input;
  BaseImage                *_Output;
  InterpolateImageFunction *_Interpolator;
  mutable int               _T;

  // -------------------------------------------------------------------------
  /// Default constructor
  ResamplingRun()
  :
    _Input(NULL),
    _Output(NULL),
    _Interpolator(NULL),
    _T(0)
  {}

  // -------------------------------------------------------------------------
  /// Copy constructor
  ResamplingRun(const ResamplingRun &other)
  :
    _Input(other._Input),
    _Output(other._Output),
    _Interpolator(other._Interpolator),
    _T(other._T)
  {}

  // -------------------------------------------------------------------------
  /// Resamples the image within given output region for frame _T
  void operator() (const blocked_range3d<int> &r) const
  {
    double x, y, z;
    for (int k = r.pages().begin(); k != r.pages().end(); ++k)
    for (int j = r.rows ().begin(); j != r.rows ().end(); ++j)
    for (int i = r.cols ().begin(); i != r.cols ().end(); ++i) {
      x = i, y = j, z = k;
      _Output->ImageToWorld(x, y, z);
      _Input ->WorldToImage(x, y, z);
      _Output->PutAsDouble(i, j, k, _T, _Interpolator->Evaluate(x, y, z, _T));
    }
  }

  // -------------------------------------------------------------------------
  /// Resamples the image for the given range of output frames
  void operator() (const blocked_range<int> &r) const
  {
    blocked_range3d<int> voxels(0, _Output->Z(),
                                0, _Output->Y(),
                                0, _Output->X());

    for (_T = r.begin(); _T != r.end(); _T++) {
      parallel_for(voxels, ResamplingRun(*this));
    }
  }

};

// ---------------------------------------------------------------------------
template <class VoxelType>
void Resampling<VoxelType>::Run()
{
  MIRTK_START_TIMING();

  this->Initialize();

  ResamplingRun body;
  body._Input        = this->_Input;
  body._Output       = this->_Output;
  body._Interpolator = this->_Interpolator;
  body(blocked_range<int>(0, this->_Output->T()));

  this->Finalize();

  MIRTK_DEBUG_TIMING(5, this->NameOfClass());
}

// ---------------------------------------------------------------------------
// Explicit template instantiations
template class Resampling<char>;
template class Resampling<unsigned char>;
template class Resampling<short>;
template class Resampling<unsigned short>;
template class Resampling<int>;
template class Resampling<float>;
template class Resampling<double>;


} // namespace mirtk
