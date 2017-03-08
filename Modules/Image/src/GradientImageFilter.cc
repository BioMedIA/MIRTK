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

#include "mirtk/GradientImageFilter.h"

#include "mirtk/Math.h"
#include "mirtk/Parallel.h"
#include "mirtk/VoxelFunction.h"

#include <limits>
#include <utility>


namespace mirtk {

// =============================================================================
// Auxiliary functors
// =============================================================================

namespace GradientImageFilterUtils {

// -----------------------------------------------------------------------------
/// Evaluate image gradient at a specific voxel using finite differences
template <class TIn, class TOut = TIn>
class EvaluateImageGradient : public VoxelFunction
{
  typedef typename GradientImageFilter<TIn>::GradientType GradientType;

  GradientType             _Type;
  const GenericImage<TIn> *_Input;
  double                   _PaddingValue;
  bool                     _UseVoxelSize;
  const Matrix            *_Orientation;
  int                      _MaxI, _MaxJ, _MaxK;
  int                      _NumberOfVoxels;

public:

  EvaluateImageGradient(GradientType type,
                        const GenericImage<TIn> *input, double padding,
                        bool use_voxel_size, const Matrix *orientation)
  :
    _Type(type),
    _Input(input),
    _PaddingValue(padding),
    _UseVoxelSize(use_voxel_size),
    _Orientation(orientation),
    _MaxI(input->X() - 1),
    _MaxJ(input->Y() - 1),
    _MaxK(input->Z() - 1),
    _NumberOfVoxels(input->NumberOfSpatialVoxels())
  {}

  void operator ()(int i, int j, int k, int l, TOut *g)
  {
    double dx = .0, dy = .0, dz = .0, v1, v2;
    const bool nan_bg = IsNaN(_PaddingValue);

    int i1 = (i == 0     ? 0     : i - 1);
    int i2 = (i == _MaxI ? _MaxI : i + 1);
    int j1 = (j == 0     ? 0     : j - 1);
    int j2 = (j == _MaxJ ? _MaxJ : j + 1);
    int k1 = (k == 0     ? 0     : k - 1);
    int k2 = (k == _MaxK ? _MaxK : k + 1);

    if (i1 < i2) {
      v1 = static_cast<double>(_Input->Get(i1, j, k, l));
      v2 = static_cast<double>(_Input->Get(i2, j, k, l));
      if (IsNaN(v1) || (!nan_bg && v1 <= _PaddingValue)) {
        if (i1 < i) {
          i1 = i;
          v1 = static_cast<double>(_Input->Get(i1, j, k, l));
          if (IsNaN(v1) || (!nan_bg && v1 <= _PaddingValue)) {
            i1 = i2;
          }
        } else {
          i1 = i2;
        }
      }
      if (IsNaN(v2) || (!nan_bg && v2 <= _PaddingValue)) {
        if (i2 > i) {
          i2 = i;
          v2 = static_cast<double>(_Input->Get(i1, j, k, l));
          if (IsNaN(v2) || (!nan_bg && v2 <= _PaddingValue)) {
            i2 = i1;
          }
        } else {
          i2 = i1;
        }
      }
      if (i1 < i2) dx = (v2 - v1) / (i2 - i1);
    }

    if (j1 < j2) {
      v1 = static_cast<double>(_Input->Get(i, j1, k, l));
      v2 = static_cast<double>(_Input->Get(i, j2, k, l));
      if (IsNaN(v1) || (!nan_bg && v1 <= _PaddingValue)) {
        if (j1 < j) {
          j1 = j;
          v1 = static_cast<double>(_Input->Get(i, j1, k, l));
          if (IsNaN(v1) || (!nan_bg && v1 <= _PaddingValue)) {
            j1 = j2;
          }
        } else {
          j1 = j2;
        }
      }
      if (IsNaN(v2) || (!nan_bg && v2 <= _PaddingValue)) {
        if (j2 > j) {
          j2 = j;
          v2 = static_cast<double>(_Input->Get(i, j2, k, l));
          if (IsNaN(v2) || (!nan_bg && v2 <= _PaddingValue)) {
            j2 = j1;
          }
        } else {
          j2 = j1;
        }
      }
      if (j1 < j2) dy = (v2 - v1) / (j2 - j1);
    }

    if (k1 < k2) {
      v1 = static_cast<double>(_Input->Get(i, j, k1, l));
      v2 = static_cast<double>(_Input->Get(i, j, k2, l));
      if (IsNaN(v1) || (!nan_bg && v1 <= _PaddingValue)) {
        if (k1 < k) {
          k1 = k;
          v1 = static_cast<double>(_Input->Get(i, j, k1, l));
          if (IsNaN(v1) || (!nan_bg && v1 <= _PaddingValue)) {
            k1 = k2;
          }
        } else {
          k1 = k2;
        }
      }
      if (IsNaN(v2) || (!nan_bg && v2 <= _PaddingValue)) {
        if (k2 > k) {
          k2 = k;
          v2 = static_cast<double>(_Input->Get(i, j, k2, l));
          if (IsNaN(v2) || (!nan_bg && v2 <= _PaddingValue)) {
            k2 = k1;
          }
        } else {
          k2 = k1;
        }
      }
      if (k1 < k2) dz = (v2 - v1) / (k2 - k1);
    }

    if (_UseVoxelSize) {
      if (_Input->GetXSize() > .0) dx /= _Input->GetXSize();
      if (_Input->GetYSize() > .0) dy /= _Input->GetYSize();
      if (_Input->GetZSize() > .0) dz /= _Input->GetZSize();
    }
    if (_Orientation) {
      // Using numerator-layout for matrix calculus.
      // http://en.wikipedia.org/wiki/Matrix_calculus#Numerator-layout_notation
      const double di = dx, dj = dy, dk = dz;
      const Matrix &R = *_Orientation;
      dx = di * R(0, 0) + dj * R(1, 0) + dk * R(2, 0);
      dy = di * R(0, 1) + dj * R(1, 1) + dk * R(2, 1);
      dz = di * R(0, 2) + dj * R(1, 2) + dk * R(2, 2);
    }

    switch (_Type) {
      case GradientType::GRADIENT_X: {
        *g = static_cast<TOut>(dx);
      } break;
      case GradientType::GRADIENT_Y: {
        *g = static_cast<TOut>(dx);
      } break;
      case GradientType::GRADIENT_Z: {
        *g = static_cast<TOut>(dz);
      } break;
      case GradientType::GRADIENT_DOT_PRODUCT: {
        *g = static_cast<TOut>(dx*dx + dy*dy + dz*dz);
      } break;
      case GradientType::GRADIENT_MAGNITUDE: {
        *g = static_cast<TOut>(sqrt(dx*dx + dy*dy + dz*dz));
      } break;
      case GradientType::GRADIENT_VECTOR: {
        g += 2 * l * _NumberOfVoxels;
        *g = static_cast<TOut>(dx), g += _NumberOfVoxels;
        *g = static_cast<TOut>(dy), g += _NumberOfVoxels;
        *g = static_cast<TOut>(dz);
      } break;
      case GradientType::NORMALISED_GRADIENT_VECTOR: {
        const double norm = sqrt(dx*dx + dy*dy + dz*dz) + 1e-10;
        g += 2 * l * _NumberOfVoxels;
        *g = static_cast<TOut>(dx / norm), g += _NumberOfVoxels;
        *g = static_cast<TOut>(dy / norm), g += _NumberOfVoxels;
        *g = static_cast<TOut>(dz / norm);
      } break;
      default: {
        *g = TOut(0);
      } break;
    }
  }
};


} // namespace GradientImageFilterUtils
using namespace GradientImageFilterUtils;

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
GradientImageFilter<VoxelType>::GradientImageFilter(GradientType type)
:
  _Type(type),
  _UseVoxelSize(true),
  _UseOrientation(false),
  _PaddingValue(-inf)
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
void GradientImageFilter<VoxelType>::Initialize()
{
  // Check inputs and outputs
  if (this->Input() == nullptr) {
    cerr << this->NameOfClass() << "::Initialize: Filter has no input" << endl;
    exit(1);
  }
  if (this->Output() == nullptr) {
    cerr << this->NameOfClass() << "::Initialize: Filter has no output" << endl;
    exit(1);
  }
  if (this->Input()->IsEmpty()) {
    cerr << this->NameOfClass() << "::Initialize: Input is empty" << endl;
    exit(1);
  }
  if (_Type < 0 || _Type >= NUM_GRADIENT_TYPES) {
    cerr << this->NameOfClass() << "::Initialize: Invalid gradient type: " << _Type << endl;
    exit(1);
  }

  // Instantiate temporary output image if necessary
  if (this->Input() == this->Output()) {
    this->Buffer(this->Output());
    this->Output(new GenericImage<VoxelType>());
  } else {
    this->Buffer(nullptr);
  }

  // Initialize output image
  if (_Type == GRADIENT_VECTOR || _Type == NORMALISED_GRADIENT_VECTOR) {
    this->Output()->Initialize(this->Input()->Attributes(), 3 * this->Input()->T());
  } else {
    this->Output()->Initialize(this->Input()->Attributes(),     this->Input()->T());
  }
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void GradientImageFilter<VoxelType>::Run()
{
  this->Initialize();

  ImageAttributes attr = this->Input()->Attributes();
  if (attr._dt == .0) attr._dt = 1.0; // Call voxel function for each scalar input

  Matrix R;
  if (_UseOrientation) R = attr.GetWorldToImageOrientation();
  EvaluateImageGradient<VoxelType> eval(
    _Type, this->Input(), _PaddingValue, _UseVoxelSize, _UseOrientation ? &R : nullptr
  );
  ParallelForEachVoxel(attr, this->Output(), eval);

  this->Finalize();
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void GradientImageFilter<VoxelType>::Finalize()
{
  if (this->Buffer()) {
    this->Buffer()->CopyFrom(this->Output()->Data());
    delete this->Output();
    this->Output(this->Buffer());
    this->Buffer(nullptr);
  }
}

// =============================================================================
// Explicit template instantiations
// =============================================================================

// -----------------------------------------------------------------------------
template class GradientImageFilter<unsigned char>;
template class GradientImageFilter<short>;
template class GradientImageFilter<float>;
template class GradientImageFilter<double>;


} // namespace mirtk
