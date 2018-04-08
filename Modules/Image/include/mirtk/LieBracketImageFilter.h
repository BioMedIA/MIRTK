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

#ifndef MIRTK_LieBracketImageFilter_H
#define MIRTK_LieBracketImageFilter_H

#include "mirtk/ImageToImage.h"

#include "mirtk/ImageAttributes.h"
#include "mirtk/BaseImage.h"
#include "mirtk/GenericImage.h"


namespace mirtk {


/**
 * Base class for image filters which compute the Lie bracket of two vector fields.
 */
template <class TVoxel>
class LieBracketImageFilter : public ImageToImage<TVoxel>
{
  mirtkAbstractImageFilterMacro(LieBracketImageFilter, TVoxel);

protected:

  using Baseclass::Input;

  /// Second input vector field
  const ImageType *_Input2;

  /// Constructor
  LieBracketImageFilter();

  /// Initialize filter
  virtual void Initialize();

public:

  /// Construct Lie bracket filter for given image domain
  static LieBracketImageFilter *New(const ImageAttributes &, bool = true);

  /// Construct Lie bracket filter for given input vector field
  static LieBracketImageFilter *New(const BaseImage *, bool = true);

  /// Destructor
  virtual ~LieBracketImageFilter();

  /// Set n-th input
  virtual void Input(int, const ImageType *);

  /// Get n-th input
  virtual const ImageType *Input(int) const;

  /// Get n-th input
  ///
  /// This overload can be used to disambiguate the getter Input(0) from
  /// the setter Input(NULL). Alternatively, use Input(int(0)).
  virtual const ImageType *GetInput(int) const;

};

///////////////////////////////////////////////////////////////////////////////
// Inline definitions
///////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
template <class VoxelType>
LieBracketImageFilter<VoxelType>::LieBracketImageFilter()
:
  _Input2(NULL)
{
}

// ----------------------------------------------------------------------------
template <class VoxelType>
LieBracketImageFilter<VoxelType>::~LieBracketImageFilter()
{
}

// --------------------------------------------------------------------------
template <class VoxelType>
void LieBracketImageFilter<VoxelType>::Input(int i, const ImageType *image)
{
  if      (i == 0) Baseclass::Input(image);
  else if (i == 1) _Input2 = image;
  else {
    cerr << this->NameOfClass() << "::Input: Input index out of range: " << i << endl;
    exit(1);
  }
}

// --------------------------------------------------------------------------
template <class VoxelType>
const typename LieBracketImageFilter<VoxelType>::ImageType *
LieBracketImageFilter<VoxelType>::Input(int i) const
{
  if      (i == 0) return Baseclass::Input();
  else if (i == 1) return _Input2;
  else {
    cerr << this->NameOfClass() << "::Input: Input index out of range: " << i << endl;
    exit(1);
  }
}

// --------------------------------------------------------------------------
template <class VoxelType>
const typename LieBracketImageFilter<VoxelType>::ImageType *
LieBracketImageFilter<VoxelType>::GetInput(int i) const
{
  return Input(i); // with explicit "int" type
}

// --------------------------------------------------------------------------
template <class VoxelType>
void LieBracketImageFilter<VoxelType>::Initialize()
{
  // Initialize base class
  Baseclass::Initialize();

  // Check second input
  if (_Input2 == NULL) {
    cerr << this->NameOfClass() << "::Initialize: Filter has no second input" << endl;
    exit(1);
  }
  if (!_Input2->HasSpatialAttributesOf(this->_Input) || _Input2->T() != this->_Input->T()) {
    cerr << this->NameOfClass() << "::Initialize: Attributes of input images do not match" << endl;
    exit(1);
  }
}


} // namespace mirtk

#endif // MIRTK_LieBracketImageFilter_H

///////////////////////////////////////////////////////////////////////////////
// Instantiation of filter implementation
///////////////////////////////////////////////////////////////////////////////

#ifndef MIRTK_LieBracketImageFilterNew_H
#define MIRTK_LieBracketImageFilterNew_H

#include "mirtk/Memory.h"
#include "mirtk/GenericImage.h"
#include "mirtk/LieBracketImageFilter2D.h"
#include "mirtk/LieBracketImageFilter3D.h"
#include "mirtk/DifferenceOfCompositionLieBracketImageFilter3D.h"


namespace mirtk {


// ----------------------------------------------------------------------------
template <class VoxelType>
LieBracketImageFilter<VoxelType> *
LieBracketImageFilter<VoxelType>::New(const ImageAttributes &attr, bool usejac)
{
  if (attr._z > 1) {
    if (usejac) return new LieBracketImageFilter3D<VoxelType>;
    else        return new DifferenceOfCompositionLieBracketImageFilter3D<VoxelType>;
  } else {
    if (usejac) return new LieBracketImageFilter2D<VoxelType>;
    else {
      cerr << NameOfType() << "::New: DifferenceOfCompositionLieBracketImageFilter2D not implemented" << endl;
      exit(1);
    }
  }
}

// ----------------------------------------------------------------------------
template <class VoxelType>
LieBracketImageFilter<VoxelType> *
LieBracketImageFilter<VoxelType>::New(const Image *image, bool usejac)
{
  return LieBracketImageFilter::New(image->Attributes(), usejac);
}

// ----------------------------------------------------------------------------
template <class VoxelType>
void liebracket(GenericImage<VoxelType> *ov,
                const GenericImage<VoxelType> *lv,
                const GenericImage<VoxelType> *rv, bool usejac = true)
{
  typedef LieBracketImageFilter<VoxelType> LieBracketFilter;
  UniquePtr<LieBracketFilter> filter(LieBracketFilter::New(lv, usejac));
  filter->Input(0, lv);
  filter->Input(1, rv);
  filter->Output(ov);
  filter->Run();
}


} // namespace mirtk

#endif // MIRTK_LieBracketImageFilterNew_H
