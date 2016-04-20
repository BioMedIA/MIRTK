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

#include "mirtk/ImageGradientFunction.hxx"


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
ImageGradientFunction::ImageGradientFunction()
:
  _WrtWorld          (false),
  _DefaultValue      (.0),
  _NumberOfDimensions(0),
  _Input             (NULL),
  _Orientation       (3, 3),
  _InfiniteInput     (NULL),
  _InfiniteInputOwner(false),
  _x1(.0), _y1(.0), _z1(.0), _t1(.0),
  _x2(.0), _y2(.0), _z2(.0), _t2(.0)
{
}

// -----------------------------------------------------------------------------
ImageGradientFunction::~ImageGradientFunction()
{
  if (_InfiniteInputOwner) delete _InfiniteInput;
}

// -----------------------------------------------------------------------------
template <class TImage>
ImageGradientFunction *NewInterpolator(InterpolationMode mode, int dim = 0)
{
  mode = InterpolationWithoutPadding(mode);
  switch (dim) {
    case 2: {
      switch (mode) {
        case Interpolation_Linear:     return new GenericLinearImageGradientFunction2D<TImage>();
        case Interpolation_FastLinear: return new GenericFastLinearImageGradientFunction2D<TImage>();
        default:                       return NULL;
      }
    }
    case 3: {
      switch (mode) {
        case Interpolation_Linear:     return new GenericLinearImageGradientFunction3D<TImage>();
        case Interpolation_FastLinear: return new GenericFastLinearImageGradientFunction3D<TImage>();
        default:                       return NULL;
      }
    }
    default: {
      switch (mode) {
        case Interpolation_Linear:     return new GenericLinearImageGradientFunction<TImage>();
        case Interpolation_FastLinear: return new GenericFastLinearImageGradientFunction<TImage>();
        default:                       return NULL;
      }
    }
  }
}

// -----------------------------------------------------------------------------
ImageGradientFunction *
ImageGradientFunction::New(enum InterpolationMode mode, const BaseImage *image)
{
  ImageGradientFunction *p = NULL;

  // Dimensionality of image to interpolate
  int dim  = 0;
  if (image) {
    if      (image->Z() == 1)                            dim = 2;
    else if (image->T() == 1 || image->GetTSize() == .0) dim = 3;
    else                                                 dim = 4;
  }

  // Instantiate interpolator for generic image (i.e., instance of GenericImage)
  if (!p && image) {
    typedef GenericImage<char>           CharImage;
    typedef GenericImage<unsigned char>  UCharImage;
    typedef GenericImage<short>          ShortImage;
    typedef GenericImage<unsigned short> UShortImage;
    typedef GenericImage<int>            IntImage;
    typedef GenericImage<unsigned int>   UIntImage;
    typedef GenericImage<float>          FloatImage;
    typedef GenericImage<double>         DoubleImage;

    if      (dynamic_cast<const CharImage   *>(image)) p = NewInterpolator<CharImage>  (mode, dim);
    else if (dynamic_cast<const UCharImage  *>(image)) p = NewInterpolator<UCharImage> (mode, dim);
    else if (dynamic_cast<const ShortImage  *>(image)) p = NewInterpolator<ShortImage> (mode, dim);
    else if (dynamic_cast<const UShortImage *>(image)) p = NewInterpolator<UShortImage>(mode, dim);
    else if (dynamic_cast<const IntImage    *>(image)) p = NewInterpolator<IntImage>   (mode, dim);
    else if (dynamic_cast<const UIntImage   *>(image)) p = NewInterpolator<UIntImage>  (mode, dim);
    else if (dynamic_cast<const FloatImage  *>(image)) p = NewInterpolator<FloatImage> (mode, dim);
    else if (dynamic_cast<const DoubleImage *>(image)) p = NewInterpolator<DoubleImage>(mode, dim);
  }
  // Instantiate interpolator for general image (i.e., subclass of BaseImage)
  if (!p) p = NewInterpolator<BaseImage>(mode, dim);
  // Initialize interpolator
  if (p) {
    p->NumberOfDimensions(dim);
    p->Input(image);
  // Throw error if no suitable interpolator available
  } else {
    cerr << "ImageGradientFunction::New: Interpolation mode (" << mode;
    cerr << ") not supported for " << (dim ? ToString(dim) : "N") << "D images" << endl;
    exit(1);
  }

  return p;
}

// -----------------------------------------------------------------------------
ExtrapolateImageFunction *
ImageGradientFunction::New(enum ExtrapolationMode mode, const BaseImage *image)
{
  return ExtrapolateImageFunction::New(mode, image);
}

// -----------------------------------------------------------------------------
ImageGradientFunction *
ImageGradientFunction::New(enum InterpolationMode imode,
                           enum ExtrapolationMode emode, const BaseImage *image)
{
  ImageGradientFunction *p = ImageGradientFunction::New(imode, image);
  if (emode != Extrapolation_Default) p->Extrapolator(p->New(emode, image), true);
  return p;
}

// -----------------------------------------------------------------------------
void ImageGradientFunction::Initialize(bool)
{
  // Check if input is a valid image
  if (!_Input) {
    cerr << this->NameOfClass() << "::Initialize: Input image not set" << endl;
    exit(1);
  }
  if (_Input->IsEmpty()) {
    cerr << this->NameOfClass() << "::Initialize: Input image has zero extent" << endl;
    exit(1);
  }

  // Image resolution and orientation matrix
  _VoxelSize._x = (_Input->GetXSize() == .0 ? 1.0 : _Input->GetXSize());
  _VoxelSize._y = (_Input->GetYSize() == .0 ? 1.0 : _Input->GetYSize());
  _VoxelSize._z = (_Input->GetZSize() == .0 ? 1.0 : _Input->GetZSize());
  _Orientation  = _Input->Attributes().GetWorldToImageOrientation();

  // Determine dimensionality of input (if not specified by subclass/New)
  if (_NumberOfDimensions == 0) {
    if (_Input->Z() > 1) {
      if (_Input->T() > 1 && _Input->GetTSize() != .0) _NumberOfDimensions = 4;
      else                                             _NumberOfDimensions = 3;
    } else                                             _NumberOfDimensions = 2;
  }

  // Default domain within which interpolation can be performed is assumed
  // to be identical to the entire finite image domain
  _x1 = .0;
  _y1 = .0;
  _z1 = .0;
  _t1 = .0;
  _x2 = _Input->X() - 1;
  _y2 = _Input->Y() - 1;
  _z2 = _Input->Z() - 1;
  _t2 = _Input->T() - 1;

  // Initialize extrapolator, i.e., infinite discrete image
  if (_InfiniteInput) {
    _InfiniteInput->Input(_Input);
    _InfiniteInput->Initialize();
  }
}


} // namespace mirtk

////////////////////////////////////////////////////////////////////////////////
// Explicit instantiations
////////////////////////////////////////////////////////////////////////////////

// ND
#include "mirtk/LinearImageGradientFunction.hxx"
#include "mirtk/FastLinearImageGradientFunction.hxx"

// 2D
#include "mirtk/LinearImageGradientFunction2D.hxx"
#include "mirtk/FastLinearImageGradientFunction2D.hxx"

// 3D
#include "mirtk/LinearImageGradientFunction3D.hxx"
#include "mirtk/FastLinearImageGradientFunction3D.hxx"


namespace mirtk {


// Base class
mirtkGradientInterpolatorInstantiations(GenericImageGradientFunction);

// ND
mirtkGradientInterpolatorInstantiations(GenericLinearImageGradientFunction);
mirtkGradientInterpolatorInstantiations(GenericFastLinearImageGradientFunction);

// 2D
mirtkGradientInterpolatorInstantiations(GenericLinearImageGradientFunction2D);
mirtkGradientInterpolatorInstantiations(GenericFastLinearImageGradientFunction2D);

// 3D
mirtkGradientInterpolatorInstantiations(GenericLinearImageGradientFunction3D);
mirtkGradientInterpolatorInstantiations(GenericFastLinearImageGradientFunction3D);


} // namespace mirtk
