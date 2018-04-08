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

#include "mirtk/InterpolateImageFunction.hxx"

#include "mirtk/GenericImage.h"


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
InterpolateImageFunction::InterpolateImageFunction()
:
  _NumberOfDimensions(0),
  _InfiniteInput     (NULL),
  _InfiniteInputOwner(false),
  _x1(.0), _y1(.0), _z1(.0), _t1(.0),
  _x2(.0), _y2(.0), _z2(.0), _t2(.0)
{
}

// -----------------------------------------------------------------------------
InterpolateImageFunction::~InterpolateImageFunction()
{
  if (_InfiniteInputOwner) delete _InfiniteInput;
}

// -----------------------------------------------------------------------------
template <class TImage>
InterpolateImageFunction *NewInterpolator(enum InterpolationMode mode, int dim = 0)
{
  if (mode == Interpolation_Default) {
    mode = DefaultInterpolationMode();
  } else {
    mode = InterpolationWithoutPadding(mode);
  }
  switch (dim) {
    case 2: {
      switch (mode) {
        case Interpolation_NN:               return new GenericNearestNeighborInterpolateImageFunction<TImage>();
        case Interpolation_Linear:           return new GenericLinearInterpolateImageFunction2D<TImage>();
        case Interpolation_FastLinear:       return new GenericLinearInterpolateImageFunction2D<TImage>();
        case Interpolation_BSpline:          return new GenericBSplineInterpolateImageFunction2D<TImage>();
        case Interpolation_CubicBSpline:     return new GenericCubicBSplineInterpolateImageFunction2D<TImage>();
        case Interpolation_FastCubicBSpline: return new GenericFastCubicBSplineInterpolateImageFunction2D<TImage>();
        case Interpolation_CSpline:          return new GenericCSplineInterpolateImageFunction2D<TImage>();
        case Interpolation_Gaussian:         return new GenericGaussianInterpolateImageFunction2D<TImage>();
        case Interpolation_Sinc:             return new GenericSincInterpolateImageFunction2D<TImage>();
        default:                             return NULL;
      }
    }
    case 3: {
      switch (mode) {
        case Interpolation_NN:               return new GenericNearestNeighborInterpolateImageFunction<TImage>();
        case Interpolation_Linear:           return new GenericLinearInterpolateImageFunction3D<TImage>();
        case Interpolation_FastLinear:       return new GenericLinearInterpolateImageFunction3D<TImage>();
        case Interpolation_BSpline:          return new GenericBSplineInterpolateImageFunction3D<TImage>();
        case Interpolation_CubicBSpline:     return new GenericCubicBSplineInterpolateImageFunction3D<TImage>();
        case Interpolation_FastCubicBSpline: return new GenericFastCubicBSplineInterpolateImageFunction3D<TImage>();
        case Interpolation_CSpline:          return new GenericCSplineInterpolateImageFunction3D<TImage>();
        case Interpolation_Gaussian:         return new GenericGaussianInterpolateImageFunction3D<TImage>();
        case Interpolation_Sinc:             return new GenericSincInterpolateImageFunction3D<TImage>();
        default:                             return NULL;
      }
    }
    case 4: {
      switch (mode) {
        case Interpolation_NN:               return new GenericNearestNeighborInterpolateImageFunction<TImage>();
        case Interpolation_Linear:           return new GenericLinearInterpolateImageFunction4D<TImage>();
        case Interpolation_FastLinear:       return new GenericLinearInterpolateImageFunction4D<TImage>();
        case Interpolation_BSpline:          return new GenericBSplineInterpolateImageFunction4D<TImage>();
        case Interpolation_CubicBSpline:     return new GenericCubicBSplineInterpolateImageFunction4D<TImage>();
        case Interpolation_FastCubicBSpline: return new GenericFastCubicBSplineInterpolateImageFunction4D<TImage>();
        case Interpolation_CSpline:          return new GenericCSplineInterpolateImageFunction4D<TImage>();
        case Interpolation_Gaussian:         return new GenericGaussianInterpolateImageFunction4D<TImage>();
        case Interpolation_Sinc:             return new GenericSincInterpolateImageFunction4D<TImage>();
        default:                             return NULL;
      }
    }
    default: {
      switch (mode) {
        case Interpolation_NN:               return new GenericNearestNeighborInterpolateImageFunction<TImage>();
        case Interpolation_Linear:           return new GenericLinearInterpolateImageFunction<TImage>();
        case Interpolation_FastLinear:       return new GenericLinearInterpolateImageFunction<TImage>();
        case Interpolation_BSpline:          return new GenericBSplineInterpolateImageFunction<TImage>();
        case Interpolation_CubicBSpline:     return new GenericCubicBSplineInterpolateImageFunction<TImage>();
        case Interpolation_FastCubicBSpline: return new GenericFastCubicBSplineInterpolateImageFunction<TImage>();
        case Interpolation_CSpline:          return new GenericCSplineInterpolateImageFunction<TImage>();
        case Interpolation_Gaussian:         return new GenericGaussianInterpolateImageFunction<TImage>();
        case Interpolation_Sinc:             return new GenericSincInterpolateImageFunction<TImage>();
        default:                             return NULL;
      }
    }
  }
}

// -----------------------------------------------------------------------------
InterpolateImageFunction *
InterpolateImageFunction::New(enum InterpolationMode mode, const BaseImage *image)
{
  InterpolateImageFunction *p = NULL;
  if (mode == Interpolation_Default) {
    mode = DefaultInterpolationMode();
  } else {
    mode = InterpolationWithoutPadding(mode);
  }

  // Dimensionality of image to interpolate
  int dim  = 0;
  if (image) {
    if      (image->Z() == 1)                            dim = 2;
    else if (image->T() == 1 || image->GetTSize() == .0) dim = 3;
    else                                                 dim = 4;
  }
  // Instantiate special purpose interpolators
  if (mode == Interpolation_SBased) {
    // Only implemented for 3D scalar images
    if (dim == 3 && (!image || image->N() == 1)) {
      p = new ShapeBasedInterpolateImageFunction();
    }
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
    cerr << "InterpolateImageFunction::New: Interpolation mode (" << mode;
    cerr << ") not supported for " << (dim ? ToString(dim) : "N") << "D images" << endl;
    exit(1);
  }

  return p;
}

// -----------------------------------------------------------------------------
ExtrapolateImageFunction *
InterpolateImageFunction::New(enum ExtrapolationMode mode, const BaseImage *image)
{
  return ExtrapolateImageFunction::New(mode, image);
}

// -----------------------------------------------------------------------------
InterpolateImageFunction *
InterpolateImageFunction::New(enum InterpolationMode imode,
                                  enum ExtrapolationMode emode, const BaseImage *image)
{
  InterpolateImageFunction *p = InterpolateImageFunction::New(imode, image);
  if (emode != Extrapolation_Default) p->Extrapolator(p->New(emode, image), true);
  return p;
}

// -----------------------------------------------------------------------------
void InterpolateImageFunction::Initialize(bool)
{
  // Initialize image function
  ImageFunction::Initialize();

  // Check if input is a valid image
  if (Input()->IsEmpty()) {
    cerr << this->NameOfClass() << "::Initialize: Input image has zero extent" << endl;
    exit(1);
  }

  // Determine dimensionality of input (if not specified by subclass/New)
  if (_NumberOfDimensions == 0) {
    if (Input()->Z() > 1) {
      if (Input()->T() > 1 && Input()->GetTSize() != .0) _NumberOfDimensions = 4;
      else                                               _NumberOfDimensions = 3;
    } else                                               _NumberOfDimensions = 2;
  }

  // Default domain within which interpolation can be performed is assumed
  // to be identical to the entire finite image domain
  _x1 = .0;
  _y1 = .0;
  _z1 = .0;
  _t1 = .0;
  _x2 = Input()->X() - 1;
  _y2 = Input()->Y() - 1;
  _z2 = Input()->Z() - 1;
  _t2 = Input()->T() - 1;

  // Initialize extrapolator, i.e., infinite discrete image
  if (_InfiniteInput) {
    _InfiniteInput->Input(this->Input());
    _InfiniteInput->Initialize();
  }
}


} // namespace mirtk

////////////////////////////////////////////////////////////////////////////////
// Explicit instantiations
////////////////////////////////////////////////////////////////////////////////

// ND
#include "mirtk/NearestNeighborInterpolateImageFunction.hxx"
#include "mirtk/LinearInterpolateImageFunction.hxx"
#include "mirtk/BSplineInterpolateImageFunction.hxx"
#include "mirtk/CubicBSplineInterpolateImageFunction.hxx"
#include "mirtk/FastCubicBSplineInterpolateImageFunction.hxx"
#include "mirtk/CSplineInterpolateImageFunction.hxx"
#include "mirtk/GaussianInterpolateImageFunction.hxx"
#include "mirtk/SincInterpolateImageFunction.hxx"

// 2D
#include "mirtk/LinearInterpolateImageFunction2D.hxx"
#include "mirtk/BSplineInterpolateImageFunction2D.hxx"
#include "mirtk/CubicBSplineInterpolateImageFunction2D.hxx"
#include "mirtk/FastCubicBSplineInterpolateImageFunction2D.hxx"
#include "mirtk/CSplineInterpolateImageFunction2D.hxx"
#include "mirtk/GaussianInterpolateImageFunction2D.hxx"
#include "mirtk/SincInterpolateImageFunction2D.hxx"

// 3D
#include "mirtk/LinearInterpolateImageFunction3D.hxx"
#include "mirtk/BSplineInterpolateImageFunction3D.hxx"
#include "mirtk/CubicBSplineInterpolateImageFunction3D.hxx"
#include "mirtk/FastCubicBSplineInterpolateImageFunction3D.hxx"
#include "mirtk/CSplineInterpolateImageFunction3D.hxx"
#include "mirtk/GaussianInterpolateImageFunction3D.hxx"
#include "mirtk/SincInterpolateImageFunction3D.hxx"

// 4D
#include "mirtk/LinearInterpolateImageFunction4D.hxx"
#include "mirtk/BSplineInterpolateImageFunction4D.hxx"
#include "mirtk/CubicBSplineInterpolateImageFunction4D.hxx"
#include "mirtk/FastCubicBSplineInterpolateImageFunction4D.hxx"
#include "mirtk/CSplineInterpolateImageFunction4D.hxx"
#include "mirtk/GaussianInterpolateImageFunction4D.hxx"
#include "mirtk/SincInterpolateImageFunction4D.hxx"


namespace mirtk {


// Base class
template class GenericInterpolateImageFunction<ByteImage>;
mirtkInterpolatorInstantiations(GenericInterpolateImageFunction);

// ND
template class GenericNearestNeighborInterpolateImageFunction<ByteImage>;
template class GenericLinearInterpolateImageFunction<ByteImage>;
mirtkInterpolatorInstantiations(GenericNearestNeighborInterpolateImageFunction);
mirtkInterpolatorInstantiations(GenericLinearInterpolateImageFunction);
mirtkInterpolatorInstantiations(GenericBSplineInterpolateImageFunction);
mirtkInterpolatorInstantiations(GenericCubicBSplineInterpolateImageFunction);
mirtkInterpolatorInstantiations(GenericFastCubicBSplineInterpolateImageFunction);
mirtkInterpolatorInstantiations(GenericCSplineInterpolateImageFunction);
mirtkInterpolatorInstantiations(GenericGaussianInterpolateImageFunction);
mirtkInterpolatorInstantiations(GenericSincInterpolateImageFunction);

// 2D
template class GenericLinearInterpolateImageFunction2D<ByteImage>;
mirtkInterpolatorInstantiations(GenericLinearInterpolateImageFunction2D);
mirtkInterpolatorInstantiations(GenericBSplineInterpolateImageFunction2D);
mirtkInterpolatorInstantiations(GenericCubicBSplineInterpolateImageFunction2D);
mirtkInterpolatorInstantiations(GenericFastCubicBSplineInterpolateImageFunction2D);
mirtkInterpolatorInstantiations(GenericCSplineInterpolateImageFunction2D);
mirtkInterpolatorInstantiations(GenericGaussianInterpolateImageFunction2D);
mirtkInterpolatorInstantiations(GenericSincInterpolateImageFunction2D);

// 3D
template class GenericLinearInterpolateImageFunction3D<ByteImage>;
mirtkInterpolatorInstantiations(GenericLinearInterpolateImageFunction3D);
mirtkInterpolatorInstantiations(GenericBSplineInterpolateImageFunction3D);
mirtkInterpolatorInstantiations(GenericCubicBSplineInterpolateImageFunction3D);
mirtkInterpolatorInstantiations(GenericFastCubicBSplineInterpolateImageFunction3D);
mirtkInterpolatorInstantiations(GenericCSplineInterpolateImageFunction3D);
mirtkInterpolatorInstantiations(GenericGaussianInterpolateImageFunction3D);
mirtkInterpolatorInstantiations(GenericSincInterpolateImageFunction3D);

// 4D
template class GenericLinearInterpolateImageFunction4D<ByteImage>;
mirtkInterpolatorInstantiations(GenericLinearInterpolateImageFunction4D);
mirtkInterpolatorInstantiations(GenericBSplineInterpolateImageFunction4D);
mirtkInterpolatorInstantiations(GenericCubicBSplineInterpolateImageFunction4D);
mirtkInterpolatorInstantiations(GenericFastCubicBSplineInterpolateImageFunction4D);
mirtkInterpolatorInstantiations(GenericCSplineInterpolateImageFunction4D);
mirtkInterpolatorInstantiations(GenericGaussianInterpolateImageFunction4D);
mirtkInterpolatorInstantiations(GenericSincInterpolateImageFunction4D);


} // namespace mirtk
