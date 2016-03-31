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

#ifndef MIRTK_ImageGradientFunction_HXX
#define MIRTK_ImageGradientFunction_HXX

#include "mirtk/ImageGradientFunction.h"


namespace mirtk {


////////////////////////////////////////////////////////////////////////////////
// GenericInterpolateImageFunction
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
template <class TImage>
inline GenericImageGradientFunction<TImage>
::GenericImageGradientFunction()
{
}

// -----------------------------------------------------------------------------
template <class TImage>
inline GenericImageGradientFunction<TImage>
::~GenericImageGradientFunction()
{
}

// -----------------------------------------------------------------------------
template <class TImage>
ExtrapolateImageFunction *
GenericImageGradientFunction<TImage>
::New(enum ExtrapolationMode mode, const BaseImage *image)
{
  const TImage *img = NULL;
  if (image) {
    img = dynamic_cast<const TImage *>(image);
    if (!img) {
      cerr << this->NameOfClass() << "::New(irtkExtrapolationMode): Invalid input image type" << endl;
      exit(1);
    }
  }
  return ExtrapolatorType::New(mode, img);
}

// =============================================================================
// Initialization
// =============================================================================

// -----------------------------------------------------------------------------
template <class TImage>
inline void GenericImageGradientFunction<TImage>::Input(const BaseImage *input)
{
  ImageGradientFunction::Input(dynamic_cast<const TImage *>(input));
  if (input && !this->_Input) {
    cerr << this->NameOfClass() << "::Input: Invalid input image type" << endl;
    exit(1);
  }
}

// -----------------------------------------------------------------------------
template <class TImage>
inline void GenericImageGradientFunction<TImage>::Initialize(bool coeff)
{
  // Ensure that input has the right type
  if (!this->_Input) {
    cerr << this->NameOfClass() << "::Initialize: Missing input image" << endl;
    exit(1);
  } else if (!dynamic_cast<const TImage *>(this->_Input)) {
    cerr << this->NameOfClass() << "::Initialize: Invalid input image type" << endl;
    exit(1);
  }
  // Initialize extrapolator, i.e., infinite discrete image
  if (this->_InfiniteInput) {
    if (!dynamic_cast<const ExtrapolatorType *>(this->_InfiniteInput)) {
      cerr << this->NameOfClass() << "::Initialize: Invalid extrapolator type" << endl;
      exit(1);
    }
  }
  // Initialize base class
  ImageGradientFunction::Initialize(coeff);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline void GenericImageGradientFunction<TImage>
::Extrapolator(ExtrapolateImageFunction *input, bool owner)
{
  ImageGradientFunction::Extrapolator(input, owner);
  if (input && !dynamic_cast<const ExtrapolatorType *>(this->_InfiniteInput)) {
    cerr << this->NameOfClass() << "::Extrapolator: Invalid extrapolator type" << endl;
    exit(1);
  }
}

////////////////////////////////////////////////////////////////////////////////
// Explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
#define mirtkGradientInterpolatorInstantiations(clsname)                       \
  template class clsname<mirtk::BaseImage>;                                    \
  template class clsname<mirtk::GenericImage<mirtk::GreyPixel> >;              \
  template class clsname<mirtk::GenericImage<float> >;                         \
  template class clsname<mirtk::GenericImage<double> >


} // namespace mirtk

#endif // MIRTK_ImageGradientFunction_HXX

////////////////////////////////////////////////////////////////////////////////
// Instantiation
////////////////////////////////////////////////////////////////////////////////

#ifndef MIRTK_ImageGradientFunctionNew_HXX
#define MIRTK_ImageGradientFunctionNew_HXX

#include "mirtk/LinearImageGradientFunction.h"
#include "mirtk/LinearImageGradientFunction2D.h"
#include "mirtk/LinearImageGradientFunction3D.h"

#include "mirtk/FastLinearImageGradientFunction.h"
#include "mirtk/FastLinearImageGradientFunction2D.h"
#include "mirtk/FastLinearImageGradientFunction3D.h"


namespace mirtk {


// -----------------------------------------------------------------------------
template <class TImage>
GenericImageGradientFunction<TImage> *
GenericImageGradientFunction<TImage>
::New(enum InterpolationMode mode, const TImage *image)
{
  mode = InterpolationWithoutPadding(mode);

  GenericImageGradientFunction<TImage> *p = NULL;

  int dim = 0;
  if (image) {
    if      (image->Z() == 1)                            dim = 2;
    else if (image->T() == 1 || image->GetTSize() == .0) dim = 3;
    else                                                 dim = 4;
  }

  switch (dim) {
    case 2: {
      switch (mode) {
        case Interpolation_Linear:
          p = new GenericLinearImageGradientFunction2D<TImage>();
          break;
        case Interpolation_FastLinear:
          p = new GenericFastLinearImageGradientFunction2D<TImage>();
          break;
        default:
          p = NULL;
      }
    }
    case 3: {
      switch (mode) {
        case Interpolation_Linear:
          p = new GenericLinearImageGradientFunction3D<TImage>();
          break;
        case Interpolation_FastLinear:
          p = new GenericFastLinearImageGradientFunction3D<TImage>();
          break;
        default:
          p = NULL;
      }
    }
    default: {
      switch (mode) {
        case Interpolation_Linear:
          p = new GenericLinearImageGradientFunction<TImage>();
          break;
        case Interpolation_FastLinear:
          p = new GenericFastLinearImageGradientFunction<TImage>();
          break;
        default:
          p = NULL;
      }
    }
  }

  // Initialize interpolator
  if (p) {
    p->NumberOfDimensions(dim);
    p->Input(image);
  // Throw error if no suitable interpolator available
  } else {
    cerr << "GenericImageGradientFunction::New: Interpolation mode (" << mode;
    cerr << ") not supported for " << (dim ? ToString(dim) : "N") << "D images" << endl;
    exit(1);
  }

  return p;
}

// -----------------------------------------------------------------------------
template <class TImage>
GenericImageGradientFunction<TImage> *
GenericImageGradientFunction<TImage>
::New(enum InterpolationMode imode, enum ExtrapolationMode emode, const TImage *image)
{
  GenericImageGradientFunction<TImage> *p;
  p = GenericImageGradientFunction<TImage>::New(imode, image);
  if (emode != Extrapolation_Default) p->Extrapolator(p->New(emode, image), true);
  return p;
}


} // namespace mirtk

#endif // MIRTK_ImageGradientFunctionNew_HXX
