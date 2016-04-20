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

#ifndef MIRTK_InterpolateImageFunction_HXX
#define MIRTK_InterpolateImageFunction_HXX

#include "mirtk/InterpolateImageFunction.h"

#include "mirtk/Voxel.h"


namespace mirtk {


////////////////////////////////////////////////////////////////////////////////
// GenericInterpolateImageFunction
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
template <class TImage>
inline GenericInterpolateImageFunction<TImage>
::GenericInterpolateImageFunction()
{
}

// -----------------------------------------------------------------------------
template <class TImage>
inline GenericInterpolateImageFunction<TImage>
::~GenericInterpolateImageFunction()
{
}

// -----------------------------------------------------------------------------
template <class TImage>
ExtrapolateImageFunction *
GenericInterpolateImageFunction<TImage>
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
inline void GenericInterpolateImageFunction<TImage>::Input(const BaseImage *input)
{
  InterpolateImageFunction::Input(dynamic_cast<const TImage *>(input));
  if (input && !this->_Input) {
    cerr << this->NameOfClass() << "::Input: Invalid input image type" << endl;
    exit(1);
  }
}

// -----------------------------------------------------------------------------
template <class TImage>
inline void GenericInterpolateImageFunction<TImage>::Initialize(bool coeff)
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
  InterpolateImageFunction::Initialize(coeff);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline void GenericInterpolateImageFunction<TImage>
::Extrapolator(ExtrapolateImageFunction *input, bool owner)
{
  InterpolateImageFunction::Extrapolator(input, owner);
  if (input && !dynamic_cast<const ExtrapolatorType *>(this->_InfiniteInput)) {
    cerr << this->NameOfClass() << "::Extrapolator: Invalid extrapolator type" << endl;
    exit(1);
  }
}

////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation macro
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
#define mirtkInterpolatorInstantiations(clsname)                               \
  template class clsname<mirtk::BaseImage>;                                    \
  template class clsname<mirtk::GenericImage<mirtk::GreyPixel> >;              \
  template class clsname<mirtk::GenericImage<float> >;                         \
  template class clsname<mirtk::GenericImage<float2> >;                        \
  template class clsname<mirtk::GenericImage<float3> >;                        \
  template class clsname<mirtk::GenericImage<mirtk::Float3> >;                 \
  template class clsname<mirtk::GenericImage<float3x3> >;                      \
  template class clsname<mirtk::GenericImage<double> >;                        \
  template class clsname<mirtk::GenericImage<double2> >;                       \
  template class clsname<mirtk::GenericImage<double3> >;                       \
  template class clsname<mirtk::GenericImage<mirtk::Double3> >;                \
  template class clsname<mirtk::GenericImage<double3x3> >


} // namespace mirtk

#endif // MIRTK_InterpolateImageFunction_HXX

////////////////////////////////////////////////////////////////////////////////
// Instantiation
////////////////////////////////////////////////////////////////////////////////

#ifndef MIRTK_InterpolateImageFunctionNew_HXX
#define MIRTK_InterpolateImageFunctionNew_HXX

// ND
#include "mirtk/NearestNeighborInterpolateImageFunction.h"
#include "mirtk/LinearInterpolateImageFunction.h"
#include "mirtk/BSplineInterpolateImageFunction.h"
#include "mirtk/CubicBSplineInterpolateImageFunction.h"
#include "mirtk/FastCubicBSplineInterpolateImageFunction.h"
#include "mirtk/CSplineInterpolateImageFunction.h"
#include "mirtk/GaussianInterpolateImageFunction.h"
#include "mirtk/SincInterpolateImageFunction.h"

// 2D
#include "mirtk/LinearInterpolateImageFunction2D.h"
#include "mirtk/BSplineInterpolateImageFunction2D.h"
#include "mirtk/CubicBSplineInterpolateImageFunction2D.h"
#include "mirtk/FastCubicBSplineInterpolateImageFunction2D.h"
#include "mirtk/CSplineInterpolateImageFunction2D.h"
#include "mirtk/GaussianInterpolateImageFunction2D.h"
#include "mirtk/SincInterpolateImageFunction2D.h"

// 3D
#include "mirtk/LinearInterpolateImageFunction3D.h"
#include "mirtk/BSplineInterpolateImageFunction3D.h"
#include "mirtk/CubicBSplineInterpolateImageFunction3D.h"
#include "mirtk/FastCubicBSplineInterpolateImageFunction3D.h"
#include "mirtk/CSplineInterpolateImageFunction3D.h"
#include "mirtk/GaussianInterpolateImageFunction3D.h"
#include "mirtk/SincInterpolateImageFunction3D.h"

#include "mirtk/ShapeBasedInterpolateImageFunction.h" // 3D scalar image only

// 4D
#include "mirtk/LinearInterpolateImageFunction4D.h"
#include "mirtk/BSplineInterpolateImageFunction4D.h"
#include "mirtk/CubicBSplineInterpolateImageFunction4D.h"
#include "mirtk/FastCubicBSplineInterpolateImageFunction4D.h"
#include "mirtk/CSplineInterpolateImageFunction4D.h"
#include "mirtk/GaussianInterpolateImageFunction4D.h"
#include "mirtk/SincInterpolateImageFunction4D.h"


namespace mirtk {


// -----------------------------------------------------------------------------
template <class TImage>
GenericInterpolateImageFunction<TImage> *
GenericInterpolateImageFunction<TImage>
::New(enum InterpolationMode mode, const TImage *image)
{
  mode = InterpolationWithoutPadding(mode);

  GenericInterpolateImageFunction<TImage> *p = NULL;

  int dim = 0;
  if (image) {
    if      (image->Z() == 1)                            dim = 2;
    else if (image->T() == 1 || image->GetTSize() == .0) dim = 3;
    else                                                 dim = 4;
  }

  switch (dim) {
    case 2: {
      switch (mode) {
        case Interpolation_NN:
          p = new GenericNearestNeighborInterpolateImageFunction<TImage>();
          break;
        case Interpolation_Linear:
        case Interpolation_FastLinear:
          p = new GenericLinearInterpolateImageFunction2D<TImage>();
          break;
        case Interpolation_BSpline:
          p = new GenericBSplineInterpolateImageFunction2D<TImage>();
          break;
        case Interpolation_CubicBSpline:
          p = new GenericCubicBSplineInterpolateImageFunction2D<TImage>();
          break;
        case Interpolation_FastCubicBSpline:
          p = new GenericFastCubicBSplineInterpolateImageFunction2D<TImage>();
          break;
        case Interpolation_CSpline:
          p = new GenericCSplineInterpolateImageFunction2D<TImage>();
          break;
        case Interpolation_Gaussian:
          p = new GenericGaussianInterpolateImageFunction2D<TImage>();
          break;
        case Interpolation_Sinc:
          p = new GenericSincInterpolateImageFunction2D<TImage>();
          break;
        default:
          p = NULL;
      }
    }
    case 3: {
      switch (mode) {
        case Interpolation_NN:
          p = new GenericNearestNeighborInterpolateImageFunction<TImage>();
          break;
        case Interpolation_Linear:
        case Interpolation_FastLinear:
          p = new GenericLinearInterpolateImageFunction3D<TImage>();
          break;
        case Interpolation_BSpline:
          p = new GenericBSplineInterpolateImageFunction3D<TImage>();
          break;
        case Interpolation_CubicBSpline:
          p = new GenericCubicBSplineInterpolateImageFunction3D<TImage>();
          break;
        case Interpolation_FastCubicBSpline:
          p = new GenericFastCubicBSplineInterpolateImageFunction3D<TImage>();
          break;
        case Interpolation_CSpline:
          p = new GenericCSplineInterpolateImageFunction3D<TImage>();
          break;
        case Interpolation_Gaussian:
          p = new GenericGaussianInterpolateImageFunction3D<TImage>();
          break;
        case Interpolation_Sinc:
          p = new GenericSincInterpolateImageFunction3D<TImage>();
          break;
        default:
          p = NULL;
      }
    }
    case 4: {
      switch (mode) {
        case Interpolation_NN:
          p = new GenericNearestNeighborInterpolateImageFunction<TImage>();
          break;
        case Interpolation_Linear:
        case Interpolation_FastLinear:
          p = new GenericLinearInterpolateImageFunction4D<TImage>();
          break;
        case Interpolation_BSpline:
          p = new GenericBSplineInterpolateImageFunction4D<TImage>();
          break;
        case Interpolation_CubicBSpline:
          p = new GenericCubicBSplineInterpolateImageFunction4D<TImage>();
          break;
        case Interpolation_FastCubicBSpline:
          p = new GenericFastCubicBSplineInterpolateImageFunction4D<TImage>();
          break;
        case Interpolation_CSpline:
          p = new GenericCSplineInterpolateImageFunction4D<TImage>();
          break;
        case Interpolation_Gaussian:
          p = new GenericGaussianInterpolateImageFunction4D<TImage>();
          break;
        case Interpolation_Sinc:
          p = new GenericSincInterpolateImageFunction4D<TImage>();
          break;
        default:
          p = NULL;
      }
    }
    default: {
      switch (mode) {
        case Interpolation_NN:
          p = new GenericNearestNeighborInterpolateImageFunction<TImage>();
          break;
        case Interpolation_Linear:
        case Interpolation_FastLinear:
          p = new GenericLinearInterpolateImageFunction<TImage>();
          break;
        case Interpolation_BSpline:
          p = new GenericBSplineInterpolateImageFunction<TImage>();
          break;
        case Interpolation_CubicBSpline:
          p = new GenericCubicBSplineInterpolateImageFunction<TImage>();
          break;
        case Interpolation_FastCubicBSpline:
          p = new GenericFastCubicBSplineInterpolateImageFunction<TImage>();
          break;
        case Interpolation_CSpline:
          p = new GenericCSplineInterpolateImageFunction<TImage>();
          break;
        case Interpolation_Gaussian:
          p = new GenericGaussianInterpolateImageFunction<TImage>();
          break;
        case Interpolation_Sinc:
          p = new GenericSincInterpolateImageFunction<TImage>();
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
    cerr << "GenericInterpolateImageFunction::New: Interpolation mode (" << mode;
    cerr << ") not supported for " << (dim ? ToString(dim) : "N") << "D images" << endl;
    exit(1);
  }

  return p;
}

// -----------------------------------------------------------------------------
template <class TImage>
GenericInterpolateImageFunction<TImage> *
GenericInterpolateImageFunction<TImage>
::New(enum InterpolationMode imode, enum ExtrapolationMode emode, const TImage *image)
{
  GenericInterpolateImageFunction<TImage> *p;
  p = GenericInterpolateImageFunction<TImage>::New(imode, image);
  if (emode != Extrapolation_Default) p->Extrapolator(p->New(emode, image), true);
  return p;
}


} // namespace mirtk

#endif // MIRTK_InterpolateImageFunctionNew_HXX
