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

#include "mirtk/RegisteredImage.h"

#include "mirtk/Assert.h"
#include "mirtk/Math.h"
#include "mirtk/Memory.h"
#include "mirtk/Parallel.h"
#include "mirtk/Profiling.h"
#include "mirtk/MultiLevelTransformation.h"
#include "mirtk/GaussianBlurring.h"
#include "mirtk/GaussianBlurringWithPadding.h"
#include "mirtk/GradientImageFilter.h"
#include "mirtk/HessianImageFilter.h"
#include "mirtk/VoxelFunction.h"
#include "mirtk/ImageGradientFunction.h"
#include "mirtk/FluidFreeFormTransformation.h"
#include "mirtk/Vector3D.h"

#include "mirtk/LinearInterpolateImageFunction.hxx"  // incl. inline definitions
#include "mirtk/FastLinearImageGradientFunction.hxx" // incl. inline definitions


namespace mirtk {


// -----------------------------------------------------------------------------
RegisteredImage::RegisteredImage()
:
  _InputImage            (NULL),
  _InputGradient         (NULL),
  _InputHessian          (NULL),
  _Transformation        (NULL),
  _InterpolationMode     (Interpolation_Default),
  _ExtrapolationMode     (Extrapolation_Default),
  _WorldCoordinates      (NULL),
  _ImageToWorld          (NULL),
  _ExternalDisplacement  (NULL),
  _FixedDisplacement     (NULL),
  _Displacement          (NULL),
  _CacheWorldCoordinates (true),  // FIXME: MUST be true to also cache anything else...
  _CacheFixedDisplacement(false), // by default, only if required by transformation
  _CacheDisplacement     (false), // (c.f. Transformation::RequiresCachingOfDisplacements)
  _SelfUpdate            (true),
  _MinIntensity          (NaN),
  _MaxIntensity          (NaN),
  _GradientSigma         (.0),
  _HessianSigma          (.0),
  _PrecomputeDerivatives (false),
  _NumberOfActiveLevels  (0),
  _NumberOfPassiveLevels (0)
{
  for (int i = 0; i < 13; ++i) _Offset[i] = -1;
}

// -----------------------------------------------------------------------------
RegisteredImage::RegisteredImage(const RegisteredImage &other)
:
  GenericImage<double>(other),
  _InputImage            (other._InputImage),
  _InputGradient         (other._InputGradient  ? new GradientImageType(*other._InputGradient)  : NULL),
  _InputHessian          (other._InputHessian   ? new GradientImageType(*other._InputHessian)   : NULL),
  _Transformation        (other._Transformation),
  _InterpolationMode     (other._InterpolationMode),
  _ExtrapolationMode     (other._ExtrapolationMode),
  _WorldCoordinates      (other._WorldCoordinates),
  _ImageToWorld          (other._ImageToWorld      ? new WorldCoordsImage (*other._ImageToWorld)      : NULL),
  _ExternalDisplacement  (other._ExternalDisplacement),
  _FixedDisplacement     (other._FixedDisplacement ? new DisplacementImageType(*other._FixedDisplacement) : NULL),
  _Displacement          (other._Displacement      ? new DisplacementImageType(*other._Displacement)      : NULL),
  _CacheWorldCoordinates (other._CacheWorldCoordinates),
  _CacheFixedDisplacement(other._CacheFixedDisplacement),
  _CacheDisplacement     (other._CacheDisplacement),
  _SelfUpdate            (other._SelfUpdate),
  _MinIntensity          (other._MinIntensity),
  _MaxIntensity          (other._MaxIntensity),
  _GradientSigma         (other._GradientSigma),
  _HessianSigma          (other._HessianSigma),
  _PrecomputeDerivatives (other._PrecomputeDerivatives),
  _NumberOfActiveLevels  (other._NumberOfActiveLevels),
  _NumberOfPassiveLevels (other._NumberOfPassiveLevels)
{
  memcpy(_Offset, other._Offset, 13 * sizeof(int));
}

// -----------------------------------------------------------------------------
RegisteredImage &RegisteredImage::operator =(const RegisteredImage &other)
{
  GenericImage<double>::operator =(other);
  _InputImage             = other._InputImage;
  _InputGradient          = other._InputGradient  ? new GradientImageType(*other._InputGradient)  : NULL;
  _InputHessian           = other._InputHessian   ? new GradientImageType(*other._InputHessian)   : NULL;
  _Transformation         = other._Transformation;
  _InterpolationMode      = other._InterpolationMode;
  _ExtrapolationMode      = other._ExtrapolationMode;
  _WorldCoordinates       = other._WorldCoordinates;
  _ImageToWorld           = other._ImageToWorld      ? new WorldCoordsImage (*other._ImageToWorld)      : NULL;
  _ExternalDisplacement   = other._ExternalDisplacement;
  _FixedDisplacement      = other._FixedDisplacement ? new DisplacementImageType(*other._FixedDisplacement) : NULL;
  _Displacement           = other._Displacement      ? new DisplacementImageType(*other._Displacement)      : NULL;
  _CacheWorldCoordinates  = other._CacheWorldCoordinates;
  _CacheFixedDisplacement = other._CacheFixedDisplacement;
  _CacheDisplacement      = other._CacheDisplacement;
  _SelfUpdate             = other._SelfUpdate;
  _MinIntensity           = other._MinIntensity;
  _MaxIntensity           = other._MaxIntensity;
  _GradientSigma          = other._GradientSigma;
  _HessianSigma           = other._HessianSigma;
  _PrecomputeDerivatives  = other._PrecomputeDerivatives;
  _NumberOfActiveLevels   = other._NumberOfActiveLevels;
  _NumberOfPassiveLevels  = other._NumberOfPassiveLevels;
  memcpy(_Offset, other._Offset, 13 * sizeof(int));
  return *this;
}

// -----------------------------------------------------------------------------
RegisteredImage::~RegisteredImage()
{
  if (_ImageToWorld != _WorldCoordinates) delete _ImageToWorld;
  delete _FixedDisplacement;
  delete _Displacement;
  if (_InputGradient != _InputImage) delete _InputGradient;
  delete _InputHessian;
}
// =============================================================================
// Initialization
// =============================================================================

// -----------------------------------------------------------------------------
void RegisteredImage::Initialize(const ImageAttributes &attr, int t)
{
  MIRTK_START_TIMING();

  // Destroy previous interpolators
  _Evaluator.reset();

  // Clear possibly previously allocated displacement cache
  if (!_Transformation) Delete(_Displacement);

  // Check if input image is set
  if (!_InputImage) {
    cerr << "RegisteredImage::Initialize: Missing input image" << endl;
    exit(1);
  }

  // Initialize base class
  if (attr._t > 1) {
    cerr << "RegisteredImage::Initialize: Split multi-channel image/temporal sequence up into separate 2D/3D images" << endl;
    exit(1);
  }
  if (t == 0) t = 1;
  if (t != 1 && t != 4 && t != 10 && t != 13) {
    cerr << "RegisteredImage::Initialize: Number of registered image channels must be either 1, 4, 10 or 13" << endl;
    exit(1);
  }
  GenericImage<double>::Initialize(attr, t);

  // Set background value/foreground mask
  if (_InputImage->HasBackgroundValue()) {
    this->PutBackgroundValueAsDouble(_InputImage->GetBackgroundValueAsDouble());
  } else {
    double min, max;
    _InputImage->GetMinMaxAsDouble(&min, &max);
    this->PutBackgroundValueAsDouble(min - 1.);
  }

  // Pre-compute world coordinates
  if (_WorldCoordinates) {
    if (_ImageToWorld != _WorldCoordinates) {
      delete _ImageToWorld;
      _ImageToWorld = _WorldCoordinates;
    }
  } else if (_CacheWorldCoordinates) {
    if (!_ImageToWorld) _ImageToWorld = new WorldCoordsImage();
    this->ImageToWorld(*_ImageToWorld, true /* i.e., always 3D vectors */);
  } else {
    Delete(_ImageToWorld);
  }

  // Determine number of active (changing) and passive (fixed) levels
  bool cache_fixed = !_ExternalDisplacement && _CacheFixedDisplacement;
  const MultiLevelTransformation *mffd = NULL;
  if ((mffd = dynamic_cast<const MultiLevelTransformation *>(_Transformation))) {
    _NumberOfPassiveLevels = 0;
    for (int l = 0; l < mffd->NumberOfLevels(); ++l) {
      if (mffd->LocalTransformationIsActive(l)) break;
      if (mffd->GetLocalTransformation(l)->RequiresCachingOfDisplacements()) cache_fixed = true;
      ++_NumberOfPassiveLevels;
    }
    _NumberOfActiveLevels = mffd->NumberOfLevels() - _NumberOfPassiveLevels;
    if (_NumberOfPassiveLevels == 0 && mffd->GetGlobalTransformation()->IsIdentity()) {
      _NumberOfPassiveLevels = -1;
    }
  } else if (_Transformation) {
    _NumberOfActiveLevels  =  1;
    _NumberOfPassiveLevels = -1;
  } else {
    _NumberOfActiveLevels  =  0;
    _NumberOfPassiveLevels = -1;
  }

  // Pre-compute fixed displacements
  if (cache_fixed && _NumberOfPassiveLevels >= 0) {
    if (!_FixedDisplacement) _FixedDisplacement = new DisplacementImageType();
    _FixedDisplacement->Initialize(attr, 3);
    mffd->Displacement(-1, _NumberOfPassiveLevels, *_FixedDisplacement, _InputImage->GetTOrigin(), _ImageToWorld);
  } else {
    Delete(_FixedDisplacement);
  }

  // Pre-compute input derivatives
  if (!_PrecomputeDerivatives) {
    enum InterpolationMode interp;
    if (_InterpolationMode == Interpolation_Default) {
      interp = DefaultInterpolationMode();
    } else {
      interp = InterpolationWithoutPadding(_InterpolationMode);
    }
    if (interp != Interpolation_Linear && interp != Interpolation_FastLinear) {
      _PrecomputeDerivatives = true;
    }
  }
  if (t > 1) ComputeInputGradient(_GradientSigma);
  if (t > 4) ComputeInputHessian (_HessianSigma );

  // Initialize offsets of registered image channels
  _Offset[0] = 0;
  _Offset[1] = this->NumberOfVoxels();
  for (int c = 2; c < 13; ++c) {
    _Offset[c] = _Offset[c-1] + _Offset[1];
  }

  // Attention: Initialization of actual image content must be forced upon first
  //            Update call. This is initiated by the ImageSimilarity::Update
  //            function which in turn is called before the first energy gradient
  //            evaluation (see GradientDescent::Gradient). Applications can
  //            call the Recompute function after Initialize to complete the
  //            initialization of the image.

  MIRTK_DEBUG_TIMING(4, "initialization of " << (_Transformation ? "moving" : "fixed") << " image");
}

// -----------------------------------------------------------------------------
void RegisteredImage::ComputeInputGradient(double sigma)
{
  MIRTK_START_TIMING();
  // Smooth input image
  InputImageType *blurred_image = _InputImage;
  if (sigma > .0) {
    blurred_image = new InputImageType;
    if (this->HasBackgroundValue()) {
      blurred_image->PutBackgroundValueAsDouble(this->GetBackgroundValueAsDouble());
      GaussianBlurringWithPadding<double> blurring(sigma * _InputImage->GetXSize(),
                                                   sigma * _InputImage->GetYSize(),
                                                   sigma * _InputImage->GetZSize(),
                                                   this->GetBackgroundValueAsDouble());
      blurring.Input (_InputImage);
      blurring.Output(blurred_image);
      blurring.Run();
    } else {
      GaussianBlurring<double> blurring(sigma * _InputImage->GetXSize(),
                                        sigma * _InputImage->GetYSize(),
                                        sigma * _InputImage->GetZSize());
      blurring.Input (_InputImage);
      blurring.Output(blurred_image);
      blurring.Run();
    }
  }
  if (_PrecomputeDerivatives) {
    // Compute image gradient using finite differences
    typedef GradientImageFilter<GradientImageType::VoxelType> FilterType;
    FilterType filter(FilterType::GRADIENT_VECTOR);
    filter.Input (blurred_image);
    filter.Output(_InputGradient ? _InputGradient : new GradientImageType);
    // Note that even though the original IRTK nreg2 implementation did divide
    // the image gradient initially by the voxel size, the similarity gradient
    // was reoriented then by ImageRegistration2::EvaluateGradient using the
    // upper 3x3 image to world matrix. This effectively multiplied by the voxel
    // size again which is equivalent to only reorienting the image gradient
    // computed w.r.t. the voxel coordinates, i.e., leaving the magnitude of the
    // gradient in voxel units rather than world units (i.e., mm's).
    filter.UseVoxelSize  (true);
    filter.UseOrientation(true);
    if (this->HasBackgroundValue()) {
      filter.PaddingValue(this->GetBackgroundValueAsDouble());
    }
    filter.Run();
    _InputGradient = filter.Output();
    _InputGradient->PutTSize(.0);
    _InputGradient->PutBackgroundValueAsDouble(.0);
    if (blurred_image != _InputImage) delete blurred_image;
    MIRTK_DEBUG_TIMING(5, "computation of 1st order image derivatives");
  } else {
    delete _InputGradient;
    _InputGradient = blurred_image;
    MIRTK_DEBUG_TIMING(5, "low-pass filtering of image for 1st order derivatives");
  }
}

// -----------------------------------------------------------------------------
void RegisteredImage::ComputeInputHessian(double sigma)
{
  MIRTK_START_TIMING();
  // Smooth input image
  InputImageType *blurred_image = _InputImage;
  if (sigma > .0) {
    blurred_image = new InputImageType;
    if (this->HasBackgroundValue()) {
      blurred_image->PutBackgroundValueAsDouble(this->GetBackgroundValueAsDouble());
      GaussianBlurringWithPadding<double> blurring(sigma * _InputImage->GetXSize(),
                                                   sigma * _InputImage->GetYSize(),
                                                   sigma * _InputImage->GetZSize(),
                                                   this->GetBackgroundValueAsDouble());
      blurring.Input (_InputImage);
      blurring.Output(blurred_image);
      blurring.Run();
    } else {
      GaussianBlurring<double> blurring(sigma * _InputImage->GetXSize(),
                                        sigma * _InputImage->GetYSize(),
                                        sigma * _InputImage->GetZSize());
      blurring.Input (_InputImage);
      blurring.Output(blurred_image);
      blurring.Run();
    }
  }
  // Compute 2nd order image derivatives using finite differences
  typedef HessianImageFilter<HessianImageType::VoxelType> FilterType;
  FilterType filter(FilterType::HESSIAN_MATRIX);
  filter.Input(blurred_image);
  filter.Output(_InputHessian ? _InputHessian : new HessianImageType);
  filter.UseVoxelSize  (true);
  filter.UseOrientation(true);
  if (this->HasBackgroundValue()) {
    filter.PaddingValue(this->GetBackgroundValueAsDouble());
  }
  filter.Run();
  _InputHessian = filter.Output();
  _InputHessian->PutTSize(.0);
  _InputHessian->PutBackgroundValueAsDouble(.0);
  if (blurred_image != _InputImage) delete blurred_image;
  MIRTK_DEBUG_TIMING(5, "computation of 2nd order image derivatives");
}

// =============================================================================
// Internally used voxel-wise update function
// =============================================================================

// -----------------------------------------------------------------------------
// Base class of voxel transformation functors
struct Transformer
{
  typedef WorldCoordsImage::VoxelType CoordType;

  /// Constructor
  Transformer()
  :
    _Input(nullptr),
    _Transformation(nullptr),
    _Output(nullptr),
    _y(0), _z(0)
  {}

  /// Initialize data members
  void Initialize(const BaseImage *i, const Transformation *t, RegisteredImage *o)
  {
    _Input          = i;
    _t              = i->GetTOrigin();
    _t0             = o->GetTOrigin();
    _Transformation = t;
    _Output         = o;

    _y = o->X() * o->Y() * o->Z();
    _z = 2 * _y;
  }

  /// Transform output voxel
  void operator ()(double &x, double &y, double &z)
  {
    _Output->ImageToWorld(x, y, z);
    _Transformation->Transform(x, y, z, _t, _t0);
    _Input->WorldToImage(x, y, z);
  }

  /// Transform output voxel using pre-computed world coordinates
  void operator ()(double &x, double &y, double &z, const CoordType *wc)
  {
    x = wc[_x], y = wc[_y], z = wc[_z];
    _Transformation->Transform(x, y, z, _t, _t0);
    _Input->WorldToImage(x, y, z);
  }

  /// Transform output voxel using pre-computed world coordinates and displacements
  void operator ()(double &x, double &y, double &z, const CoordType *wc, const double *dx)
  {
    x = wc[_x] + dx[_x];
    y = wc[_y] + dx[_y];
    z = wc[_z] + dx[_z];
    _Input->WorldToImage(x, y, z);
  }

protected:

  const BaseImage      *_Input;
  const Transformation *_Transformation;
  RegisteredImage      *_Output;

  static const int _x = 0; ///< Offset of x component
  int              _y;     ///< Offset of y component
  int              _z;     ///< Offset of z component

  double _t;  ///< Time point
  double _t0; ///< Time point of target
};

// -----------------------------------------------------------------------------
// Transformer used when no fixed transformation is cached
struct DefaultTransformer : public Transformer
{
  using Transformer::operator();

  /// As this transformer is only used when no fixed transformation is cached,
  /// this overloaded operator should never be invoked
  void operator ()(double &, double &, double &, const CoordType *, const double *, const double *)
  {
    cerr << "RegisteredImage::DefaultTransformer used even though _FixedDisplacement assumed to be NULL ?!?" << endl;
    exit(1);
  }
};

// -----------------------------------------------------------------------------
// Transform output voxel using additive composition of displacements
struct AdditiveTransformer : public Transformer
{
  using Transformer::operator();

  /// Transform output voxel using pre-computed world coordinates and displacements
  void operator ()(double &x, double &y, double &z, const CoordType *wc, const double *d1, const double *d2)
  {
    x = wc[_x] + d1[_x] + d2[_x];
    y = wc[_y] + d1[_y] + d2[_y];
    z = wc[_z] + d1[_z] + d2[_z];
    _Input->WorldToImage(x, y, z);
  }
};

// -----------------------------------------------------------------------------
// Transform output voxel using fluid composition of displacements
struct FluidTransformer : public Transformer
{
  using Transformer::operator();

  /// Transform output voxel using pre-computed world coordinates and displacements
  ///
  /// Because fluid composition of displacement fields would require interpolation,
  /// let Transformation::Displacement handle the fluid composition already when
  /// computing the second displacement field.
  void operator ()(double &x, double &y, double &z, const CoordType *wc, const double *, const double *dx)
  {
    x = wc[_x] + dx[_x];
    y = wc[_y] + dx[_y];
    z = wc[_z] + dx[_z];
    _Input->WorldToImage(x, y, z);
  }
};

// -----------------------------------------------------------------------------
// Transformer used when no transformation is set or custom displacement field given
struct FixedTransformer : public Transformer
{
  // Visual Studio 2013 has troubles resolving
  //   void operator()(double&, double&, double&)
  // if a using Transformer::operator() statement is used because it is also
  // defined by this subclass. Instead, just re-implement the only other
  // overloaded version as well.

  /// Transform output voxel
  void operator ()(double &x, double &y, double &z)
  {
    _Output->ImageToWorld(x, y, z);
    _Input ->WorldToImage(x, y, z);
  }

  /// Transform output voxel using pre-computed world coordinates
  void operator ()(double &x, double &y, double &z, const CoordType *wc)
  {
    x = wc[_x], y = wc[_y], z = wc[_z];
    _Input->WorldToImage(x, y, z);
  }

  /// Transform output voxel using pre-computed world coordinates and displacements
  void operator ()(double &x, double &y, double &z, const CoordType *wc, const double *dx)
  {
    Transformer::operator()(x, y, z, wc, dx);
  }

  /// As this transformer is only used when no transformation is set,
  /// this overloaded operator should never be invoked
  void operator ()(double &, double &, double &, const CoordType *, const double *, const double *)
  {
    cerr << "RegisteredImage::FixedTransformer(..., d1, d2) used even though _Transformation assumed to be NULL ?!?" << endl;
    exit(1);
  }
};

// -----------------------------------------------------------------------------
// Auxiliary function to allocate and initialize image interpolate function
template <class ImageFunction>
void New(
  ImageFunction *&f, const BaseImage *image,
#ifndef NDEBUG
  InterpolationMode interp,
#else
  InterpolationMode,
#endif
  ExtrapolationMode extrap, double padding, double default_value)
{
  if (image) {
    f = new ImageFunction();
#ifndef NDEBUG
    if (interp == Interpolation_Default) {
      interp = DefaultInterpolationMode();
    } else {
      interp = InterpolationWithoutPadding(interp);
    }
    if (interp == Interpolation_FastLinear) interp = Interpolation_Linear;
    InterpolationMode mode = f->InterpolationMode();
    if (mode == Interpolation_FastLinear) mode = Interpolation_Linear;
    if (mode != interp) {
      cout << endl;
      cerr << __FILE__ << ":" << __LINE__ << ": Mismatch of interpolation mode: expected \""
                << ToString(interp) << "\", but got \"" << ToString(mode) << "\"" << endl;
      exit(1);
    }
#endif
    ImageGradientFunction *g = dynamic_cast<ImageGradientFunction *>(f);
    if (g) g->WrtWorld(true);
    f->Input(const_cast<BaseImage *>(image));
    if (extrap != Extrapolation_Default) {
      f->Extrapolator(f->New(extrap, image), true);
    }
    f->DefaultValue(default_value);
    f->Initialize();
    if (f->Extrapolator()) f->Extrapolator()->DefaultValue(padding);
  }
}

template <>
void New(InterpolateImageFunction *&f, const BaseImage *image,
         InterpolationMode interp, ExtrapolationMode extrap,
         double padding, double default_value)
{
  if (image) {
    f = InterpolateImageFunction::New(interp, const_cast<BaseImage *>(image));
    f->Input(const_cast<BaseImage *>(image));
    if (extrap != Extrapolation_Default) {
      f->Extrapolator(f->New(extrap, image), true);
    }
    f->DefaultValue(default_value);
    f->Initialize();
    if (f->Extrapolator()) f->Extrapolator()->DefaultValue(padding);
  }
}

template <>
void New(ImageGradientFunction *&f, const BaseImage *image,
         InterpolationMode interp, ExtrapolationMode extrap,
         double padding, double default_value)
{
  if (image) {
    f = ImageGradientFunction::New(interp, const_cast<BaseImage *>(image));
    f->WrtWorld(true);
    f->Input(const_cast<BaseImage *>(image));
    if (extrap != Extrapolation_Default) {
      f->Extrapolator(f->New(extrap, image), true);
    }
    f->DefaultValue(default_value);
    f->Initialize();
    if (f->Extrapolator()) f->Extrapolator()->DefaultValue(padding);
  }
}

// TODO: Add template specialization for ImageHessianFunction

// -----------------------------------------------------------------------------
// Auxiliary evaluator functions
template <class TIntensityFunction, class TGradientFunction, class THessianFunction>
class Interpolator : public Object
{
  mirtkObjectMacro(RegisteredImageInterpolator);

public:

  typedef TIntensityFunction  IntensityFunction;
  typedef TGradientFunction   GradientFunction;
  typedef THessianFunction    HessianFunction;

protected:

  IntensityFunction *_IntensityFunction;
  GradientFunction  *_GradientFunction;
  HessianFunction   *_HessianFunction;
  bool               _InterpolateWithPadding;
  bool               _PrecomputedDerivatives;
  double             _PaddingValue;
  double             _MinIntensity;
  double             _MaxIntensity;
  double             _RescaleSlope;
  double             _RescaleIntercept;
  int                _NumberOfVoxels;
  int                _NumberOfChannels;
  Vector3D<int>      _InputSize;

public:

  /// Constructor
  Interpolator(RegisteredImage *o = nullptr)
  :
    _IntensityFunction     (nullptr),
    _GradientFunction      (nullptr),
    _HessianFunction       (nullptr),
    _InterpolateWithPadding(false),
    _PrecomputedDerivatives(false),
    _PaddingValue          (-1.),
    _MinIntensity          (NaN),
    _MaxIntensity          (NaN),
    _RescaleSlope          (1.),
    _RescaleIntercept      (.0),
    _NumberOfVoxels        (0),
    _NumberOfChannels      (0)
  {
    if (o) Initialize(o);
  }

  /// Copy constructor
  Interpolator(const Interpolator &other)
  :
    _IntensityFunction     (nullptr),
    _GradientFunction      (nullptr),
    _HessianFunction       (nullptr),
    _InterpolateWithPadding(other._InterpolateWithPadding),
    _PrecomputedDerivatives(other._PrecomputedDerivatives),
    _PaddingValue          (other._PaddingValue),
    _MinIntensity          (other._MinIntensity),
    _MaxIntensity          (other._MaxIntensity),
    _RescaleSlope          (other._RescaleSlope),
    _RescaleIntercept      (other._RescaleIntercept),
    _NumberOfVoxels        (other._NumberOfVoxels),
    _NumberOfChannels      (other._NumberOfChannels),
    _InputSize             (other._InputSize)
  {
    if (other._IntensityFunction) {
      const BaseImage *f = other._IntensityFunction->Input();
      const double f_bg = (f->HasBackgroundValue() ? f->GetBackgroundValueAsDouble() : MIN_GREY);
      New<IntensityFunction>(_IntensityFunction, f,
                             other._IntensityFunction->InterpolationMode(),
                             other._IntensityFunction->ExtrapolationMode(),
                             f_bg, f_bg);
    }
    if (other._GradientFunction) {
      const BaseImage *g = other._GradientFunction->Input();
      const double g_bg = (g->HasBackgroundValue() ? g->GetBackgroundValueAsDouble() : .0);
      New<GradientFunction>(_GradientFunction, g,
                            other._GradientFunction->InterpolationMode(),
                            other._GradientFunction->ExtrapolationMode(),
                            g_bg, .0);
    }
    if (other._HessianFunction) {
      const BaseImage *h = other._HessianFunction->Input();
      const double h_bg = (h->HasBackgroundValue() ? h->GetBackgroundValueAsDouble() : .0);
      New<HessianFunction>(_HessianFunction, h,
                           other._HessianFunction->InterpolationMode(),
                           other._HessianFunction->ExtrapolationMode(),
                           h_bg, .0);
    }
  }

  /// Destructor
  ~Interpolator()
  {
    Delete(_IntensityFunction);
    Delete(_GradientFunction);
    Delete(_HessianFunction);
  }

  /// Initialize interpolators
  void Initialize(RegisteredImage *o)
  {
    InterpolationMode interp;
    if (o->InterpolationMode() == Interpolation_Default) {
      interp = DefaultInterpolationMode();
    } else {
      interp = o->InterpolationMode();
    }
    _PrecomputedDerivatives = o->PrecomputeDerivatives();
    _InterpolateWithPadding = (interp != InterpolationWithoutPadding(interp));
    if (o->HasBackgroundValue()) _PaddingValue = o->GetBackgroundValueAsDouble();
    _NumberOfVoxels   = o->X() * o->Y() * o->Z();
    _NumberOfChannels = o->T();
    const BaseImage * const f = o->InputImage();
    const BaseImage * const g = o->InputGradient();
    const BaseImage * const h = o->InputHessian();
    const double f_bg = (f->HasBackgroundValue() ? f->GetBackgroundValueAsDouble() : MIN_GREY);
    const double g_bg = (g && g->HasBackgroundValue() ? g->GetBackgroundValueAsDouble() : .0);
    const double h_bg = (h && h->HasBackgroundValue() ? h->GetBackgroundValueAsDouble() : .0);
    New<IntensityFunction>(_IntensityFunction, f, interp, o->ExtrapolationMode(), f_bg, f_bg);
    New<GradientFunction >(_GradientFunction,  g, interp, Extrapolation_Default,  g_bg, .0);
    New<HessianFunction  >(_HessianFunction,   h, interp, Extrapolation_Default,  h_bg, .0);
    f->GetMinMaxAsDouble(_MinIntensity, _MaxIntensity);
    double omin = o->MinIntensity();
    double omax = o->MaxIntensity();
    if (!IsNaN(omin) || !IsNaN(omax)) {
      if (IsNaN(omin)) omin = _MinIntensity;
      if (IsNaN(omax)) omax = _MaxIntensity;
      _RescaleSlope     = (omax - omin) / (_MaxIntensity - _MinIntensity);
      _RescaleIntercept = omin - _RescaleSlope * _MinIntensity;
      _MinIntensity     = omin;
      _MaxIntensity     = omax;
    } else {
      _RescaleSlope     = 1.0;
      _RescaleIntercept = 0.0;
      if (f->HasBackgroundValue() && _MinIntensity > f->GetBackgroundValueAsDouble()) {
        _MinIntensity = f->GetBackgroundValueAsDouble();
      }
    }
    _InputSize = Vector3D<int>(f->X(), f->Y(), f->Z());
  }

  /// Get input image
  const BaseImage *Input() const
  {
    return _IntensityFunction->Input();
  }

  /// Determine interpolation mode at given location
  ///
  /// \retval  1 Output channels can be interpolated without boundary checks.
  /// \retval  0 Output channels must be interpolated with boundary checks.
  /// \retval -1 Output channels should be padded.
  int Mode(double x, double y, double z, bool check_value = true) const
  {
    bool inside;
    if (_InputSize._z == 1) {
      inside = (abs(z) < 1e-3) && _IntensityFunction->IsInside(x, y) &&
               (!_GradientFunction || _GradientFunction->IsInside(x, y)) &&
               (!_HessianFunction  || _HessianFunction ->IsInside(x, y));
    } else {
      inside = _IntensityFunction->IsInside(x, y, z) &&
               (!_GradientFunction || _GradientFunction->IsInside(x, y, z)) &&
               (!_HessianFunction  || _HessianFunction ->IsInside(x, y, z));
    }
    if (inside && check_value && _InterpolateWithPadding) {
      double value = _IntensityFunction->EvaluateWithPadding(x, y, z);
      if (value == _IntensityFunction->DefaultValue()) inside = false;
    }
    if (inside) return 1;
    if (0. <= x && x <= _InputSize._x - 1. &&
        0. <= y && y <= _InputSize._y - 1. &&
        0. <= z && z <= _InputSize._z - 1.) return 0;
    return -1;
  }

  /// Interpolate input intensity function
  ///
  /// \return The interpolation mode, i.e., result of inside/outside domain check.
  int Intensity(double x, double y, double z, double *o) const
  {
    // Check if location is inside image domain
    int mode = Mode(x, y, z, false);
    if (mode != -1) {
      // Either interpolate using the input padding value to exclude background
      if (_InterpolateWithPadding) {
        if (mode == 1) {
          *o = _IntensityFunction->EvaluateWithPaddingInside(x, y, z);
        } else {
          *o = _IntensityFunction->EvaluateWithPaddingOutside(x, y, z);
        }
      // or simply ignore the input background value as done by nreg2
      } else {
        if (mode == 1) {
          *o = _IntensityFunction->EvaluateInside(x, y, z);
        } else {
          *o = _IntensityFunction->EvaluateOutside(x, y, z);
        }
      }
      // Set background to output padding value
      if (*o == _IntensityFunction->DefaultValue()) {
        *o = _PaddingValue;
        if (_InterpolateWithPadding) return -1;
      // Rescale/clamp foreground to [min, max] range
      } else {
        *o = (*o) * _RescaleSlope + _RescaleIntercept;
        if      (*o < _MinIntensity) *o = _MinIntensity;
        else if (*o > _MaxIntensity) *o = _MaxIntensity;
      }
    // Otherwise, set output intensity to outside value
    } else {
      *o = _PaddingValue;
    }
    // Pass inside/outside check result on to derivative interpolation
    // functions such that these boundary checks are only done once.
    // This requires the same interpolation mode for all channels.
    return mode;
  }

  /// Interpolate 1st order derivatives of input intensity function
  void Gradient(double x, double y, double z, double *o, int mode = 0) const
  {
    o += _NumberOfVoxels;
    switch (mode) {
      // Inside
      case 1:
        if (_InterpolateWithPadding) {
          _GradientFunction->EvaluateWithPaddingInside(o, x, y, z, _NumberOfVoxels);
        } else {
          _GradientFunction->EvaluateInside(o, x, y, z, _NumberOfVoxels);
        }
        break;
      // Boundary
      case 0:
        if (_PrecomputedDerivatives) {
          if (_InterpolateWithPadding) {
            _GradientFunction->EvaluateWithPaddingOutside(o, x, y, z, _NumberOfVoxels);
          } else {
            _GradientFunction->EvaluateOutside(o, x, y, z, _NumberOfVoxels);
          }
          break;
        } // otherwise continue with Outside case
      // Outside
      default: {
        for (int c = 1; c <= 3; ++c, o += _NumberOfVoxels) *o = .0;
      } break;
    }
  }

  /// Interpolate 2nd order derivatives of input intensity function
  void Hessian(double x, double y, double z, double *o, int mode = 0) const
  {
    o += 4 * _NumberOfVoxels;
    switch (mode) {
      // Inside
      case 1:
        if (_InterpolateWithPadding) {
          _HessianFunction->EvaluateWithPaddingInside(o, x, y, z, _NumberOfVoxels);
        } else {
          _HessianFunction->EvaluateInside(o, x, y, z, _NumberOfVoxels);
        }
        break;
      // Boundary
      case 0:
        if (_PrecomputedDerivatives) {
          if (_InterpolateWithPadding) {
            _HessianFunction->EvaluateWithPaddingOutside(o, x, y, z, _NumberOfVoxels);
          } else {
            _HessianFunction->EvaluateOutside(o, x, y, z, _NumberOfVoxels);
          }
          break;
        } // otherwise continue with Outside case
      // Outside
      default: {
        for (int c = 4; c < _NumberOfChannels; ++c, o += _NumberOfVoxels) *o = .0;
      } break;
    }
  }
};

// -----------------------------------------------------------------------------
// Interpolates intensity
template <class IntensityFunction, class GradientFunction, class HessianFunction>
struct IntensityInterpolator
{
  typedef Interpolator<IntensityFunction, GradientFunction, HessianFunction> Evaluator;

  const Evaluator *_Evaluate;
  IntensityInterpolator(const Evaluator *eval = nullptr)
  :
    _Evaluate(eval)
  {}

  void operator()(double x, double y, double z, double *o)
  {
    _Evaluate->Intensity(x, y, z, o);
  }
};

// -----------------------------------------------------------------------------
// Interpolates 1st order derivatives
template <class IntensityFunction, class GradientFunction, class HessianFunction>
struct GradientInterpolator
{
  typedef Interpolator<IntensityFunction, GradientFunction, HessianFunction> Evaluator;

  const Evaluator *_Evaluate;
  GradientInterpolator(const Evaluator *eval = nullptr)
  :
    _Evaluate(eval)
  {}

  void operator()(double x, double y, double z, double *o)
  {
    int mode = _Evaluate->Mode(x, y, z);
    _Evaluate->Gradient(x, y, z, o, mode);
  }
};

// -----------------------------------------------------------------------------
// Interpolates 2nd order derivatives
template <class IntensityFunction, class GradientFunction, class HessianFunction>
struct HessianInterpolator
{
  typedef Interpolator<IntensityFunction, GradientFunction, HessianFunction> Evaluator;

  const Evaluator *_Evaluate;
  HessianInterpolator(const Evaluator *eval = nullptr)
  :
    _Evaluate(eval)
  {}

  void operator()(double x, double y, double z, double *o)
  {
    int mode = _Evaluate->Mode(x, y, z);
    _Evaluate->Hessian(x, y, z, o, mode);
  }
};

// -----------------------------------------------------------------------------
// Interpolates intensity and 1st order derivatives
template <class IntensityFunction, class GradientFunction, class HessianFunction>
struct IntensityAndGradientInterpolator
{
  typedef Interpolator<IntensityFunction, GradientFunction, HessianFunction> Evaluator;

  const Evaluator *_Evaluate;
  IntensityAndGradientInterpolator(const Evaluator *eval = nullptr)
  :
    _Evaluate(eval)
  {}

  void operator()(double x, double y, double z, double *o)
  {
    int mode = _Evaluate->Intensity(x, y, z, o);
    _Evaluate->Gradient(x, y, z, o, mode);
  }
};

// -----------------------------------------------------------------------------
// Interpolates intensity and 2nd order derivatives
template <class IntensityFunction, class GradientFunction, class HessianFunction>
struct IntensityAndHessianInterpolator
{
  typedef Interpolator<IntensityFunction, GradientFunction, HessianFunction> Evaluator;

  const Evaluator *_Evaluate;
  IntensityAndHessianInterpolator(const Evaluator *eval = nullptr)
  :
    _Evaluate(eval)
  {}

  void operator()(double x, double y, double z, double *o)
  {
    int mode = _Evaluate->Intensity(x, y, z, o);
    _Evaluate->Hessian(x, y, z, o, mode);
  }
};

// -----------------------------------------------------------------------------
// Interpolates 1st and 2nd order derivatives
template <class IntensityFunction, class GradientFunction, class HessianFunction>
struct GradientAndHessianInterpolator
{
  typedef Interpolator<IntensityFunction, GradientFunction, HessianFunction> Evaluator;

  const Evaluator *_Evaluate;
  GradientAndHessianInterpolator(const Evaluator *eval = nullptr)
  :
    _Evaluate(eval)
  {}

  void operator()(double x, double y, double z, double *o)
  {
    int mode = _Evaluate->Mode(x, y, z);
    _Evaluate->Gradient(x, y, z, o, mode);
    _Evaluate->Hessian(x, y, z, o, mode);
  }
};

// -----------------------------------------------------------------------------
// Interpolates intensity, 1st and 2nd order derivatives
template <class IntensityFunction, class GradientFunction, class HessianFunction>
struct IntensityAndGradientAndHessianInterpolator
{
  typedef Interpolator<IntensityFunction, GradientFunction, HessianFunction> Evaluator;

  const Evaluator *_Evaluate;
  IntensityAndGradientAndHessianInterpolator(const Evaluator *eval = nullptr)
  :
    _Evaluate(eval)
  {}

  void operator()(double x, double y, double z, double *o)
  {
    int mode = _Evaluate->Intensity(x, y, z, o);
    _Evaluate->Gradient(x, y, z, o, mode);
    _Evaluate->Hessian(x, y, z, o, mode);
  }
};

// -----------------------------------------------------------------------------
// Voxel update function
template <class Transformer, class Interpolator>
struct UpdateVoxelFunction : public VoxelFunction
{
private:

  typedef typename Transformer::CoordType  CoordType;
  typedef typename Interpolator::Evaluator Evaluator;

  Transformer  _Transform;
  Interpolator _Interpolate;

public:

  /// Constructor
  UpdateVoxelFunction(const Evaluator *eval, const Transformation *t, RegisteredImage *o)
  :
    _Interpolate(eval)
  {
    _Transform.Initialize(eval->Input(), t, o);
  }

  /// Resample input without pre-computed maps
  void operator ()(int i, int j, int k, int, double *o)
  {
    double x = i, y = j, z = k;
    _Transform  (x, y, z);
    _Interpolate(x, y, z, o);
  }

  /// Resample input using pre-computed world coordinates
  void operator ()(int i, int j, int k, int, const CoordType *wc, double *o)
  {
    double x = i, y = j, z = k;
    _Transform  (x, y, z, wc);
    _Interpolate(x, y, z, o);
  }

  /// Resample input using pre-computed world coordinates and displacements
  void operator ()(int i, int j, int k, int, const CoordType *wc, const double *dx, double *o)
  {
    double x = i, y = j, z = k;
    _Transform  (x, y, z, wc, dx);
    _Interpolate(x, y, z, o);
  }

  /// Resample input using pre-computed world coordinates and additive displacements
  void operator ()(int i, int j, int k, int, const CoordType *wc, const double *d1, const double *d2, double *o)
  {
    double x = i, y = j, z = k;
    _Transform  (x, y, z, wc, d1, d2);
    _Interpolate(x, y, z, o);
  }
};

// =============================================================================
// Update
// =============================================================================

// -----------------------------------------------------------------------------
template <class Transformer, class Interpolator>
void RegisteredImage::Update3(const blocked_range3d<int> &region, bool, bool, bool)
{
  typedef typename Interpolator::Evaluator Evaluator;
  typedef UpdateVoxelFunction<Transformer, Interpolator> UpdateFunction;
  Evaluator *eval = dynamic_cast<Evaluator *>(_Evaluator.get());
  if (eval == nullptr) {
    eval = new Evaluator(this);
    _Evaluator.reset(eval);
  }
  UpdateFunction func(eval, _Transformation, this);
  if (_ImageToWorld) {
    if (_ExternalDisplacement) {
      ParallelForEachVoxel(region, _ImageToWorld, _ExternalDisplacement, this, func);
    } else if (_FixedDisplacement && _Displacement) {
      ParallelForEachVoxel(region, _ImageToWorld, _FixedDisplacement, _Displacement, this, func);
    } else if (_Displacement) {
      ParallelForEachVoxel(region, _ImageToWorld, _Displacement, this, func);
    } else if (_FixedDisplacement) {
      ParallelForEachVoxel(region, _ImageToWorld, _FixedDisplacement, this, func);
    } else {
      ParallelForEachVoxel(region, _ImageToWorld, this, func);
    }
  } else {
    ParallelForEachVoxel(region, this, func);
  }
}

// -----------------------------------------------------------------------------
template <class Transformer, class IntensityFunction, class GradientFunction, class HessianFunction>
void RegisteredImage::Update2(const blocked_range3d<int> &region,
                              bool intensity, bool gradient, bool hessian)
{
  // Auxiliary macro -- undefined again at the end of this body
  #define _update_using(Interpolator)                                          \
    Update3<                                                                   \
        Transformer,                                                           \
        Interpolator<IntensityFunction, GradientFunction, HessianFunction>     \
      >(region, intensity, gradient, hessian)
  if (intensity) {
    if (gradient) {
      if (hessian) {
        _update_using(IntensityAndGradientAndHessianInterpolator);
      } else {
        _update_using(IntensityAndGradientInterpolator);
      }
    } else {
      if (hessian) {
        _update_using(IntensityAndHessianInterpolator);
      } else {
        _update_using(IntensityInterpolator);
      }
    }
  } else {
    if (gradient) {
      if (hessian) {
        _update_using(GradientAndHessianInterpolator);
      } else {
        _update_using(GradientInterpolator);
      }
    } else {
      if (hessian) {
        _update_using(HessianInterpolator);
      } else {
        cerr << "RegisteredImage::Update: At least one output channel should be updated" << endl;
        exit(1);
      }
    }
  }
  #undef _update_using
}

// -----------------------------------------------------------------------------
template <class Transformer>
void RegisteredImage::Update1(const blocked_range3d<int> &region,
                              bool intensity, bool gradient, bool hessian)
{
  enum InterpolationMode interp;
  if (_InterpolationMode == Interpolation_Default) {
    interp = DefaultInterpolationMode();
  } else {
    interp = InterpolationWithoutPadding(_InterpolationMode);
  }

  if (_PrecomputeDerivatives) {
    // Instantiate image functions for commonly used interpolation methods
    // to allow the compiler to generate optimized code for these
    if (interp == Interpolation_Linear || interp == Interpolation_FastLinear) {
      // Auxiliary macro -- undefined again at the end of this body
      #define _update_using(InterpolatorType)                                  \
        Update2<Transformer, InterpolatorType<InputImageType>,                 \
                             InterpolatorType<GradientImageType>,              \
                             InterpolatorType<HessianImageType> >              \
            (region, intensity, gradient, hessian)
      if (this->GetZ() == 1) {
        _update_using(GenericLinearInterpolateImageFunction2D);
      } else {
        _update_using(GenericLinearInterpolateImageFunction3D);
      }
      #undef _update_using
    // Otherwise use generic interpolate image function interface
    } else {
      Update2<Transformer, InterpolateImageFunction,
                           InterpolateImageFunction,
                           InterpolateImageFunction>
          (region, intensity, gradient, hessian);
    }
  } else {
    // Auxiliary macro -- undefined again at the end of this body
    // TODO: Use also some HessianInterpolatorType
    #define _update_using(InterpolatorType, GradientInterpolatorType)          \
      Update2<Transformer, InterpolatorType<InputImageType>,                   \
                           GradientInterpolatorType<InputImageType>,           \
                           InterpolatorType<HessianImageType> >                \
          (region, intensity, gradient, hessian)
    // Instantiate image functions for commonly used interpolation methods
    // to allow the compiler to generate optimized code for these
    if (interp == Interpolation_Linear) {

      if (this->GetZ() == 1) {
        _update_using(GenericLinearInterpolateImageFunction2D,
                      GenericLinearImageGradientFunction2D);
      } else {
        _update_using(GenericLinearInterpolateImageFunction3D,
                      GenericLinearImageGradientFunction3D);
      }
    } else if (interp == Interpolation_FastLinear) {
      if (this->GetZ() == 1) {
        _update_using(GenericLinearInterpolateImageFunction2D,
                      GenericFastLinearImageGradientFunction2D);
      } else {
        _update_using(GenericLinearInterpolateImageFunction3D,
                      GenericFastLinearImageGradientFunction3D);
      }
    // Otherwise use generic interpolate image function interface
    // TODO: Implement and use ImageHessianFunction
    } else {
      Update2<Transformer, InterpolateImageFunction,
                           ImageGradientFunction,
                           InterpolateImageFunction>
          (region, intensity, gradient, hessian);
    }
    #undef _update_using
  }
}

// -----------------------------------------------------------------------------
template <class TOut, class TIn> inline
void CopyChannels(GenericImage<TOut> *tgt, int l, const GenericImage<TIn> *src)
{
  mirtkAssert(tgt->X() == src->X(), "images have identical size");
  mirtkAssert(tgt->Y() == src->Y(), "images have identical size");
  mirtkAssert(tgt->Z() == src->Z(), "images have identical size");
  mirtkAssert(tgt->T() >= l + src->T(), "target has at sufficient number of channels");
  TOut      *out  = tgt->GetPointerToVoxels(0, 0, 0, l);
  const TIn *in   = src->GetPointerToVoxels();
  const int  nvox = src->GetNumberOfVoxels();
  for (int idx = 0; idx < nvox; ++idx) {
    (*out++) = static_cast<TOut>(*in++);
  }
}

// -----------------------------------------------------------------------------
template <> inline
void CopyChannels(GenericImage<RegisteredImage::VoxelType> *tgt, int l,
                  const GenericImage<RegisteredImage::VoxelType> *src)
{
  mirtkAssert(tgt->X() == src->X(), "images have identical size");
  mirtkAssert(tgt->Y() == src->Y(), "images have identical size");
  mirtkAssert(tgt->Z() == src->Z(), "images have identical size");
  mirtkAssert(tgt->T() >= l + src->T(), "target has at sufficient number of channels");
  memcpy(tgt->GetPointerToVoxels(0, 0, 0, l),
         src->GetPointerToVoxels(),
         src->GetNumberOfVoxels() * sizeof(RegisteredImage::VoxelType));
}

// -----------------------------------------------------------------------------
void RegisteredImage::Update(const blocked_range3d<int> &region,
                             bool intensity, bool gradient, bool hessian,
                             bool force)
{
  // Update only channels that were initialized even if requested
  gradient = gradient && this->T() >=  4;
  hessian  = hessian  && this->T() >= 10;

  // Do nothing if no output should be updated
  if (!intensity && !gradient && !hessian) return;

  // Do nothing if no changing transformation is set or self-update is disabled
  // (i.e., external process is responsible for update of registered image)
  if (!force && (!(_Transformation && _NumberOfActiveLevels > 0) || !_SelfUpdate)) return;

  MIRTK_START_TIMING();

  if (_ExternalDisplacement &&
      region.cols ().begin() == 0 && region.cols ().end() == _ExternalDisplacement->X() &&
      region.rows ().begin() == 0 && region.rows ().end() == _ExternalDisplacement->Y() &&
      region.pages().begin() == 0 && region.pages().end() == _ExternalDisplacement->Z()) {

    // Always use provided externally updated displacement field if given
    Update1<DefaultTransformer>(region, intensity, gradient, hessian);

  } else {

    // End time point of deformation and initial time for velocity-based
    // transformations, i.e., time point of initial condition of ODE
    const double t  = _InputImage->GetTOrigin();
    const double t0 = this       ->GetTOrigin();

    // -------------------------------------------------------------------------
    // Update moving image (i.e., constantly changing transformation is set)
    if (_Transformation && _NumberOfActiveLevels > 0) {

      // For some transformations, it is faster to compute the displacements
      // all at once such as those which are represented by velocity fields.
      const bool cache = _CacheDisplacement || _Transformation->RequiresCachingOfDisplacements();
      if (cache && !_Displacement) _Displacement = new DisplacementImageType();

      // If we pre-computed the fixed displacement of the passive MFFD levels
      const MultiLevelTransformation *mffd;
      if (_FixedDisplacement && (mffd = dynamic_cast<const MultiLevelTransformation *>(_Transformation))) {

        if (dynamic_cast<const FluidFreeFormTransformation *>(mffd)) {

          if (_Displacement) {
            *_Displacement = *_FixedDisplacement;
            mffd->Displacement(_NumberOfPassiveLevels, -1,
                               *_Displacement, t, t0, _ImageToWorld);
          }
          Update1<FluidTransformer>(region, intensity, gradient, hessian);

        } else {

          if (_Displacement) {
            _Displacement->Initialize(_attr, 3);
            mffd->Displacement(_NumberOfPassiveLevels, -1,
                               *_Displacement, t, t0, _ImageToWorld);
          }
          Update1<AdditiveTransformer>(region, intensity, gradient, hessian);

        }

      // Otherwise, simply let the (non-)MFFD compute the total transformation
      } else {

        if (_Displacement) {
          _Displacement->Initialize(_attr, 3);
          _Transformation->Displacement(*_Displacement, t, t0, _ImageToWorld);
        }
        Update1<DefaultTransformer>(region, intensity, gradient, hessian);

      }

    // -------------------------------------------------------------------------
    // Update fixed image (i.e. no transformation or no active levels)
    } else if (_Transformation) {

      // Pre-compute fixed transformation displacements if needed
      if (!_FixedDisplacement) {
        if (_CacheFixedDisplacement || _Transformation->RequiresCachingOfDisplacements()) {
          _FixedDisplacement = new DisplacementImageType(_attr, 3);
          _Transformation->Displacement(*_FixedDisplacement, t, t0, _ImageToWorld);
        }
      }

      // Apply fixed transformation
      Update1<FixedTransformer>(region, intensity, gradient, hessian);

      // Discard chached displacements when caching is disabled
      if (!_CacheFixedDisplacement) Delete(_FixedDisplacement);

    // Copy input images when no (fixed) transformation is set and the
    // attributes of input image grid matches those of the output image grid
    } else if (this->HasSpatialAttributesOf(_InputImage)) {

      // Copy intensities
      if (intensity) {
        // Rescale foreground intensities to [_MinIntensity, _MaxIntensity]
        if (!IsNaN(_MinIntensity) || !IsNaN(_MaxIntensity)) {
          const int nvox = NumberOfVoxels();
          if (nvox > 0) {
            InputImageType::VoxelType *iptr = _InputImage->Data();
            InputImageType::VoxelType  imin;
            InputImageType::VoxelType  imax;
            imin = voxel_limits<InputImageType::VoxelType>::max();
            imax = voxel_limits<InputImageType::VoxelType>::min();
            for (int idx = 0; idx < nvox; ++idx, ++iptr) {
              if (_InputImage->IsForeground(idx)) {
                if (*iptr < imin) imin = *iptr;
                if (*iptr > imax) imax = *iptr;
              }
            }
            if (imin <= imax) {
              double omin = _MinIntensity;
              double omax = _MaxIntensity;
              if (IsNaN(omin)) omin = imin;
              if (IsNaN(omax)) omax = imax;
              const double slope = (omax - omin) / static_cast<double>(imax - imin);
              const double inter = omin - slope * static_cast<double>(imin);
              iptr = _InputImage->Data();
              VoxelType *optr = this->Data();
              const VoxelType bg = voxel_cast<VoxelType>(this->_bg);
              for (int idx = 0; idx < nvox; ++idx, ++iptr, ++optr) {
                if (_InputImage->IsForeground(idx)) {
                  *optr = voxel_cast<VoxelType>(inter + slope * static_cast<double>(*iptr));
                  if      (*optr < _MinIntensity) *optr = _MinIntensity;
                  else if (*optr > _MaxIntensity) *optr = _MaxIntensity;
                } else {
                  *optr = bg;
                }
              }
            }
          }
        } else {
          CopyChannels(this, 0, _InputImage);
        }
      }

      // Copy derivatives
      if (gradient) CopyChannels(this, 1, _InputGradient);
      if (hessian)  CopyChannels(this, 4, _InputHessian);

      // Copy background mask (if set)
      this->PutMask(_InputImage->GetMask());

    // Otherwise, resample input images on output image grid
    } else {
      Update1<FixedTransformer>(region, intensity, gradient, hessian);
    }
  }

  MIRTK_DEBUG_TIMING(4, "update of " << (_Transformation ? "moving" : "fixed") << " image"
                     << " (intensity=" << (intensity ? "yes" : "no")
                     << ", gradient="  << (gradient  ? "yes" : "no")
                     << ", hessian="   << (hessian   ? "yes" : "no") << ")");
}

// -----------------------------------------------------------------------------
void RegisteredImage::Update(bool intensity, bool gradient, bool hessian, bool force)
{
  blocked_range3d<int> region(0, Z(), 0, Y(), 0, X());
  this->Update(region, intensity, gradient, hessian, force);
}

// -----------------------------------------------------------------------------
void RegisteredImage::Update(const blocked_range3d<int>  &region,
                             const DisplacementImageType *disp,
                             bool intensity, bool gradient, bool hessian)
{
  // Update only channels that were initialized even if requested
  gradient = gradient && this->T() >=  4;
  hessian  = hessian  && this->T() >= 10;

  // Do nothing if no output should be updated
  if (!intensity && !gradient && !hessian) return;

  // Image to world map required by Update3
  if (!_ImageToWorld) {
    _ImageToWorld = new WorldCoordsImage();
    this->ImageToWorld(*_ImageToWorld, true /* i.e., always 3D vectors */);
  }

  // Keep pointers to own displacement fields
  DisplacementImageType * const _disp  = _Displacement;
  DisplacementImageType * const _fixed = _FixedDisplacement;

  // Replace displacement fields by user arguments
  _Displacement      = const_cast<DisplacementImageType *>(disp);
  _FixedDisplacement = NULL;

  // Interpolate within specified region using fixed transfomer
  Update1<FixedTransformer>(region, intensity, gradient, hessian);

  // Reset pointers to own displacement fields
  _Displacement      = _disp;
  _FixedDisplacement = _fixed;
}

// -----------------------------------------------------------------------------
void RegisteredImage::Recompute(const blocked_range3d<int> &region)
{
  const bool update_intensity = true;
  const bool update_gradient  = true;
  const bool update_hessian   = true;
  const bool force_update     = true;
  this->Update(region, update_intensity, update_gradient, update_hessian, force_update);
}

// -----------------------------------------------------------------------------
void RegisteredImage::Recompute()
{
  this->Recompute(blocked_range3d<int>(0, Z(), 0, Y(), 0, X()));
}


} // namespace mirtk
