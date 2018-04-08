/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2017 Imperial College London
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

#include "mirtk/GenericRegistrationFilter.h"

#include "mirtk/Config.h" // WINDOWS

#include "mirtk/Array.h"
#include "mirtk/Utils.h"
#include "mirtk/Math.h"
#include "mirtk/Memory.h"
#include "mirtk/Version.h"
#include "mirtk/Matrix.h"
#include "mirtk/Parallel.h"
#include "mirtk/Profiling.h"
#include "mirtk/Vector3D.h"
#include "mirtk/VoxelFunction.h"

#include "mirtk/Downsampling.h"
#include "mirtk/Resampling.h"
#include "mirtk/ResamplingWithPadding.h"
#include "mirtk/GaussianBlurring.h"
#include "mirtk/GaussianBlurringWithPadding.h"
#include "mirtk/GaussianPyramidFilter.h"
#include "mirtk/LinearInterpolateImageFunction.h"
#include "mirtk/NearestNeighborInterpolateImageFunction.h"

#include "mirtk/InverseAffineTransformation.h"
#include "mirtk/PartialAffineTransformation.h"
#include "mirtk/PartialBSplineFreeFormTransformationSV.h"
#include "mirtk/FreeFormTransformation.h"
#include "mirtk/FluidFreeFormTransformation.h"
#include "mirtk/MultiLevelFreeFormTransformation.h"
#include "mirtk/MultiLevelStationaryVelocityTransformation.h"
#include "mirtk/PartialMultiLevelStationaryVelocityTransformation.h"

#include "mirtk/ImageSimilarity.h"
#include "mirtk/TransformationConstraint.h"
#include "mirtk/MeanSquaredDisplacementError.h"

#if MIRTK_Registration_WITH_PointSet
#  include "mirtk/RegisteredPointSet.h"
#  include "mirtk/PointSetUtils.h"
#  include "mirtk/PointSetDistance.h"
#  include "mirtk/SurfaceRemeshing.h"
#  if MIRTK_Registration_WITH_Deformable
#    include "mirtk/InternalForce.h"
#  endif
#endif

#include "RegistrationEnergyParser.h"


namespace mirtk {


// =============================================================================
// Auxiliary functions/functors, constants, and types
// =============================================================================

namespace GenericRegistrationFilterUtils {

// -----------------------------------------------------------------------------
// Constants
// -----------------------------------------------------------------------------

// Names of configuration parameters
static const string MINSTEP     = "Minimum length of steps";
static const string MAXSTEP     = "Maximum length of steps";
static const string MAXLINEITER = "Maximum no. of line search iterations";
static const string STRICTRANGE = "Strict step length range";
static const string REUSESTEP   = "Reuse previous step length";
static const string MAXREJECTED = "Maximum streak of rejected steps";

// Tolerance used for voxel size comparisons
static const double TOL = 1.0e-6;

// -----------------------------------------------------------------------------
// Types
// -----------------------------------------------------------------------------

// Global type redefinitions for auxiliary non-class member functions
typedef GenericRegistrationFilter::ResampledImageType ResampledImageType;
typedef GenericRegistrationFilter::ResampledImageList ResampledImageList;
typedef GenericRegistrationFilter::VoxelType          VoxelType;

// -----------------------------------------------------------------------------
// Auxiliaries used by Set parameter function
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
int ParseValues(const char *str, double *x, double *y, double *z)
{
#ifdef WINDOWS
  return sscanf_s(str, "%lf %lf %lf", x, y, z);
#else
  return sscanf(str, "%lf %lf %lf", x, y, z);
#endif
}

// -----------------------------------------------------------------------------
// Auxiliaries used by GuessParameter and InitializePyramid
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
/// Minimum image margin required by Gaussian blurring filter
Vector3D<int> BackgroundMargin(const ImageAttributes &attr, const Vector3D<double> &res, const Vector3D<double> &sigma)
{
  Vector3D<int> margin(0);
  if (res._x > 0. && sigma._x > 0.) {
    int wx = GaussianBlurring<double>::KernelSize(sigma._x / res._x);
    margin._x = iceil(static_cast<double>(wx / 2) * res._x / attr._dx);
  }
  if (res._y > 0. && sigma._y > 0.) {
    auto wy = GaussianBlurring<double>::KernelSize(sigma._y / res._y);
    margin._y = iceil(static_cast<double>(wy / 2) * res._y / attr._dy);
  }
  if (res._z > 0. && sigma._z > 0.) {
    auto wz = GaussianBlurring<double>::KernelSize(sigma._z / res._z);
    margin._z = iceil(static_cast<double>(wz / 2) * res._z / attr._dz);
  }
  return margin;
}

// -----------------------------------------------------------------------------
/// Get foreground image domain with background margin required by Gaussian blurring
ImageAttributes ForegroundDomain(const BaseImage *image, double padding,
                                 const Vector3D<double> &res, const Vector3D<double> &sigma,
                                 bool orthogonal = true)
{
  ImageAttributes attr = image->ForegroundDomain(padding, orthogonal);
  if (sigma > 0.) {
    const auto margin = BackgroundMargin(attr, res, sigma);
    attr._x += 2 * margin._x;
    attr._y += 2 * margin._y;
    attr._z += 2 * margin._z;
  }
  return attr;
}

// -----------------------------------------------------------------------------
/// Get foreground image domain with background margin required by Gaussian blurring
ImageAttributes ForegroundDomain(const BaseImage *image, double padding,
                                 const Vector3D<double> &res, double sigma,
                                 bool orthogonal = true)
{
  return ForegroundDomain(image, padding, res, Vector3D<double>(sigma), orthogonal);
}

// -----------------------------------------------------------------------------
/// Choose default background value if undefined by user
class SetDefaultBackgroundValue
{
  const Array<const BaseImage *> &_Input;
  Array<double>                  &_Background;
  double                          _Default;

public:

  SetDefaultBackgroundValue(const Array<const BaseImage *> &input,
                            Array<double>                  &background,
                            double                          default_bg = NaN)
  :
    _Input(input), _Background(background), _Default(default_bg)
  {}

  void operator()(const blocked_range<int> &re) const
  {
    for (int n = re.begin(); n != re.end(); ++n) {
      if (IsNaN(_Background[n])) {
        if (IsNaN(_Default)) {
          double min_intensity = _Input[n]->GetAsDouble(0);
          const int nvox = _Input[n]->NumberOfVoxels();
          for (int idx = 1; idx < nvox; ++idx) {
            double value = _Input[n]->GetAsDouble(idx);
            if (value < min_intensity) min_intensity = value;
          }
          _Background[n] = min(min_intensity - 1., 0.);
        } else {
          _Background[n] = _Default;
        }
      }
    }
  }
};

// -----------------------------------------------------------------------------
/// Clamp intensities below background without change of image size
class PadImages
{
  const Array<const BaseImage *> &_Input;
  Array<double>                  &_Background;
  ResampledImageList             &_Image;

public:

  PadImages(const Array<const BaseImage *> &input,
            Array<double>                  &background,
            ResampledImageList             &image)
  :
    _Input(input), _Background(background), _Image(image)
  {}

  void operator()(const blocked_range<int> &re) const
  {
    for (int n = re.begin(); n != re.end(); ++n) {
      _Image[n] = *_Input[n];
      _Image[n].PutBackgroundValueAsDouble(_Background[n], true);
    }
  }
};

// -----------------------------------------------------------------------------
/// Clamp intensities below background and crop image
class CropImages
{
  Array<const BaseImage *>  _Input;
  Array<double>            &_Background;
  Array<double>            &_Outside;
  Array<Vector3D<double> > &_Resolution;
  Array<double>            &_Blurring;
  ResampledImageList       &_Output;

public:

  CropImages(const Array<const BaseImage *> &input,
             Array<double>                  &background,
             Array<double>                  &outside,
             Array<Vector3D<double> >       &resolution,
             Array<double>                  &blurring,
             ResampledImageList             &output)
  :
    _Input(input),
    _Background(background),
    _Outside(outside),
    _Resolution(resolution),
    _Blurring(blurring),
    _Output(output)
  {}

  CropImages(const ResampledImageList &input,
             Array<double>            &background,
             Array<double>            &outside,
             Array<Vector3D<double> > &resolution,
             Array<double>            &blurring,
             ResampledImageList       &output)
  :
    _Input(input.size()),
    _Background(background),
    _Outside(outside),
    _Resolution(resolution),
    _Blurring(blurring),
    _Output(output)
  {
    for (size_t n = 0; n < input.size(); ++n) {
      _Input[n] = &input[n];
    }
  }

  void operator()(const blocked_range<int> &re) const
  {
    for (int n = re.begin(); n != re.end(); ++n) {
      // Determine minimal lattice containing foreground voxels
      auto attr = ForegroundDomain(_Input[n], _Background[n], _Resolution[n], _Blurring[n], false);
      // Add extra margin for derivatives computation
      const int margin = 1;
      if (attr._x > 1) attr._x += 2 * margin;
      if (attr._y > 1) attr._y += 2 * margin;
      if (attr._z > 1) attr._z += 2 * margin;
      // Resample input image on foreground lattice
      // (interpolation not required as voxel centers should still match)
      _Output[n].Initialize(attr, 1);
      auto *value = _Output[n].Data();
      for (int k2 = 0; k2 < attr._z; ++k2)
      for (int j2 = 0; j2 < attr._y; ++j2)
      for (int i2 = 0; i2 < attr._x; ++i2) {
        double x = i2, y = j2, z = k2;
        _Output[n]. ImageToWorld(x, y, z);
        _Input [n]->WorldToImage(x, y, z);
        const int i1 = iround(x), j1 = iround(y), k1 = iround(z);
        if (_Input[n]->IsInside(i1, j1, k1)) {
          *value = _Input[n]->GetAsDouble(i1, j1, k1);
          if (*value < _Background[n]) *value = _Background[n];
        } else {
          *value = _Outside[n];
        }
        ++value;
      }
      _Output[n].PutBackgroundValueAsDouble(_Background[n]);
    }
  }
};

// -----------------------------------------------------------------------------
/// Rescale foreground intensities to [min, max]
class Rescale
{
  ResampledImageList &_Image;
  VoxelType           _Min;
  VoxelType           _Max;

public:

  Rescale(ResampledImageList &image, VoxelType min_value, VoxelType max_value)
  :
    _Image(image), _Min(min_value), _Max(max_value)
  {}

  void operator()(const blocked_range<int> &re) const
  {
    for (int n = re.begin(); n != re.end(); ++n) {
      _Image[n].PutMinMax(_Min, _Max);
    }
  }
};

// -----------------------------------------------------------------------------
/// Crop mask to minimum bounding box enclosing foreground
BinaryImage *CropMask(const BinaryImage *input)
{
  bool all_non_zero = true;
  for (int vox = 0; vox < input->NumberOfVoxels(); ++vox) {
    if (input->Get(vox) != 0) {
      all_non_zero = false;
      break;
    }
  }
  // Determine minimal lattice containing foreground voxels
  ImageAttributes attr = input->ForegroundDomain(0., false);
  if (attr._x > 1 && attr._dx > 0.) attr._x += 4;
  if (attr._y > 1 && attr._dy > 0.) attr._y += 4;
  if (attr._z > 1 && attr._dz > 0.) attr._z += 4;
  // Resample input image on foreground lattice
  // (interpolation not required as voxel centers should still match)
  BinaryImage *output = new BinaryImage(attr, 1);
  if (all_non_zero) {
    *output = 1;
  } else {
    BinaryPixel *value  = output->Data();
    for (int k2 = 0; k2 < attr._z; ++k2)
    for (int j2 = 0; j2 < attr._y; ++j2)
    for (int i2 = 0; i2 < attr._x; ++i2) {
      double x = i2, y = j2, z = k2;
      output->ImageToWorld(x, y, z);
      input ->WorldToImage(x, y, z);
      int i1 = iround(x), j1 = iround(y), k1 = iround(z);
      if (input->IsInside(i1, j1, k1)) {
        *value = input->Get(i1, j1, k1);
        if (*value != 0) *value = 1;
      }
      ++value;
    }
  }
  return output;
}

// -----------------------------------------------------------------------------
/// Copy images from first level without cropping
class CopyImages
{
  const ResampledImageList &_Input;
  ResampledImageList       &_Output;

public:

  CopyImages(const ResampledImageList &input, ResampledImageList &output)
  :
    _Input(input), _Output(output)
  {}

  void operator()(const blocked_range<int> &re) const
  {
    for (int n = re.begin(); n != re.end(); ++n) _Output[n] = _Input[n];
  }
};

// -----------------------------------------------------------------------------
/// Blur images either after downsampling using Gaussian resolution pyramid
/// or before applying the resampling filter
class BlurImages
{
  Array<ResampledImageList> &_Image;
  const Array<double>       *_Sigma;
  const Array<double>       *_Padding;

public:

  BlurImages(Array<ResampledImageList> &image,
             const Array<double>       *sigma,
             const Array<double>       *padding = nullptr)
  :
    _Image(image), _Sigma(sigma), _Padding(padding)
  {}

  void operator()(const blocked_range2d<int> &re) const
  {
    for (int l = re.rows().begin(); l != re.rows().end(); ++l)
    for (int n = re.cols().begin(); n != re.cols().end(); ++n) {
      if (_Sigma[l][n] > .0) {
        if (_Padding) {
          GaussianBlurringWithPadding<VoxelType> blurring(_Sigma[l][n], (*_Padding)[n]);
          blurring.Input (&_Image[l][n]);
          blurring.Output(&_Image[l][n]);
          blurring.Run();
        } else {
          GaussianBlurring<VoxelType> blurring(_Sigma[l][n]);
          blurring.Input (&_Image[l][n]);
          blurring.Output(&_Image[l][n]);
          blurring.Run();
        }
      }
    }
  }
};

// -----------------------------------------------------------------------------
/// Resample images (after user defined blurring)
class ResampleImages
{
  Array<ResampledImageList>      &_Image;
  const Array<Vector3D<double> > *_Resolution;
  const Array<double>            *_Background;
  const Array<double>            *_Padding;

public:

  ResampleImages(Array<ResampledImageList>      &image,
                 const Array<Vector3D<double> >  res[],
                 const Array<double>            &background,
                 const Array<double>            *padding = nullptr)
  :
    _Image(image), _Resolution(res), _Background(&background), _Padding(padding)
  {}

  void operator()(const blocked_range2d<int> &re) const
  {
    GenericLinearInterpolateImageFunction<ResampledImageType> f;
    for (int l = re.rows().begin(); l != re.rows().end(); ++l)
    for (int n = re.cols().begin(); n != re.cols().end(); ++n) {
      const Vector3D<double> &res = _Resolution[l][n];
      if (res._x > 0 && res._y > 0 && res._z > 0) {
        double dx, dy, dz;
        _Image[l][n].GetPixelSize(&dx, &dy, &dz);
        if (!fequal(res._x, dx, TOL) ||
            !fequal(res._y, dy, TOL) ||
            !fequal(res._z, dz, TOL)) {
          f.DefaultValue((*_Background)[n]);
          if (_Padding) {
            ResamplingWithPadding<VoxelType> resample(res._x, res._y, res._z, (*_Padding)[n]);
            resample.Interpolator(&f);
            resample.Input (&_Image[l][n]);
            resample.Output(&_Image[l][n]);
            resample.Run();
          } else {
            Resampling<VoxelType> resample(res._x, res._y, res._z);
            resample.Interpolator(&f);
            resample.Input (&_Image[l][n]);
            resample.Output(&_Image[l][n]);
            resample.Run();
          }
        }
      }
    }
  }
};

// -----------------------------------------------------------------------------
/// Instead of blurring and resampling to a custom user defined resolution,
/// downsample images using a standard Gaussian image pyramid approach
class DownsampleImages
{
  Array<ResampledImageList> &_Image;
  const Array<double>       *_Background;
  const Array<double>       *_Outside;
  const Array<double>       *_Padding;
  const Array<double>       *_Blurring;
  bool                       _CropPad;
  int                        _Level;

public:

  DownsampleImages(Array<ResampledImageList> &image, int l,
                   const Array<double> *background,
                   const Array<double> *outside  = nullptr,
                   const Array<double> *padding  = nullptr,
                   const Array<double> *blurring = nullptr,
                   bool crop_pad = false)
  :
    _Image(image),
    _Background(background),
    _Outside(outside ? outside : background),
    _Padding(padding),
    _Blurring(blurring),
    _CropPad(crop_pad),
    _Level(l)
  {}

  void operator()(const blocked_range<int> &re) const
  {
    GenericLinearInterpolateImageFunction<ResampledImageType> f;
    for (int n = re.begin(); n != re.end(); ++n) {
      // Outside padding value
      f.DefaultValue((*_Outside)[n]);
      // Determine spacing and blurring sigma for this level
      double dx = (_Image[1][n].X() > 1 ? _Image[1][n].XSize() : .0);
      double dy = (_Image[1][n].Y() > 1 ? _Image[1][n].YSize() : .0);
      double dz = (_Image[1][n].Z() > 1 ? _Image[1][n].ZSize() : .0);
      Vector3D<double> var(.0);
      for (int l = 2; l <= _Level; ++l) {
        var._x +=  pow(0.7355 * dx, 2);
        var._y +=  pow(0.7355 * dy, 2);
        var._z +=  pow(0.7355 * dz, 2);
        dx *= 2., dy *= 2., dz *= 2.;
      }
      Vector3D<double> sigma(sqrt(var._x), sqrt(var._y), sqrt(var._z));
      // Determine minimal lattice containing foreground voxels
      ImageAttributes attr;
      if (_CropPad) {
        // Crop image with extra margin for Gaussian filtering
        const auto res = Vector3D<double>(dx, dy, dz);
        auto blurring = sigma;
        if (_Blurring) {
          blurring._x = max(blurring._x, (*_Blurring)[n]);
          blurring._y = max(blurring._y, (*_Blurring)[n]);
          blurring._z = max(blurring._z, (*_Blurring)[n]);
        }
        attr = ForegroundDomain(&_Image[1][n], (*_Background)[n], res, blurring, false);
      } else {
        attr = _Image[1][n].Attributes();
      }
      // Ensure lattice can be divided by 2, given that dx = 2^n * attr._dx
      if (attr._x % 2 != 0) attr._x += 1, attr._xorigin += .5 * attr._dx;
      if (attr._y % 2 != 0) attr._y += 1, attr._yorigin += .5 * attr._dy;
      if (attr._z % 2 != 0) attr._z += 1, attr._zorigin += .5 * attr._dz;
      // Calculate number of voxels preserving the image size
      for (int l = 2; l <= _Level; ++l) {
        if (attr._x < 64) {
          dx *= .5;
        } else {
          attr._x /= 2;
        }
        if (attr._y < 64) {
          dy *= .5;
        } else {
          attr._y /= 2;
        }
        if (attr._z < 64) {
          dz *= .5;
        } else {
          attr._z /= 2;
        }
      }
      // Add extra margin for derivatives computation
      if (_CropPad) {
        const int margin = 1;
        if (attr._x > 1) attr._x += 2 * margin;
        if (attr._y > 1) attr._y += 2 * margin;
        if (attr._z > 1) attr._z += 2 * margin;
      }
      // If background value set, consider foreground only
      if (_Padding) {
        typedef GaussianBlurringWithPadding<VoxelType> BlurFilter;
        typedef ResamplingWithPadding      <VoxelType> ResampleFilter;
        BlurFilter blurring(sigma._x, sigma._y, sigma._z, (*_Padding)[n]);
        blurring.Input (&_Image[1][n]);
        blurring.Output(&_Image[_Level][n]);
        blurring.Run();
        ResampleFilter resample(attr._x, attr._y, attr._z, dx, dy, dz, (*_Padding)[n]);
        resample.Interpolator(&f);
        resample.Input (&_Image[_Level][n]);
        resample.Output(&_Image[_Level][n]);
        resample.Run();
      // Otherwise, downsample using all image intensities
      } else if (_CropPad) {
        typedef GaussianBlurring<VoxelType> BlurFilter;
        typedef Resampling      <VoxelType> ResampleFilter;
        // Ensure that blurring creates smooth edges all around even when
        // the image (foreground) is very close to the image boundary
        int wx = BlurFilter::KernelSize(sigma._x / _Image[1][n].XSize());
        int wy = BlurFilter::KernelSize(sigma._y / _Image[1][n].YSize());
        int wz = BlurFilter::KernelSize(sigma._z / _Image[1][n].ZSize());
        ResampleFilter resize(_Image[1][n].X() + wx,
                              _Image[1][n].Y() + wy,
                              _Image[1][n].Z() + wz);
        resize.Interpolator(&f);
        resize.Input (&_Image[1][n]);
        resize.Output(&_Image[_Level][n]);
        resize.Run();
        BlurFilter blurring(sigma._x, sigma._y, sigma._z);
        blurring.Input (&_Image[_Level][n]);
        blurring.Output(&_Image[_Level][n]);
        blurring.Run();
        ResampleFilter resample(attr._x, attr._y, attr._z, dx, dy, dz);
        resample.Interpolator(&f);
        resample.Input (&_Image[_Level][n]);
        resample.Output(&_Image[_Level][n]);
        resample.Run();
      } else {
        typedef GaussianBlurring<VoxelType> BlurFilter;
        typedef Resampling      <VoxelType> ResampleFilter;
        BlurFilter blurring(sigma._x, sigma._y, sigma._z);
        blurring.Input (&_Image[1][n]);
        blurring.Output(&_Image[_Level][n]);
        blurring.Run();
        ResampleFilter resample(attr._x, attr._y, attr._z, dx, dy, dz);
        resample.Interpolator(&f);
        resample.Input (&_Image[_Level][n]);
        resample.Output(&_Image[_Level][n]);
        resample.Run();
      }
    }
  }
};

// -----------------------------------------------------------------------------
/// Resample provided foreground mask
class ResampleMask
{
  BinaryImage                    *_Domain;
  Array<BinaryImage *>           &_Mask;
  const Array<Vector3D<double> > &_Resolution;
  Vector3D<double>                _VoxelSize;

public:

  ResampleMask(BinaryImage                    *domain,
               Array<BinaryImage *>           &mask,
               const Array<Vector3D<double> > &res)
  :
    _Domain(domain), _Mask(mask), _Resolution(res)
  {
    _Domain->GetPixelSize(&_VoxelSize._x, &_VoxelSize._y, &_VoxelSize._z);
  }

  void operator()(const blocked_range<int> &re) const
  {
    GenericLinearInterpolateImageFunction<BinaryImage> interpolator;
    interpolator.Input(_Domain);
    interpolator.Initialize();
    interpolator.DefaultValue(0.);
    for (int l = re.begin(); l != re.end(); ++l) {
      if (_Mask[l] != _Domain) Delete(_Mask[l]);
      _Mask[l] = _Domain;
      const double dx = _Resolution[l]._x;
      const double dy = _Resolution[l]._y;
      const double dz = _Resolution[l]._z;
      if ((_Domain->X() > 1 && dx > .0) ||
          (_Domain->Y() > 1 && dy > .0) ||
          (_Domain->Z() > 1 && dz > .0)) {
        if (!fequal(dx, _VoxelSize._x, TOL) ||
            !fequal(dy, _VoxelSize._y, TOL) ||
            !fequal(dz, _VoxelSize._z, TOL)) {
          // Allocate output mask
          _Mask[l] = new BinaryImage();
          // Calculate number of voxels preserving the image size
          int nx = (dx > 0. ? iceil(_Domain->X() * _Domain->XSize() / dx) : 1);
          int ny = (dy > 0. ? iceil(_Domain->Y() * _Domain->YSize() / dy) : 1);
          int nz = (dz > 0. ? iceil(_Domain->Z() * _Domain->ZSize() / dz) : 1);
          // Resample mask using linear interpolation, voxels with interpolated
          // value < 0.5 are cast to zero and other are cast to 1
          Resampling<BinaryPixel> resample(nx, ny, nz, dx, dy, dz);
          resample.Input (_Domain);
          resample.Output(_Mask[l]);
          resample.Interpolator(&interpolator);
          resample.Run();
        }
      }
    }
  }
};

// -----------------------------------------------------------------------------
// Functor types used by InitializePointSets
// -----------------------------------------------------------------------------

#if MIRTK_Registration_WITH_PointSet
// -----------------------------------------------------------------------------
/// Remesh input surface meshes
class RemeshSurfaces
{
  int                                           _Level;
  const Array<vtkSmartPointer<vtkPointSet> >   *_Input;
  Array<Array<vtkSmartPointer<vtkPointSet> > > *_Output;
  const Array<double>                          *_MinEdgeLength;
  const Array<double>                          *_MaxEdgeLength;

public:

  RemeshSurfaces(int level, const Array<vtkSmartPointer<vtkPointSet> > &input,
                 Array<Array<vtkSmartPointer<vtkPointSet> > > &output,
                 const Array<double> *dmin, const Array<double> *dmax)
  :
    _Level(level),
    _Input(&input),
    _Output(&output),
    _MinEdgeLength(dmin),
    _MaxEdgeLength(dmax)
  {}

  void operator ()(const blocked_range<int> &re) const
  {
    for (int n = re.begin(); n != re.end(); ++n) {
      if (IsSurfaceMesh((*_Input)[n]) &&
          (_MinEdgeLength[_Level][n] > .0 ||
           _MaxEdgeLength[_Level][n] < inf)) {
        SurfaceRemeshing remesher;
        remesher.Input(vtkPolyData::SafeDownCast((*_Output)[_Level-1][n]));
        remesher.MinEdgeLength(_MinEdgeLength[_Level][n]);
        remesher.MaxEdgeLength(_MaxEdgeLength[_Level][n]);
        remesher.MeltingOrder(SurfaceRemeshing::AREA);
        remesher.MeltNodesOn();
        remesher.MeltTrianglesOn();
        remesher.Run();
        (*_Output)[_Level][n] = remesher.Output();
      }
    }
  }
};
#endif // MIRTK_Registration_WITH_PointSet

// -----------------------------------------------------------------------------
// Functor types used by InitializeStatus
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
class InitializeCPStatus
{
  const ResampledImageList &_Image;
  const Array<bool>        &_IsTargetImage;
  const int                 _NumberOfImages;
  FreeFormTransformation   *_FFD;

  #if MIRTK_Registration_WITH_PointSet
    const Array<vtkSmartPointer<vtkPointSet> > &_PointSetInput;
    const Array<bool>                          &_IsMovingPointSet;
    const int                                   _NumberOfPointSets;
  #endif // MIRTK_Registration_WITH_PointSet

  bool _RegisterX, _RegisterY, _RegisterZ;

public:

  InitializeCPStatus(const ResampledImageList &image,
                     const Array<bool>        &is_target_image,
                     FreeFormTransformation   *ffd,
                     #if MIRTK_Registration_WITH_PointSet
                       const Array<vtkSmartPointer<vtkPointSet> > &pointsets,
                       const Array<bool> &is_moving_pointset,
                     #endif // MIRTK_Registration_WITH_PointSet
                     bool regx, bool regy, bool regz)
  :
    _Image           (image),
    _IsTargetImage   (is_target_image),
    _NumberOfImages  (static_cast<int>(image.size())),
    _FFD             (ffd),
    #if MIRTK_Registration_WITH_PointSet
      _PointSetInput    (pointsets),
      _IsMovingPointSet (is_moving_pointset),
      _NumberOfPointSets(static_cast<int>(pointsets.size())),
    #endif // MIRTK_Registration_WITH_PointSet
    _RegisterX(regx), _RegisterY(regy), _RegisterZ(regz)
  {
    #if MIRTK_Registration_WITH_PointSet
      // vtkPointSet::GetBounds is only thread-safe if bounds are precomputed
      for (size_t i = 0; i < pointsets.size(); ++i) pointsets[i]->ComputeBounds();
    #endif // MIRTK_Registration_WITH_PointSet
  }

  void operator()(const blocked_range<int> &re) const
  {
    Transformation::DOFStatus sx, sy, sz;
    double x1, y1, z1, x2, y2, z2;
    int    ci, cj, ck, cl, i1, j1, k1, l1, i2, j2, k2, l2;
    bool   fg;

    for (int cp = re.begin(); cp != re.end(); ++cp) {
      fg = false;
      // Bounding box of control point
      _FFD->BoundingBox(cp, x1, y1, z1, x2, y2, z2);
      // Check if any non-padded input image voxel is influenced by control point
      for (int n = 0; n < _NumberOfImages; ++n) {
        if (_IsTargetImage[n]) {
          if (_FFD->BoundingBox(&_Image[n], cp, i1, j1, k1, l1, i2, j2, k2, l2)) {
            for (int l = l1; l <= l2; ++l)
            for (int k = k1; k <= k2; ++k)
            for (int j = j1; j <= j2; ++j)
            for (int i = i1; i <= i2; ++i) {
              if (_Image[n].IsForeground(i, j, k, l)) {
                fg = true;
                j = j2, k = k2, l = l2, n = _NumberOfImages; // Break out of all loops
                break;
              }
            }
          }
        }
      }
      // Otherwise, check if any point is influenced by control point
      #if MIRTK_Registration_WITH_PointSet
        for (int n = 0; !fg && n < _NumberOfPointSets; ++n) {
          if (_IsMovingPointSet[n]) {
            double b[6];
            _PointSetInput[n]->GetBounds(b);
            fg = (x1 <= b[1] && x2 >= b[0] &&
                  y1 <= b[3] && y2 >= b[2] &&
                  z1 <= b[5] && z2 >= b[4]);
          }
        }
      #endif // MIRTK_Registration_WITH_PointSet
      // Set status of unused DoFs to passive
      _FFD->IndexToLattice(cp, ci, cj, ck, cl);
      if (fg) {
        _FFD->GetStatus(ci, cj, ck, cl, sx, sy, sz);
        if (!_RegisterX) sx = Passive;
        if (!_RegisterY) sy = Passive;
        if (!_RegisterZ) sz = Passive;
      } else {
        sx = sy = sz = Passive;
      }
      _FFD->PutStatus(ci, cj, ck, cl, sx, sy, sz);
    }
  }
};


// -----------------------------------------------------------------------------
class InitializeCPStatusGivenDomainMask
{
  const BinaryImage        *_Mask;
  FreeFormTransformation   *_FFD;
  bool _RegisterX, _RegisterY, _RegisterZ;

public:

  InitializeCPStatusGivenDomainMask(const BinaryImage      *mask,
                                    FreeFormTransformation *ffd,
                                    bool regx, bool regy, bool regz)
  :
    _Mask(mask), _FFD(ffd), _RegisterX(regx), _RegisterY(regy), _RegisterZ(regz)
  {}

  void operator()(const blocked_range<int> &re) const
  {
    Transformation::DOFStatus sx, sy, sz;
    double x1, y1, z1, x2, y2, z2;
    int    ci, cj, ck, cl, i1, j1, k1, l1, i2, j2, k2, l2;
    bool   fg;

    for (int cp = re.begin(); cp != re.end(); ++cp) {
      // Bounding box of control point
      _FFD->BoundingBox(cp, x1, y1, z1, x2, y2, z2);
      // Check if control point is in vicinity of domain mask
      fg = false;
      if (_FFD->BoundingBox(_Mask, cp, i1, j1, k1, l1, i2, j2, k2, l2)) {
        for (int l = l1; l <= l2; ++l)
        for (int k = k1; k <= k2; ++k)
        for (int j = j1; j <= j2; ++j)
        for (int i = i1; i <= i2; ++i) {
          if (_Mask->Get(i, j, k, l)) {
            fg = true;
            j = j2, k = k2, l = l2; // Break out of all loops
            break;
          }
        }
      }
      // Set status of unused DoFs to passive
      _FFD->IndexToLattice(cp, ci, cj, ck, cl);
      if (fg) {
        _FFD->GetStatus(ci, cj, ck, cl, sx, sy, sz);
        if (!_RegisterX) sx = Passive;
        if (!_RegisterY) sy = Passive;
        if (!_RegisterZ) sz = Passive;
      } else {
        sx = sy = sz = Passive;
      }
      _FFD->PutStatus(ci, cj, ck, cl, sx, sy, sz);
    }
  }
};


} // namespace GenericRegistrationFilterUtils
using namespace GenericRegistrationFilterUtils;

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
void GenericRegistrationFilter::Reset()
{
  _CurrentModel = TM_Unknown;
  _CurrentLevel = 0;
  // Free memory allocated upon previous Run
  _Image .clear();
  _Energy.Clear();
  #if MIRTK_Registration_WITH_PointSet
    for (size_t i = 0; i < _PointSetOutput.size(); ++i) {
      delete _PointSetOutput[i];
    }
  #endif // MIRTK_Registration_WITH_PointSet
  _PointSet.clear();
  _PointSetOutput.clear();
  _PointSetOutputInfo.clear();
  for (size_t i = 0; i < _TransformationInstance.size(); ++i) {
    Delete(_TransformationInstance[i]);
  }
  _TransformationInfo    .clear();
  _TransformationInstance.clear();
  for (size_t i = 0; i < _DisplacementField.size(); ++i) {
    Delete(_DisplacementField[i]);
  }
  _DisplacementInfo.clear();
  _DisplacementField.clear();
  Delete(_Transformation);
  Delete(_Optimizer);
  if (!_Mask.empty()) {
    for (size_t l = 1; l < _Mask.size(); ++l) {
      if (_Mask[l] != _Mask[0]) Delete(_Mask[l]);
    }
    if (_Mask[0] != _Domain) Delete(_Mask[0]);
    _Mask.clear();
  }
  // Invalidate all settings such that GuessParameter knows which have been
  // set by the user and which have not been properly initialized
  InterpolationMode(Interpolation_Default);
  ExtrapolationMode(Extrapolation_Default);
  _TransformationModel.clear();
  _TargetTransformationErrorWeight     = 0.;
  _TargetTransformationErrorName       = "MSDE";
  _NumberOfLevels                      = -1;
  _FinalLevel                          = 1;
  _MultiLevelMode                      = MFFD_Default;
  _MergeGlobalAndLocalTransformation   = false;
  _DirichletBoundaryCondition          = false;
  _PrecomputeDerivatives               = true;
  _SimilarityMeasure                   = SIM_NMI;
  _PointSetDistanceMeasure             = PDM_FRE;
  _OptimizationMethod                  = OM_ConjugateGradientDescent;
  _RegisterX = _RegisterY = _RegisterZ = true;
  _DownsampleWithPadding               = true;
  _CropPadImages                       = true;
  _CropPadFFD                          = -1;
  _NormalizeWeights                    = true;
  _AdaptiveRemeshing                   = false;
  _TargetOffset = _SourceOffset = Point();
  _EnergyFormula.clear();
  _ImageSimilarityInfo.clear();
  _PointSetDistanceInfo.clear();
  _PointSetConstraintInfo.clear();
  _Background.clear();
  _DefaultBackground = NaN;
  _MaxRescaledIntensity = inf;
  memset(_MinControlPointSpacing, 0, 4 * MAX_NO_RESOLUTIONS * sizeof(double));
  memset(_MaxControlPointSpacing, 0, 4 * MAX_NO_RESOLUTIONS * sizeof(double));
  for (int level = 0; level < MAX_NO_RESOLUTIONS; ++level) {
    _Centering [level] = -1;
    _Parameter [level].clear();
    _Resolution[level].clear();
    _Blurring  [level].clear();
    _MinEdgeLength[level].clear();
    _MaxEdgeLength[level].clear();
    _Subdivide [level][0] = true;
    _Subdivide [level][1] = true;
    _Subdivide [level][2] = true;
    _Subdivide [level][3] = false;
  }
  _UseGaussianResolutionPyramid = -1;
  _TargetOffset = _SourceOffset = Point(.0, .0, .0);
}

// -----------------------------------------------------------------------------
void GenericRegistrationFilter::Clear()
{
  // Reset settings
  this->Reset();
  // Clear inputs
  _Input.clear();
  _InitialGuess = NULL;
}

// -----------------------------------------------------------------------------
GenericRegistrationFilter::GenericRegistrationFilter()
:
  _InitialGuess(nullptr),
  _TargetTransformation(nullptr),
  _Domain(nullptr),
  _Transformation(nullptr),
  _Optimizer(nullptr)
{
  // Bind broadcast method to optimizer events (excl. Start/EndEvent!)
  _EventDelegate.Bind(IterationEvent,                MakeDelegate(this, &Observable::Broadcast));
  _EventDelegate.Bind(IterationStartEvent,           MakeDelegate(this, &Observable::Broadcast));
  _EventDelegate.Bind(IterationEndEvent,             MakeDelegate(this, &Observable::Broadcast));
  _EventDelegate.Bind(LineSearchStartEvent,          MakeDelegate(this, &Observable::Broadcast));
  _EventDelegate.Bind(LineSearchIterationStartEvent, MakeDelegate(this, &Observable::Broadcast));
  _EventDelegate.Bind(LineSearchIterationEndEvent,   MakeDelegate(this, &Observable::Broadcast));
  _EventDelegate.Bind(LineSearchEndEvent,            MakeDelegate(this, &Observable::Broadcast));
  _EventDelegate.Bind(AcceptedStepEvent,             MakeDelegate(this, &Observable::Broadcast));
  _EventDelegate.Bind(RejectedStepEvent,             MakeDelegate(this, &Observable::Broadcast));
  _EventDelegate.Bind(RestartEvent,                  MakeDelegate(this, &Observable::Broadcast));
  _EventDelegate.Bind(LogEvent,                      MakeDelegate(this, &Observable::Broadcast));
  // Bind pre-update callback function
  _PreUpdateDelegate = MakeDelegate(this, &GenericRegistrationFilter::PreUpdateCallback);
  // Invalidate all settings
  Reset();
}

// -----------------------------------------------------------------------------
GenericRegistrationFilter::~GenericRegistrationFilter()
{
  // Stop forwarding events
  _Energy.DeleteObserver(_EventDelegate);
  // Clear allocated memory and input
  Clear();
}

// =============================================================================
// Input images
// =============================================================================

// -----------------------------------------------------------------------------
void GenericRegistrationFilter::AddInput(const BaseImage *image)
{
  _Input.push_back(image);
}

// -----------------------------------------------------------------------------
void GenericRegistrationFilter::Input(const BaseImage *target, const BaseImage *source)
{
  _Input.clear();
  AddInput(target);
  AddInput(source);
}

// -----------------------------------------------------------------------------
void GenericRegistrationFilter::Input(int num, const BaseImage **image)
{
  _Input.clear();
  for (int n = 0; n < num; ++n) AddInput(image[n]);
}

// -----------------------------------------------------------------------------
int GenericRegistrationFilter::NumberOfImages() const
{
  return static_cast<int>(_Input.size());
}

// -----------------------------------------------------------------------------
int GenericRegistrationFilter::NumberOfRequiredImages() const
{
  int n = 0;
  Array<ImageSimilarityInfo>::const_iterator sim;
  for (sim = _ImageSimilarityInfo.begin(); sim != _ImageSimilarityInfo.end(); ++sim) {
    if (sim->_TargetIndex + 1 > n) n = sim->_TargetIndex + 1;
    if (sim->_SourceIndex + 1 > n) n = sim->_SourceIndex + 1;
  }
  Array<PointSetConstraintInfo>::const_iterator cst;
  for (cst = _PointSetConstraintInfo.begin(); cst != _PointSetConstraintInfo.end(); ++cst) {
    if (cst->_RefImageIndex + 1 > n) n = cst->_RefImageIndex + 1;
  }
  return n;
}

// -----------------------------------------------------------------------------
bool GenericRegistrationFilter::IsTargetImage(int n) const
{
  // Note: An image can be both a source and a target (cf. inverse consistent energy)
  Array<ImageSimilarityInfo>::const_iterator it;
  for (it = _ImageSimilarityInfo.begin(); it != _ImageSimilarityInfo.end(); ++it) {
    if (it->_TargetIndex == n && !it->_TargetTransformation.IsForwardTransformation()) return true;
    if (it->_SourceIndex == n && !it->_SourceTransformation.IsForwardTransformation()) return true;
  }
  return false;
}

// -----------------------------------------------------------------------------
bool GenericRegistrationFilter::IsSourceImage(int n) const
{
  // Note: An image can be both a source and a target (cf. inverse consistent energy)
  Array<ImageSimilarityInfo>::const_iterator it;
  for (it = _ImageSimilarityInfo.begin(); it != _ImageSimilarityInfo.end(); ++it) {
    if (it->_TargetIndex == n && it->_TargetTransformation.IsForwardTransformation()) return true;
    if (it->_SourceIndex == n && it->_SourceTransformation.IsForwardTransformation()) return true;
  }
  return false;
}

// -----------------------------------------------------------------------------
bool GenericRegistrationFilter::IsFixedImage(int n) const
{
  return !IsMovingImage(n);
}

// -----------------------------------------------------------------------------
bool GenericRegistrationFilter::IsMovingImage(int n) const
{
  // Note: Image is moving if it is being transformed by any similarity term
  Array<ImageSimilarityInfo>::const_iterator it;
  for (it = _ImageSimilarityInfo.begin(); it != _ImageSimilarityInfo.end(); ++it) {
    // Note: For longitudinal registration, target image deformed by FFD at t=0 may be compared to itself
    if (it->_TargetIndex == n && it->_SourceIndex != n && !it->_TargetTransformation.IsIdentity()) return true;
    if (it->_TargetIndex != n && it->_SourceIndex == n && !it->_SourceTransformation.IsIdentity()) return true;
  }
  return false;
}

// -----------------------------------------------------------------------------
void GenericRegistrationFilter::InterpolationMode(enum InterpolationMode mode)
{
  _InterpolationMode.clear();
  _DefaultInterpolationMode = mode;
}

// -----------------------------------------------------------------------------
void GenericRegistrationFilter::InterpolationMode(int n, enum InterpolationMode mode)
{
  if (static_cast<size_t>(n) >= _InterpolationMode.size()) {
    _InterpolationMode.resize(n + 1, Interpolation_Default);
  }
  _InterpolationMode[n] = mode;
}

// -----------------------------------------------------------------------------
enum InterpolationMode GenericRegistrationFilter::InterpolationMode(int n) const
{
  if (n < 0 || static_cast<size_t>(n) >= _InterpolationMode.size() || _InterpolationMode[n] == Interpolation_Default) {
    return _DefaultInterpolationMode;
  } else {
    return _InterpolationMode[n];
  }
}

// -----------------------------------------------------------------------------
void GenericRegistrationFilter::ExtrapolationMode(enum ExtrapolationMode mode)
{
  _ExtrapolationMode.clear();
  _DefaultExtrapolationMode = mode;
}

// -----------------------------------------------------------------------------
void GenericRegistrationFilter::ExtrapolationMode(int n, enum ExtrapolationMode mode)
{
  if (static_cast<size_t>(n) >= _ExtrapolationMode.size()) {
    _ExtrapolationMode.resize(n + 1, Extrapolation_Default);
  }
  _ExtrapolationMode[n] = mode;
}

// -----------------------------------------------------------------------------
enum ExtrapolationMode GenericRegistrationFilter::ExtrapolationMode(int n) const
{
  if (n < 0 || static_cast<size_t>(n) >= _ExtrapolationMode.size() || _ExtrapolationMode[n] == Extrapolation_Default) {
    return _DefaultExtrapolationMode;
  } else {
    return _ExtrapolationMode[n];
  }
}

// -----------------------------------------------------------------------------
double GenericRegistrationFilter::BackgroundValue(int n) const
{
  if (n < 0 || static_cast<size_t>(n) >= _Background.size()) {
    return _DefaultBackground;
  } else {
    return _Background[n];
  }
}

// =============================================================================
// Input points, lines, and/or surfaces
// =============================================================================
#if MIRTK_Registration_WITH_PointSet

// -----------------------------------------------------------------------------
void GenericRegistrationFilter::AddInput(vtkPointSet *data, double t)
{
  _PointSetInput.push_back(data);
  _PointSetTime .push_back(t);
}

// -----------------------------------------------------------------------------
void GenericRegistrationFilter::Input(vtkPointSet *target, vtkPointSet *source,
                                      double t, double t0)
{
  _PointSetInput.clear();
  _PointSetTime .clear();
  AddInput(target, t);
  AddInput(source, t0);
}

// -----------------------------------------------------------------------------
void GenericRegistrationFilter::Input(int num, vtkPointSet **points, double *t)
{
  _PointSetInput.clear();
  _PointSetTime .clear();
  for (int n = 0; n < num; ++n) AddInput(points[n], t ? t[n] : .0);
}

#endif // MIRTK_Registration_WITH_PointSet

// -----------------------------------------------------------------------------
int GenericRegistrationFilter::NumberOfPointSets() const
{
  return static_cast<int>(_PointSetInput.size());
}

// -----------------------------------------------------------------------------
int GenericRegistrationFilter::NumberOfRequiredPointSets() const
{
  int n = 0;
  Array<PointSetDistanceInfo>::const_iterator dist;
  for (dist = _PointSetDistanceInfo.begin(); dist != _PointSetDistanceInfo.end(); ++dist) {
    if (dist->_TargetIndex + 1 > n) n = dist->_TargetIndex + 1;
    if (dist->_SourceIndex + 1 > n) n = dist->_SourceIndex + 1;
  }
  Array<PointSetConstraintInfo>::const_iterator cst;
  for (cst = _PointSetConstraintInfo.begin(); cst != _PointSetConstraintInfo.end(); ++cst) {
    if (cst->_PointSetIndex    + 1 > n) n = cst->_PointSetIndex    + 1;
    if (cst->_RefPointSetIndex + 1 > n) n = cst->_RefPointSetIndex + 1;
  }
  return n;
}

// -----------------------------------------------------------------------------
bool GenericRegistrationFilter::IsTargetPointSet(int n) const
{
  // Note: Target points are transformed to source space, not vice versa!
  //       A point set can be both a source and a target (cf. inverse consistent energy).
  Array<PointSetDistanceInfo>::const_iterator it;
  for (it = _PointSetDistanceInfo.begin(); it != _PointSetDistanceInfo.end(); ++it) {
    if (it->_TargetIndex == n && it->_TargetTransformation.IsForwardTransformation()) return true;
    if (it->_SourceIndex == n && it->_SourceTransformation.IsForwardTransformation()) return true;
  }
  return false;
}

// -----------------------------------------------------------------------------
bool GenericRegistrationFilter::IsSourcePointSet(int n) const
{
  // Note: Target points are transformed to source space, not vice versa!
  //       A point set can be both a source and a target (cf. inverse consistent energy).
  Array<PointSetDistanceInfo>::const_iterator it;
  for (it = _PointSetDistanceInfo.begin(); it != _PointSetDistanceInfo.end(); ++it) {
    if (it->_TargetIndex == n && !it->_TargetTransformation.IsForwardTransformation()) return true;
    if (it->_SourceIndex == n && !it->_SourceTransformation.IsForwardTransformation()) return true;
  }
  return false;
}

// -----------------------------------------------------------------------------
bool GenericRegistrationFilter::IsFixedPointSet(int n) const
{
  return !IsMovingPointSet(n);
}

// -----------------------------------------------------------------------------
bool GenericRegistrationFilter::IsMovingPointSet(int n) const
{
  // Note: Point set is moving if it is being transformed by any similarity term
  Array<PointSetDistanceInfo>::const_iterator it;
  for (it = _PointSetDistanceInfo.begin(); it != _PointSetDistanceInfo.end(); ++it) {
    // Note: For longitudinal registration, target point set deformed by FFD at t=0 may be compared to itself
    if (it->_TargetIndex == n && it->_SourceIndex != n && !it->_TargetTransformation.IsIdentity()) return true;
    if (it->_TargetIndex != n && it->_SourceIndex == n && !it->_SourceTransformation.IsIdentity()) return true;
  }
  return false;
}

// =============================================================================
// Parameter
// =============================================================================

// -----------------------------------------------------------------------------
void GenericRegistrationFilter::TransformationModel(enum TransformationModel model)
{
  _TransformationModel.resize(1, model);
}

// -----------------------------------------------------------------------------
bool GenericRegistrationFilter::AtInitialLevel() const
{
  return (_CurrentLevel == NumberOfLevels());
}

// -----------------------------------------------------------------------------
bool GenericRegistrationFilter::AtFinalLevel() const
{
  return (_CurrentLevel == _FinalLevel);
}

// -----------------------------------------------------------------------------
inline void rtrim(char *s)
{
  size_t l = strlen(s);
  if (l == 0) return;
  char *p = s + l - 1;
  while (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n') {
    *p = '\0';
    if (p == s) break;
    --p;
  }
}

// -----------------------------------------------------------------------------
static bool read_next(istream &in, char *buffer, int n, char *&name, char *&value, int &no)
{
  char *p = nullptr;
  // skip comments and blank lines
  do {
    if (!in) return false;
    in.getline(buffer, n);
    const size_t len = strlen(buffer);
    if (len > 0u && buffer[len-1] == '\r') buffer[len-1] = '\0';
    ++no;
    // discard leading whitespace characters
    name = buffer;
    while (*name == ' ' || *name == '\t') ++name;
    // discard comment and trailing whitespace characters at end of line
    if ((p = strstr(name, "#")) != nullptr) {
      if (p == name) {
        name[0] = '\0';
      } else if (p[-1] == ' ' || p[-1] == '\t') {
        *p = '\0';
        rtrim(name);
      }
    }
  } while (name[0] == '\0' || name[0] == '\n' || name[0] == '\r');
  // parse configuration section header
  if (name[0] == '[') {
    const size_t len = strlen(name);
    if (name[len-1] != ']') return false;
    // discard [ and leading whitespace characters
    name[0] = '\0';
    value = name + 1;
    while (*value == ' ' || *value == '\t') ++value;
    // discard ] and trailing whitespace characters
    name[len-1] = '\0';
    rtrim(value);
    // identify empty section name
    if (value[0] == '\0') return false;
    return true;
  }
  // find '=' character
  if ((p = strchr(name, '=')) == nullptr) return false;
  // skip leading whitespace characters of parameter value
  value = p;
  do { ++value; } while (*value == ' ' || *value == '\t');
  // truncate parameter name, skipping trailing whitespace characters
  *p = '\0';
  rtrim(name);
  // concatenate multi-line energy function value
  if (strcmp(name, "Energy function") == 0) {
    rtrim(value);
    while (true) {
      size_t len = strlen(value);
      if (len > 3 && strcmp(value + len - 3, "...") == 0) {
        // separate multi-line statements by single space (i.e. replace first '.' by ' ')
        value[len - 3] = ' ';
        // read next line into buffer at the end of current 'value'
        char *next_value = value + len - 2;
        in.getline(next_value, n - int(next_value - buffer));
        if (!in) return false;
        rtrim(next_value);
        ++no;
        // discard leading whitespace characters
        p = next_value;
        while (*p == ' ' || *p == '\t') ++p;
        if (p != next_value) {
          char *q = next_value;
          while (*p) {
            *q = *p;
            ++q, ++p;
          }
          *q = '\0';
        }
      } else break;
    }
  }
  return true;
}

// -----------------------------------------------------------------------------
bool GenericRegistrationFilter::Read(istream &from, bool echo)
{
  const size_t sz         = 4096;
  char         buffer[sz] = {0};
  char        *name       = NULL;
  char        *value      = NULL;
  int          no         = 0;
  int          level      = 0;
  bool         ok         = true;
  string       prefix, param;

  while (!from.eof()) {
    if (read_next(from, buffer, sz, name, value, no)) {
      if (name[0] == '\0') { // section header, value contains section name
        if (echo) cout << (no > 0 ? "\n" : "") << "[ " << value << " ]" << endl;
        level = -1;
        prefix.clear();
        const char c = toupper(value[0]);
        if      (c == 'R' && strncmp(value+1, "esolution level ", 16) == 0) value += 17;
        else if (c == 'L' && strncmp(value+1, "evel ",             5) == 0) value +=  6;
        else level = 0; // reset resolution level in case of other section
        if (level == -1) {
          // ignore any further text in section name following the level number
          char *p = strchr(value, ' ');
          if (p != NULL) *p = '\0';
          // parse resolution level number
          if (!FromString(value, level) || level < 0) {
            cerr << "Error in configuration at line: " << no << ": " << name << " = " << value << endl;
            ok = false;
          }
        } else {
          // section name ending with '...' indicates common prefix of following parameters
          const size_t len = strlen(value);
          if (len > 3u && strcmp(value + len - 3u, "...") == 0) {
            value[len - 3u] = '\0';
            prefix  = c;
            prefix += value + 1u;
          }
        }
      } else if (strcmp(name, "Level") == 0 || strcmp(name, "Resolution level") == 0) {
        if (!FromString(value, level)) {
          cerr << "Error in configuration at line: " << no << ": " << name << " = " << value << endl;
          ok = false;
        }
        if (level < 0) level = 0;
        if (echo) PrintParameter(cout, "\nResolution level", level);
      } else {
        if (prefix.empty()) {
          param = name;
        } else {
          param  = prefix;
          param += ' ';
          param += tolower(name[0]);
          param += name + 1;
        }
        if (!this->Set(param.c_str(), value, level)) {
          cerr << "Error in configuration at line: " << no << ": " << name << " = " << value << endl;
          ok = false;
        }
        if (echo) PrintParameter(cout, name, value);
      }
      if (level >= MAX_NO_RESOLUTIONS) {
        cerr << "Error in configuration at line: " << no << ": " << name << " = " << value << endl;
        cerr << "Maximum number of resolution levels is " << MAX_NO_RESOLUTIONS << endl;
        return false;
      }
    }
  }

  return ok;
}

// -----------------------------------------------------------------------------
// Note: Only set the specified level, which can be the zero-th level in case of
//       common settings for all levels. The settings for the remaining levels
//       are filled in by GuessParameter.
bool GenericRegistrationFilter::Set(const char *param, const char *value, int level)
{
  if (level < 0 || level >= MAX_NO_RESOLUTIONS) {
    cerr << this->NameOfType() << "::Set: Level index out-of-bounds" << endl;
    exit(1);
  }

  string name, type = ParameterUnits(param, &name);

  // Version
  if (name == "Version") {
    if (!FromString(value, version) || version > current_version) return false;
    if (!version) version = current_version;
    return true;

  // (Default) Similarity measure
  } else if (name == "Image (dis-)similarity measure" ||
             name == "Image dissimilarity measure" ||
             name == "Image similarity measure" ||
             name == "(Dis-)similarity measure" ||
             name == "Dissimilarity measure" ||
             name == "Similarity measure" ||
             name == "SIM") {
    return FromString(value, _SimilarityMeasure);
  } else if (name == "Precompute image derivatives") {
    return FromString(value, _PrecomputeDerivatives);

  // (Default) Point set distance measure
  } else if (name == "Point set distance measure" ||
             name == "Polydata distance measure" || // legacy
             name == "PDM") {
    return FromString(value, _PointSetDistanceMeasure);

  // Whether to remesh surfaces adaptively
  } else if (name == "Adaptive remeshing" ||
             name == "Adaptive surface remeshing") {
    return FromString(value, _AdaptiveRemeshing);

  // Transformation model
  } else if (name == "Transformation model") {
    _TransformationModel.clear();
    for (auto val : Split(value, '+')) {
      enum TransformationModel model;
      if (FromString(Trim(val), model)) {
        _TransformationModel.push_back(model);
      } else {
        _TransformationModel.clear();
        return false;
      }
    }
    return true;

  // Restrict deformation along coordinate axis
  // Note: E.g., useful for EPI distortion correction
  } else if (name == "Allow deformation in X") {
    return FromString(value, _RegisterX);
  } else if (name == "Allow deformation in Y") {
    return FromString(value, _RegisterY);
  } else if (name == "Allow deformation in Z") {
    return FromString(value, _RegisterZ);

  // Multi-level transformation model
  } else if (name == "Multi-level transformation" ||
             name == "Multi level transformation" ||
             name == "Multilevel transformation") {
    return FromString(value, _MultiLevelMode);

  // Energy function
  } else if (name == "Energy function" ||
             name == "Registration energy") {
    _EnergyFormula = value;
    return true;

  } else if (name == "Normalize weights of energy terms" ||
             name == "Normalise weights of energy terms") {
    return FromString(value, _NormalizeWeights);

  // Optimization method
  } else if (name == "Optimization method" ||
             name == "Optimisation method") {
    return FromString(value, _OptimizationMethod);

  // Number of resolution levels
  } else if (name == "No. of levels" ||
             name == "Number of levels" ||
             name == "No. of resolution levels" ||
             name == "Number of resolution levels" ||
             name == "First level" ||
             name == "First resolution level") {
    return FromString(value, _NumberOfLevels) && (_NumberOfLevels >= 1 || _NumberOfLevels == -1);

  } else if (name == "Final level" ||
             name == "Final resolution level" ||
             name == "Last level" ||
             name == "Last resolution level") {
    return FromString(value, _FinalLevel);

  // Whether to use Gaussian resolution pyramid
  } else if (name == "Use Gaussian resolution pyramid" ||
             name == "Use Gaussian image resolution pyramid") {
    bool use_pyramid;
    if (!FromString(value, use_pyramid)) return false;
    _UseGaussianResolutionPyramid = use_pyramid;
    return true;

  // Image resolution
  } else if (name.compare(0, 10, "Resolution", 10) == 0) {
    // Make parameter and value units consistent or return false in case of mismatch
    string number, units = ValueUnits(value, &number, type.c_str());
    if      (units == "rel") units = "vox";
    else if (units == "abs") units = "mm";
    if (type.empty()) {
      if (units.empty()) type = units = "signed";
      else type = units;
    } else {
      if      (type == "rel") type = "vox";
      else if (type == "abs") type = "mm";
    }
    if (type != units) {
      cerr << "Mismatching units specification for value of paramter '" << name << "'!" << endl;
      return false;
    }
    // Parse values
    double dx = .0, dy = .0, dz = .0;
    int n = ParseValues(number.c_str(), &dx, &dy, &dz);
    if (n == 0) return false;
    if (n == 1) dz = dy = dx;
    // Convert to signed "units"
    if (units != "signed") {
      if (dx < .0 || dy < .0 || dz < .0) {
        cerr << "Value of paramter '" << name << "' with units " << units << " must be positive!" << endl;
        return false;
      }
      if (units == "vox") {
        dx = - dx, dy = - dy, dz = - dz;
      } else if (units == "%") {
        dx = - 100.0 / dx;
        dy = - 100.0 / dy;
        dz = - 100.0 / dz;
      } else if (units != "mm") {
        cerr << "Units of parameter '" << name << "' must be either 'vox', 'mm', '%', or 'signed' (neg: 'vox', pos: 'mm')!" << endl;
        return false;
      }
    }
    // Assign parameter values
    n = 0; // used for image index next
    if (name.compare(10, 10, " of image ", 10) == 0) {
      if (!FromString(name.substr(20), n) || n < 1) return false;
      if (_Resolution[level].size() < static_cast<size_t>(n)) {
        _Resolution[level].resize(n, NaN);
      }
      --n;
    } else {
      _Resolution[level].resize(1);
    }
    _Resolution[level][n]._x = dx;
    _Resolution[level][n]._y = dy;
    _Resolution[level][n]._z = dz;
    return true;

  // Image blurring
  } else if (name.compare(0, 8, "Blurring", 8) == 0) {
    // Make parameter and value units consistent or return false in case of mismatch
    string number, units = ValueUnits(value, &number, type.c_str());
    if      (units == "rel") units = "vox";
    else if (units == "abs") units = "mm";
    if (type.empty()) {
      if (units.empty()) type = units = "signed";
      else type = units;
    } else {
      if      (type == "rel") type = "vox";
      else if (type == "abs") type = "mm";
    }
    if (type != units) {
      cerr << "Mismatching units specification for value of paramter '" << name << "'!" << endl;
      return false;
    }
    // Parse values
    double sigma = .0;
    if (!FromString(number, sigma)) return false;
    // Convert to signed "units"
    if (units != "signed") {
      if (sigma < .0) {
        cerr << "Value of parameter '" << name << "' with units " << units << " must be positive!" << endl;
        return false;
      }
      if (units == "vox") {
        sigma = - sigma;
      } else if (units == "%") {
        sigma = - sigma / 100.0;
      } else if (units != "mm") {
        cerr << "Units of parameter '" << name << "' must be either 'vox', 'mm', '%', or 'signed' (neg: 'vox', pos: 'mm')!" << endl;
        return false;
      }
    }
    // Assign value
    int n = 0;
    if (name.compare(8, 10, " of image ", 10) == 0) {
      if (!FromString(name.substr(18), n) || n < 1) return false;
      if (_Blurring[level].size() < static_cast<size_t>(n)) {
        _Blurring[level].resize(n, NaN);
      }
      --n;
    } else {
      _Blurring[level].resize(1);
    }
    _Blurring[level][n] = sigma;
    return true;

  // Image background
  } else if (name.compare(0, 16, "Background value", 16) == 0) {
    if (name.compare(16, 10, " of image ", 10) == 0) {
      int n = 0;
      if (!FromString(name.substr(26), n) || n < 1) return false;
      if (_Background.size() < static_cast<size_t>(n)) {
        _Background.resize(n, NaN);
      }
      return FromString(value, _Background[n-1]);
    } else {
      return FromString(value, _DefaultBackground);
    }

  // Image padding (deprecated option)
  } else if (name.compare(0, 13, "Padding value", 13) == 0) {
    if (name.compare(13, 10, " of image ", 10) == 0) {
      int n = 0;
      if (!FromString(name.substr(23), n) || n < 1) return false;
      if (_Background.size() < static_cast<size_t>(n)) {
        _Background.resize(n, NaN);
      }
      return FromString(value, _Background[n-1]);
    } else {
      return FromString(value, _DefaultBackground);
    }

  // Image interpolation
  } else if (name == "Image interpolation" || name == "Image interpolation mode" || name == "Interpolation mode") {
    return FromString(value, _DefaultInterpolationMode);
  } else if (name.compare(0, 23, "Interpolation of image ", 23) == 0) {
    int n = -1;
    if (!FromString(name.substr(23), n) || n < 1) return false;
    enum InterpolationMode mode;
    if (!FromString(value, mode)) return false;
    InterpolationMode(n, mode);
    return true;
  } else if (name.compare(0, 28, "Interpolation mode of image ", 28) == 0) {
    int n = -1;
    if (!FromString(name.substr(28), n) || n < 1) return false;
    enum InterpolationMode mode;
    if (!FromString(value, mode)) return false;
    InterpolationMode(n, mode);
    return true;

  // Image extrapolation
  } else if (name == "Image extrapolation" || name == "Image extrapolation mode" || name == "Extrapolation mode") {
    return FromString(value, _DefaultExtrapolationMode);
  } else if (name.compare(0, 23, "Extrapolation of image ", 23) == 0) {
    int n = -1;
    if (!FromString(name.substr(23), n) || n < 1) return false;
    enum ExtrapolationMode mode;
    if (!FromString(value, mode)) return false;
    ExtrapolationMode(n, mode);
    return true;
  } else if (name.compare(0, 28, "Extrapolation mode of image ", 28) == 0) {
    int n = -1;
    if (!FromString(name.substr(28), n) || n < 1) return false;
    enum ExtrapolationMode mode;
    if (!FromString(value, mode)) return false;
    ExtrapolationMode(n, mode);
    return true;

  // Image centering
  } else if (name == "Foreground-centric global transformation") {
    bool center;
    if (!FromString(value, center)) return false;
    _Centering[level] = center;
    return true;

  // Surface resolution
  } else if (name.compare(0, 11, "Edge length", 11) == 0 ||
             name.compare(0, 19, "Minimum edge length", 19) == 0 ||
             name.compare(0, 19, "Maximum edge length", 19) == 0) {

    const size_t nskip  = (name.compare(0, 11, "Edge length", 11) == 0 ? 11 : 19);
    const bool   setmin = (name.compare(0, 3, "Min") == 0);
    const bool   setmax = (name.compare(0, 3, "Max") == 0);

    string number, units = ValueUnits(value, &number, type.c_str());
    if (units.empty()) {
      units = "mm";
    } else if (units != type) {
      cerr << "Mismatching units specification for value of parameter '" << name << "'!" << endl;
      return false;
    }
    if (units != "mm") {
      cerr << "Value of parameter '" << name << "' must be in mm units!" << endl;
      return false;
    }

    double length;
    if (!FromString(number, length)) return false;
    if (length < .0) {
      cerr << "Value of parameter '" << name << "' must be positive!" << endl;
      return false;
    }

    int n = 0;
    if (name.compare(nskip, 12, " of surface ", 12) == 0) {
      if (!FromString(name.substr(nskip + 12), n) || n < 1) return false;
      if (setmin && _MinEdgeLength[level].size() < static_cast<size_t>(n)) {
        _MinEdgeLength[level].resize(n, -1.0);
      }
      if (setmax && _MaxEdgeLength[level].size() < static_cast<size_t>(n)) {
        _MaxEdgeLength[level].resize(n, -1.0);
      }
      --n;
    } else {
      if (setmin) _MinEdgeLength[level].resize(1);
      if (setmax) _MaxEdgeLength[level].resize(1);
    }

    if (setmin) _MinEdgeLength[level][n] = length;
    if (setmax) _MaxEdgeLength[level][n] = length;
    return true;

  // Merge global input transformation into FFD
  } else if (name == "Merge global and local transformation" ||
             name == "Merge global and local transformations") {
    return FromString(value, _MergeGlobalAndLocalTransformation);

  } else if (name == "Dirichlet boundary condition") {
    return FromString(value, _DirichletBoundaryCondition);

  // FFD control point spacing
  } else if (name.compare(0, 21, "Control point spacing")         == 0 ||
             name.compare(0, 29, "Minimum control point spacing") == 0 ||
             name.compare(0, 29, "Maximum control point spacing") == 0) {
    string number, units = ValueUnits(value, &number, type.c_str());
    if (units.empty()) {
      units = "mm";
    } else if (units != type) {
      cerr << "Mismatching units specification for value of parameter '" << name << "'!" << endl;
      return false;
    }
    if (units != "mm") {
      cerr << "Value of parameter '" << name << "' must be in mm units!" << endl;
      return false;
    }
    double ds = .0;
    if (!FromString(number, ds)) return false;
    if (ds < .0) {
      cerr << "Value of parameter '" << name << "' must be positive!" << endl;
      return false;
    }
    if        (name == "Control point spacing in X") {
      _MaxControlPointSpacing[level][0] = _MinControlPointSpacing[level][0] = ds;
    } else if (name == "Control point spacing in Y") {
      _MaxControlPointSpacing[level][1] = _MinControlPointSpacing[level][1] = ds;
    } else if (name == "Control point spacing in Z") {
      _MaxControlPointSpacing[level][2] = _MinControlPointSpacing[level][2] = ds;
    } else if (name == "Control point spacing in T") {
      _MaxControlPointSpacing[level][3] = _MinControlPointSpacing[level][3] = ds;
    } else if (name == "Control point spacing") {
      // Set only spatial dimensions, temporal resolution usually differs
      _MinControlPointSpacing[level][0] = _MaxControlPointSpacing[level][0] = ds;
      _MinControlPointSpacing[level][1] = _MaxControlPointSpacing[level][1] = ds;
      _MinControlPointSpacing[level][2] = _MaxControlPointSpacing[level][2] = ds;
    } else if (name == "Minimum control point spacing in X") {
      _MinControlPointSpacing[level][0] = ds;
    } else if (name == "Minimum control point spacing in Y") {
      _MinControlPointSpacing[level][1] = ds;
    } else if (name == "Minimum control point spacing in Z") {
      _MinControlPointSpacing[level][2] = ds;
    } else if (name == "Minimum control point spacing in T") {
      _MinControlPointSpacing[level][3] = ds;
    } else if (name == "Minimum control point spacing") {
      // Set only spatial dimensions, temporal resolution usually differs
      _MinControlPointSpacing[level][0] = ds;
      _MinControlPointSpacing[level][1] = ds;
      _MinControlPointSpacing[level][2] = ds;
    } else if (name == "Maximum control point spacing in X") {
      _MaxControlPointSpacing[level][0] = ds;
    } else if (name == "Maximum control point spacing in Y") {
      _MaxControlPointSpacing[level][1] = ds;
    } else if (name == "Maximum control point spacing in Z") {
      _MaxControlPointSpacing[level][2] = ds;
    } else if (name == "Maximum control point spacing in T") {
      _MaxControlPointSpacing[level][3] = ds;
    } else if (name == "Maximum control point spacing") {
      // Set only spatial dimensions, temporal resolution usually differs
      _MaxControlPointSpacing[level][0] = ds;
      _MaxControlPointSpacing[level][1] = ds;
      _MaxControlPointSpacing[level][2] = ds;
    } else {
      return false;
    }
    return true;

  // FFD subdivision
  } else if (name == "Subdivision dimension") {

    int dim = 0;
    if (FromString(value, dim)) {
      if        (dim == 0) { // no subdivision
        _Subdivide[level][0] = false;
        _Subdivide[level][1] = false;
        _Subdivide[level][2] = false;
        _Subdivide[level][3] = false;
      } else if   (dim == 2) { // 2D
        _Subdivide[level][0] = true;
        _Subdivide[level][1] = true;
        _Subdivide[level][2] = false;
        _Subdivide[level][3] = false;
      } else if (dim == 3) { // 3D
        _Subdivide[level][0] = true;
        _Subdivide[level][1] = true;
        _Subdivide[level][2] = true;
        _Subdivide[level][3] = false;
      } else if (dim == 4) { // 4D
        _Subdivide[level][0] = true;
        _Subdivide[level][1] = true;
        _Subdivide[level][2] = true;
        _Subdivide[level][3] = true;
      } else {
        return false;
      }
    } else { // e.g., xy, x+Y z, yxt,...

      _Subdivide[level][0] = false;
      _Subdivide[level][1] = false;
      _Subdivide[level][2] = false;
      _Subdivide[level][3] = false;
      for (const char *p = value; *p; ++p) {
        if      (*p == 'x' || *p == 'X') _Subdivide[level][0] = true;
        else if (*p == 'y' || *p == 'Y') _Subdivide[level][0] = true;
        else if (*p == 'z' || *p == 'Z') _Subdivide[level][0] = true;
        else if (*p == 't' || *p == 'T') _Subdivide[level][0] = true;
        else if (*p != '+' && *p != ' ' && *p != '\t') return false;
      }

    }

    // Because bool cannot be set to an "invalid" value such as most
    // other double parameters (which usually may not be zero), the
    // common subdivision setting of the "zero-th level" has to be
    // copied to the other levels already here. This requires that the
    // user sets/specifies the zero-th level settings before those of
    // the actual registration levels... this was anyway already required
    // by former MIRTK registration packages as well.
    if (level == 0) {
      for (int lvl = 1; lvl < MAX_NO_RESOLUTIONS; ++lvl) {
        _Subdivide[lvl][0] = _Subdivide[0][0];
        _Subdivide[lvl][1] = _Subdivide[0][1];
        _Subdivide[lvl][2] = _Subdivide[0][2];
        _Subdivide[lvl][3] = _Subdivide[0][3];
      }
    }
    return true;

  } else if (name == "Maximum rescaled intensity") {
    return FromString(value, _MaxRescaledIntensity);

  } else if (name == "Downsample images with padding") {
    return FromString(value, _DownsampleWithPadding);

  } else if (name == "Crop/pad images") {
    return FromString(value, _CropPadImages);

  } else if (name == "Crop/pad FFD lattice" ||
             name == "Crop/pad lattice") {
    bool do_crop_pad;
    if (!FromString(value, do_crop_pad)) return false;
    _CropPadFFD = (do_crop_pad ? 1 : 0);
    return true;

  // Sub-module parameter - store for later
  } else {
    // TODO: How can be checked which parameters are accepted by sub-modules
    //       before having them instantiated yet? Shall unknown parameters indeed
    //       simply be ignored? At least check once all modules are instantiated
    //       that at least one accepts each given parameter. -as12312
    Insert(_Parameter[level], param, value);
    return true;
  }

  return false;
}

// -----------------------------------------------------------------------------
bool GenericRegistrationFilter::Set(const char *name, const char *value)
{
  return this->Set(name, value, 0);
}

// -----------------------------------------------------------------------------
ParameterList GenericRegistrationFilter::Parameter(int level) const
{
  if (level < 0 || level > _NumberOfLevels) {
    cerr << "GenericRegistrationFilter::Parameter: Invalid resolution level: " << level << endl;
    exit(1);
  }
  ParameterList params = _Parameter[level];
  if (level == 0) {
    Remove(params, MINSTEP);
    Remove(params, MAXSTEP);
    if (version && version != current_version) Insert(params, "Version", version.ToString());
    string model;
    for (size_t i = 0; i < _TransformationModel.size(); ++i) {
      if (i > 0) model += '+';
      model += ToString(_TransformationModel[i]);
    }
    Insert(params, "Transformation model",                  model);
    Insert(params, "Multi-level transformation",            _MultiLevelMode);
    Insert(params, "Merge global and local transformation", _MergeGlobalAndLocalTransformation);
    Insert(params, "Dirichlet boundary condition", _DirichletBoundaryCondition);
    Insert(params, "Optimization method",                   _OptimizationMethod);
    Insert(params, "No. of resolution levels",              _NumberOfLevels);
    Insert(params, "Final level",                           _FinalLevel);
    Insert(params, "Precompute image derivatives",          _PrecomputeDerivatives);
    Insert(params, "Maximum rescaled intensities",          _MaxRescaledIntensity);
    Insert(params, "Normalize weights of energy terms",     _NormalizeWeights);
    Insert(params, "Downsample images with padding",        _DownsampleWithPadding);
    Insert(params, "Crop/pad images",                       _CropPadImages);
    if (_CropPadFFD != -1) {
      Insert(params, "Crop/pad lattice", _CropPadFFD != 0 ? true : false);
    }
    Insert(params, "Adaptive surface remeshing", _AdaptiveRemeshing);
    if (!_EnergyFormula.empty()) Insert(params, "Energy function", _EnergyFormula);
    if (_EnergyFormula.find("SIM") != string::npos) {
      Insert(params, "Image (dis-)similarity measure", _SimilarityMeasure);
    }
    if (_EnergyFormula.find("PDM") != string::npos) {
      Insert(params, "Point set distance measure", _PointSetDistanceMeasure);
    }
    if (NumberOfImages() > 0) {
      int n = 1;
      double bg = _Background[0];
      while (n < NumberOfImages() && AreEqualOrNaN(bg, _Background[n])) {
        ++n;
      }
      if (n == NumberOfImages()) {
        Insert(params, "Background value", bg);
      } else {
        for (n = 0; n < NumberOfImages(); ++n) {
          Insert(params, string("Background value of image ") + ToString(n + 1), _Background[n]);
        }
      }
      Insert(params, "Image interpolation", _DefaultInterpolationMode);
      Insert(params, "Image extrapolation", _DefaultExtrapolationMode);
      for (n = 0; n < NumberOfImages(); ++n) {
        if (static_cast<size_t>(n) < _InterpolationMode.size() && _InterpolationMode[n] != Interpolation_Default) {
          Insert(params, string("Interpolation of image ") + ToString(n + 1), _InterpolationMode[n]);
        }
      }
      for (n = 0; n < NumberOfImages(); ++n) {
        if (static_cast<size_t>(n) < _ExtrapolationMode.size() && _ExtrapolationMode[n] != Extrapolation_Default) {
          Insert(params, string("Extrapolation of image ") + ToString(n + 1), _ExtrapolationMode[n]);
        }
      }
      Insert(params, "Use Gaussian image resolution pyramid", _UseGaussianResolutionPyramid != 0 ? true : false);
    }
  }
  // Image pyramid
  if (level > 0 && NumberOfImages() > 0) {
    int n;
    // Resolution
    n = 1;
    Vector3D<double> res = _Resolution[level][0];
    while (n < NumberOfImages() && res == _Resolution[level][n]) ++n;
    if (n == NumberOfImages()) {
      Insert(params, "Resolution [mm]", ToString(res._x) + " " +
                                        ToString(res._y) + " " +
                                        ToString(res._z));
    } else {
      char name[64];
      for (int n = 0; n < NumberOfImages(); ++n) {
        snprintf(name, 64, "Resolution of image %d [mm]", n+1);
        Insert(params, name, ToString(_Resolution[level][n]._x) + " " +
                             ToString(_Resolution[level][n]._y) + " " +
                             ToString(_Resolution[level][n]._z));
      }
    }
    // Blurring
    n = 1;
    double sigma = _Blurring[level][0];
    while (n < NumberOfImages() && sigma == _Blurring[level][n]) ++n;
    if (n == NumberOfImages()) {
      Insert(params, "Blurring [mm]", ToString(sigma));
    } else {
      char name[64];
      for (int n = 0; n < NumberOfImages(); ++n) {
        snprintf(name, 64, "Blurring of image %d [mm]", n+1);
        Insert(params, name, ToString(_Blurring[level][n]));
      }
    }
  }
  // Control point spacing
  if (!IsLinear(_TransformationModel)) {
    const bool dim[4] = {
      _RegistrationDomain._x > 1,
      _RegistrationDomain._y > 1,
      _RegistrationDomain._z > 1,
      _RegistrationDomain._t > 1
    };
    double minds = .0;
    for (int d = 0; d < 4; ++d) {
      if (!dim[d]) continue;
      if (minds == .0) minds = _MinControlPointSpacing[level][d];
      else if (minds != _MinControlPointSpacing[level][d]) {
        minds = NaN;
        break;
      }
    }
    double maxds = .0;
    for (int d = 0; d < 4; ++d) {
      if (!dim[d]) continue;
      if (maxds == .0) maxds = _MaxControlPointSpacing[level][d];
      else if (maxds != _MaxControlPointSpacing[level][d]) {
        maxds = NaN;
        break;
      }
    }
    if (minds > .0 && maxds > .0) {
      if (minds == maxds) {
        Insert(params, "Control point spacing [mm]", ToString(minds));
      } else {
        string name;
        if (!IsNaN(minds)) {
          Insert(params, "Minimum control point spacing [mm]", ToString(minds));
        } else {
          for (int d = 0; d < 4; ++d) {
            if (!dim[d]) continue;
            name = "Minimum control point spacing in ";
            if      (d == 0) name += 'X';
            else if (d == 1) name += 'Y';
            else if (d == 2) name += 'Z';
            else if (d == 3) name += 'T';
            name += " [mm]";
            Insert(params, name, ToString(_MinControlPointSpacing[level][d]));
          }
        }
        if (!IsNaN(maxds)) {
          Insert(params, "Maximum control point spacing [mm]", ToString(maxds));
        } else {
          for (int d = 0; d < 4; ++d) {
            if (!dim[d]) continue;
            name = "Maximum control point spacing in ";
            if      (d == 0) name += 'X';
            else if (d == 1) name += 'Y';
            else if (d == 2) name += 'Z';
            else if (d == 3) name += 'T';
            name += " [mm]";
            Insert(params, name, ToString(_MaxControlPointSpacing[level][d]));
          }
        }
      }
    }
  }
  return params;
}

// -----------------------------------------------------------------------------
ParameterList GenericRegistrationFilter::Parameter() const
{
  return this->Parameter(0);
}

// -----------------------------------------------------------------------------
// Determine temporal attributes of set of 2D/3D input images/polydata using
// STL container which stores an ordered list of unique temporal coordinates
int GenericRegistrationFilter::NumberOfFrames(double *mint, double *maxt, double *avgdt) const
{
  if (NumberOfImages() > 0 || NumberOfPointSets() > 0) {
    OrderedSet<double> t;
    for (size_t n = 0; n < _Input.size(); ++n) {
      t.insert(_Input[n]->GetTOrigin());
    }
    for (size_t n = 0; n < _PointSetTime.size(); ++n) {
      t.insert(_PointSetTime[n]);
    }
    if (mint)  *mint  = (*t. begin());
    if (maxt)  *maxt  = (*t.rbegin());
    if (avgdt) *avgdt = AverageInterval(t);
    return static_cast<int>(t.size());
  }
  return 1;
}

// -----------------------------------------------------------------------------
Vector3D<double> GenericRegistrationFilter::AverageOutputResolution(int level) const
{
  if (level < 0) level = _CurrentLevel;
  Vector3D<double> res(.0, .0, .0);
  int              num = 0;

  for (size_t i = 0; i < _ImageSimilarityInfo.size(); ++i) {
    const int &t = _ImageSimilarityInfo[i]._TargetIndex;
    const int &s = _ImageSimilarityInfo[i]._SourceIndex;
    if (_ImageSimilarityInfo[i]._TargetTransformation) res += _Resolution[level][t], num += 1;
    if (_ImageSimilarityInfo[i]._SourceTransformation) res += _Resolution[level][s], num += 1;
  }
  if (num == 0) {
    if (_Domain) {
      res = _Resolution[level][0];
      num = 1;
    } else {
      const FreeFormTransformation *ffd = nullptr;
      if (_TargetTransformation) {
        if (!IsLinear(_TransformationModel)) {
          res._x = _MinControlPointSpacing[level][0];
          res._y = _MinControlPointSpacing[level][1];
          res._z = _MinControlPointSpacing[level][2];
        } else {
          const MultiLevelTransformation *mffd;
          mffd = dynamic_cast<const MultiLevelTransformation *>(_TargetTransformation);
          if (mffd) {
            if (mffd->NumberOfLevels() > 0) {
              ffd = mffd->GetLocalTransformation(-1);
            }
          } else {
            ffd = dynamic_cast<const FreeFormTransformation *>(_TargetTransformation);
          }
        }
      }
      if (ffd) {
        res._x = ffd->GetXSpacing();
        res._y = ffd->GetYSpacing();
        res._z = ffd->GetZSpacing();
      } else {
        double d = 1.0;
        if (level > 0) d *= pow(2.0, level - 1);
        if (_RegistrationDomain._x > 1) res._x = d;
        if (_RegistrationDomain._y > 1) res._y = d;
        if (_RegistrationDomain._z > 1) res._z = d;
      }
      num = 1;
    }
  }

  if (num > 1) res /= num;
  return res;
}

// -----------------------------------------------------------------------------
void GenericRegistrationFilter::ParseEnergyFormula(int nimages, int npsets, int nframes)
{
  string formula = _EnergyFormula;

  if (nimages < 0) nimages = NumberOfImages();
  if (npsets  < 0) npsets  = NumberOfPointSets();
  if (nframes < 1) nframes = NumberOfFrames();

  // Default registration energy function, i.e.,
  // either cross-sectional (N == 2) or longitudinal (N > 2) registration
  // of source image (frames) to first target image (frame)
  // with non-linear transformation constraints for regularization
  //
  // TODO: Automatically detect and group input channels with equal _torigin
  if (formula.empty()) {

    if (nimages > 0) {
      if (!formula.empty()) formula += " + ";
      if (nimages > 2) {
        if (nframes > 2 || (nframes == 2 && IsSpatioTemporal(_TransformationModel))) {
          formula += "SIM[Sim of I{s}](I(1), I(1:end) o T)";
        } else {
          formula += "SIM[Sim of I{s}](I(1), I(2:end) o T)";
        }
      } else {
        formula += "SIM[Image dissimilarity](I(1), I(2) o T)";
      }
    }
    if (npsets > 0) {
      if (!formula.empty()) formula += " + ";
      if (npsets > 2) {
        if (nframes > 2 || (nframes == 2 && IsSpatioTemporal(_TransformationModel))) {
          formula += "PDM[Dist of P{s}](T o P(1), P(1:end))";
        } else {
          formula += "PDM[Dist of P{s}](T o P(1), P(2:end))";
        }
      } else {
        formula += "PDM[Point set distance](T o P(1), P(2))";
      }
    }
    if (IsNonLinear(_TransformationModel)) {
      if (!formula.empty()) formula += " + ";
      const double be_w = ((nimages >= 2) ? 0.001 : .0);
      formula += ToString(be_w);
      formula +=     " BE[Bending energy](T)"
                 " + 0 LE[Linear energy](T)"
                 " + 0 TP[Topology preservation](T)"
                 " + 0 VP[Volume preservation](T)"
                 " + 0 LogJac[LogJac penalty](T)"
                 " + 0 NegJac[NegJac penalty](T)"
                 " + 0 Sparsity(T)";
    }
    formula += " + 0 MSDE[Displacement error](T)";
  }

  // Parse registration energy function
  RegistrationEnergyParser parser(this);
  parser.ParseEnergyFormula(formula, nimages, npsets);
}

// -----------------------------------------------------------------------------
void GenericRegistrationFilter::GuessParameter()
{
  const int _NumberOfImages    = NumberOfImages();
  const int _NumberOfPointSets = NumberOfPointSets();

  if (_NumberOfImages == 0 && _NumberOfPointSets == 0 && _TargetTransformation == nullptr) {
    cerr << "GenericRegistrationFilter::GuessParameter: Filter has no input data" << endl;
    exit(1);
  }
  for (int n = 0; n < _NumberOfImages; ++n) {
    if (_Input[n]->IsEmpty()) {
      cerr << "GenericRegistrationFilter::GuessParameter: Input image " << (n+1) << " is empty" << endl;
      exit(1);
    }
    if (_Input[n]->GetT() > 1) {
      cerr << "GenericRegistrationFilter::GuessParameter: Input image " << (n+1) << " has fourth dimension (_t > 1)." << endl;
      cerr << "  This registration filter only supports 2D/3D input images. Split multi-channel and/or" << endl;
      cerr << "  temporal image sequences into separate images if supported by this filter." << endl;
      exit(1);
    }
  }
  if (_NumberOfLevels < 1) _NumberOfLevels = (_NumberOfImages > 0 ? 4 : 1);
  if (_NumberOfLevels >= MAX_NO_RESOLUTIONS) {
    cerr << "GenericRegistrationFilter::Run: Maximum number of levels ("
         // Note that the "zero-th" level is only used to store default settings,
         // but the registration goes from level _NumberOfLevels to level 1.
         << (MAX_NO_RESOLUTIONS-1) << ") exceeded" << endl;
    exit(1);
  }
  if (_FinalLevel <= 0) _FinalLevel = 1;
  if (_FinalLevel > _NumberOfLevels) {
    Throw(ERR_InvalidArgument, __FUNCTION__, "Final level cannot be greater than total number of levels!");
  }

  // Default transformation model(s)
  if (_TransformationModel.empty()) {
    _TransformationModel.resize(3);
    _TransformationModel[0] = TM_Rigid;
    _TransformationModel[1] = TM_Affine;
    _TransformationModel[2] = TM_BSplineFFD;
  }

  // Determine temporal attributes of set of 2D/3D input images/point sets
  double mint, maxt, avgdt;
  int numt = NumberOfFrames(&mint, &maxt, &avgdt);

  // (Re-)parse energy formula (considering available input images)
  this->ParseEnergyFormula();

  // Set background values if undefined
  _Background.resize(_NumberOfImages, NaN);
  if (_NumberOfImages > 0) {
    SetDefaultBackgroundValue setbg(_Input, _Background, _DefaultBackground);
    parallel_for(blocked_range<int>(0,  _NumberOfImages), setbg);
  }

  // Initialize resolution pyramid
  const int nimages = max(_NumberOfImages, (_Domain ? 1 : 0));
  for (int level = 0; level <= _NumberOfLevels; ++level) {
    if (_Centering[level] == -1) {
      _Centering[level] = ((level == 0) ? true : _Centering[0]);
    }
    if (_Resolution[level].size() == 1) {
      _Resolution[level].resize(nimages, _Resolution[level][0]);
    } else {
      _Resolution[level].resize(nimages, Vector3D<double>(NaN));
    }
    if (_Blurring[level].size() == 1) {
      _Blurring[level].resize(nimages, _Blurring[level][0]);
    } else {
      _Blurring[level].resize(nimages, NaN);
    }
  }
  if (_UseGaussianResolutionPyramid == -1) {
    _UseGaussianResolutionPyramid = true;
    for (int level = 1; level <= _NumberOfLevels; ++level)
    for (int n     = 0; n     <  nimages;         ++n    ) {
      if ((!IsNaN(_Resolution[level][n]._x) && _Resolution[level][n]._x != .0) ||
          (!IsNaN(_Resolution[level][n]._y) && _Resolution[level][n]._y != .0) ||
          (!IsNaN(_Resolution[level][n]._z) && _Resolution[level][n]._x != .0)) {
        _UseGaussianResolutionPyramid = false;
      }
    }
  }
  for (int level = 0; level <= _NumberOfLevels; ++level)
  for (int n     = 0; n     <  nimages;         ++n) {
    const BaseImage * const image = (_Domain ? _Domain : _Input[n]);
    // Set initial image resolution
    if (level == 0) {
      // Resolution
      if (IsNaN(_Resolution[0][n]._x) || _Resolution[0][n]._x == .0) _Resolution[0][n]._x = -1.0;
      if (IsNaN(_Resolution[0][n]._y) || _Resolution[0][n]._y == .0) _Resolution[0][n]._y = -1.0;
      if (IsNaN(_Resolution[0][n]._z) || _Resolution[0][n]._z == .0) _Resolution[0][n]._z = -1.0;
      if (_Resolution[0][n]._x < .0) _Resolution[0][n]._x = abs(_Resolution[0][n]._x) * image->XSize();
      if (_Resolution[0][n]._y < .0) _Resolution[0][n]._y = abs(_Resolution[0][n]._y) * image->YSize();
      if (_Resolution[0][n]._z < .0) _Resolution[0][n]._z = abs(_Resolution[0][n]._z) * image->ZSize();
      if (image->Z() == 1) _Resolution[0][n]._z = .0;
    // Initialize lower resolution levels
    } else {
      // Resolution
      const double s = (level == 1 ? 1.0 : 2.0);
      if (_UseGaussianResolutionPyramid) {
        _Resolution[level][n]._x = s * _Resolution[level-1][n]._x;
        _Resolution[level][n]._y = s * _Resolution[level-1][n]._y;
        _Resolution[level][n]._z = s * _Resolution[level-1][n]._z;
      } else {
        if (_Resolution[level][n]._x < .0) {
          _Resolution[level][n]._x = abs(_Resolution[level][n]._x) * image->XSize();
        } else if (IsNaN(_Resolution[level][n]._x) || _Resolution[level][n]._x == .0) {
          _Resolution[level][n]._x = s * _Resolution[level-1][n]._x;
        }
        if (_Resolution[level][n]._y < .0) {
          _Resolution[level][n]._y = abs(_Resolution[level][n]._y) * image->YSize();
        } else if (IsNaN(_Resolution[level][n]._y) || _Resolution[level][n]._y == .0) {
          _Resolution[level][n]._y = s * _Resolution[level-1][n]._y;
        }
        if (_Resolution[level][n]._z < .0) {
          _Resolution[level][n]._z = abs(_Resolution[level][n]._z) * image->ZSize();
        } else if (IsNaN(_Resolution[level][n]._z) || _Resolution[level][n]._z == .0) {
          _Resolution[level][n]._z = s * _Resolution[level-1][n]._z;
        }
      }
      if (image->Z() == 1) {
        _Resolution[level][n]._z = _Resolution[0][n]._z;
      }
      // Blurring
      const double ds = max(max(_Resolution[level][n]._x,
                                _Resolution[level][n]._y),
                                _Resolution[level][n]._z);
      if (IsNaN(_Blurring[level][n])) {
        if (IsNaN(_Blurring[0][n])) {
          if (_UseGaussianResolutionPyramid && level > 1) {
            _Blurring[level][n] = .0;
          } else {
            _Blurring[level][n] = -.5;
          }
        } else {
          _Blurring[level][n] = _Blurring[0][n];
        }
      }
      if (_Blurring[level][n] < .0) {
        _Blurring[level][n] = abs(_Blurring[level][n]) * ds;
      }
    }
  }

  // Compute domain on which transformation is applied
  if (_Domain) {
    if (_CropPadImages) {
      _RegistrationDomain = _Domain->ForegroundDomain(0., true);
    } else {
      _RegistrationDomain = _Domain->Attributes();
    }
  } else {
    // Collect attributes of all "target" regions. A target region is either the
    // finite discrete domain of one input of the image similarity measure which
    // is not transformed at all, or the union of the finite discrete image domain
    // of both input images which are compared by the similarity measure in case
    // of a symmetric transformation model. Note that the image domain may further
    // be restricted to the minimal bounding region of the foreground only.
    Array<struct ImageAttributes> attrs;
    for (size_t i = 0; i < _ImageSimilarityInfo.size(); ++i) {
      const int &t = _ImageSimilarityInfo[i]._TargetIndex;
      const int &s = _ImageSimilarityInfo[i]._SourceIndex;
      if (_ImageSimilarityInfo[i].IsSymmetric() || !_ImageSimilarityInfo[i]._TargetTransformation) {
        if (!IsNaN(_Background[t]) && _CropPadImages) {
          Vector3D<double> sigma(0.);
          Vector3D<double> res = _Resolution[0][t];
          for (int l = 1; l <= _NumberOfLevels; ++l) {
            res._x = max(res._x, _Resolution[l][t]._x);
            res._y = max(res._y, _Resolution[l][t]._y);
            res._z = max(res._z, _Resolution[l][t]._z);
            sigma._x = max(sigma._x, _Blurring[l][t]);
            sigma._y = max(sigma._y, _Blurring[l][t]);
            sigma._z = max(sigma._z, _Blurring[l][t]);
          }
          if (_UseGaussianResolutionPyramid) {
            sigma._x = max(sigma._x, .5 * res._x);
            sigma._y = max(sigma._y, .5 * res._y);
            sigma._z = max(sigma._z, .5 * res._z);
          }
          attrs.push_back(ForegroundDomain(_Input[t], _Background[t], res, sigma, true));
        } else {
          attrs.push_back(OrthogonalFieldOfView(_Input[t]->Attributes()));
        }
      }
      if (_ImageSimilarityInfo[i].IsSymmetric() || !_ImageSimilarityInfo[i]._SourceTransformation) {
        if (!IsNaN(_Background[s]) && _CropPadImages) {
          Vector3D<double> sigma(0.);
          Vector3D<double> res = _Resolution[0][s];
          for (int l = 1; l <= _NumberOfLevels; ++l) {
            res._x = max(res._x, _Resolution[l][s]._x);
            res._y = max(res._y, _Resolution[l][s]._y);
            res._z = max(res._z, _Resolution[l][s]._z);
            sigma._x = max(sigma._x, _Blurring[l][s]);
            sigma._y = max(sigma._y, _Blurring[l][s]);
            sigma._z = max(sigma._z, _Blurring[l][s]);
          }
          if (_UseGaussianResolutionPyramid) {
            sigma._x = max(sigma._x, .5 * res._x);
            sigma._y = max(sigma._y, .5 * res._y);
            sigma._z = max(sigma._z, .5 * res._z);
          }
          attrs.push_back(ForegroundDomain(_Input[s], _Background[s], res, sigma, true));
        } else {
          attrs.push_back(OrthogonalFieldOfView(_Input[s]->Attributes()));
        }
      }
    }
    #if MIRTK_Registration_WITH_PointSet
      if (attrs.empty()) {
        // FIXME: PointSetDomain should estimate a more reasonable grid spacing
        double dx = 1.0, dy = 1.0, dz = 1.0;
        for (size_t i = 0; i < _PointSetDistanceInfo.size(); ++i) {
          const int &t = _PointSetDistanceInfo[i]._TargetIndex;
          const int &s = _PointSetDistanceInfo[i]._SourceIndex;
          if (_PointSetDistanceInfo[i]._TargetTransformation) {
            attrs.push_back(PointSetDomain(_PointSetInput[t], dx, dy, dz));
          }
          if (_PointSetDistanceInfo[i]._SourceTransformation) {
            attrs.push_back(PointSetDomain(_PointSetInput[s], dx, dy, dz));
          }
        }
      }
      #if MIRTK_Registration_WITH_Deformable
        if (attrs.empty()) {
          // FIXME: PointSetDomain should estimate a more reasonable grid spacing
          double dx = 1.0, dy = 1.0, dz = 1.0;
          for (size_t i = 0; i < _PointSetConstraintInfo.size(); ++i) {
            if (_PointSetConstraintInfo[i]._Transformation) {
              const int &t = _PointSetConstraintInfo[i]._PointSetIndex;
              attrs.push_back(PointSetDomain(_PointSetInput[t], dx, dy, dz));
            }
          }
        }
      #endif // MIRTK_Registration_WITH_Deformable
    #endif // MIRTK_Registration_WITH_PointSet
    // Remove any invalid domains and see if any valid ones are remaining
    Array<struct ImageAttributes>::iterator it;
    for (it = attrs.begin(); it != attrs.end(); ++it) {
      if (it->_x < 1 || it->_y < 1 || it->_z < 1) {
        Array<struct ImageAttributes>::iterator pos = it; --it;
        attrs.erase(pos);
      }
    }
    if (attrs.empty()) {
      if (_TargetTransformation) {
        const FreeFormTransformation *ffd = nullptr;
        const MultiLevelTransformation *mffd;
        mffd = dynamic_cast<const MultiLevelTransformation *>(_TargetTransformation);
        if (mffd) {
          if (mffd->NumberOfLevels() > 0) {
            ffd = mffd->GetLocalTransformation(-1);
          }
        } else {
          ffd = dynamic_cast<const FreeFormTransformation *>(_TargetTransformation);
        }
        if (ffd) {
          attrs.push_back(ffd->Attributes());
        }
      }
    }
    if (attrs.empty()) {
      cerr << "GenericRegistrationFilter::GuessParameter:";
      cerr << " Cannot determine domain of target input data!";
      cerr << " Is there any input given? If yes, try providing a mask image.";
      cerr << endl;
      exit(1);
    }
    // Now compute "minimal" finite grid which fully contains all "target" regions
    // This is the (minimal) domain for which the transformation has to be defined
    _RegistrationDomain = OverallFieldOfView(attrs);
  }
  // Set z spacing to zero for 2D domain
  if (_RegistrationDomain._z == 1) {
    _RegistrationDomain._dz = .0;
  }
  // Adjust temporal attributes
  if (_RegistrationDomain._t == 1) {
    if (numt > 2 || (numt == 2 && IsSpatioTemporal(_TransformationModel))) {
      _RegistrationDomain._t       = numt;
      _RegistrationDomain._torigin = mint;
      _RegistrationDomain._dt      = (maxt - mint) / (numt - 1);
    } else {
      _RegistrationDomain._dt = .0;
    }
  }

  // Average edge length range for surface meshes
  const int npointsets = NumberOfPointSets();
  if (_MinEdgeLength[0].size() == 1) {
    _MinEdgeLength[0].resize(npointsets, _MinEdgeLength[0][0]);
  } else {
    _MinEdgeLength[0].resize(npointsets, -1.);
  }
  if (_MaxEdgeLength[0].size() == 1) {
    _MaxEdgeLength[0].resize(npointsets, _MaxEdgeLength[0][0]);
  } else {
    _MaxEdgeLength[0].resize(npointsets, -1.);
  }
  for (int level = 1; level <= _NumberOfLevels; ++level) {
    const Vector3D<double> avgres = this->AverageOutputResolution(level);
    const double dmin = min(min(avgres._x, avgres._y), avgres._z);
    if (_MinEdgeLength[level].size() == 1) {
      _MinEdgeLength[level].resize(npointsets, _MinEdgeLength[level][0]);
    } else {
      _MinEdgeLength[level].resize(npointsets, -1.);
    }
    if (_MaxEdgeLength[level].size() == 1) {
      _MaxEdgeLength[level].resize(npointsets, _MaxEdgeLength[level][0]);
    } else {
      _MaxEdgeLength[level].resize(npointsets, -1.);
    }
    for (int n = 0; n < npointsets; ++n) {
      if (_MinEdgeLength[level][n] < 0.) {
        if (_MinEdgeLength[0][n] < 0.) {
          _MinEdgeLength[level][n] = (level == 1 ? 0. : dmin);
        } else {
          _MinEdgeLength[level][n] = _MinEdgeLength[0][n];
        }
      }
      if (_MaxEdgeLength[level][n] < 0.) {
        if (_MaxEdgeLength[0][n] < 0.) {
          _MaxEdgeLength[level][n] = 2. * sqrt(3.) * max(dmin, _MinEdgeLength[level][n]);
        } else {
          _MaxEdgeLength[level][n] = _MaxEdgeLength[0][n];
        }
      }
    }
  }
  for (int n = 0; n < npointsets; ++n) {
    const Vector3D<double> avgres = this->AverageOutputResolution(0);
    const double dmin = min(min(avgres._x, avgres._y), avgres._z);
    if (_MinEdgeLength[0][n] < 0.) {
      _MinEdgeLength[0][n] = dmin;
    }
    if (_MaxEdgeLength[0][n] < 0.) {
      _MaxEdgeLength[0][n] = 2. * sqrt(3.) * max(dmin, _MinEdgeLength[0][n]);
    }
  }

  // Final spatial FFD control point spacing
  const Vector3D<double> avgd = this->AverageOutputResolution(1);
  const double avgres[3] = { avgd._x, avgd._y, avgd._z };
  const double relres = (IsLinearFFD(_TransformationModel) ? 1.0 : 4.0);

  for (int d = 0; d < 3; ++d) {
    if (!_MinControlPointSpacing[0][d] && !_MaxControlPointSpacing[0][d]) {
      _MinControlPointSpacing[0][d] =     _MaxControlPointSpacing[0][d] = avgres[d] * relres;
    } else if (!_MinControlPointSpacing[0][d]) {
      _MinControlPointSpacing[0][d] = min(_MaxControlPointSpacing[0][d],  avgres[d]);
    } else if (!_MaxControlPointSpacing[0][d]) {
      _MaxControlPointSpacing[0][d] = max(_MinControlPointSpacing[0][d],  avgres[d] * relres);
    }
  }

  // Final temporal FFD control point spacing
  if (!_MinControlPointSpacing[0][3] && !_MaxControlPointSpacing[0][3]) {
    _MinControlPointSpacing[0][3] = _MaxControlPointSpacing[0][3] = avgdt;
  } else if (!_MinControlPointSpacing[0][3]) {
    _MinControlPointSpacing[0][3] = min(_MaxControlPointSpacing[0][3], avgdt);
  } else if (!_MaxControlPointSpacing[0][3]) {
    _MaxControlPointSpacing[0][3] = max(_MinControlPointSpacing[0][3], avgdt);
  }

  // By default, crop/pad FFD lattice only if domain not explicitly specified
  // and none of the transformation models is parameterised by velocities
  if (_CropPadFFD == -1) {
    _CropPadFFD = (_Domain == nullptr && !IsDiffeo(_TransformationModel) ? 1 : 0);
  }

  // Compute centers of foreground mass (if needed)
  bool centering = false;
  for (int l = 1; !centering && l <= _NumberOfLevels; ++l) {
    centering = (_Centering[l] != 0 ? true : false);
  }
  if (centering) {
    centering = (NumberOfImages() > 0);
    for (int n = 0; centering && n < NumberOfImages(); ++n) {
      centering = _Input[n]->GetAffineMatrix().IsIdentity();
    }
  }
  if (NumberOfImages() > 0 && (!_InitialGuess || centering)) {
    if (_Centroid.empty()) {
      Broadcast(LogEvent, "Computing centroids .....");
      _Centroid.resize(NumberOfImages());
      for (int n = 0; n < NumberOfImages(); ++n) {
        if (_Input[n]->CenterOfForeground(_Centroid[n], _Background[n]) == 0) {
          Broadcast(LogEvent, " failed\n");
          Throw(ERR_InvalidArgument, "Input image ", n + 1, " contains background only!");
        }
      }
      Broadcast(LogEvent, " done\n");
    }
  } else {
    _Centroid.clear();
  }

  // Optimization parameters
  double value;

  if (!Contains(_Parameter[0], REUSESTEP)) {
    // Start next line search with previously found optimal step length
    Insert(_Parameter[0], REUSESTEP, "Yes");
  }
  if (!Contains(_Parameter[0], STRICTRANGE)) {
    // Do not allow line search to reuse previously accepted step length
    // for next iteration when it is greater than the maximum allowed step.
    // This would reduce the running time of the optimization, but may not
    // find as good a solution as with smaller steps and more re-evaluations
    // of the gradient of the energy function.
    Insert(_Parameter[0], STRICTRANGE, "Yes");
  }
  // Default maximum no. of line search iterations, use the same defaults
  // as the previous rreg2, areg2, and nreg2 implementations have used
  const bool is_linear_model = IsLinear(_TransformationModel);
  int nlineiter = is_linear_model ? 20 : 12;
  if (Contains(_Parameter[0], MAXLINEITER)) {
    FromString(Get(_Parameter[0], MAXLINEITER), nlineiter);
  } else {
    // Insert at start of list as it may be overridden by other alternative
    // names for this setting such as, e.g., "No. of line iterations"...
    _Parameter[0].insert(_Parameter[0].begin(), make_pair(MAXLINEITER, ToString(nlineiter)));
  }
  if (!Contains(_Parameter[0], MAXREJECTED)) {
    // Interrupt line search sooner in case of deformable registration in order
    // to save function evaluations which with high chance are rejected anyway
    Insert(_Parameter[0], MAXREJECTED, max(nlineiter / 4, 2));
  }
  if (!Contains(_Parameter[0], MINSTEP) && !Contains(_Parameter[0], MAXSTEP)) {
    if (is_linear_model) {
      Insert(_Parameter[0], MINSTEP, .01);
      Insert(_Parameter[0], MAXSTEP, 1.);
    } else {
      // By default, limit step length to one (target) voxel unit
      const double maxres = max(avgres[0], max(avgres[1], avgres[2]));
      Insert(_Parameter[0], MINSTEP, maxres / 100.);
      Insert(_Parameter[0], MAXSTEP, maxres);
    }
  } else if (!Contains(_Parameter[0], MINSTEP)) {
    if (!FromString(Get(_Parameter[0], MAXSTEP).c_str(), value)) {
      cerr << "GenericRegistrationFilter::GuessParameter: Invalid '"
           << MAXSTEP << "' argument: " << Get(_Parameter[0], MAXSTEP) << endl;
      exit(1);
    }
    Insert(_Parameter[0], MINSTEP, value / 100.);
  } else if (!Contains(_Parameter[0], MAXSTEP)) {
    if (!FromString(Get(_Parameter[0], MINSTEP).c_str(), value)) {
      cerr << "GenericRegistrationFilter::GuessParameter: Invalid '"
           << MINSTEP << "' argument: " << Get(_Parameter[0], MINSTEP) << endl;
      exit(1);
    }
    Insert(_Parameter[0], MAXSTEP, value * 100.);
  }
  for (int level = 1; level <= _NumberOfLevels; ++level) {
    if (!Contains(_Parameter[level], MINSTEP) && !Contains(_Parameter[level], MAXSTEP)) {
      const auto   avgres = this->AverageOutputResolution(level);
      const double maxres = max(avgres._x, max(avgres._x, avgres._x));
      if (!FromString(Get(_Parameter[level-1], MINSTEP).c_str(), value)) {
        cerr << "GenericRegistrationFilter::GuessParameter: Invalid '"
             << MINSTEP << "' argument: " << Get(_Parameter[level-1], MINSTEP) << endl;
        exit(1);
      }
      Insert(_Parameter[level], MINSTEP, ToString(value));
      if (!FromString(Get(_Parameter[level-1], MAXSTEP).c_str(), value)) {
        cerr << "GenericRegistrationFilter::GuessParameter: Invalid '"
             << MAXSTEP << "' argument: " << Get(_Parameter[level-1], MAXSTEP) << endl;
        exit(1);
      }
      if (!is_linear_model && level > 1 && (value - maxres) < 1e-3) {
        value = 2. * value;
      }
      Insert(_Parameter[level], MAXSTEP, ToString(value));
    } else if (!Contains(_Parameter[level], MINSTEP)) {
      if (!FromString(Get(_Parameter[level], MAXSTEP).c_str(), value)) {
        cerr << "GenericRegistrationFilter::GuessParameter: Invalid '"
             << MAXSTEP << "' argument: " << Get(_Parameter[level], MAXSTEP) << endl;
        exit(1);
      }
      Insert(_Parameter[level], MINSTEP, value / 100.);
    } else if (!Contains(_Parameter[level], MAXSTEP)) {
      if (!FromString(Get(_Parameter[level], MINSTEP).c_str(), value)) {
        cerr << "GenericRegistrationFilter::GuessParameter: Invalid '"
             << MINSTEP << "' argument: " << Get(_Parameter[level], MINSTEP) << endl;
        exit(1);
      }
      Insert(_Parameter[level], MAXSTEP, value * 100.);
    }
  }
}

// -----------------------------------------------------------------------------
void GenericRegistrationFilter::Write(const char *fname) const
{
  ofstream to(fname);
  if (!to) {
    cerr << "Could not write configuration file " << fname << endl;
    exit(1);
  }
  PrintVersion(to, "## Version");
  for (int level = 0; level <= _NumberOfLevels; ++level) {
    to << "\n";
    if (level == 0) to << "[default]\n";
    else to << "[level " << level << "]\n";
    ParameterList params = this->Parameter(level);
    for (ParameterConstIterator it = params.begin(); it != params.end(); ++it) {
      PrintParameter(to, it->first, it->second);
    }
  }
  to.close();
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void GenericRegistrationFilter::Run()
{
  MIRTK_START_TIMING();

  // Guess parameters not specified by user
  this->GuessParameter();

  // Initialize image resolution pyramid
  this->InitializePyramid();
  this->InitializePointSets();

  // Make initial guess of transformation if none provided
  const Transformation * const dofin = _InitialGuess;
  if (!_InitialGuess && IsLinear(_TransformationModel.front())) {
    _InitialGuess = this->MakeInitialGuess();
  }

  MIRTK_DEBUG_TIMING(1, "initialization of registration");

  // For each transformation model (usually increasing number of DoFs)...
  Iteration model(0, static_cast<int>(_TransformationModel.size()));
  while (!model.End()) {
    _CurrentModel = _TransformationModel[model.Iter()];

    // Broadcast status message
    if (_TransformationModel.size() > 1) {
      string msg = "\n\nRegistration with ";
      msg += ToPrettyString(_CurrentModel);
      msg += " model\n";
      Broadcast(StatusEvent, msg.c_str());
    }

    // Run multi-resolution registration
    // (memorizing and restoring settings that will possibly be modified)
    bool merge = _MergeGlobalAndLocalTransformation;
    this->MultiResolutionOptimization();
    _MergeGlobalAndLocalTransformation = merge;

    // Delete previous transformation
    if (_InitialGuess != dofin) delete _InitialGuess;
    _InitialGuess   = _Transformation;
    _Transformation = NULL;

    // Continue with next transformation model
    ++model;
  }

  // Restore initial user guess
  _InitialGuess = dofin;

  MIRTK_DEBUG_TIMING(1, "registration");
}

// -----------------------------------------------------------------------------
void GenericRegistrationFilter::MultiResolutionOptimization()
{
  // For each resolution level (coarse to fine)...
  // Note: Zero-th level is used to store global registration settings
  Iteration level(_NumberOfLevels, _FinalLevel - 1);
  while (!level.End()) {
    _CurrentLevel = level.Iter();
    MIRTK_START_TIMING();

    // Initialize registration at current resolution
    Broadcast(InitEvent, &level);
    this->Initialize();

    // Solve registration problem by optimizing energy function
    Broadcast(StartEvent, &level);
    _Optimizer->Run();
    Broadcast(EndEvent, &level);

    // Finalize registration at current resolution
    this->Finalize();
    Broadcast(FinishEvent, &level);

    MIRTK_DEBUG_TIMING(2, "registration at level " << level.Iter());
    ++level;
  }
}

// -----------------------------------------------------------------------------
void GenericRegistrationFilter::InitializePyramid()
{
  MIRTK_START_TIMING();

  // Note: Level indices are in the range [1, N]
  const blocked_range  <int> images (0,  NumberOfImages());
  const blocked_range  <int> levels (1, _NumberOfLevels + 1);
  const blocked_range2d<int> pyramid(1, _NumberOfLevels + 1, 0, NumberOfImages());
  const blocked_range  <int> &level = images;

  // Allocate image list for each level even if empty
  _Image.resize(_NumberOfLevels + 1);

  if (NumberOfImages() > 0) {

    // Instantiate resolution pyramid
    for (int l = 1; l <= _NumberOfLevels; ++l) {
      _Image[l].resize(NumberOfImages());
    }

    // Use minimum intensity value to pad image unless user specified
    // a background value above the minimum intensity value. The background
    // value is not used here such that when user set no background value,
    // the blurring extends into the extended image margin.
    //
    // Note: Outside value always greater or equal background value.
    Array<double> outside(NumberOfImages());
    for (int n = 0; n < NumberOfImages(); ++n) {
      outside[n] = +inf;
      const int nvox = _Input[n]->NumberOfVoxels();
      for (int vox = 0; vox < nvox; ++vox) {
        outside[n] = min(outside[n], _Input[n]->GetAsDouble(vox));
      }
      if (IsInf(outside[n]) || outside[n] < _Background[n]) {
        outside[n] = _Background[n];
      }
    }

    // Copy/cast foreground of input images
    if (_CropPadImages) {
      Broadcast(LogEvent, "Crop/pad images .........");
      CropImages crop(_Input, _Background, outside, _Resolution[1], _Blurring[1], _Image[1]);
      parallel_for(images, crop);
    } else {
      Broadcast(LogEvent, "Padding images ..........");
      PadImages pad(_Input, _Background, _Image[1]);
      parallel_for(images, pad);
    }
    Broadcast(LogEvent, " done\n");

    MIRTK_DEBUG_TIMING(1, (_CropPadImages ? "cropping" : "padding") << " of images");
    MIRTK_RESET_TIMING();

    // Rescale input intensities from [min, max] to [eps, 1] such that image
    // gradient vectors used in transformation gradient computation have
    // comparable magnitude. This is only useful for symmetric and inverse
    // consistent registration energy functions where neither image should
    // have a stronger influence simply because of different intensity range.
    //
    // Note: The background is mapped to 0, hence the use of a small eps.
    const double _MinRescaledIntensity = 1e-3;
    if (_MaxRescaledIntensity > _MinRescaledIntensity && !IsInf(_MaxRescaledIntensity)) {
      Broadcast(LogEvent, "Rescaling images ........");
      for (int n = 0; n < NumberOfImages(); ++n) {
        _Image[1][n].ResetBackgroundValueAsDouble(NaN);
      }
      Rescale rescale(_Image[1], VoxelType(_MinRescaledIntensity), VoxelType(_MaxRescaledIntensity));
      parallel_for(images, rescale);
      for (int n = 0; n < NumberOfImages(); ++n) {
        _Image[1][n].ResetBackgroundValueAsDouble(0.);
        _Background[n] = outside[n] = 0.;
      }
      Broadcast(LogEvent, " done\n");
      MIRTK_DEBUG_TIMING(1, "rescaling of images");
      MIRTK_RESET_TIMING();
    }

    // Resample images to current resolution
    //
    // The images are in fact resampled after each step by the respective
    // RegisteredImage instance using the current transformation estimate.
    // As this class yet computes the image derivatives using finite differences
    // with step size equal to one voxel, the input images must still be
    // downsampled beforehand. Alternatively, if only downsampling by a factor
    // of two is applied at each resolution level, the finite differences step
    // size could be increased by this downsampling factor instead to save the
    // intermediate resampling and interpolation. Actually downsampling the images
    // beforehand, however, does not have the potential of surprise when inspecting
    // the debug output images and might result in buggy implementations which
    // assume the images would be downsampled during the initialization as done
    // also by previous registration packages.
    //
    // Downsampling with padding ensures that no background is "smeared" into
    // the foreground and is faster. However, in particular for low resolutions,
    // the edges of the images downsampled without considering the background
    // are smoother and may therefore better guide the optimization.
    //
    // On the other hand, at the finest resolution, sharp foreground boundaries
    // can result in high gradients which may lead to undesired results. Either
    // exclude background from gradient computation or smooth before finite
    // difference compuation (or use derivative of Gaussian filter).
    const Array<double> *padding = (_DownsampleWithPadding ? &_Background : nullptr);

    // Downsample (blur and resample by factor 2) if Gaussian pyramid is used.
    // Otherwise just copy first level to remaining levels which will be
    // blurred and resampled by the following processing steps. Optionally,
    // the images are cropped to minimal foreground size while being copied.
    // The recursive computation of the Gaussian pyramid is more efficient
    // than blurring and resampling the input image for each level separately.
    if (_UseGaussianResolutionPyramid && _NumberOfLevels > 1) {
      Broadcast(LogEvent, "Downsample images .......");
      if (debug_time) Broadcast(LogEvent, "\n");
    }
    for (int l = 2; l <= _NumberOfLevels; ++l) {
      if (_UseGaussianResolutionPyramid) {
        DownsampleImages downsample(_Image, l, &_Background, &outside, padding, &_Blurring[l], _CropPadImages);
        parallel_for(level, downsample);
      } else if (_CropPadImages) {
        CropImages crop(_Image[1], _Background, outside, _Resolution[l], _Blurring[l], _Image[l]);
        parallel_for(level, crop);
      } else {
        CopyImages copy(_Image[1], _Image[l]);
        parallel_for(level, copy);
      }
    }
    if (_UseGaussianResolutionPyramid && _NumberOfLevels > 1) {
      if (debug_time) Broadcast(LogEvent, "Downsample images .......");
      Broadcast(LogEvent, " done\n");
    }

    // Blur images (by default only if no Gaussian pyramid with implicit blurring is used)
    bool anything_to_blur = false;
    for (int l = 1; l <= _NumberOfLevels;   ++l)
    for (int n = 0; n <   NumberOfImages(); ++n) {
      if (_Blurring[l][n] > .0) anything_to_blur = true;
    }
    if (anything_to_blur) {
      Broadcast(LogEvent, "Blurring images .........");
      if (debug_time) Broadcast(LogEvent, "\n");
      BlurImages blur(_Image, _Blurring, padding);
      parallel_for(pyramid, blur);
      if (debug_time) Broadcast(LogEvent, "Blurring images .........");
      Broadcast(LogEvent, " done\n");
    }

    // Resample images after blurring if no Gaussian pyramid is used
    if (_UseGaussianResolutionPyramid) {
      for (int l = 1; l <= _NumberOfLevels;   ++l)
      for (int n = 0; n <   NumberOfImages(); ++n) {
        _Resolution[l][n]._x = _Image[l][n].XSize();
        _Resolution[l][n]._y = _Image[l][n].YSize();
        _Resolution[l][n]._z = _Image[l][n].ZSize();
      }
    } else {
      Broadcast(LogEvent, "Resample images .........");
      if (debug_time) Broadcast(LogEvent, "\n");
      ResampleImages resample(_Image, _Resolution, outside, padding);
      parallel_for(pyramid, resample);
      if (debug_time) Broadcast(LogEvent, "Resample images .........");
      Broadcast(LogEvent, " done\n");
    }

    // Set background value to be considered by RegisteredImage for
    // image gradient computation and image interpolation (resampling)
    for (int l = 1; l <= _NumberOfLevels;   ++l)
    for (int n = 0; n <   NumberOfImages(); ++n) {
      if (padding) {
        _Image[l][n].PutBackgroundValueAsDouble((*padding)[n]);
      } else {
        _Image[l][n].ClearBackgroundValue();
      }
    }
  } // if (NumberOfImages() > 0)

  // Resample domain mask
  _Mask.resize(_NumberOfLevels + 1, nullptr);
  if (_Domain) {
    if (_CropPadImages) {
      Broadcast(LogEvent, "Cropping mask ...........");
      if (debug_time) Broadcast(LogEvent, "\n");
      _Mask[0] = CropMask(_Domain);
      Broadcast(LogEvent, " done\n");
    } else {
      _Mask[0] = new BinaryImage(*_Domain);
      for (int vox = 0; vox < _Mask[0]->NumberOfVoxels(); ++vox) {
        if (_Mask[0]->Get(vox) != 0) _Mask[0]->Put(vox, 1);
      }
    }
    Broadcast(LogEvent, "Resample mask ...........");
    if (debug_time) Broadcast(LogEvent, "\n");
    Array<Vector3D<double> > res(_NumberOfLevels + 1);
    for (int l = 1; l <= _NumberOfLevels; ++l) {
      res[l] = this->AverageOutputResolution(l);
    }
    ResampleMask resample(_Mask[0], _Mask, res);
    parallel_for(levels, resample);
    if (debug_time) Broadcast(LogEvent, "Resample mask ...........");
    Broadcast(LogEvent, " done\n");
  }

  MIRTK_DEBUG_TIMING(1, "downsampling of images");
}

// -----------------------------------------------------------------------------
void GenericRegistrationFilter::InitializePointSets()
{
  _PointSet.resize(_NumberOfLevels + 1); // also when unused!
  if (NumberOfPointSets() == 0) return;

#if MIRTK_Registration_WITH_PointSet
  // Copy input point sets to every level (pointers only)
  for (int l = 0; l <= _NumberOfLevels; ++l) _PointSet[l] = _PointSetInput;

  // Remesh input surfaces
  for (int l = 1; l <= _NumberOfLevels; ++l) {
    RemeshSurfaces remesh(l, _PointSetInput, _PointSet, _MinEdgeLength, _MaxEdgeLength);
    parallel_for(blocked_range<int>(0, NumberOfPointSets()), remesh);
  }
#endif // MIRTK_Registration_WITH_PointSet
}

// -----------------------------------------------------------------------------
void GenericRegistrationFilter::Initialize()
{
  MIRTK_START_TIMING();

  // Initialize output of sub-registration at current resolution
  this->InitializeOutput();

  // Initialize registration energy function
  this->InitializeEnergy();

  // Initialize optimizer of registration energy
  this->InitializeOptimizer();

  MIRTK_DEBUG_TIMING(2, "initialization of level " << _CurrentLevel);
}

// -----------------------------------------------------------------------------
void GenericRegistrationFilter::InitializeStatus(HomogeneousTransformation *lin)
{
  const int ndofs = lin->NumberOfDOFs();
  if (!_RegisterX) {
    _Transformation->PutStatus(TX, Passive);
    _Transformation->PutStatus(RY, Passive);
    _Transformation->PutStatus(RZ, Passive);
    if (SX  < ndofs) _Transformation->PutStatus(SX,  Passive);
    if (SXY < ndofs) _Transformation->PutStatus(SXY, Passive);
    if (SXZ < ndofs) _Transformation->PutStatus(SXZ, Passive);
  }
  if (!_RegisterY) {
    _Transformation->PutStatus(TY, Passive);
    _Transformation->PutStatus(RX, Passive);
    _Transformation->PutStatus(RZ, Passive);
    if (SY  < ndofs) _Transformation->PutStatus(SY,  Passive);
    if (SXY < ndofs) _Transformation->PutStatus(SXY, Passive);
    if (SYZ < ndofs) _Transformation->PutStatus(SYZ, Passive);
  }
  if (!_RegisterZ) {
    _Transformation->PutStatus(TZ, Passive);
    _Transformation->PutStatus(RX, Passive);
    _Transformation->PutStatus(RY, Passive);
    if (SZ  < ndofs) _Transformation->PutStatus(SZ,  Passive);
    if (SXZ < ndofs) _Transformation->PutStatus(SXZ, Passive);
    if (SYZ < ndofs) _Transformation->PutStatus(SYZ, Passive);
  }
}

// -----------------------------------------------------------------------------
void GenericRegistrationFilter::InitializeStatus(FreeFormTransformation *ffd)
{
  const bool diffeo = IsDiffeo(_CurrentModel);

  if (_CropPadFFD || !diffeo) {
    if (_Mask[_CurrentLevel]) {

      // Initialize status of control points
      InitializeCPStatusGivenDomainMask init_status(_Mask[_CurrentLevel], ffd,
                                                    _RegisterX, _RegisterY, _RegisterZ);
      blocked_range<int> cps(0, ffd->NumberOfCPs());
      parallel_for(cps, init_status);

    } else {

      // Determine target data sets
      Array<bool> is_target_image(NumberOfImages());
      for (int n = 0; n < NumberOfImages(); ++n) {
        is_target_image[n] = !IsMovingImage(n);
      }
      #if MIRTK_Registration_WITH_PointSet
        Array<bool> is_moving_pointset(NumberOfPointSets());
        for (int n = 0; n < NumberOfPointSets(); ++n) {
          is_moving_pointset[n] = IsMovingPointSet(n);
        }
      #endif // MIRTK_Registration_WITH_PointSet

      // In case of fluid multi-level transformation, apply the global transformation
      // to the target images because the FFDs are defined on this transformed lattice
      UniquePtr<Matrix[]> smat;
      ResampledImageType *image;
      const FluidFreeFormTransformation *fluid;
      if ((fluid = dynamic_cast<const FluidFreeFormTransformation *>(_Transformation))) {
        smat.reset(new Matrix[NumberOfImages()]);
        for (int n = 0; n < NumberOfImages(); ++n) {
          if (is_target_image[n]) {
            image   = const_cast<ResampledImageType *>(&_Image[_CurrentLevel][n]);
            smat[n] = image->GetAffineMatrix();
            image->PutAffineMatrix(fluid->GetGlobalTransformation()->GetMatrix());
          }
        }
      }

      // Initialize status of control points
      InitializeCPStatus init_status(_Image[_CurrentLevel], is_target_image, ffd,
                                     #if MIRTK_Registration_WITH_PointSet
                                       _PointSet[_CurrentLevel], is_moving_pointset,
                                     #endif // MIRTK_Registration_WITH_PointSet
                                     _RegisterX, _RegisterY, _RegisterZ);
      blocked_range<int> cps(0, ffd->NumberOfCPs());
      parallel_for(cps, init_status);

      // Restore affine transformation matrices of input images
      if (smat) {
        for (int n = 0; n < NumberOfImages(); ++n) {
          if (is_target_image[n]) {
            image = const_cast<ResampledImageType *>(&_Image[_CurrentLevel][n]);
            image->PutAffineMatrix(smat[n]);
          }
        }
      }

    }

    // Discard passive DoFs to reduce memory/disk use and speed up computations
    if (_CropPadFFD) {
      const int margin = (diffeo ? 2 : 1) * ffd->KernelSize();
      ffd->CropPadPassiveCPs(margin, margin, margin, 0, true);
    }
  }

  // In case of a transformation parameterized by a (stationary) velocity
  // field, all control points are active as each of them may influence a
  // trajectory that starts (and ends) within the foreground image region.
  if (diffeo) {
    Transformation::DOFStatus sx, sy, sz;
    for (int cl = 0; cl < ffd->T(); ++cl)
    for (int ck = 0; ck < ffd->Z(); ++ck)
    for (int cj = 0; cj < ffd->Y(); ++cj)
    for (int ci = 0; ci < ffd->X(); ++ci) {
      sx = _RegisterX ? Active : Passive;
      sy = _RegisterY ? Active : Passive;
      sz = _RegisterZ ? Active : Passive;
      ffd->PutStatus(ci, cj, ck, cl, sx, sy, sz);
    }
  }

  // Enforce Dirichlet boundary condition by setting parameters at FFD boundary
  // to zero and forcing the CPs to be passive such that they are not modified
  //
  // Note: Condition not applied in temporal dimension as we usually only
  //       have at most one extra volume of CPs outside the period of the
  //       FFD time domain. Moreover, extrapolation could be periodic in time.
  if (_DirichletBoundaryCondition) {
    const int border_width = ffd->KernelRadius();
    for (int cl = 0; cl < ffd->T(); ++cl) {
      for (int m = 0; m < border_width; ++m) {
        for (int ck = m; ck < ffd->Z(); ck += ffd->Z() - m - 1)
        for (int cj = m; cj < ffd->Y(); cj += ffd->Y() - m - 1)
        for (int ci = m; ci < ffd->X(); ci += ffd->X() - m - 1) {
          ffd->Put(ci, cj, ck, cl, 0., 0., 0.);
          ffd->PutStatus(ci, cj, ck, cl, Passive, Passive, Passive);
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------
void GenericRegistrationFilter::InitializeStatus()
{
  MIRTK_START_TIMING();

  FreeFormTransformation    *ffd = NULL;
  MultiLevelTransformation *mffd = NULL;
  HomogeneousTransformation *lin = NULL;

  ( ffd = dynamic_cast<FreeFormTransformation    *>(_Transformation)) ||
  (mffd = dynamic_cast<MultiLevelTransformation  *>(_Transformation)) ||
  ( lin = dynamic_cast<HomogeneousTransformation *>(_Transformation));

  if (mffd) {
    Vector3D<double> res = this->AverageOutputResolution();
    for (int lvl = mffd->NumberOfLevels() - 1; lvl >= 0; --lvl) {
      if (!mffd->LocalTransformationIsActive(lvl)) continue;
      FreeFormTransformation *ffd = mffd->GetLocalTransformation(lvl);
      if (mffd->NumberOfActiveLevels() > 2 &&                // keep at least 2 active levels
          ((ffd->X() > 1 && ffd->GetXSpacing() < res._x) ||  // to be able to distinguish
           (ffd->Y() > 1 && ffd->GetYSpacing() < res._y) ||  // simultaneous optimization of
           (ffd->Z() > 1 && ffd->GetZSpacing() < res._z))) { // multiple levels from sequential
        mffd->LocalTransformationStatus(lvl, Passive);       // optimization
      } else {
        this->InitializeStatus(ffd);
      }
    }
  }
  else if (ffd) this->InitializeStatus(ffd);
  else if (lin) this->InitializeStatus(lin);

  MIRTK_DEBUG_TIMING(4, "initialization of " << (lin ? "DoF" : "control point") << " status");
}

// -----------------------------------------------------------------------------
TransformationType GenericRegistrationFilter::TransformationType()
{
  enum TransformationType type = ToTransformationType(_CurrentModel, _RegistrationDomain);
  if (type == TRANSFORMATION_UNKNOWN) {
    cerr << "GenericRegistrationFilter::TransformationType: Unknown transformation model: " << _CurrentModel << endl;
    exit(1);
  }
  return type;
}

// -----------------------------------------------------------------------------
void GenericRegistrationFilter::InitializeTransformation()
{
  MIRTK_START_TIMING();

  // ---------------------------------------------------------------------------
  // Determine if global initial guess must be merged into local transformation.
  // This is currently required by an inverse consistent registration which is
  // not based on a stationary velocity field transformation model.
  //
  // See also ApplyInitialGuess
  if (_CurrentModel != TM_BSplineSVFFD || _MultiLevelMode != MFFD_LogSum) {
    if (!_MergeGlobalAndLocalTransformation) {
      Array<ImageSimilarityInfo>::const_iterator sim;
      for (sim = _ImageSimilarityInfo.begin(); sim != _ImageSimilarityInfo.end(); ++sim) {
        if (sim->_SourceTransformation.IsBackwardTransformation() ||
            sim->_TargetTransformation.IsBackwardTransformation()) {
          _MergeGlobalAndLocalTransformation = true;
          break;
        }
      }
    }
    if (!_MergeGlobalAndLocalTransformation) {
      Array<PointSetDistanceInfo>::const_iterator dist;
      for (dist = _PointSetDistanceInfo.begin(); dist != _PointSetDistanceInfo.end(); ++dist) {
        if (dist->_SourceTransformation.IsBackwardTransformation() ||
            dist->_TargetTransformation.IsBackwardTransformation()) {
          _MergeGlobalAndLocalTransformation = true;
          break;
        }
      }
    }
  }

  // ---------------------------------------------------------------------------
  // Instantiate new transformation of selected model/type
  _Transformation = Transformation::New(this->TransformationType());

  // Set non-DoF parameters
  _Transformation->Parameter(_Parameter[0]);
  _Transformation->Parameter(_Parameter[_CurrentLevel]);

  FreeFormTransformation *ffd;
  ffd = dynamic_cast<FreeFormTransformation *>(_Transformation);
  if (!ffd) return;

  // ---------------------------------------------------------------------------
  // Instantiate multi-level free-form deformation of selected model/type
  MultiLevelTransformation *mffd = NULL;

  switch (_MultiLevelMode) {
    case MFFD_None: // discarded below again
    case MFFD_Default:
    case MFFD_Sum:
      mffd = new MultiLevelFreeFormTransformation();
      break;
    case MFFD_Fluid:
      mffd = new FluidFreeFormTransformation();
      break;
    case MFFD_LogSum:
      if (dynamic_cast<BSplineFreeFormTransformationSV *>(_Transformation) == NULL) {
        cout << endl;
        cerr << "Multi-level mode " << ToString(MFFD_LogSum)
        << " requires the BSplineSVFFD transformation model!" << endl;
        exit(1);
      }
      mffd = new MultiLevelStationaryVelocityTransformation();
      break;
    default:
      cout << endl;
      cerr << "The " << ToString(_MultiLevelMode) << " multi-level transformation mode is not supported" << endl;
      exit(1);
  }

  // Set non-DoF parameters
  mffd->Parameter(_Parameter[0]);
  mffd->Parameter(_Parameter[_CurrentLevel]);

  // ---------------------------------------------------------------------------
  // Compute domain on which free-form deformation must be defined
  struct ImageAttributes domain = _RegistrationDomain;

  if (_InitialGuess && !_InitialGuess->IsIdentity()) {
    const HomogeneousTransformation *ilin;
    ilin = dynamic_cast<const HomogeneousTransformation *>(_InitialGuess);
    if (!ilin) {
      const MultiLevelFreeFormTransformation *iffd;
      iffd = dynamic_cast<const MultiLevelFreeFormTransformation *>(_InitialGuess);
      if (iffd) ilin = iffd->GetGlobalTransformation();
    }
    if (ilin) {
      // In case the global transformation is to be merged into the local one,
      // make sure that the domain on which the local transformation is defined
      // is large enough to avoid unpleasant boundary effects
      if (_MergeGlobalAndLocalTransformation) {
        domain = ffd->ApproximationDomain(domain, _InitialGuess);
      // In case of fluid composition of global and local transformation,
      // i.e., T = T_local o T_global, linearly transform attributes such that
      // local transformation is defined for the mapped image region
      } else if (_MultiLevelMode == MFFD_Fluid) {
        domain.PutAffineMatrix(ilin->GetMatrix());
      }
    }
  }

  // ---------------------------------------------------------------------------
  // Initialize levels of multi-level transformation

  // Note: If minimum control point spacing is actually greater than the
  //       maximum control point spacing, the spacing in this dimension is
  //       decreased from level to level instead of increased. By swapping
  //       the order of these input parameters, the user can reverse the order
  //       of the levels in the multi-level free-form deformation stack.
  //       This makes no difference for the classical MFFD, but for the
  //       fluid transformation where the order of composition matters.

  const double * const min_ds = _MinControlPointSpacing[_CurrentLevel];
  const double * const max_ds = _MaxControlPointSpacing[_CurrentLevel];

  const MultiLevelTransformation *mffdin = NULL;
  const FreeFormTransformation   *ffdin  = NULL;
  (mffdin = dynamic_cast<const MultiLevelTransformation *>(_InitialGuess)) ||
  (ffdin  = dynamic_cast<const FreeFormTransformation   *>(_InitialGuess));
  if (mffdin && mffdin->NumberOfLevels() == 1) {
    ffdin = mffdin->GetLocalTransformation(0);
  }

  // If...
  if (// ...levels of FFD are not optimized simultaneously
      fequal(min_ds[0], max_ds[0]) && fequal(min_ds[1], max_ds[1]) && fequal(min_ds[2], max_ds[2]) &&
      // ...input FFD model is identical to output FFD model
      ffdin && strcmp(ffdin->NameOfClass(), ffd->NameOfClass()) == 0 &&
      // ...spacing of FFD models is identical
      fequal(ffdin->GetXSpacing(), max_ds[0]) &&
      fequal(ffdin->GetYSpacing(), max_ds[1]) &&
      fequal(ffdin->GetZSpacing(), max_ds[2]) &&
      // ...input MFFD model is identical to output MFFD model (if any used)
      (!mffdin || (_MultiLevelMode != MFFD_None && strcmp(mffdin->NameOfClass(), mffd->NameOfClass()) == 0))) {
    // then copy input transformation and continue optimizing it
    if (_MultiLevelMode == MFFD_None) {
      delete mffd;
      _Transformation = Transformation::New(ffdin);
      _Transformation->Parameter(_Parameter[0]);
      _Transformation->Parameter(_Parameter[_CurrentLevel]);
    } else {
      Transformation *copy = Transformation::New(ffdin);
      copy->Parameter(_Parameter[0]);
      copy->Parameter(_Parameter[_CurrentLevel]);
      mffd->PushLocalTransformation(dynamic_cast<FreeFormTransformation *>(copy));
      mffd->LocalTransformationStatus(mffd->NumberOfLevels()-1, Active);
      if (mffdin) {
        mffd->GetGlobalTransformation()->PutMatrix(mffdin->GetGlobalTransformation()->GetMatrix());
      }
      _Transformation = mffd;
    }
    MIRTK_DEBUG_TIMING(4, "copy of input transformation");
    return;
  }

  double ds  [4] = {.0, .0, .0, .0};
  double prev[4] = {.0, .0, .0, .0};
  double next[4] = {max_ds[0], max_ds[1], max_ds[2], max_ds[3]};
  bool   done = false, skip = false;

  // Extend domain to add a layer of (active) control points at the boundary
  if (_CropPadImages) {
    const double margin = (!_CropPadFFD && IsDiffeo(_CurrentModel) ? 2. : 1.) * ffd->KernelRadius();
    if (domain._x > 1) domain._dx += margin * max(min_ds[0], max_ds[0]) / (domain._x - 1);
    if (domain._y > 1) domain._dy += margin * max(min_ds[1], max_ds[1]) / (domain._y - 1);
    if (domain._z > 1) domain._dz += margin * max(min_ds[2], max_ds[2]) / (domain._z - 1);
  }
  // Add layer of (active) control points at the outermost time points
  // of a 4D FFD if it is not extended periodically by the extrapolator
  enum ExtrapolationMode m = ffd->ExtrapolationMode();
  const bool periodic = (m == ExtrapolationWithPeriodicTime(m));
  if (!periodic && _CropPadFFD) {
    if (domain._t > 1) {
      // Note: Image/FFD lattices are not centered in the temporal domain
      const double margin = ffd->KernelRadius();
      domain._dt      += margin * max(min_ds[3], max_ds[3]) / (domain._t - 1);
      domain._torigin -= margin * max(min_ds[3], max_ds[3]) / 2.0;
    }
  }

  // At least two control points in each (used) dimension
  const double threshold[4] = {
    (domain._x - 1) * domain._dx,
    (domain._y - 1) * domain._dy,
    (domain._z - 1) * domain._dz,
    (domain._t - 1) * domain._dt
  };

  do {
    skip = true;
    for (int d = 0; d < 4; ++d) {
      ds[d] = next[d];
      if (ds[d] > threshold[d]) ds[d] = threshold[d];
      if (ds[d] != prev[d]) skip  = false;
    }

    if (!skip) {
      // Instantiate free-form transformation
      // Note: At first iteration, use the previously instantiated object.
      if (mffd->NumberOfLevels() > 0) {
        _Transformation = Transformation::New(this->TransformationType());
        _Transformation->Parameter(_Parameter[0]);
        _Transformation->Parameter(_Parameter[_CurrentLevel]);
        // Type is the same as before, therefore no need for dynamic_cast
        ffd = reinterpret_cast<FreeFormTransformation *>(_Transformation);
      }
      // Initialize FFD with desired control point spacing
      ffd->Initialize(ffd->DefaultAttributes(domain, ds[0], ds[1], ds[2], ds[3]));
      memcpy(prev, ds, 4 * sizeof(double));
      // Push FFD onto MFFD stack
      mffd->PushLocalTransformation(ffd);
    }

    // Control point spacing of next level
    done = true;
    for (int d = 0; d < 4; ++d) {
      if (min_ds[d] < max_ds[d]) {
        next[d] = next[d] / 2.0;
        if (next[d] >= min_ds[d]) done = false;
      } else if (max_ds[d] < min_ds[d]) {
        next[d] = next[d] * 2.0;
        if (next[d] <= min_ds[d]) done = false;
      }
    }
  } while (!done);

  if (mffd->NumberOfLevels() == 0) {
    cout << endl;
    cerr << "Invalid registration domain! Try using a foreground mask image." << endl;
    exit(1);
  }

  // Make all levels active, status may be set to passive by InitializeStatus
  for (int l = 0; l < mffd->NumberOfLevels(); ++l) {
    mffd->LocalTransformationStatus(l, Active);
  }

  // ---------------------------------------------------------------------------
  // Discard multi-level transformation if multi-level mode is None
  if (_MultiLevelMode == MFFD_None) {
    if (mffd->NumberOfLevels() != 1) {
      cout << endl;
      cerr << "Multi-level transformation mode cannot be None if multiple levels are optimized simultaneously" << endl;
      exit(1);
    }
    _Transformation = mffd->PopLocalTransformation();
    delete mffd;
  } else {
    _Transformation = mffd;
  }

  MIRTK_DEBUG_TIMING(4, "instantiation of transformation");
}

// -----------------------------------------------------------------------------
Transformation *GenericRegistrationFilter::MakeInitialGuess()
{
  UniquePtr<Transformation> dofin;

  // Use identity transformation as initial guess for longitudinal registration
  if (_RegistrationDomain._t > 2 || (_RegistrationDomain._t == 2 && IsSpatioTemporal(_TransformationModel))) {
    return nullptr;
  }

  double tx = .0, ty = .0, tz = .0;
  int    t, s, n = 0;

  // ---------------------------------------------------------------------------
  Array<ImageSimilarityInfo>::const_iterator sim;
  for (sim = _ImageSimilarityInfo.begin(); sim != _ImageSimilarityInfo.end(); ++sim) {
    // Determine which of the two input images is the target
    t = -1;
    if (!sim->_TargetTransformation.IsForwardTransformation()) t = sim->_TargetIndex;
    if (!sim->_SourceTransformation.IsForwardTransformation()) t = (t == -1) ? sim->_SourceIndex : -1;
    // Skip if the energy formulation is ambiguous
    if (t == -1) continue;
    s = (t == sim->_TargetIndex) ? sim->_SourceIndex : sim->_TargetIndex;
    // Add displacement of image centroids
    tx += _Centroid[s]._x - _Centroid[t]._x;
    ty += _Centroid[s]._y - _Centroid[t]._y;
    tz += _Centroid[s]._z - _Centroid[t]._z;
    ++n;
  }

  // ---------------------------------------------------------------------------
  #if MIRTK_Registration_WITH_PointSet
    double tc[3], sc[3];
    Array<PointSetDistanceInfo>::const_iterator pdm;
    for (pdm = _PointSetDistanceInfo.begin(); pdm != _PointSetDistanceInfo.end(); ++pdm) {
      // Determine which of the two input data sets is the target
      t = -1;
      if (pdm->_TargetTransformation.IsForwardTransformation()) t = pdm->_TargetIndex;
      if (pdm->_SourceTransformation.IsForwardTransformation()) t = (t == -1) ? pdm->_SourceIndex : -1;
      // Skip if the energy formulation is ambiguous
      if (t == -1) continue;
      s = (t == pdm->_TargetIndex) ? pdm->_SourceIndex : pdm->_TargetIndex;
      // Add displacement of bounding box centers
      _PointSetInput[t]->GetCenter(tc);
      _PointSetInput[s]->GetCenter(sc);
      tx += sc[0] - tc[0];
      ty += sc[1] - tc[1];
      tz += sc[2] - tc[2];
      ++n;
    }
  #endif // MIRTK_Registration_WITH_PointSet

  // ---------------------------------------------------------------------------
  // Set initial translation to average displacement of input data centroids
  if (n > 0) {
    Broadcast(LogEvent, "Make initial guess ......");
    UniquePtr<RigidTransformation> rigid(new RigidTransformation());
    rigid->PutTranslationX(tx / n);
    rigid->PutTranslationY(ty / n);
    rigid->PutTranslationZ(tz / n);
    dofin.reset(rigid.release());
    Broadcast(LogEvent, " done\n");
  }

  return dofin.release();
}

// -----------------------------------------------------------------------------
void GenericRegistrationFilter::ApplyInitialGuess()
{
  // Just copy parameters whenever possible
  if (_Transformation->CopyFrom(_InitialGuess)) return;

  // Input...
  const HomogeneousTransformation        *ilin = NULL; // ...linear transformation
  const FreeFormTransformation           *iffd = NULL; // or non-linear  FFD
  const MultiLevelFreeFormTransformation *isum = NULL; // or multi-level FFD

  (ilin = dynamic_cast<const HomogeneousTransformation        *>(_InitialGuess)) ||
  (iffd = dynamic_cast<const FreeFormTransformation           *>(_InitialGuess)) ||
  (isum = dynamic_cast<const MultiLevelFreeFormTransformation *>(_InitialGuess));

  if (isum) {
    if (isum->NumberOfLevels() == 0) {
      ilin = isum->GetGlobalTransformation();
      isum = nullptr;
    } else if (isum->NumberOfLevels() == 1 && isum->GetGlobalTransformation()->IsIdentity()) {
      iffd = isum->GetLocalTransformation(0);
      isum = nullptr;
    }
  }

  // Output...
  HomogeneousTransformation        *olin  = NULL; // ...linear transformation
  FreeFormTransformation           *offd  = NULL; // or non-linear FFD
  MultiLevelTransformation         *omffd = NULL; // or multi-level FFD
  MultiLevelFreeFormTransformation *osum  = NULL; // (i.e., additive MFFD)

  ( olin = dynamic_cast<HomogeneousTransformation *>(_Transformation)) ||
  ( offd = dynamic_cast<FreeFormTransformation    *>(_Transformation)) ||
  (omffd = dynamic_cast<MultiLevelTransformation  *>(_Transformation));

  if (omffd) {
    const int nactive = omffd->NumberOfActiveLevels();
    if (nactive == 0) {
      cerr << "GenericRegistrationFilter::ApplyInitialGuess:"
              " Expected output MFFD to have at least one active level!" << endl;
      exit(1);
    } else if (nactive == 1) {
      for (int l = omffd->NumberOfLevels(); l >= 0; --l) {
        if (!omffd->LocalTransformationIsActive(l)) continue;
        offd = omffd->GetLocalTransformation(l);
      }
    }
    osum = dynamic_cast<MultiLevelFreeFormTransformation *>(omffd);
  }

  // Copy global transformation
  if (ilin) {
    if (olin) {
      olin->CopyFrom(ilin);
      return;
    } else if (omffd && !_MergeGlobalAndLocalTransformation) {
      omffd->GetGlobalTransformation()->CopyFrom(ilin);
      return;
    }
  // Copy local transformation
  } else if (iffd && offd) {
    if (offd->CopyFrom(iffd)) return;
  // Copy global and local transformation (additive MFFD only!)
  } else if (isum && isum->NumberOfLevels() == 1 &&
             osum && offd && !_MergeGlobalAndLocalTransformation) {
    osum->GetGlobalTransformation()->CopyFrom(isum->GetGlobalTransformation());
    if (offd->CopyFrom(isum->GetLocalTransformation(0))) return;
  }

  // Determine common attributes of (downsampled) target images
  struct ImageAttributes domain;
  if (_Mask[_CurrentLevel]) {
    domain = _Mask[_CurrentLevel]->Attributes();
  } else {
    Array<struct ImageAttributes> attrs;
    Array<ImageSimilarityInfo>::const_iterator sim;
    for (sim = _ImageSimilarityInfo.begin(); sim != _ImageSimilarityInfo.end(); ++sim) {
      if (sim->IsSymmetric() || sim->_TargetTransformation.IsIdentity()) {
        attrs.push_back(OrthogonalFieldOfView(this->ImageAttributes(sim->_TargetIndex)));
      }
      if (sim->IsSymmetric() || sim->_SourceTransformation.IsIdentity()) {
        attrs.push_back(OrthogonalFieldOfView(this->ImageAttributes(sim->_SourceIndex)));
      }
    }
    if (attrs.empty()) domain = this->RegistrationDomain();
    else               domain = OverallFieldOfView(attrs);
  }

  // Otherwise, approximate the initial guess
  double error = _Transformation->ApproximateAsNew(domain, _InitialGuess);
  ostringstream msg;
  msg << endl << "RMS error of initial guess approximation = " << error << endl;
  Broadcast(LogEvent, msg.str().c_str());
}

// -----------------------------------------------------------------------------
void GenericRegistrationFilter::InitializeOutput()
{
  MIRTK_START_TIMING();

  // Non-degenerated dimensions of registration domain
  const bool dim[4] = { _RegistrationDomain._x > 1,
                        _RegistrationDomain._y > 1,
                        _RegistrationDomain._z > 1,
                        _RegistrationDomain._t > 1 };

  // Control point spacing of final and current level
  double * const fmin = _MinControlPointSpacing[0];
  double * const fmax = _MaxControlPointSpacing[0];
  double * const cmin = _MinControlPointSpacing[_CurrentLevel];
  double * const cmax = _MaxControlPointSpacing[_CurrentLevel];

  // Tolerance for equality check of control point spacings
  double tol = 1e-6;
  for (int d = 0; d < 4; ++d) {
    if (!dim[d]) continue;
    tol = min(tol, 1e-6 * min(fmin[d], fmax[d]));
  }

  // Initialize transformation for first resolution level
  if (AtInitialLevel()) {

    // Leave no potential memory leaks...
    if (_Transformation != Output()) delete _Transformation;
    _Transformation = NULL;

    // Initial control point spacing
    for (int d = 0; d < 4; ++d) {
      if (!dim[d]) continue;
      double scale = 1.0;
      if (fequal(fmin[d], fmax[d], tol)) {
        for (int l = _CurrentLevel; l > 1; --l) {
          scale *= (_Subdivide[l][d] ? 2.0 : 1.0);
        }
      }
      if (!cmin[d]) cmin[d] = scale * fmin[d];
      if (!cmax[d]) cmax[d] = scale * fmax[d];
    }

    // Initialize transformation
    this->InitializeTransformation();

    // Apply original initial guess
    if (_InitialGuess && !_InitialGuess->IsIdentity()) {
      MIRTK_START_TIMING();
      this->ApplyInitialGuess();
      MIRTK_DEBUG_TIMING(4, "applying initial guess");
    }

  // Initialize transformation for consecutive levels
  } else {

    // Set non-DoF parameters of this level
    _Transformation->Parameter(_Parameter[_CurrentLevel]);

    // Memorize original initial guess and by default make previous
    // transformation the initial guess of the current level
    const Transformation * const initial_guess = _InitialGuess;
    _InitialGuess = _Transformation;

    // Determine type of previous transformation
    FreeFormTransformation   *ffd  = NULL;
    MultiLevelTransformation *mffd = NULL;

    (mffd = dynamic_cast<MultiLevelTransformation *>(_Transformation)) ||
    (ffd  = dynamic_cast<FreeFormTransformation   *>(_Transformation));

    // In case of a non-linear multi-level transformation...
    if (mffd) {
      // ...activate all levels again if simultaneously optimized
      if (mffd->NumberOfActiveLevels() > 1) { // c.f. InitializeStatus
        for (int l = 0; l < mffd->NumberOfLevels(); ++l) {
          mffd->LocalTransformationStatus(l, Active);
        }
      }
      // ...get last active local transformation level
      for (int l = mffd->NumberOfLevels() - 1; l >= 0; --l) {
        if (mffd->LocalTransformationIsActive(l)) ffd = mffd->GetLocalTransformation(l);
      }
    }

    // In case of a non-linear (multi-level) transformation
    if (ffd) {

      // Control point spacing of previous level
      const int      p    = _CurrentLevel + 1;
      double * const pmin = _MinControlPointSpacing[p];
      double * const pmax = _MaxControlPointSpacing[p];

      bool reuse  = false; // Whether to re-use    FFD
      bool subdiv = false; // Whether to subdivide FFD

      const int nactive = (mffd ? mffd->NumberOfActiveLevels() : 1);

      // If same control point configuration explicitly requested for this
      // resolution level, re-use previous transformation
      reuse = true;
      for (int d = 0; d < 4; ++d) {
        if (!dim[d]) continue;
        reuse = reuse && fequal(cmin[d], pmin[d], tol)
                      && fequal(cmax[d], pmax[d], tol);
      }

      // If no control point configuration was explicitly requested for this
      // level, and on the previous resolution level a transformation with
      // multiple active levels has been optimized, continue optimizing it
      if (!reuse) {
        reuse = (nactive > 1);
        for (int d = 0; d < 4; ++d) {
          if (!dim[d]) continue;
          reuse = reuse && !cmin[d] && !cmax[d];
        }
      }

      // If subdivision of a single (active) level (M)FFD is requested,
      // subdivide this level or keep it intact whenever possible
      //
      // Note: Uses subdivision settings of previous level.
      if (!reuse) {
        const double ds[4] = { _Subdivide[p][0] ? ffd->GetXSpacingAfterSubdivision() : ffd->GetXSpacing(),
                               _Subdivide[p][1] ? ffd->GetYSpacingAfterSubdivision() : ffd->GetYSpacing(),
                               _Subdivide[p][2] ? ffd->GetZSpacingAfterSubdivision() : ffd->GetZSpacing(),
                               _Subdivide[p][3] ? ffd->GetTSpacingAfterSubdivision() : ffd->GetTSpacing() };
        subdiv = (nactive == 1);
        for (int d = 0; d < 4; ++d) {
          if (!dim[d]) continue;
          subdiv = subdiv && (!cmin[d] || fequal(cmin[d], ds[d], tol))
                          && (!cmax[d] || fequal(cmax[d], ds[d], tol));
        }
      }

      // Subdivide FFD estimate of previous level whenever possible to
      // obtain initial transformation for current resolution level
      if (subdiv) {

        // Keeps FFD intact if all arguments are false or ignored (e.g., no subdivision in t for 3D FFD...)
        MIRTK_START_TIMING();
        ffd->Subdivide(_Subdivide[p][0], _Subdivide[p][1], _Subdivide[p][2], _Subdivide[p][3]);
        MIRTK_DEBUG_TIMING(4, "subdivision of FFD");

      // Otherwise, create new (M)FFD with desired control point spacing
      } else if (!reuse) {

        // Keep control point spacing of previous level if none specified for this level
        for (int d = 0; d < 4; ++d) {
          if (!cmin[d]) cmin[d] = pmin[d];
          if (!cmax[d]) cmax[d] = pmax[d];
        }

        // Initialize new (M)FFD with desired control point spacing
        _Transformation = NULL; // we still have the (m)ffd pointer to it
        this->InitializeTransformation();

        FreeFormTransformation   *cffd;
        MultiLevelTransformation *cmffd;
        cffd  = dynamic_cast<FreeFormTransformation   *>(_Transformation);
        cmffd = dynamic_cast<MultiLevelTransformation *>(_Transformation);

        // Push new FFD onto existing MFFD stack if single active level is optimized
        if (mffd) {
          for (int l = 0; l < mffd->NumberOfLevels(); ++l) {
            mffd->LocalTransformationStatus(l, Passive);
          }
          if (cmffd) {
            if (cmffd->NumberOfLevels() == 1) {
              mffd->PushLocalTransformation(cmffd->PopLocalTransformation());
              mffd->LocalTransformationStatus(-1, Active);
              ffd             = mffd->GetLocalTransformation(-1);
              _Transformation = mffd;
              delete cmffd;
            } else {
              mffd = cmffd;
              ffd  = NULL;
            }
          } else if (cffd) {
            mffd->PushLocalTransformation(cffd);
            mffd->LocalTransformationStatus(-1, Active);
            ffd             = cffd;
            _Transformation = mffd;
          } else {
            cout << endl;
            cerr << "GenericRegistrationFilter::InitializeOutput:"
                    " Expected new transformation to be (multi-level) FFD!" << endl;
            exit(1);
          }
        }
      }

      // Memorize actually used control point spacing at this resolution level
      if (mffd) {
        for (int l = mffd->NumberOfLevels() - 1; l >= 0; --l) {
          if (mffd->LocalTransformationIsActive(l)) ffd = mffd->GetLocalTransformation(l);
        }
      }
      if (!ffd) {
        cout << endl;
        cerr << "GenericRegistrationFilter::InitializeOutput:"
                " Expected at least one active FFD!" << endl;
        exit(1);
      }
      cmin[0] = ffd->GetXSpacing();
      cmin[1] = ffd->GetYSpacing();
      cmin[2] = ffd->GetZSpacing();
      cmin[3] = ffd->GetTSpacing();
      if (mffd) {
        for (int l = 0; l < mffd->NumberOfLevels(); ++l) {
          if (mffd->LocalTransformationIsActive(l)) ffd = mffd->GetLocalTransformation(l);
        }
      }
      cmax[0] = ffd->GetXSpacing();
      cmax[1] = ffd->GetYSpacing();
      cmax[2] = ffd->GetZSpacing();
      cmax[3] = ffd->GetTSpacing();
    }

    // Apply initial guess for the current resolution level and destroy
    // transformation of previous resolution level if not reused
    if (_InitialGuess != _Transformation) {
      if (!_InitialGuess->IsIdentity()) {
        MIRTK_START_TIMING();
        this->ApplyInitialGuess();
        MIRTK_DEBUG_TIMING(4, "applying initial guess");
      }
      Delete(_InitialGuess);
    }

    // Restore the pointer to the original initial guess
    _InitialGuess = initial_guess;
  }

  // In case of a linear transformation model, the centroids of the images are
  // moved to the origin of the world coordinate system such that the linear
  // operations (i.e., rotation, scaling, and shearing) occur around this point.
  HomogeneousTransformation *lin;
  _TargetOffset = _SourceOffset = Point(.0, .0, .0);
  if ((lin = dynamic_cast<HomogeneousTransformation *>(_Transformation))) {
    // FIXME: The centering as done in the following is only correct when
    //        the input images include no implicit affine (scanner) to
    //        other anatomy transformation. Once this can be applied also in
    //        such cases, then make sure to precompute the foreground centroids
    //        of the input images in InitializePyramid also in these cases.
    bool centering = (_Centering[_CurrentLevel] && !_Centroid.empty());
    for (int n = 0; centering && n < NumberOfImages(); ++n) {
      centering = _Input[n]->GetAffineMatrix().IsIdentity();
    }
    if (centering) {
      // Compute common centroid of target and source channel(s)
      int tn = 0, sn = 0;
      for (int n = 0; n < NumberOfImages(); ++n) {
        if (IsTargetImage(n)) _TargetOffset += _Centroid[n], ++tn;
        else                  _SourceOffset += _Centroid[n], ++sn;
      }
      if (tn > 0) _TargetOffset /= tn;
      if (sn > 0) _SourceOffset /= sn;

      // No translation in dimensions where this is disabled
      if (lin->GetStatus(TX) == Passive) {
        _TargetOffset._x = _SourceOffset._x = 0.;
      }
      if (lin->GetStatus(TY) == Passive) {
        _TargetOffset._y = _SourceOffset._y = 0.;
      }
      if (lin->GetStatus(TZ) == Passive) {
        _TargetOffset._z = _SourceOffset._z = 0.;
      }
      if (AreEqual(_TargetOffset.Distance(), 0.) &&
          AreEqual(_SourceOffset.Distance(), 0.)) {
        _TargetOffset = _SourceOffset = 0.;
      } else {

        // Translate origin of images such that
        for (int n = 0; n < NumberOfImages(); ++n) {
          Point origin = _Image[_CurrentLevel][n].GetOrigin();
          if (IsTargetImage(n)) {
            _Image[_CurrentLevel][n].PutOrigin(origin - _TargetOffset);
          } else {
            _Image[_CurrentLevel][n].PutOrigin(origin - _SourceOffset);
          }
        }
        if (_Mask[_CurrentLevel]) {
          Point origin = _Mask[_CurrentLevel]->GetOrigin();
          _Mask[_CurrentLevel]->PutOrigin(origin - _TargetOffset);
        }

        // Translates points of point sets
        #if MIRTK_Registration_WITH_PointSet
          for (int n = 0; n < NumberOfPointSets(); ++n) {
            double p[3], offset[3];
            if (IsTargetPointSet(n)) {
              offset[0] = _TargetOffset._x;
              offset[1] = _TargetOffset._y;
              offset[2] = _TargetOffset._z;
            } else {
              offset[0] = _SourceOffset._x;
              offset[1] = _SourceOffset._y;
              offset[2] = _SourceOffset._z;
            }
            vtkPoints *points = _PointSet[_CurrentLevel][n]->GetPoints();
            for (vtkIdType ptId = 0; ptId < points->GetNumberOfPoints(); ++ptId) {
              points->GetPoint(ptId, p);
              p[0] -= offset[0];
              p[1] -= offset[1];
              p[2] -= offset[2];
              points->SetPoint(ptId, p);
            }
          }
        #endif // MIRTK_Registration_WITH_PointSet

        // Adjust linear transformation
        Matrix pre(4, 4);
        pre.Ident();
        pre(0, 3)  = + _TargetOffset._x;
        pre(1, 3)  = + _TargetOffset._y;
        pre(2, 3)  = + _TargetOffset._z;

        Matrix post(4, 4);
        post.Ident();
        post(0, 3) = - _SourceOffset._x;
        post(1, 3) = - _SourceOffset._y;
        post(2, 3) = - _SourceOffset._z;

        lin->PutMatrix(post * lin->GetMatrix() * pre);
      }
    }
  }

  // Set status of transformation parameters (DoFs)
  this->InitializeStatus();

  // Verify transformation is well-initialised
  _Transformation->Verify();

  // Memorize (all) used (non-DoF) transformation settings for -parout file
  Insert(_Parameter[_CurrentLevel], _Transformation->Parameter());

  MIRTK_DEBUG_TIMING(3, "initialization of output");
}

// -----------------------------------------------------------------------------
ImageAttributes GenericRegistrationFilter::ImageAttributes(int n, int level) const
{
  // Get desired resolution of registered image at given level
  if (level < 0) level = _CurrentLevel;
  const Vector3D<double> &res = _Resolution[level][n];
  // Use attributes of previously resampled image if possible
  if (res._x == _Image[level][n].GetXSize() &&
      res._y == _Image[level][n].GetYSize() &&
      res._z == _Image[level][n].GetZSize()) {
    return _Image[level][n].Attributes();
  // Derive attributes with desired resolution from attributes of final level
  } else {
    struct ImageAttributes attr = _Image[1][n].Attributes();

    attr._x = iround(attr._x * attr._dx / res._x);
    attr._y = iround(attr._y * attr._dy / res._y);
    attr._z = iround(attr._z * attr._dz / res._z);

    if (attr._x < 1) attr._x  = 1;
    else             attr._dx = res._x;
    if (attr._y < 1) attr._y  = 1;
    else             attr._dy = res._y;
    if (attr._z < 1) attr._z  = 1;
    else             attr._dz = res._z;

    return attr;
  }
}

// -----------------------------------------------------------------------------
ImageAttributes GenericRegistrationFilter::RegistrationDomain(int level) const
{
  if (level < 0) level = _CurrentLevel;

  struct ImageAttributes attr = _RegistrationDomain;
  Vector3D<double> res = this->AverageOutputResolution(level);

  attr._x = iround(attr._x * attr._dx / res._x);
  attr._y = iround(attr._y * attr._dy / res._y);
  attr._z = iround(attr._z * attr._dz / res._z);

  if (attr._x < 1) attr._x  = 1;
  else             attr._dx = res._x;
  if (attr._y < 1) attr._y  = 1;
  else             attr._dy = res._y;
  if (attr._z < 1) attr._z  = 1;
  else             attr._dz = res._z;

  return attr;
}

// -----------------------------------------------------------------------------
Transformation *GenericRegistrationFilter::OutputTransformation(TransformationInfo ti)
{
  // Identity, i.e., no transformation: T^0
  if (ti == TransformationInfo::Identity()) return NULL;

  // Forward transformation: T, T^1, T^+1
  if (ti == TransformationInfo::Full()) return _Transformation;

  // Return previously created instance if one exists
  for (size_t i = 0; i < _TransformationInstance.size(); ++i) {
    if (_TransformationInfo[i] == ti) return _TransformationInstance[i];
  }

  // Otherwise, instantiate new transformation object
  HomogeneousTransformation                  *lin    = NULL;
  AffineTransformation                       *affine = NULL;
  BSplineFreeFormTransformationSV            *svffd  = NULL;
  MultiLevelStationaryVelocityTransformation *msvffd = NULL;
  Transformation                             *dof    = NULL;

  // Inverse transformation: T^-1
  if (ti == TransformationInfo::Inverse()) {
    lin    = dynamic_cast<HomogeneousTransformation                  *>(_Transformation);
    affine = dynamic_cast<AffineTransformation                       *>(_Transformation);
    svffd  = dynamic_cast<BSplineFreeFormTransformationSV            *>(_Transformation);
    msvffd = dynamic_cast<MultiLevelStationaryVelocityTransformation *>(_Transformation);
    if (affine) {
      // More efficient than PartialAffineTransformation
      dof = new InverseAffineTransformation(affine);
    } else if (lin) {
      dof = new PartialAffineTransformation(lin, -1.0);
    } else if (msvffd) {
      dof = new PartialMultiLevelStationaryVelocityTransformation(msvffd, -1.0);
    } else if (svffd) {
      dof = new PartialBSplineFreeFormTransformationSV(svffd, -1.0);
    } else {
      cerr << "Cannot currently use inverse transformation of selected transformation model during optimization." << endl;
      if (_MultiLevelMode != MFFD_None) {
        cerr << "Try again with parameter \"Multi-level transformation = None\"." << endl;
      }
      exit(1);
    }
  // Partial transformation: e.g., T^-0.5, T^0.5
  } else {
    lin    = dynamic_cast<HomogeneousTransformation                  *>(_Transformation);
    svffd  = dynamic_cast<BSplineFreeFormTransformationSV            *>(_Transformation);
    msvffd = dynamic_cast<MultiLevelStationaryVelocityTransformation *>(_Transformation);
    if (lin) {
      dof = new PartialAffineTransformation(lin, ti._Exponent);
    } else if (msvffd) {
      dof = new PartialMultiLevelStationaryVelocityTransformation(msvffd, ti._Exponent);
    } else if (svffd) {
      dof = new PartialBSplineFreeFormTransformationSV(svffd, ti._Exponent);
    } else {
      cerr << "Cannot currently use partial (inverse) transformation of selected transformation model during optimization." << endl;
      if (_MultiLevelMode != MFFD_None) {
        cerr << "Try again with parameter \"Multi-level transformation = None\"." << endl;
      }
      exit(1);
    }
  }

  // Keep pointer to newly instantiated auxiliary object which will be deleted
  // again during finalization of the registration at the current level (cf. Finalize)
  _TransformationInfo    .push_back(ti);
  _TransformationInstance.push_back(dof);
  return dof;
}

// -----------------------------------------------------------------------------
void GenericRegistrationFilter::SetInputOf(RegisteredImage *output, const struct ImageAttributes &domain, int n, TransformationInfo ti)
{
  mirtkAssert(0 <= n && size_t(n) < _Image[_CurrentLevel].size(), "input image index is within bounds");

  // Set input of registered image
  output->InputImage           (&_Image[_CurrentLevel][n]);
  output->InterpolationMode    (InterpolationMode(n));
  output->ExtrapolationMode    (ExtrapolationMode(n));
  output->PrecomputeDerivatives(_PrecomputeDerivatives);
  output->Transformation       (this->OutputTransformation(ti));
  output->PutBackgroundValueAsDouble(_Background[n]);

  // Add/Amend displacement cache entry
  if (output->Transformation() && output->Transformation()->RequiresCachingOfDisplacements()) {
    const double t = _Image[_CurrentLevel][n].GetTOrigin();
    Array<DisplacementInfo>::iterator i;
    for (i = _DisplacementInfo.begin(); i != _DisplacementInfo.end(); ++i) {
      if (i->_Transformation == output->Transformation() &&
          i->_Domain.EqualInSpace(domain) &&
          fequal(i->_Domain._torigin, output->GetTOrigin(), 1e-9) &&
          (fequal(i->_InputTime, t, 1e-9) || IsNaN(i->_InputTime))) {
        i->_InputTime = t;
        break;
      }
    }
    if (i == _DisplacementInfo.end()) {
      _DisplacementInfo.resize(_DisplacementInfo.size() + 1);
      i = _DisplacementInfo.end() - 1;
      i->_InputTime      = t;
      i->_DispIndex      = static_cast<int>(_DisplacementField.size());
      i->_Domain         = domain;
      i->_Transformation = output->Transformation();
      _DisplacementField.push_back(new DisplacementImageType(domain, 3));
    }
    output->ExternalDisplacement(_DisplacementField[i->_DispIndex]);
  }
}

#if MIRTK_Registration_WITH_PointSet

// -----------------------------------------------------------------------------
static bool ContainsPoints(const ImageAttributes &domain, vtkPoints *points)
{
  double p[3];
  for (vtkIdType id = 0; id < points->GetNumberOfPoints(); ++id) {
    points->GetPoint(id, p);
    domain.WorldToLattice(p[0], p[1], p[2]);
    if (p[0] <= -.5 || p[0] >= domain._x - .5 ||
        p[1] <= -.5 || p[1] >= domain._y - .5 ||
        p[2] <= -.5 || p[2] >= domain._z - .5) {
      return false;
    }
  }
  return true;
}

// -----------------------------------------------------------------------------
RegisteredPointSet *GenericRegistrationFilter::OutputPointSet(int n, double t, TransformationInfo ti)
{
  mirtkAssert(0 <= n && size_t(n) < _PointSetInput.size(), "input point set index is within bounds");

  // Return existing output point set
  for (size_t i = 0; i < _PointSetOutput.size(); ++i) {
    if (_PointSetOutputInfo[i]._InputIndex     == n &&
        _PointSetOutputInfo[i]._Transformation == ti &&
        fequal(_PointSetOutput[i]->Time(), t, 1e-9)) {
      return _PointSetOutput[i];
    }
  }

  // Otherwise, instantiate new output point set
  RegisteredPointSet *output = new RegisteredPointSet();
  output->InputPointSet (_PointSet[_CurrentLevel][n]);
  output->InputTime     (_PointSetTime[n]);
  output->Time          (t);
  output->Transformation(this->OutputTransformation(ti));
  output->SelfUpdate    (false);

  if (output->Transformation() && output->Transformation()->RequiresCachingOfDisplacements()) {

    // Add/Amend displacement cache entry
    Array<DisplacementInfo>::iterator i;
    for (i = _DisplacementInfo.begin(); i != _DisplacementInfo.end(); ++i) {
      if (i->_Transformation == output->Transformation() &&
          ContainsPoints(i->_Domain, _PointSetInput[n]->GetPoints()) &&
           fequal(i->_Domain._torigin, output->Time(),      1e-9) &&
          (fequal(i->_InputTime,       output->InputTime(), 1e-9) || IsNaN(i->_InputTime))) {
        break;
      }
    }
    if (i == _DisplacementInfo.end()) {
      Vector3D<double> res = this->AverageOutputResolution();
      if (ContainsPoints(_RegistrationDomain, _PointSetInput[n]->GetPoints())) {
        struct ImageAttributes domain = _RegistrationDomain;
        domain._x  = (res._x ? int(floor(_RegistrationDomain._x * _RegistrationDomain._dx / res._x)) + 1 : 1);
        domain._y  = (res._y ? int(floor(_RegistrationDomain._y * _RegistrationDomain._dy / res._y)) + 1 : 1);
        domain._z  = (res._z ? int(floor(_RegistrationDomain._z * _RegistrationDomain._dz / res._z)) + 1 : 1);
        domain._dx = (res._x ? (_RegistrationDomain._x * _RegistrationDomain._dx / domain._x) : .0);
        domain._dy = (res._y ? (_RegistrationDomain._y * _RegistrationDomain._dy / domain._y) : .0);
        domain._dz = (res._z ? (_RegistrationDomain._z * _RegistrationDomain._dz / domain._z) : .0);
        output->Domain(domain);
      } else {
        output->Domain(PointSetDomain(_PointSetInput[n], res));
        if (output->Domain()._x > 1) output->Domain()._x += 8;
        if (output->Domain()._y > 1) output->Domain()._y += 8;
        if (output->Domain()._z > 1) output->Domain()._z += 8;
      }
      if (output->Domain()) {
        output->Domain()._torigin = output->Time();
        _DisplacementInfo.resize(_DisplacementInfo.size() + 1);
        i = _DisplacementInfo.end() - 1;
        i->_InputTime      = output->InputTime();
        i->_DispIndex      = static_cast<int>(_DisplacementField.size());
        i->_Domain         = output->Domain();
        i->_Transformation = output->Transformation();
        _DisplacementField.push_back(new DisplacementImageType(output->Domain(), 3));
      }
    } else {
      output->Domain(i->_Domain);
    }
    if (output->Domain()) {
      output->ExternalDisplacement(_DisplacementField[i->_DispIndex]);
    }
  }

  // Initialize output
  output->Initialize();

  // Keep pointer to newly instantiated auxiliary object which will be deleted
  // again during destruction of the registration filter. The output objects
  // are further updated by PreUpdateCallback when requested by the optimizer.
  PointSetOutputInfo info;
  info._InputIndex     = n;
  info._Transformation = ti;
  info._InitialUpdate  = true;
  _PointSetOutputInfo.push_back(info);
  _PointSetOutput    .push_back(output);
  return output;
}

#endif // MIRTK_Registration_WITH_PointSet

// -----------------------------------------------------------------------------
void GenericRegistrationFilter::AddImageSimilarityTerm()
{
  ImageSimilarity *similarity;
  struct ImageAttributes  attr;
  Array<ImageSimilarityInfo>::const_iterator sim;
  for (sim = _ImageSimilarityInfo.begin(); sim != _ImageSimilarityInfo.end(); ++sim) {
    // Determine common attributes of output images
    if (_Mask[_CurrentLevel]) {
      attr = _Mask[_CurrentLevel]->Attributes();
    } else {
      Array<struct ImageAttributes> attrs;
      if (sim->IsSymmetric() || sim->_TargetTransformation.IsIdentity()) {
        attrs.push_back(OrthogonalFieldOfView(this->ImageAttributes(sim->_TargetIndex)));
      }
      if (sim->IsSymmetric() || sim->_SourceTransformation.IsIdentity()) {
        attrs.push_back(OrthogonalFieldOfView(this->ImageAttributes(sim->_SourceIndex)));
      }
      if (attrs.empty()) attr = this->RegistrationDomain();
      else               attr = OverallFieldOfView(attrs);
    }
    // Instantiate new similarity measure term
    similarity = ImageSimilarity::New(sim->_Measure);
    // Use sign of default weight if similarity measure was not explicitly
    // named in the energy formula in which case the weight for the similarity
    // term is such that SIM in the energy formula is to be minimized.
    // The default constructor of the actual similarity measure initializes
    // its weight such that the sign of the resulting similarity measure
    // corresponds to a minimization problem as well.
    if (sim->_DefaultSign) similarity->Weight(copysign(sim->_Weight, similarity->Weight()));
    else                   similarity->Weight(sim->_Weight);
    similarity->Name(sim->_Name);
    // Domain on which to evaluate similarity
    similarity->Mask  (_Mask[_CurrentLevel]);
    similarity->Domain(attr);
    // Set input of registered images
    this->SetInputOf(similarity->Target(), attr, sim->_TargetIndex, sim->_TargetTransformation);
    this->SetInputOf(similarity->Source(), attr, sim->_SourceIndex, sim->_SourceTransformation);
    // Add similarity term to energy function
    _Energy.Add(similarity);
  }
}

// -----------------------------------------------------------------------------
void GenericRegistrationFilter::AddPointSetDistanceTerm()
{
  #if MIRTK_Registration_WITH_PointSet
    PointSetDistance *dist = NULL;
    Array<PointSetDistanceInfo>::const_iterator d;
    for (d = _PointSetDistanceInfo.begin(); d != _PointSetDistanceInfo.end(); ++d) {
      // Instantiate new point set distance term
      dist = PointSetDistance::New(d->_Measure);
      // Use sign of default weight if distance measure was not explicitly
      // named in the energy formula in which case the weight for the distance
      // term is such that PDM in the energy formula is to be minimized.
      // The default constructor of the actual distance measure initializes
      // its weight such that the sign of the resulting distance measure
      // corresponds to a minimization problem as well.
      if (d->_DefaultSign) dist->Weight(copysign(d->_Weight, dist->Weight()));
      else                 dist->Weight(d->_Weight);
      dist->Name(d->_Name);
      // Set registered polydata objects whose distance is to be minimized
      double tt = _PointSetTime[d->_TargetTransformation ? d->_SourceIndex : d->_TargetIndex];
      double ts = _PointSetTime[d->_SourceTransformation ? d->_TargetIndex : d->_SourceIndex];
      dist->Target(this->OutputPointSet(d->_TargetIndex, tt, d->_TargetTransformation));
      dist->Source(this->OutputPointSet(d->_SourceIndex, ts, d->_SourceTransformation));
      // Add fiducial registration error term to energy function
      _Energy.Add(dist);
    }
  #else // MIRTK_Registration_WITH_PointSet
    if (!_PointSetDistanceInfo.empty()) {
      cerr << "GenericRegistrationFilter: Cannot register point sets when MIRTK was built without PointSet module!" << endl;
      exit(1);
    }
  #endif // MIRTK_Registration_WITH_PointSet
}

// -----------------------------------------------------------------------------
#if MIRTK_Registration_WITH_PointSet && MIRTK_Registration_WITH_Deformable
inline void AddNewPointSetConstraintTerm(RegistrationEnergy &energy,
                                         const GenericRegistrationFilter::PointSetConstraintInfo &info,
                                         RegisteredPointSet *output, int i = 0)
{
  string name = RegistrationEnergyParser::Substitute(info._Name, "{i}", ++i);
  // Instantiate new point set constraint term
  InternalForce *constraint;
  constraint = InternalForce::New(static_cast<InternalForceTerm>(info._Measure));
  constraint->PointSet(output);
  constraint->Weight  (info._Weight);
  constraint->Name    (name);
  // Add constraint term to energy function
  energy.Add(constraint);
}
#endif // MIRTK_Registration_WITH_PointSet && MIRTK_Registration_WITH_Deformable

// -----------------------------------------------------------------------------
void GenericRegistrationFilter::AddPointSetConstraintTerm()
{
  #if MIRTK_Registration_WITH_PointSet && MIRTK_Registration_WITH_Deformable
    RegisteredPointSet *output;
    Array<PointSetConstraintInfo>::const_iterator cst;
    for (cst = _PointSetConstraintInfo.begin(); cst != _PointSetConstraintInfo.end(); ++cst) {
      if (cst->_RefPointSetIndex > -1) {
        const int    &n = cst->_PointSetIndex;
        const double t  = _PointSetTime[cst->_RefPointSetIndex];
        output = this->OutputPointSet(n, t, cst->_Transformation);
        AddNewPointSetConstraintTerm(_Energy, *cst, output);
      } else if (cst->_RefImageIndex > -1) {
        const int    &n = cst->_PointSetIndex;
        const double t  = _Input[cst->_RefImageIndex]->GetTOrigin();
        output = this->OutputPointSet(n, t, cst->_Transformation);
        AddNewPointSetConstraintTerm(_Energy, *cst, output);
      } else {
        int i = 0;
        for (size_t n = 0; n < _PointSetOutput.size(); ++n) {
          if (_PointSetOutputInfo[n]._InputIndex     == cst->_PointSetIndex &&
              _PointSetOutputInfo[n]._Transformation == cst->_Transformation) {
            AddNewPointSetConstraintTerm(_Energy, *cst, _PointSetOutput[n], ++i);
          }
        }
        if (i == 0) {
          const int    &n = cst->_PointSetIndex;
          const double t0 = _PointSetTime[n];
          output = this->OutputPointSet(n, t0, cst->_Transformation);
          AddNewPointSetConstraintTerm(_Energy, *cst, output, ++i);
        }
      }
    }
  #else // MIRTK_Registration_WITH_PointSet && MIRTK_Registration_WITH_Deformable
    if (!_PointSetConstraintInfo.empty()) {
      cerr << "GenericRegistrationFilter: Cannot constrain point set when MIRTK was"
              " built without PointSet or Deformable module!" << endl;
      exit(1);
    }
  #endif // MIRTK_Registration_WITH_PointSet && MIRTK_Registration_WITH_Deformable
}

// -----------------------------------------------------------------------------
void GenericRegistrationFilter::AddPenaltyTerm()
{
  // Determine common attributes of (downsampled) target images
  struct ImageAttributes attr;
  if (_Mask[_CurrentLevel]) {
    attr = _Mask[_CurrentLevel]->Attributes();
  } else {
    Array<struct ImageAttributes> attrs;
    Array<ImageSimilarityInfo>::const_iterator sim;
    for (sim = _ImageSimilarityInfo.begin(); sim != _ImageSimilarityInfo.end(); ++sim) {
      if (sim->IsSymmetric() || sim->_TargetTransformation.IsIdentity()) {
        attrs.push_back(OrthogonalFieldOfView(this->ImageAttributes(sim->_TargetIndex)));
      }
      if (sim->IsSymmetric() || sim->_SourceTransformation.IsIdentity()) {
        attrs.push_back(OrthogonalFieldOfView(this->ImageAttributes(sim->_SourceIndex)));
      }
    }
    if (attrs.empty()) attr = this->RegistrationDomain();
    else               attr = OverallFieldOfView(attrs);
  }
  // Add constraint terms
  TransformationConstraint *constraint;
  Array<ConstraintInfo>::const_iterator cst;
  for (cst = _ConstraintInfo.begin(); cst != _ConstraintInfo.end(); ++cst) {
    // Instantiate new constraint measure term
    constraint = TransformationConstraint::New(cst->_Measure);
    constraint->Weight(cst->_Weight);
    constraint->Name  (cst->_Name);
    constraint->Domain(attr);
    // Add constraint term to energy function
    _Energy.Add(constraint);
  }
  if (_TargetTransformationErrorWeight > 0.) {
    UniquePtr<MeanSquaredDisplacementError> msde;
    msde.reset(new MeanSquaredDisplacementError());
    msde->Name  (_TargetTransformationErrorName);
    msde->Weight(_TargetTransformationErrorWeight);
    msde->Domain(attr);
    msde->TargetTransformation(_TargetTransformation);
    _Energy.Add(msde.release());
  }
}

// -----------------------------------------------------------------------------
void GenericRegistrationFilter::PreUpdateCallback(bool gradient)
{
  // Update cached displacements
  if (_Transformation->Changed() || gradient) {
    MIRTK_START_TIMING();
    Array<DisplacementInfo>::const_iterator i;
    for (i = _DisplacementInfo.begin(); i != _DisplacementInfo.end(); ++i) {
      if (IsNaN(i->_InputTime)) {
        i->_Transformation->Displacement(*_DisplacementField[i->_DispIndex]);
      } else {
        i->_Transformation->Displacement(*_DisplacementField[i->_DispIndex], i->_InputTime);
      }
    }
    MIRTK_DEBUG_TIMING(2, "caching of displacements");
  }

  // Update transformed point sets
#if MIRTK_Registration_WITH_PointSet
  if (_Transformation->Changed() || gradient) {
    MIRTK_START_TIMING();
    bool reinit_pointset_terms = false;
    Array<bool> remeshed(_PointSetOutput.size(), false);
    for (size_t i = 0; i < _PointSetOutput.size(); ++i) {
      if (_PointSetOutputInfo[i]._InitialUpdate) {
        _PointSetOutput[i]->Update(true);
        _PointSetOutputInfo[i]._InitialUpdate = false;
      } else if (_PointSetOutput[i]->Transformation()) {
        if (gradient) {
          // Remesh deformed surface after each gradient step
          //
          // The original (downsampled) input surface is remeshed instead of the
          // deformed surface mesh such that for each point in the remeshed
          // surface we know the untransformed coordinates needed for gradient
          // computation in Transformation::ParametricGradient.
          const int    j    = _PointSetOutputInfo[i]._InputIndex;
          const double dmin = _MinEdgeLength[_CurrentLevel][j];
          const double dmax = _MaxEdgeLength[_CurrentLevel][j];
          if (_AdaptiveRemeshing && IsSurfaceMesh(_PointSetInput[j]) && (dmin > .0 || !IsInf(dmax))) {
            MIRTK_START_TIMING();
            SurfaceRemeshing remesher;
            remesher.Input(vtkPolyData::SafeDownCast(_PointSetOutput[i]->InputPointSet()));
            remesher.Transformation(_PointSetOutput[i]->Transformation());
            remesher.MeltingOrder(SurfaceRemeshing::AREA);
            remesher.MeltNodesOff();
            remesher.MeltTrianglesOn();
            remesher.MinEdgeLength(dmin);
            remesher.MaxEdgeLength(dmax);
            remesher.Run();
            _PointSetOutput[i]->InputPointSet(remesher.Output());
            _PointSetOutput[i]->Initialize();
            MIRTK_DEBUG_TIMING(7, "remeshing moving surface");
            remeshed[i] = true;
            reinit_pointset_terms = true;
          }
        }
        _PointSetOutput[i]->Update(true);
        _PointSetOutputInfo[i]._InitialUpdate = false;
      }
    }
    if (reinit_pointset_terms) {
      for (int i = 0; i < _Energy.NumberOfTerms(); ++i) {
        EnergyTerm *term = _Energy.Term(i);
        PointSetDistance *pdm = dynamic_cast<PointSetDistance *>(term);
        if (pdm) {
          for (size_t j = 0; j < _PointSetOutput.size(); ++j) {
            if ((pdm->Target() == _PointSetOutput[j] ||
                 pdm->Source() == _PointSetOutput[j]) && remeshed[j]) {
              pdm->Reinitialize();
              break;
            }
          }
        }
        #if MIRTK_Registration_WITH_Deformable
          InternalForce *pcm = dynamic_cast<InternalForce *>(term);
          if (pcm) {
            for (size_t j = 0; j < _PointSetOutput.size(); ++j) {
              if (pcm->PointSet() == _PointSetOutput[j] && remeshed[j]) {
                pcm->Reinitialize();
                break;
              }
            }
          }
        #endif // MIRTK_Registration_WITH_Deformable
      }
    }
    MIRTK_DEBUG_TIMING(2, "update of moving point sets");
  }
#endif // MIRTK_Registration_WITH_PointSet
}

// -----------------------------------------------------------------------------
void GenericRegistrationFilter::InitializeEnergy()
{
  MIRTK_START_TIMING();

  // Construct registration energy function
  this->AddImageSimilarityTerm();
  this->AddPointSetDistanceTerm();
  this->AddPointSetConstraintTerm();
  this->AddPenaltyTerm();

  // Forward energy function events
  _Energy.AddObserver(_EventDelegate);

  // Set pre-update function for update of cached displacements
  if (!_PointSetOutput.empty() || !_DisplacementInfo.empty()) {
    _Energy.PreUpdateFunction(_PreUpdateDelegate);
  }

  // Set current output transformation
  _Energy.Transformation(_Transformation);

  // Set parameters of energy function such as penalty weights
  _Energy.Parameter(_Parameter[0]);
  _Energy.Parameter(_Parameter[_CurrentLevel]);

  // Normalize weights of data fidelity terms
  if (_NormalizeWeights) {
    double W = .0;
    for (int i = 0; i < _Energy.NumberOfTerms(); ++i) {
      if (_Energy.IsDataTerm(i)) {
        W += abs(_Energy.Term(i)->Weight());
      }
    }
    for (int i = 0; i < _Energy.NumberOfTerms(); ++i) {
      if (_Energy.IsDataTerm(i)) {
        _Energy.Term(i)->Weight(_Energy.Term(i)->Weight() / W);
      }
    }
  }

  // Initialize energy terms
  _Energy.Initialize();

  // Memorize (all) used (non-DoF) energy settings for -parout file
  Insert(_Parameter[_CurrentLevel], _Energy.Parameter());

  MIRTK_DEBUG_TIMING(3, "initialization of energy function");
}

// -----------------------------------------------------------------------------
void GenericRegistrationFilter::InitializeOptimizer()
{
  MIRTK_START_TIMING();

  if (AtInitialLevel()) {
    // Instantiate optimizer
    _Optimizer = LocalOptimizer::New(_OptimizationMethod, &_Energy);
    // Enable forwarding of optimization events
    _Optimizer->AddObserver(_EventDelegate);
  }

  // Set optimization parameters (global first)
  _Optimizer->Parameter(_Parameter[0]);
  _Optimizer->Parameter(_Parameter[_CurrentLevel]);

  // Initialize optimizer
  _Optimizer->Initialize();

  // Memorize (all) used (non-DoF) optimizer settings for -parout file
  Insert(_Parameter[_CurrentLevel], _Optimizer->Parameter());

  MIRTK_DEBUG_TIMING(3, "initialization of optimizer");
}

// -----------------------------------------------------------------------------
void GenericRegistrationFilter::Finalize()
{
  MIRTK_START_TIMING();

  // Destruct energy function and related auxiliary objects
  _Energy.Clear();
  #if MIRTK_Registration_WITH_PointSet
    for (size_t i = 0; i < _PointSetOutput.size(); ++i) {
      Delete(_PointSetOutput[i]);
    }
  #endif // MIRTK_Registration_WITH_PointSet
  _PointSetOutput.clear();
  _PointSetOutputInfo.clear();
  for (size_t i = 0; i < _TransformationInstance.size(); ++i) {
    Delete(_TransformationInstance[i]);
  }
  _TransformationInstance.clear();
  _TransformationInfo    .clear();
  for (size_t i = 0; i < _DisplacementField.size(); ++i) {
    Delete(_DisplacementField[i]);
  }
  _DisplacementField.clear();
  _DisplacementInfo .clear();

  // Include centering transformations in final linear transformation
  HomogeneousTransformation *lin = NULL;
  if ((lin = dynamic_cast<HomogeneousTransformation *>(_Transformation))) {
    const Matrix mat = lin->GetMatrix();
    Matrix pre(4, 4);
    pre.Ident();
    pre(0, 3)  = - _TargetOffset._x;
    pre(1, 3)  = - _TargetOffset._y;
    pre(2, 3)  = - _TargetOffset._z;
    Matrix post(4, 4);
    post.Ident();
    post(0, 3) = + _SourceOffset._x;
    post(1, 3) = + _SourceOffset._y;
    post(2, 3) = + _SourceOffset._z;
    lin->PutMatrix(post * mat * pre);
  }

  // Restore origin of images in resolution pyramid
  // (as a non-rigid registration might follow this initial linear alignment)
  for (int n = 0; n <  NumberOfImages(); ++n) {
    Point origin = _Image[_CurrentLevel][n].GetOrigin();
    if (IsTargetImage(n)) {
      _Image[_CurrentLevel][n].PutOrigin(origin + _TargetOffset);
    } else {
      _Image[_CurrentLevel][n].PutOrigin(origin + _SourceOffset);
    }
  }
  if (_Mask[_CurrentLevel]) {
    Point origin = _Mask[_CurrentLevel]->GetOrigin();
    _Mask[_CurrentLevel]->PutOrigin(origin + _TargetOffset);
  }

  // Restore positions of point sets
  #if MIRTK_Registration_WITH_PointSet
    for (int n = 0; n < NumberOfPointSets(); ++n) {
      double p[3], offset[3];
      if (IsTargetPointSet(n)) {
        offset[0] = _TargetOffset._x;
        offset[1] = _TargetOffset._y;
        offset[2] = _TargetOffset._z;
      } else {
        offset[0] = _SourceOffset._x;
        offset[1] = _SourceOffset._y;
        offset[2] = _SourceOffset._z;
      }
      vtkPoints *points = _PointSet[_CurrentLevel][n]->GetPoints();
      for (vtkIdType ptId = 0; ptId < points->GetNumberOfPoints(); ++ptId) {
        points->GetPoint(ptId, p);
        p[0] += offset[0];
        p[1] += offset[1];
        p[2] += offset[2];
        points->SetPoint(ptId, p);
      }
    }
  #endif // MIRTK_Registration_WITH_PointSet

  // Reset stored offsets
  _TargetOffset = _SourceOffset = Point(.0, .0, .0);

  // Destroy optimizer
  if (AtFinalLevel()) Delete(_Optimizer);

  // Update output transformation
  Output(_Transformation);

  MIRTK_DEBUG_TIMING(2, "finalization of level " << _CurrentLevel);
}


} // namespace mirtk
