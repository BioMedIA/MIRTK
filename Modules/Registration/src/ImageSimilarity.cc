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

#include "mirtk/ImageSimilarity.h"

#include "mirtk/Assert.h"
#include "mirtk/Math.h"
#include "mirtk/Memory.h"
#include "mirtk/Parallel.h"
#include "mirtk/Profiling.h"
#include "mirtk/VoxelFunction.h"
#include "mirtk/MultiLevelTransformation.h"


namespace mirtk {


// =============================================================================
// Factory
// =============================================================================

// -----------------------------------------------------------------------------
ImageSimilarity *ImageSimilarity::New(SimilarityMeasure sim, const char *name, double w)
{
  enum EnergyMeasure em = static_cast<enum EnergyMeasure>(sim);
  if (SIM_Begin < em && em < SIM_End) {
    EnergyTerm *term = EnergyTerm::TryNew(em, name, w);
    if (term) return dynamic_cast<ImageSimilarity *>(term);
    cerr << NameOfType() << "::New: Image similarity measure not available: ";
  } else {
    cerr << NameOfType() << "::New: Energy term is not an image similarity measure: ";
  }
  cerr << ToString(em) << " (" << em << ")" << endl;
  exit(1);
  return NULL;
}

// =============================================================================
// Auxiliary functor classes for parallel execution
// =============================================================================

namespace ImageSimilarityUtils {


// -----------------------------------------------------------------------------
/// Voxel function used by MultiplyByImageGradient to post-multiply similarity
/// gradient by transformed image gradient according to the chain rule.
class MultiplySimilarityGradientByImageGradient : public VoxelFunction
{
private:

  double _Scale;
  int _NumberOfVoxels;
  int _dx, _dy, _dz;

public:

  MultiplySimilarityGradientByImageGradient(const RegisteredImage *image, double scale = 1.)
  :
    _Scale(scale), _NumberOfVoxels(image->NumberOfVoxels()),
    _dx            (image->Offset(RegisteredImage::Dx)),
    _dy            (image->Offset(RegisteredImage::Dy)),
    _dz            (image->Offset(RegisteredImage::Dz))
  {}

  template <class Image, class TScalar, class TReal>
  void operator ()(const Image &, int, const TScalar *dI, TReal *gradient)
  {
    (*gradient) *= static_cast<TReal>(_Scale * dI[_dx]); gradient += _NumberOfVoxels;
    (*gradient) *= static_cast<TReal>(_Scale * dI[_dy]); gradient += _NumberOfVoxels;
    (*gradient) *= static_cast<TReal>(_Scale * dI[_dz]);
  }
};

// -----------------------------------------------------------------------------
/// Determine maximum norm of voxel-wise image similarity gradient
class MaxVoxelWiseSimilarityGradient : public VoxelReduction
{
private:

  double    _norm;
  const int _x, _y, _z;

public:

  /// Constructor
  MaxVoxelWiseSimilarityGradient(const BaseImage *img)
  :
    _norm(.0), _x(0), _y(img->X() * img->Y() * img->Z()), _z(2 * _y)
  {}

  /// Split "constructor"
  void split(const MaxVoxelWiseSimilarityGradient &other)
  {
    _norm = other._norm;
  }

  /// Join results
  void join(const MaxVoxelWiseSimilarityGradient &other)
  {
    if (other._norm > _norm) _norm = other._norm;
  }

  /// Get maximum norm
  double Norm() const { return sqrt(_norm); }

  /// Compute norm of similarity gradient
  template <class Image, class TReal>
  void operator ()(const Image &, int, TReal *gradient)
  {
    double gx = static_cast<double>(gradient[_x]);
    double gy = static_cast<double>(gradient[_y]);
    double gz = static_cast<double>(gradient[_z]);
    double norm = gx * gx + gy * gy + gz * gz;
    if (norm > _norm) _norm = norm;
  }
};

// -----------------------------------------------------------------------------
/// Determine maximum norm of node-based image similarity gradient
class MaxNodeBasedSimilarityGradient
{
private:

  const FreeFormTransformation *_FFD;
  const double                 *_Gradient;
  double                        _MaxNorm;

public:

  /// Constructor
  MaxNodeBasedSimilarityGradient(const FreeFormTransformation *ffd,
                                 const double                 *gradient)
  :
    _FFD(ffd), _Gradient(gradient), _MaxNorm(.0)
  {}

  /// Copy constructor
  MaxNodeBasedSimilarityGradient(const MaxNodeBasedSimilarityGradient &other)
  :
    _FFD     (other._FFD),
    _Gradient(other._Gradient),
    _MaxNorm (other._MaxNorm)
  {}

  /// Split constructor
  MaxNodeBasedSimilarityGradient(const MaxNodeBasedSimilarityGradient &other, split)
  :
    _FFD     (other._FFD),
    _Gradient(other._Gradient),
    _MaxNorm (other._MaxNorm)
  {}

  /// Join results
  void join(const MaxNodeBasedSimilarityGradient &other)
  {
    if (other._MaxNorm > _MaxNorm) _MaxNorm = other._MaxNorm;
  }

  /// Maximum norm
  double Norm() const { return sqrt(_MaxNorm); }

  /// Determine maximum norm of specified control point gradients
  void operator()(const blocked_range<int> &re)
  {
    double norm;
    int    x, y, z;

    for (int cp = re.begin(); cp != re.end(); ++cp) {
      _FFD->IndexToDOFs(cp, x, y, z);
      norm = pow(_Gradient[x], 2) + pow(_Gradient[y], 2) + pow(_Gradient[z], 2);
      if (norm > _MaxNorm) _MaxNorm = norm;
    }
  }
};

// -----------------------------------------------------------------------------
/// Normalize voxel-wise image similarity gradient
class NormalizeVoxelWiseSimilarityGradient : public VoxelFunction
{
private:

  double    _Sigma;
  const int _x, _y, _z;

public:

  /// Constructor
  NormalizeVoxelWiseSimilarityGradient(const BaseImage *img, double sigma = .0)
  :
    _Sigma(sigma), _x(0), _y(img->X() * img->Y() * img->Z()), _z(2 * _y)
  {}

  /// Normalize similarity gradient
  template <class Image, class TReal>
  void operator ()(const Image &, int, TReal *gradient)
  {
    double gx = static_cast<double>(gradient[_x]);
    double gy = static_cast<double>(gradient[_y]);
    double gz = static_cast<double>(gradient[_z]);
    double norm = gx * gx + gy * gy + gz * gz;
    if (norm) {
      norm = sqrt(norm) + _Sigma;
      gradient[_x] = static_cast<TReal>(gx / norm);
      gradient[_y] = static_cast<TReal>(gy / norm);
      gradient[_z] = static_cast<TReal>(gz / norm);
    }
  }
};

// -----------------------------------------------------------------------------
/// Normalize node-based image similarity gradient
class NormalizeNodeBasedSimilarityGradient
{
private:

  const FreeFormTransformation *_FFD;
  double                       *_Gradient;
  double                        _Sigma;

public:

  /// Constructor
  NormalizeNodeBasedSimilarityGradient(const FreeFormTransformation *ffd,
                                       double                       *gradient,
                                       double                        sigma)
  :
    _FFD(ffd), _Gradient(gradient), _Sigma(sigma)
  {}

  /// Copy constructor
  NormalizeNodeBasedSimilarityGradient(const NormalizeNodeBasedSimilarityGradient &other)
  :
    _FFD     (other._FFD),
    _Gradient(other._Gradient),
    _Sigma   (other._Sigma)
  {}

  /// Normalize node-based similarity gradient
  void operator ()(const blocked_range<int> &re) const
  {
    double norm;
    int    x, y, z;

    for (int cp = re.begin(); cp != re.end(); ++cp) {
      _FFD->IndexToDOFs(cp, x, y, z);
      norm = pow(_Gradient[x], 2) + pow(_Gradient[y], 2) + pow(_Gradient[z], 2);
      if (norm) {
        norm = sqrt(norm) + _Sigma;
        _Gradient[x] /= norm;
        _Gradient[y] /= norm;
        _Gradient[z] /= norm;
      }
    }
  }
};


} // namespace ImageSimilarityUtils
using namespace ImageSimilarityUtils;

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
ImageSimilarity::ImageSimilarity(const char *name, double weight)
:
  DataFidelity(name, weight),
  _Target(new RegisteredImage()), _TargetOwner(true),
  _Source(new RegisteredImage()), _SourceOwner(true),
  _Foreground              (FG_Overlap),
  _Mask                    (nullptr),
  _GradientWrtTarget       (nullptr),
  _GradientWrtSource       (nullptr),
  _Gradient                (nullptr),
  _NumberOfVoxels          (0),
  _NormalizeImageGradient  (true),
  _UseApproximateGradient  (false),
  _VoxelWisePreconditioning(0.),
  _NodeBasedPreconditioning(0.),
  _SkipTargetInitialization(false),
  _SkipSourceInitialization(false),
  _InitialUpdate           (false)
{
  _ParameterPrefix.push_back("Image (dis-)similarity ");
  _ParameterPrefix.push_back("Image dissimilarity ");
  _ParameterPrefix.push_back("Image similarity ");
  _ParameterPrefix.push_back("(Dis-)similarity ");
  _ParameterPrefix.push_back("Dissimilarity ");
  _ParameterPrefix.push_back("Similarity ");
}

// -----------------------------------------------------------------------------
void ImageSimilarity::CopyAttributes(const ImageSimilarity &other)
{
  if (_TargetOwner) delete _Target;
  _Target      = (other._TargetOwner ? new RegisteredImage(*other._Target) : other._Target);
  _TargetOwner = other._TargetOwner;

  if (_SourceOwner) delete _Source;
  _Source      = (other._SourceOwner ? new RegisteredImage(*other._Source) : other._Source);
  _SourceOwner = other._SourceOwner;

  Delete(_GradientWrtTarget);
  Delete(_GradientWrtSource);
  Deallocate(_Gradient);

  _Domain                   = other._Domain;
  _Foreground               = other._Foreground;
  _Mask                     = other._Mask;
  _NumberOfVoxels           = other._NumberOfVoxels;
  _NormalizeImageGradient   = other._NormalizeImageGradient;
  _UseApproximateGradient   = other._UseApproximateGradient;
  _VoxelWisePreconditioning = other._VoxelWisePreconditioning;
  _NodeBasedPreconditioning = other._NodeBasedPreconditioning;
  _SkipTargetInitialization = other._SkipTargetInitialization;
  _SkipSourceInitialization = other._SkipSourceInitialization;
  _InitialUpdate            = other._InitialUpdate;
}

// -----------------------------------------------------------------------------
ImageSimilarity::ImageSimilarity(const ImageSimilarity &other)
:
  DataFidelity(other),
  _Target           (nullptr),
  _Source           (nullptr),
  _GradientWrtTarget(nullptr),
  _GradientWrtSource(nullptr),
  _Gradient         (nullptr)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
ImageSimilarity &ImageSimilarity::operator =(const ImageSimilarity &other)
{
  if (this != &other) {
    DataFidelity::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
ImageSimilarity::~ImageSimilarity()
{
  if (_TargetOwner) Delete(_Target);
  if (_SourceOwner) Delete(_Source);
  Delete(_GradientWrtTarget);
  Delete(_GradientWrtSource);
  Deallocate(_Gradient);
}

// =============================================================================
// Initialization
// =============================================================================

// -----------------------------------------------------------------------------
void ImageSimilarity::InitializeInput(const ImageAttributes &domain)
{
  if (!_SkipTargetInitialization) {
    _Target->Initialize(domain, _Target->Transformation() ? 4 : 1);
  }
  if (!_SkipSourceInitialization) {
    _Source->Initialize(domain, _Source->Transformation() ? 4 : 1);
  }
  _Domain = domain;
  _Domain._t = 1;
}

// -----------------------------------------------------------------------------
void ImageSimilarity::Initialize()
{
  // Free previously allocated memory
  Delete(_GradientWrtTarget);
  Delete(_GradientWrtSource);
  Deallocate(_Gradient);
  // Check if all inputs are set
  mirtkAssert(_Target != nullptr, "target image component cannot be nullptr");
  mirtkAssert(_Source != nullptr, "source image component cannot be nullptr");
  if (( _Mask &&  _Mask->IsEmpty()) ||
      (!_Mask && (_Domain._x == 0 || _Domain._y == 0 || _Domain._z == 0))) {
    cerr << "ImageSimilarity::Initialize: No image domain specified" << endl;
    exit(1);
  }
  // Initialize base class
  DataFidelity::Initialize();
  // Initialize registered images
  this->InitializeInput(_Mask ? _Mask->Attributes() : _Domain);
  _InitialUpdate = true; // i.e., initialize image content upon first Update
  // Allocate memory for temporary similarity gradient
  if (_NodeBasedPreconditioning > .0) {
    const class Transformation *T1 = _Target->Transformation();
    const class Transformation *T2 = _Source->Transformation();
    const int n = max(T1 ? T1->NumberOfDOFs() : 0,
                           T2 ? T2->NumberOfDOFs() : 0);
    Allocate(_Gradient, n);
  }
  // Number of voxels per registered image
  _NumberOfVoxels = _Domain.NumberOfSpatialPoints();
}

// =============================================================================
// Parameters
// =============================================================================

// -----------------------------------------------------------------------------
bool ImageSimilarity::SetWithoutPrefix(const char *param, const char *value)
{
  if (strcmp(param, "Foreground") == 0 || strcmp(param, "Foreground region") == 0) {
    return FromString(value, _Foreground);
  }
  if (strcmp(param, "Normalize image gradient") == 0 ||
      strcmp(param, "Normalise image gradient") == 0) {
    return FromString(value, _NormalizeImageGradient);
  }
  if (strcmp(param, "Approximate gradient") == 0) {
    return FromString(value, _UseApproximateGradient);
  }
  if (strcmp(param, "Preconditioning (voxel-wise)") == 0) {
    return FromString(value, _VoxelWisePreconditioning);
  }
  if (strcmp(param, "Preconditioning (node-based)") == 0 ||
      strcmp(param, "Preconditioning") == 0) {
    return FromString(value, _NodeBasedPreconditioning);
  }
  if (strcmp(param, "Blurring of 1st order image derivatives") == 0 ||
      strcmp(param, "Blurring of image jacobian") == 0 ||
      strcmp(param, "Blurring of image gradient") == 0 ||
      strcmp(param, "Source gradient blurring") == 0 ||
      strcmp(param, "Source jacobian blurring") == 0) {
    double sigma;
    if (!FromString(value, sigma)) return false;
    _Target->GradientSigma(sigma);
    _Source->GradientSigma(sigma);
    return true;
  }
  if (strcmp(param, "Blurring of 2nd order image derivatives") == 0 ||
      strcmp(param, "Blurring of image hessian") == 0 ||
      strcmp(param, "Source hessian blurring") == 0) {
    double sigma;
    if (!FromString(value, sigma)) return false;
    _Target->HessianSigma(sigma);
    _Source->HessianSigma(sigma);
    return true;
  }
  if (strcmp(param, "Source gradient threshold") == 0) {
    int pct;
    string str;
    if (ValueUnits(value, &str, "%") != "%" || !FromString(str, pct)) return false;
    _Target->MaxGradientPercentile(pct);
    _Source->MaxGradientPercentile(pct);
    return true;
  }
  if (strcmp(param, "Source gradient maximum") == 0) {
    double norm;
    if (!FromString(value, norm)) return false;
    _Target->MaxGradientMagnitude(norm);
    _Source->MaxGradientMagnitude(norm);
    return true;
  }
  return DataFidelity::SetWithoutPrefix(param, value);
}

// -----------------------------------------------------------------------------
ParameterList ImageSimilarity::Parameter() const
{
  ParameterList params = DataFidelity::Parameter();
  InsertWithPrefix(params, "Foreground region",            _Foreground);
  InsertWithPrefix(params, "Normalize image gradient",     _NormalizeImageGradient);
  InsertWithPrefix(params, "Approximate gradient",         _UseApproximateGradient);
  InsertWithPrefix(params, "Preconditioning (voxel-wise)", _VoxelWisePreconditioning);
  InsertWithPrefix(params, "Preconditioning (node-based)", _NodeBasedPreconditioning);
  InsertWithPrefix(params, "Blurring of image gradient",   _Target->GradientSigma());
  InsertWithPrefix(params, "Blurring of image hessian",    _Target->HessianSigma());
  return params;
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
void ImageSimilarity::Update(bool gradient)
{
  if (_InitialUpdate || _Target->Transformation()) {
    _Target->Update(true, gradient, false, _InitialUpdate);
  }
  if (_InitialUpdate || _Source->Transformation()) {
    _Source->Update(true, gradient, false, _InitialUpdate);
  }
  _InitialUpdate = false;
}

// -----------------------------------------------------------------------------
void ImageSimilarity::Exclude(const blocked_range3d<int> &)
{
  // By default, call Update upon Include
}

// -----------------------------------------------------------------------------
void ImageSimilarity::Include(const blocked_range3d<int> &)
{
  this->Update(false);
}

// -----------------------------------------------------------------------------
void ImageSimilarity::MultiplyByImageGradient(const RegisteredImage *image,
                                              GradientImageType     *gradient)
{
  // Copy (dSimilarity / dI) from x component also to y and z components
  const int nbytes = image->NumberOfVoxels() * sizeof(GradientType);
  GradientType *gx = gradient->GetPointerToVoxels();
  memcpy(gradient->GetPointerToVoxels(0, 0, 0, 1), gx, nbytes);
  memcpy(gradient->GetPointerToVoxels(0, 0, 0, 2), gx, nbytes);
  // Scale transformed image gradient as if it were computed from an image
  // with an intensity range of [0, 10]. This allows penalty weights of
  // comparable magnitude to be used with both NMI and LNCC (and possibly
  // other similarity measures if scaled appropriately).
  double scale = 1.;
  if (_NormalizeImageGradient) {
    double min_intensity = image->MinIntensity();
    double max_intensity = image->MaxIntensity();
    if (IsNaN(min_intensity)) min_intensity = image->MinInputIntensity();
    if (IsNaN(max_intensity)) max_intensity = image->MaxInputIntensity();
    scale = 10. / (max_intensity - min_intensity);
  }
  // Apply chain rule (dSimilarity / dy) = (dSimilarity / dI) * (dI / dy)
  // where y = T(x) to obtain the non-parametric similarity gradient
  MultiplySimilarityGradientByImageGradient times_dIdx(image, scale);
  ParallelForEachVoxel(image, gradient, times_dIdx);
}

// -----------------------------------------------------------------------------
bool ImageSimilarity::NonParametricGradient(const RegisteredImage *,
                                            GradientImageType     *)
{
  return false; // By default, ApproximateGradient using finite differences
}

// -----------------------------------------------------------------------------
void ImageSimilarity::NormalizeGradient(GradientImageType *gradient)
{
  MIRTK_START_TIMING();

  // Determine maximum norm of control point gradients
  MaxVoxelWiseSimilarityGradient maximum(gradient);
  ParallelForEachVoxel(gradient, maximum);

  // Sigma value used to suppress noise
  const double sigma = _VoxelWisePreconditioning * maximum.Norm();

  // Normalize control point gradients to be possibly similar
  NormalizeVoxelWiseSimilarityGradient norm(gradient, sigma);
  ParallelForEachVoxel(gradient, norm);

  MIRTK_DEBUG_TIMING(2, "normalization of voxel-wise (dis-)similarity gradient");
}

// -----------------------------------------------------------------------------
void ImageSimilarity::ParametricGradient(const RegisteredImage *image,
                                         GradientImageType     *np_gradient,
                                         double                *gradient,
                                         double                 weight)
{
  const class Transformation *T = image->Transformation();
  mirtkAssert(T != NULL, "image is being transformed");
  const WorldCoordsImage *i2w = image->ImageToWorld();
  const double            t0  = image->GetTOrigin();
  np_gradient->PutTOrigin(image->InputImage()->GetTOrigin());
  T->ParametricGradient(np_gradient, gradient, i2w, t0, weight);
}

// -----------------------------------------------------------------------------
void ImageSimilarity::ApproximateGradient(RegisteredImage        *image,
                                          FreeFormTransformation *ffd,
                                          double *gradient, double step,
                                          double weight)
{
  weight /= 2.0 * step;
  double a, b, value;
  int i1, j1, k1, i2, j2, k2, dof[3];
  for (int cp = 0; cp < ffd->NumberOfCPs(); ++cp) {
    if (ffd->IsActive(cp) &&
        ffd->BoundingBox(image, cp, i1, j1, k1, i2, j2, k2)) {
      blocked_range3d<int> region(k1, k2+1, j1, j2+1, i1, i2+1);
      ffd->IndexToDOFs(cp, dof[0], dof[1], dof[2]);
      for (int i = 0; i < 3; ++i) {
        if (ffd->GetStatus(dof[i]) == Active) {
          value = ffd->Get(dof[i]);

          ffd->Put(dof[i], value + step);
          this->Exclude(region);
          image->Update(region);
          this->Include(region);
          a = this->Evaluate();

          ffd->Put(dof[i], value - step);
          this->Exclude(region);
          image->Update(region);
          this->Include(region);
          b = this->Evaluate();

          ffd->Put(dof[i], value);
          this->Exclude(region);
          image->Update(region);
          this->Include(region);

          gradient[dof[i]] += weight * (a - b);
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------
void ImageSimilarity::ApproximateGradient(RegisteredImage *image,
                                          double *gradient, double step,
                                          double weight)
{
  MIRTK_START_TIMING();

  class Transformation *T = const_cast<class Transformation *>(image->Transformation());
  mirtkAssert(T != NULL, "image is being transformed");
  FreeFormTransformation   *ffd  = NULL;
  MultiLevelTransformation *mffd = NULL;
  (mffd = dynamic_cast<MultiLevelTransformation *>(T)) ||
  (ffd  = dynamic_cast<FreeFormTransformation   *>(T));

  if (mffd) {
    for (int i = 0; i < mffd->NumberOfLevels(); ++i) {
      if (mffd->LocalTransformationIsActive(i)) {
        ffd = mffd->GetLocalTransformation(i);
        this->ApproximateGradient(image, ffd, gradient, step, weight);
        gradient += ffd->NumberOfDOFs();
      }
    }
  } else if (ffd) {
    this->ApproximateGradient(image, ffd, gradient, step, weight);
  } else {
    weight /= 2.0 * step;
    double a, b, value;
    blocked_range3d<int> region(0, image->Z(), 0, image->Y(), 0, image->X());
    for (int dof = 0; dof < T->NumberOfDOFs(); ++dof) {
      if (T->GetStatus(dof) == Active) {
        value = T->Get(dof);

        T->Put(dof, value + step);
        this->Exclude(region);
        image->Update(region);
        this->Include(region);
        a = this->Evaluate();

        T->Put(dof, value - step);
        this->Exclude(region);
        image->Update(region);
        this->Include(region);
        b = this->Evaluate();

        T->Put(dof, value);
        this->Exclude(region);
        image->Update(region);
        this->Include(region);

        gradient[dof] += weight * (a - b);
      }
    }
  }

  MIRTK_DEBUG_TIMING(2, "approximation of similarity gradient");
}

// -----------------------------------------------------------------------------
void ImageSimilarity::NormalizeGradient(const RegisteredImage *image, double *gradient)
{
  const MultiLevelTransformation *mffd = NULL;
  const FreeFormTransformation   *affd = NULL;

  (mffd = dynamic_cast<const MultiLevelTransformation *>(image->Transformation())) ||
  (affd = dynamic_cast<const FreeFormTransformation   *>(image->Transformation()));

  const int nlevels = (mffd ? mffd->NumberOfLevels() : (affd ? 1 : 0));
  if (nlevels == 0) return; // Skip if transformation is not a FFD

  MIRTK_START_TIMING();

  for (int lvl = 0; lvl < nlevels; ++lvl) {
    if (mffd) {
      if (!mffd->LocalTransformationIsActive(lvl)) continue;
      affd = mffd->GetLocalTransformation(lvl);
    }

    // Range of control point indices
    blocked_range<int> cps(0, affd->NumberOfCPs());

    // Determine maximum norm of control point gradients
    MaxNodeBasedSimilarityGradient maximum(affd, gradient);
    parallel_reduce(cps, maximum);

    // Sigma value used to suppress noise
    const double sigma = _NodeBasedPreconditioning * maximum.Norm();

    // Normalize control point gradients to be possibly similar
    NormalizeNodeBasedSimilarityGradient norm(affd, gradient, sigma);
    parallel_for(cps, norm);

    // Gradient w.r.t parameters of next active level
    gradient += affd->NumberOfDOFs();
  }

  MIRTK_DEBUG_TIMING(2, "normalization of (dis-)similarity gradient");
}

// -----------------------------------------------------------------------------
void ImageSimilarity::EvaluateGradient(RegisteredImage    *image,
                                       GradientImageType *&np_gradient,
                                       double             *gradient,
                                       double step, double weight)
{
  if (!image->Transformation()) return;
  const int ndofs = image->Transformation()->NumberOfDOFs();
  // Directly update gradient if no node-based preconditioning is used
  GradientType * const tmp_gradient = _Gradient;
  if (_NodeBasedPreconditioning <= .0) _Gradient = gradient;
  // Compute parametric gradient w.r.t. given transformed image
  if (_Gradient != gradient) {
    memset(_Gradient, 0, ndofs * sizeof(double));
  }
  // Compute analytic gradient if implemented by subclass
  if (!_UseApproximateGradient) {
    if (!np_gradient) {
      np_gradient = new GradientImageType(_Domain, 3);
    }
    if (this->NonParametricGradient(image, np_gradient)) {
      // Normalize voxel-based gradient
      if (_VoxelWisePreconditioning > .0) {
        this->NormalizeGradient(np_gradient);
      }
      this->ParametricGradient(image, np_gradient, _Gradient, weight);
    // ...otherwise, use finite differences approximation
    } else {
      _UseApproximateGradient = true;
    }
  }
  // If no analytic gradient computation implemented or approximation chosen,
  // approximate the image similarity measure gradient using finite differences
  if (_UseApproximateGradient) {
    Delete(np_gradient);
    this->ApproximateGradient(image, _Gradient, step, weight);
  }
  // Normalize node-based gradient
  if (_NodeBasedPreconditioning > .0) {
    this->NormalizeGradient(image, _Gradient);
    for (int dof = 0; dof < ndofs; ++dof) {
      // Note: Normalization cancels out weight factor!
      gradient[dof] += weight * _Gradient[dof];
    }
  }
  // Restore gradient pointer
  _Gradient = tmp_gradient;
}

// -----------------------------------------------------------------------------
void ImageSimilarity::EvaluateGradient(double *gradient, double step, double weight)
{
  // Compute parametric gradient w.r.t target transformation
  this->EvaluateGradient(_Target, _GradientWrtTarget, gradient, step, weight);
  // If target and source are transformed by different transformations,
  // the gradient vector contains first the derivative values w.r.t the
  // parameters of the target transformation followed by those computed
  // w.r.t the parameters of the source transformation. Otherwise, if
  // both images are transformed by the same transformation, i.e., a
  // velocity based transformation integrated half way in both directions,
  // the derivative values are summed up instead.
  if (_Target->Transformation() && _Source->Transformation()) {
    if (!HaveSameDOFs(_Target->Transformation(), _Source->Transformation())) {
      gradient += _Target->Transformation()->NumberOfDOFs();
    }
  }
  // Compute parametric gradient w.r.t source transformation
  this->EvaluateGradient(_Source, _GradientWrtSource, gradient, step, weight);
}

// =============================================================================
// Debugging
// =============================================================================

// -----------------------------------------------------------------------------
void ImageSimilarity::Print(Indent indent) const
{
  EnergyTerm::Print(indent);
  cout << "Target image:";
  if (_Target) {
    cout << endl;
    _Target->Print(indent + 1);
  } else {
    cout << " Missing!" << endl;
  }
  cout << "Source image:";
  if (_Source) {
    cout << endl;
    _Source->Print(indent + 1);
  } else {
    cout << " Missing!" << endl;
  }
  cout << "Image domain:" << endl;
  _Domain.Print(indent + 1);
  cout << "Overlap mask:";
  if (_Mask) {
    cout << endl;
    _Mask->Print(indent + 1);
  } else {
    cout << " None" << endl;
  }
}

// -----------------------------------------------------------------------------
void ImageSimilarity::WriteDataSets(const char *p, const char *suffix, bool all) const
{
  const int   sz = 1024;
  char        fname[sz];
  string _prefix = Prefix(p);
  const char  *prefix = _prefix.c_str();

  if (_Target->Transformation() || all) {
    snprintf(fname, sz, "%starget%s", prefix, suffix);
    _Target->GetFrame(0).Write(fname);
    if (_Target->T() > 4) {
      snprintf(fname, sz, "%starget_gradient%s", prefix, suffix);
      _Target->GetFrame(1, 3).Write(fname);
    }
  }
  if (_Source->Transformation() || all) {
    snprintf(fname, sz, "%ssource%s", prefix, suffix);
    _Source->GetFrame(0).Write(fname);
    if (_Source->T() > 4) {
      snprintf(fname, sz, "%ssource_gradient%s", prefix, suffix);
      _Source->GetFrame(1, 3).Write(fname);
    }
  }
}

// -----------------------------------------------------------------------------
void ImageSimilarity::WriteGradient(const char *p, const char *suffix) const
{
  const int   sz = 1024;
  char        fname[sz];
  string _prefix = Prefix(p);
  const char  *prefix = _prefix.c_str();

  // Image derivatives needed for gradient computation
  if (_Target->T() >= 4) {
    snprintf(fname, sz, "%starget_gradient%s", prefix, suffix);
    _Target->GetFrame(1, 3).Write(fname);
  }
  if (_Target->T() > 4) {
    snprintf(fname, sz, "%starget_hessian%s", prefix, suffix);
    _Target->GetFrame(4, _Target->T()-1).Write(fname);
  }

  if (_Source->T() >= 4) {
    snprintf(fname, sz, "%ssource_gradient%s", prefix, suffix);
    _Source->GetFrame(1, 3).Write(fname);
  }
  if (_Source->T() > 4) {
    snprintf(fname, sz, "%ssource_hessian%s", prefix, suffix);
    _Source->GetFrame(4, _Source->T()-1).Write(fname);
  }

  // Computed non-parametric similarity gradient
  if (_GradientWrtTarget) {
    snprintf(fname, sz, "%sgradient_wrt_target%s", prefix, suffix);
    _GradientWrtTarget->Write(fname);
  }
  if (_GradientWrtSource) {
    snprintf(fname, sz, "%sgradient_wrt_source%s", prefix, suffix);
    _GradientWrtSource->Write(fname);
  }
}


} // namespace mirtk
