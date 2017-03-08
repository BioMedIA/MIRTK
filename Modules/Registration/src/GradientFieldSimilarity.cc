/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2015 Imperial College London
 * Copyright 2013-2015 Stefan Pszczolkowski Parraguez, Andreas Schuh
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

#include "mirtk/GradientFieldSimilarity.h"

#include "mirtk/Memory.h"
#include "mirtk/Parallel.h"
#include "mirtk/GenericImage.h"
#include "mirtk/VoxelFunction.h"
#include "mirtk/Matrix.h"

#include "mirtk/FreeFormTransformation.h"
#include "mirtk/MultiLevelFreeFormTransformation.h"


namespace mirtk {


// =============================================================================
// Auxiliary functors
// =============================================================================

namespace GradientFieldSimilarityUtils {


/// Multiply similarity gradient by transformed Hessian of image
struct MultiplyGradientByHessian : public VoxelFunction
{
  static const int _x = 0; int _y, _z;
  int _dxx, _dxy, _dxz, _dyx, _dyy, _dyz, _dzx, _dzy, _dzz;

  MultiplyGradientByHessian(const RegisteredImage *image)
  :
    _y  (image->NumberOfVoxels()), _z(2 * _y),
    _dxx(image->Offset(RegisteredImage::Dxx)),
    _dxy(image->Offset(RegisteredImage::Dxy)),
    _dxz(image->Offset(RegisteredImage::Dxz)),
    _dyx(image->Offset(RegisteredImage::Dyx)),
    _dyy(image->Offset(RegisteredImage::Dyy)),
    _dyz(image->Offset(RegisteredImage::Dyz)),
    _dzx(image->Offset(RegisteredImage::Dzx)),
    _dzy(image->Offset(RegisteredImage::Dzy)),
    _dzz(image->Offset(RegisteredImage::Dzz))
  {}

  template <class Image, class TGradient, class TReal>
  void operator()(const Image &, int, const TGradient *dI, TReal *g)
  {
    if (g[_x] != .0 || g[_y] != .0 || g[_z] != .0) {
      TReal gx = g[_x] * dI[_dxx] + g[_y] * dI[_dyx] + g[_z] * dI[_dzx];
      TReal gy = g[_x] * dI[_dxy] + g[_y] * dI[_dyy] + g[_z] * dI[_dzy];
      TReal gz = g[_x] * dI[_dxz] + g[_y] * dI[_dyz] + g[_z] * dI[_dzz];
      g[_x] = gx, g[_y] = gy, g[_z] = gz;
    }
  }
};

// -----------------------------------------------------------------------------
class AddOtherPartialDerivatives
{
public:
  const Transformation                           *_Transformation;
  const ImageSimilarity::GradientImageType       *_NPGradient;
  const GenericImage<RegisteredImage::VoxelType> *_ImageGradient;
  const WorldCoordsImage                         *_Image2World;
  double                                         *_Output;
  double                                          _Weight;

  double _t;            ///< Time corrresponding to input gradient image (in ms)
  double _t0;           ///< Second time argument for velocity-based transformations

private:

  int    _X;            ///< Number of voxels along x axis
  int    _Y;            ///< Number of voxels along y axis
  int    _NumberOfDOFs; ///< Number of transformation parameters

public:

  // ---------------------------------------------------------------------------
  /// Default constructor
  AddOtherPartialDerivatives()
  :
    _Transformation(NULL),
    _NPGradient    (NULL),
    _ImageGradient (NULL),
    _Image2World   (NULL),
    _Output        (NULL),
    _Weight        ( 1.0),
    _t             ( 0.0),
    _t0            (-1.0),
    _X             (0),
    _Y             (0),
    _NumberOfDOFs  (0)
  {}

  // ---------------------------------------------------------------------------
  /// Split constructor
  AddOtherPartialDerivatives(AddOtherPartialDerivatives &rhs, split)
  :
    _Transformation(rhs._Transformation),
    _NPGradient    (rhs._NPGradient),
    _ImageGradient (rhs._ImageGradient),
    _Image2World   (rhs._Image2World),
    _Weight        (rhs._Weight),
    _t             (rhs._t),
    _t0            (rhs._t0),
    _X             (rhs._X),
    _Y             (rhs._Y),
    _NumberOfDOFs  (rhs._NumberOfDOFs)
  {
    _Output = new double[_NumberOfDOFs];
    memset(_Output, 0, _NumberOfDOFs * sizeof(double));
  }

  // ---------------------------------------------------------------------------
  /// Destructor
  ~AddOtherPartialDerivatives() {}

  // ---------------------------------------------------------------------------
  /// Join results of right-hand body with this body
  void join(AddOtherPartialDerivatives &rhs)
  {
    for (int dof = 0; dof < _NumberOfDOFs; ++dof) {
      _Output[dof] += rhs._Output[dof];
    }
    delete[] rhs._Output;
    rhs._Output = NULL;
  }

  // ---------------------------------------------------------------------------
  /// Calculates the gradient of the similarity term w.r.t. the transformation
  /// parameters for each discrete voxel of the specified image region
  void operator ()(const blocked_range3d<int> &re) const
  {
    const int i1 = re.cols ().begin();
    const int j1 = re.rows ().begin();
    const int k1 = re.pages().begin();
    const int i2 = re.cols ().end();
    const int j2 = re.rows ().end();
    const int k2 = re.pages().end();

    Matrix dJdp; // Derivative of the Jacobian of transformation (w.r.t. world coordinates) w.r.t DoF

    const auto *gx = _NPGradient->Data(i1, j1, k1, 0);
    const auto *gy = _NPGradient->Data(i1, j1, k1, 1);
    const auto *gz = _NPGradient->Data(i1, j1, k1, 2);

    const auto *ix = _ImageGradient->Data(i1, j1, k1, 0);
    const auto *iy = _ImageGradient->Data(i1, j1, k1, 1);
    const auto *iz = _ImageGradient->Data(i1, j1, k1, 2);

    // s1=1
    const int s2 =  _X - (i2 - i1);
    const int s3 = (_Y - (j2 - j1)) * _X;

    // With pre-computed world coordinates
    if (_Image2World) {
      const double *wx = _Image2World->Data(i1, j1, k1, 0);
      const double *wy = _Image2World->Data(i1, j1, k1, 1);
      const double *wz = _Image2World->Data(i1, j1, k1, 2);
      // Loop over image voxels
      for (int k = k1; k < k2; ++k, wx += s3, wy += s3, wz += s3, gx += s3, gy += s3, gz += s3, ix += s3, iy += s3, iz += s3)
      for (int j = j1; j < j2; ++j, wx += s2, wy += s2, wz += s2, gx += s2, gy += s2, gz += s2, ix += s2, iy += s2, iz += s2)
      for (int i = i1; i < i2; ++i, wx +=  1, wy +=  1, wz +=  1, gx +=  1, gy +=  1, gz +=  1, ix +=  1, iy +=  1, iz +=  1) {
        // Check whether reference point is valid
        if (*gx || *gy || *gz) {
          // Loop over DoFs
          for (int dof = 0; dof < _NumberOfDOFs; ++dof) {
            // Check status of DoF
            if (_Transformation->GetStatus(dof) == Active) {
              // Add gradient term
              _Transformation->DeriveJacobianWrtDOF(dJdp, dof, *wx, *wy, *wz, _t, _t0);
              _Output[dof] += _Weight * ((*ix) * (dJdp(0, 0) * (*gx) + dJdp(0, 1) * (*gy) + dJdp(0, 2) * (*gz)) +
                                         (*iy) * (dJdp(1, 0) * (*gx) + dJdp(1, 1) * (*gy) + dJdp(1, 2) * (*gz)) +
                                         (*iz) * (dJdp(2, 0) * (*gx) + dJdp(2, 1) * (*gy) + dJdp(2, 2) * (*gz)));
            }
          }
        }
      }
    // Without pre-computed world coordinates
    } else {
      double x, y, z;
      // Loop over image voxels
      for (int k = k1; k < k2; ++k, gx += s3, gy += s3, gz += s3, ix += s3, iy += s3, iz += s3)
      for (int j = j1; j < j2; ++j, gx += s2, gy += s2, gz += s2, ix += s2, iy += s2, iz += s2)
      for (int i = i1; i < i2; ++i, gx +=  1, gy +=  1, gz +=  1, ix +=  1, iy +=  1, iz +=  1) {
        // Check whether reference point is valid
        if (*gx || *gy || *gz) {
          // Convert voxel to world coordinates
          x = i, y = j, z = k;
          _NPGradient->ImageToWorld(x, y, z);
          // Loop over DoFs
          for (int dof = 0; dof < _NumberOfDOFs; ++dof) {
            // Check status of DoF
            if (_Transformation->GetStatus(dof) == Active) {
              // Add gradient term
              _Transformation->DeriveJacobianWrtDOF(dJdp, dof, x, y, z, _t, _t0);
              _Output[dof] += _Weight * ((*ix) * (dJdp(0, 0) * (*gx) + dJdp(0, 1) * (*gy) + dJdp(0, 2) * (*gz)) +
                                         (*iy) * (dJdp(1, 0) * (*gx) + dJdp(1, 1) * (*gy) + dJdp(1, 2) * (*gz)) +
                                         (*iz) * (dJdp(2, 0) * (*gx) + dJdp(2, 1) * (*gy) + dJdp(2, 2) * (*gz)));
            }
          }
        }
      }
    }
  }

  // ---------------------------------------------------------------------------
  void operator ()()
  {
    // Initialize members to often accessed data
    _X            = _NPGradient->X();
    _Y            = _NPGradient->Y();
    const int _Z  = _NPGradient->Z();
    _NumberOfDOFs = _Transformation->NumberOfDOFs();
    // Check input
    if (_Image2World && (_Image2World->X() != _X ||
                         _Image2World->Y() != _Y ||
                         _Image2World->Z() != _Z ||
                         _Image2World->T() !=  3)) {
      cerr << "GradientFieldSimilarity::ParametricGradient: Invalid world coordinates map" << endl;
      exit(1);
    }
    // Calculate parametric gradient
    blocked_range3d<int> voxels(0, _Z, 0, _Y, 0, _X);
    parallel_reduce(voxels, *this);
  }

}; // AddOtherPartialDerivatives


// -----------------------------------------------------------------------------
template <class FreeFormTransformation>
class AddOtherPartialDerivativesOfFFD
{
public:
  const FreeFormTransformation                   *_Transformation;
  const ImageSimilarity::GradientImageType       *_NPGradient;
  const GenericImage<RegisteredImage::VoxelType> *_ImageGradient;
  double                                         *_Output;
  double                                          _Weight;

  double _t;            ///< Time corrresponding to input gradient image (in ms)
  double _t0;           ///< Second time argument for velocity-based transformations

  int    _NumberOfDOFs; ///< Number of transformation parameters

  // ---------------------------------------------------------------------------
  /// Default constructor
  AddOtherPartialDerivativesOfFFD()
  :
    _Transformation(NULL),
    _NPGradient    (NULL),
    _ImageGradient (NULL),
    _Output        (NULL),
    _Weight        ( 1.0),
    _t             ( 0.0),
    _t0            (-1.0),
    _NumberOfDOFs  (0)
  {}

  // ---------------------------------------------------------------------------
  /// Split constructor
  AddOtherPartialDerivativesOfFFD(AddOtherPartialDerivativesOfFFD &rhs, split)
  :
    _Transformation(rhs._Transformation),
    _NPGradient    (rhs._NPGradient),
    _ImageGradient (rhs._ImageGradient),
    _Weight        (rhs._Weight),
    _t             (rhs._t),
    _t0            (rhs._t0),
    _NumberOfDOFs  (rhs._NumberOfDOFs)
  {
    _Output = new double[_NumberOfDOFs];
    memset(_Output, 0, _NumberOfDOFs * sizeof(double));
  }

  // ---------------------------------------------------------------------------
  /// Destructor
  ~AddOtherPartialDerivativesOfFFD() {}

  // ---------------------------------------------------------------------------
  /// Join results of right-hand body with this body
  void join(AddOtherPartialDerivativesOfFFD &rhs)
  {
    for (int dof = 0; dof < _NumberOfDOFs; ++dof) {
      _Output[dof] += rhs._Output[dof];
    }
    delete[] rhs._Output;
    rhs._Output = NULL;
  }

  // ---------------------------------------------------------------------------
  /// Calculates the gradient of the similarity term w.r.t. the transformation
  /// parameters for each discrete voxel of the specified image region
  void operator ()(const blocked_range<int> &dofs) const
  {
    Matrix dJdp; // Derivative of the Jacobian of transformation (w.r.t. world coordinates) w.r.t DoF

    double x, y, z, gx, gy, gz, ix, iy, iz;
    int i1, j1, k1, i2, j2, k2;

    // Loop over transformation parameters
    for (int dof = dofs.begin(); dof < dofs.end(); ++dof) {
      if (_Transformation->GetStatus(dof) == Active) {
        // Get support region of parameter
        _Transformation->DOFBoundingBox(_NPGradient, dof, i1, j1, k1, i2, j2, k2);
        // Loop over voxels in bounding box
        for (int k = k1; k < k2; ++k)
        for (int j = j1; j < j2; ++j)
        for (int i = i1; i < i2; ++i) {
          // Gradient computed by NonParametericGradient
          gx = _NPGradient->GetAsDouble(i, j, k, 0);
          gy = _NPGradient->GetAsDouble(i, j, k, 1);
          gz = _NPGradient->GetAsDouble(i, j, k, 2);
          // Check whether reference point is valid
          if (gx || gy || gz) {
            // Convert voxel to world coordinates
            x = i, y = j, z = k;
            _NPGradient->ImageToWorld(x, y, z);
            // Obtain image gradient values
            ix = _ImageGradient->GetAsDouble(i, j, k, 0);
            iy = _ImageGradient->GetAsDouble(i, j, k, 1);
            iz = _ImageGradient->GetAsDouble(i, j, k, 2);
            // Multiply according to chain rule and add to output gradient
            _Transformation->DeriveJacobianWrtDOF(dJdp, dof, x, y, z, _t, _t0);
            _Output[dof] += _Weight * (ix * (dJdp(0, 0) * gx + dJdp(0, 1) * gy + dJdp(0, 2) * gz) +
                                       iy * (dJdp(1, 0) * gx + dJdp(1, 1) * gy + dJdp(1, 2) * gz) +
                                       iz * (dJdp(2, 0) * gx + dJdp(2, 1) * gy + dJdp(2, 2) * gz));
          }
        }
      }
    }
  }

  // ---------------------------------------------------------------------------
  void operator ()()
  {
    // Initialize members to often accessed data
    _NumberOfDOFs = _Transformation->NumberOfDOFs();
    // Calculate parametric gradient
    parallel_reduce(blocked_range<int>(0, _NumberOfDOFs), *this);
  }

}; // AddOtherPartialDerivativesOfFFD


} // namespace GradientFieldSimilarityUtils
using namespace GradientFieldSimilarityUtils;

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
GradientFieldSimilarity::GradientFieldSimilarity(const char *name, double weight)
:
  ImageSimilarity(name, weight),
  _IgnoreJacobianGradientWrtDOFs(false)
{
}

// -----------------------------------------------------------------------------
GradientFieldSimilarity::GradientFieldSimilarity(const GradientFieldSimilarity &other)
:
  ImageSimilarity(other),
  _IgnoreJacobianGradientWrtDOFs(false)
{
}

// -----------------------------------------------------------------------------
GradientFieldSimilarity::~GradientFieldSimilarity()
{
}

// =============================================================================
// Parameters
// =============================================================================

// -----------------------------------------------------------------------------
bool GradientFieldSimilarity::Set(const char *name, const char *value)
{
  if (strcmp(name, "Ignore jacobian gradient w.r.t DOFs")           == 0 ||
      strcmp(name, "Ignore jacobian gradient w.r.t. DOFs")          == 0 ||
      strcmp(name, "Ignore jacobian gradient wrt DOFs")             == 0 ||
      strcmp(name, "Ignore jacobian gradient with respect to DOFs") == 0) {
    return FromString(value, _IgnoreJacobianGradientWrtDOFs);
  }
  return ImageSimilarity::Set(name, value);
}

// -----------------------------------------------------------------------------
ParameterList GradientFieldSimilarity::Parameter() const
{
  ParameterList params = ImageSimilarity::Parameter();
  Insert(params, "Ignore jacobian gradient w.r.t. DOFs", ToString(_IgnoreJacobianGradientWrtDOFs));
  return params;
}

// =============================================================================
// Initialization/Update
// =============================================================================

// -----------------------------------------------------------------------------
/// Reorients the transformed image derivatives
///
/// The gradient vectors are reoriented according to dI(y)/dy * dy/dx,
/// where y = T(x) and dI(y)/dy is the transformed image gradient.
/// The 2nd order derivatives needed for the similarity gradient computation
/// are optionally post-multiplied by the Jacobian of the transformation
/// as well, according to d2I(y)/dydx = d2I(y)/dy2 * dy/dx.
///
/// The template argument must be either one of
/// - Transformation1
/// - Transformation2
/// - MultiLevelTransformation
/// - FluidTransformation
/// The static Jacobian member function of this template type is used to
/// compute the Jacobian of the total transformation, i.e., dT/dx.
class ReorientTransformedDerivatives : public VoxelFunction
{
protected:

  typedef WorldCoordsImage::VoxelType CoordType;

  const RegisteredImage *_Image;
  static const int       _x = 0;
  int                    _y, _z;
  int                    _dx, _dy, _dz;
  int                    _dxx, _dxy, _dxz, _dyx, _dyy, _dyz, _dzx, _dzy, _dzz;
  bool                   _UpdateHessian;

  template <class T>
  void ReorientGradient(T *dI, const Matrix &dTdx) const
  {
    double dIdx = dI[_dx] * dTdx(0, 0) + dI[_dy] * dTdx(1, 0) + dI[_dz] * dTdx(2, 0);
    double dIdy = dI[_dx] * dTdx(0, 1) + dI[_dy] * dTdx(1, 1) + dI[_dz] * dTdx(2, 1);
    double dIdz = dI[_dx] * dTdx(0, 2) + dI[_dy] * dTdx(1, 2) + dI[_dz] * dTdx(2, 2);

    dI[_dx] = dIdx, dI[_dy] = dIdy, dI[_dz] = dIdz;
  }

  template <class T>
  void ReorientHessian(T *dI, const Matrix &dTdx) const
  {
    double dIdxx = dI[_dxx] * dTdx(0, 0) + dI[_dxy] * dTdx(1, 0) + dI[_dxz] * dTdx(2, 0);
    double dIdxy = dI[_dxx] * dTdx(0, 1) + dI[_dxy] * dTdx(1, 1) + dI[_dxz] * dTdx(2, 1);
    double dIdxz = dI[_dxx] * dTdx(0, 2) + dI[_dxy] * dTdx(1, 2) + dI[_dxz] * dTdx(2, 2);

    double dIdyx = dI[_dxy] * dTdx(0, 0) + dI[_dyy] * dTdx(1, 0) + dI[_dyz] * dTdx(2, 0);
    double dIdyy = dI[_dxy] * dTdx(0, 1) + dI[_dyy] * dTdx(1, 1) + dI[_dyz] * dTdx(2, 1);
    double dIdyz = dI[_dxy] * dTdx(0, 2) + dI[_dyy] * dTdx(1, 2) + dI[_dyz] * dTdx(2, 2);

    double dIdzx = dI[_dxz] * dTdx(0, 0) + dI[_dyz] * dTdx(1, 0) + dI[_dzz] * dTdx(2, 0);
    double dIdzy = dI[_dxz] * dTdx(0, 1) + dI[_dyz] * dTdx(1, 1) + dI[_dzz] * dTdx(2, 1);
    double dIdzz = dI[_dxz] * dTdx(0, 2) + dI[_dyz] * dTdx(1, 2) + dI[_dzz] * dTdx(2, 2);

    dI[_dxx] = dIdxx, dI[_dxy] = dIdxy, dI[_dxz] = dIdxz;
    dI[_dyx] = dIdyx, dI[_dyy] = dIdyy, dI[_dyz] = dIdyz;
    dI[_dzx] = dIdzx, dI[_dzy] = dIdzy, dI[_dzz] = dIdzz;
  }

  template <class T>
  void Reorient(T *dI, double x, double y, double z) const
  {
    Matrix dTdx;
    _Image->Transformation()->Jacobian(dTdx, x, y, z);
    ReorientGradient(dI, dTdx);
    if (_UpdateHessian) ReorientHessian(dI, dTdx);
  }

public:

  ReorientTransformedDerivatives(const RegisteredImage *image, bool hessian)
  :
    _Image        (image),
    _y            (image->NumberOfVoxels()),
    _z            (2 * _y),
    _dx           (image->Offset(RegisteredImage::Dx)),
    _dy           (image->Offset(RegisteredImage::Dy)),
    _dz           (image->Offset(RegisteredImage::Dz)),
    _dxx          (image->Offset(RegisteredImage::Dxx)),
    _dxy          (image->Offset(RegisteredImage::Dxy)),
    _dxz          (image->Offset(RegisteredImage::Dxz)),
    _dyx          (image->Offset(RegisteredImage::Dyx)),
    _dyy          (image->Offset(RegisteredImage::Dyy)),
    _dyz          (image->Offset(RegisteredImage::Dyz)),
    _dzx          (image->Offset(RegisteredImage::Dzx)),
    _dzy          (image->Offset(RegisteredImage::Dzy)),
    _dzz          (image->Offset(RegisteredImage::Dzz)),
    _UpdateHessian(hessian)
  {}

  template <class T>
  void operator()(int i, int j, int k, int, T *dI)
  {
    double x = i, y = j, z = k;
    _Image->ImageToWorld(x, y, z);
    Reorient(dI, x, y, z);
  }

  template <class Image, class T>
  void operator()(const Image &, int, const CoordType *wc, T *dI)
  {
    Reorient(dI, wc[_x], wc[_y], wc[_z]);
  }
};

// -----------------------------------------------------------------------------
void GradientFieldSimilarity::ReorientGradient(RegisteredImage *image, bool hessian)
{
  if (image->Transformation()) {
    ReorientTransformedDerivatives reorient(image, hessian);
    if (image->ImageToWorld()) {
      ParallelForEachVoxel(image->ImageToWorld(), image, reorient);
    } else {
      ParallelForEachVoxel(Domain(), image, reorient);
    }
  }
}

// -----------------------------------------------------------------------------
void GradientFieldSimilarity::InitializeInput(const ImageAttributes &domain)
{
  // Force precomputation of gradients
  _Target->PrecomputeDerivatives(true);
  _Source->PrecomputeDerivatives(true);

  // Initialize registered input images
  _Target->Initialize(domain, _Target->Transformation() ? 13 : 4);
  _Source->Initialize(domain, _Source->Transformation() ? 13 : 4);
  // Memorize image domain on which to evaluate similarity
  _Domain = domain;
}

// -----------------------------------------------------------------------------
void GradientFieldSimilarity
::ParametricGradient(const RegisteredImage *image,
                     GradientImageType     *np_gradient,
                     double                *gradient,
                     double                 weight)
{
  if (_IgnoreJacobianGradientWrtDOFs) {
    // Multiply np_gradient by (reoriented) Hessian of image
    MultiplyByImageHessian(image, np_gradient);
    ImageSimilarity::ParametricGradient(image, np_gradient, gradient, weight);
  } else {
    const class Transformation *T = image->Transformation();
    mirtkAssert(T != NULL, "image is being transformed");

    // Multiply np_gradient by (reoriented) Hessian of image (NOT in place!!)
    GradientImageType tmp(*np_gradient);
    MultiplyByImageHessian(image, &tmp);

    // Multiply gradient by dy/dx
    const WorldCoordsImage *i2w = image->ImageToWorld();
    const double                t0  = image->GetTOrigin();
    np_gradient->PutTOrigin(image->InputImage()->GetTOrigin());
    T->ParametricGradient(&tmp, gradient, i2w, t0, weight);

    // Select transformed gradient of image
    const GenericImage<RegisteredImage::VoxelType> &dIdy =
        (image == Target() ? _TargetTransformedGradient : _SourceTransformedGradient);

    // Select transformation type
    const FreeFormTransformation           *ffd  = NULL;
    const MultiLevelFreeFormTransformation *mffd = NULL;

    // Add np_gradient * [d(dy/dx)/dp]^T * [dI(y)/dy]^T  to gradient (algorithm depends on the transformation type)
    if ((mffd = dynamic_cast<const MultiLevelFreeFormTransformation *>(T))) {
      AddOtherPartialDerivativesOfFFD<MultiLevelFreeFormTransformation> body;
      body._Transformation = mffd;
      body._NPGradient     = np_gradient;
      body._ImageGradient  = &dIdy;
      body._Output         = gradient;
      body._Weight         = weight;
      body._t              = np_gradient->ImageToTime(.0);
      body._t0             = t0;
      body();
    } else if ((ffd = dynamic_cast<const FreeFormTransformation *>(T))) {
      AddOtherPartialDerivativesOfFFD<FreeFormTransformation> body;
      body._Transformation = ffd;
      body._NPGradient     = np_gradient;
      body._ImageGradient  = &dIdy;
      body._Output         = gradient;
      body._Weight         = weight;
      body._t              = np_gradient->ImageToTime(.0);
      body._t0             = t0;
      body();
    } else {
      AddOtherPartialDerivatives body;
      body._Transformation = T;
      body._NPGradient     = np_gradient;
      body._ImageGradient  = &dIdy;
      body._Output         = gradient;
      body._Weight         = weight;
      body._Image2World    = i2w;
      body._t              = np_gradient->ImageToTime(.0);
      body._t0             = t0;
      body();
    }
  }
}

// -----------------------------------------------------------------------------
void GradientFieldSimilarity::Update(bool hessian)
{
  if (_InitialUpdate || _Target->Transformation()) {
    // Update target image
    _Target->Update(true, true, hessian, _InitialUpdate);

    // Update copy of the transformed gradient of the target image
    if (_Target->Transformation() && !_IgnoreJacobianGradientWrtDOFs) {
      if (_InitialUpdate) {
        _TargetTransformedGradient.Initialize(_Target->Attributes(), 3);
      }
      _TargetTransformedGradient.CopyFrom(_Target->Data(0, 0, 0, 1));
    }

    // Reorient gradient and hessian (if needed) of the target image according
    // to the jacobian of the transformation
    this->ReorientGradient(_Target, hessian);
  }
  if (_InitialUpdate || _Source->Transformation()) {
    // Update source image
    _Source->Update(true, true, hessian, _InitialUpdate);

    // Update copy of the transformed gradient of the source image
    if (_Source->Transformation() && !_IgnoreJacobianGradientWrtDOFs) {
      if (_InitialUpdate) {
        _SourceTransformedGradient.Initialize(_Source->Attributes(), 3);
      }
      _SourceTransformedGradient.CopyFrom(_Source->Data(0, 0, 0, 1));
    }

    // Reorient gradient and hessian (if needed) of the source image according
    // to the jacobian of the transformation
    this->ReorientGradient(_Source, hessian);
  }
  _InitialUpdate = false;
}

// -----------------------------------------------------------------------------
void GradientFieldSimilarity
::MultiplyByImageHessian(const RegisteredImage *image, GradientImageType *gradient)
{
  MultiplyGradientByHessian times_hessian(image);
  ParallelForEachVoxel(image, gradient, times_hessian);
}


} // namespace mirtk
