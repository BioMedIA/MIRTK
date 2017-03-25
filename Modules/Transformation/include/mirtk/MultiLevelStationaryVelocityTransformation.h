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

#ifndef MIRTK_MultiLevelStationaryVelocityTransformation_H
#define MIRTK_MultiLevelStationaryVelocityTransformation_H

#include "mirtk/MultiLevelTransformation.h"
#include "mirtk/BSplineFreeFormTransformationSV.h"

#include "mirtk/Array.h"
#include "mirtk/Matrix.h"
#include "mirtk/EventDelegate.h"
#include "mirtk/VelocityToDisplacementFieldSS.h"
#include "mirtk/VoxelFunction.h"
#include "mirtk/Profiling.h"


namespace mirtk {


/**
 * Multi-level SV FFD where global and local transformations are summed up
 * in the log space before exponentiation
 *
 * T_msvffd(x) = exp(log(T_global(x)) + sum_i log(T_local^i(x)))
 */
class MultiLevelStationaryVelocityTransformation : public MultiLevelTransformation
{
  mirtkTransformationMacro(MultiLevelStationaryVelocityTransformation);

protected:

  /// Logarithm of global transformation matrix
  mirtkAttributeMacro(Matrix, LogA);

  /// Observes changes of global transformation matrix
  EventDelegate _GlobalTransformationObserver;

  /// Update logarithm of global transformation matrix
  void UpdateLogMatrix();

public:

  /// Get n-th local SV FFD
  BSplineFreeFormTransformationSV *SVFFD(int n);

  /// Get n-th local SV FFD
  const BSplineFreeFormTransformationSV *SVFFD(int n) const;

  /// Get active SV FFD whose parameters are being optimized
  BSplineFreeFormTransformationSV *ActiveSVFFD();

  /// Get active SV FFD whose parameters are being optimized
  const BSplineFreeFormTransformationSV *ActiveSVFFD() const;

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Default constructor
  MultiLevelStationaryVelocityTransformation();

  /// Construct multi-level transformation given a rigid transformation
  MultiLevelStationaryVelocityTransformation(const RigidTransformation &);

  /// Construct multi-level transformation given an affine transformation
  MultiLevelStationaryVelocityTransformation(const AffineTransformation &);

  /// Copy constructor
  MultiLevelStationaryVelocityTransformation(const MultiLevelStationaryVelocityTransformation &);

  /// Destructor
  virtual ~MultiLevelStationaryVelocityTransformation();

  // ---------------------------------------------------------------------------
  // Transformation parameters (DoFs)

  /// Get norm of the gradient vector
  virtual double DOFGradientNorm(const double *) const;

  /// Get number of transformation parameters
  virtual int NumberOfDOFs() const;

  /// Put value of transformation parameter
  virtual void Put(int, double);

  /// Put values of transformation parameters
  virtual void Put(const DOFValue *);

  /// Add change to transformation parameters
  virtual void Add(const DOFValue *);

  /// Update transformation parameters given parametric gradient
  virtual double Update(const DOFValue *);

  /// Get value of transformation parameter
  virtual double Get(int) const;

  /// Get values of transformation parameters
  virtual void Get(DOFValue *) const;

  /// Put status of transformation parameter
  virtual void PutStatus(int, DOFStatus);

  /// Get status of transformation parameter
  virtual DOFStatus GetStatus(int) const;

  // ---------------------------------------------------------------------------
  // Levels

  /// Combine local transformations on stack
  virtual void CombineLocalTransformation();

  /// Convert the global transformation from a matrix representation to a
  /// FFD and incorporate it with any existing local transformation
  virtual void MergeGlobalIntoLocalDisplacement();

  // ---------------------------------------------------------------------------
  // Point transformation

  using MultiLevelTransformation::Transform;
  using MultiLevelTransformation::Inverse;
  using MultiLevelTransformation::Displacement;
  using MultiLevelTransformation::InverseDisplacement;

  /// Whether the caching of the transformation displacements is required
  /// (or preferred) by this transformation. For some transformations such as
  /// those parameterized by velocities, caching of the displacements for
  /// each target voxel results in better performance or is needed for example
  /// for the scaling and squaring method.
  virtual bool RequiresCachingOfDisplacements() const;

protected:

  /// Upper integration limit used given the temporal origin of both target
  /// and source images. If both images have the same temporal origin, the value
  /// of the member variable @c T is returned. Note that @c T may be set to zero
  /// in order to force no displacement between images located at the same point
  /// in time to be zero. This is especially useful when animating the
  /// deformation between two images with successively increasing upper
  /// integration limit. Otherwise, if the temporal origin of the two images
  /// differs, the signed difference between these is returned, i.e., t - t0,
  /// which corresponds to a forward integration of the target point (x, y, z)
  /// to the time point of the source image.
  double UpperIntegrationLimit(double t, double t0) const;

  /// Get number of integration steps for a given temporal integration interval
  int NumberOfStepsForIntervalLength(double) const;

  /// Set number of integration steps for a given temporal integration interval
  void NumberOfStepsForIntervalLength(double, int) const;

  /// Whether to use scaling and squaring when possible
  bool UseScalingAndSquaring() const;

  /// Whether to use fast scaling and squaring when possible
  bool FastScalingAndSquaring() const;

  /// Whether scaling and squaring can be used to obtain the displacement field
  /// for each voxel of the specified vector field
  template <class VoxelType>
  bool CanUseScalingAndSquaring(const GenericImage<VoxelType> &) const;

  /// Maximum norm of scaled velocity for scaling and squaring
  double MaxScaledVelocity() const;

  /// Get stationary velocity field as 3D+t vector field
  template <class ScalarType>
  void VelocityComponents(int, int, GenericImage<ScalarType> &, bool = false) const;

  /// Get stationary velocity field
  template <class VoxelType>
  void Velocity(int, int, GenericImage<VoxelType> &, bool = false) const;

  /// Compute group exponential map using the scaling and squaring (SS) method
  ///
  /// \attention The scaling and squaring cannot be used to obtain the displacements
  ///            generated by a 3D velocity field on a 2D slice only. It always
  ///            must be applied to the entire domain of the velocity field.
  template <class VoxelType>
  void ScalingAndSquaring(int, int, GenericImage<VoxelType> &, double,
                          const WorldCoordsImage * = NULL) const;

  /// Transform point using the forward Euler integration method
  void RKE1(int, int, double &x, double &y, double &z, double t) const;

public:


  /// Transforms a single point
  virtual void Transform(int, int, double &, double &, double &, double = 0, double = NaN) const;

  /// Transforms a single point using the inverse of the transformation
  virtual bool Inverse(int, int, double &, double &, double &, double = 0, double = NaN) const;

  /// Calculates the displacement vectors for a whole image domain
  ///
  /// \attention The displacements are computed at the positions after applying the
  ///            current displacements at each voxel. These displacements are then
  ///            added to the current displacements. Therefore, set the input
  ///            displacements to zero if only interested in the displacements of
  ///            this transformation at the voxel positions.
  virtual void Displacement(int, int, GenericImage<double> &, double, double = NaN, const WorldCoordsImage * = NULL) const;

  /// Calculates the displacement vectors for a whole image domain
  ///
  /// \attention The displacements are computed at the positions after applying the
  ///            current displacements at each voxel. These displacements are then
  ///            added to the current displacements. Therefore, set the input
  ///            displacements to zero if only interested in the displacements of
  ///            this transformation at the voxel positions.
  virtual void Displacement(int, int, GenericImage<float> &, double, double = NaN, const WorldCoordsImage * = NULL) const;

  /// Calculates the inverse displacement vectors for a whole image domain
  ///
  /// \attention The displacements are computed at the positions after applying the
  ///            current displacements at each voxel. These displacements are then
  ///            added to the current displacements. Therefore, set the input
  ///            displacements to zero if only interested in the displacements of
  ///            this transformation at the voxel positions.
  ///
  /// \returns Always zero.
  virtual int InverseDisplacement(int, int, GenericImage<double> &, double, double = NaN, const WorldCoordsImage * = NULL) const;

  /// Calculates the inverse displacement vectors for a whole image domain
  ///
  /// \attention The displacements are computed at the positions after applying the
  ///            current displacements at each voxel. These displacements are then
  ///            added to the current displacements. Therefore, set the input
  ///            displacements to zero if only interested in the displacements of
  ///            this transformation at the voxel positions.
  ///
  /// \returns Always zero.
  virtual int InverseDisplacement(int, int, GenericImage<float> &, double, double = NaN, const WorldCoordsImage * = NULL) const;

  // ---------------------------------------------------------------------------
  // Derivatives
  using Transformation::ParametricGradient;

protected:

  /// Applies the chain rule to convert spatial non-parametric gradient
  /// to a gradient w.r.t the parameters of this transformation.
  virtual void ParametricGradient(const GenericImage<double> *, double *,
                                  const WorldCoordsImage *,
                                  const WorldCoordsImage *,
                                  double, double, double) const;

public:

  /// Applies the chain rule to convert spatial non-parametric gradient
  /// to a gradient w.r.t the parameters of this transformation.
  virtual void ParametricGradient(const GenericImage<double> *, double *,
                                  const WorldCoordsImage *,
                                  const WorldCoordsImage *,
                                  double = NaN, double = 1) const;

  // ---------------------------------------------------------------------------
  // I/O

  /// Prints the parameters of the transformation
  virtual void Print(ostream &, Indent = 0) const;

  // ---------------------------------------------------------------------------
  // Others

  /// Invert transformation
  virtual void Invert();

  friend class PartialMultiLevelStationaryVelocityTransformation;
};

////////////////////////////////////////////////////////////////////////////////
// Auxiliary voxel functions for subclasses
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
/// Add parameters of B-spline SV FFD
class AddDOFsOfBSplineSVFFD : public VoxelFunction
{
public:

  AddDOFsOfBSplineSVFFD(const BSplineFreeFormTransformationSV *input,
                        BaseImage                             *output)
  :
    _Input(input), _Output(output)
  {}

//  template <class T>
//  void Add(Vector2D<T> *v, double vx, double vy, double)
//  {
//    v->_x += vx, v->_y += vy;
//  }

  template <class T>
  void Add(Vector3D<T> *v, double vx, double vy, double vz)
  {
    v->_x += vx, v->_y += vy, v->_z += vz;
  }

  template <class VectorType>
  void operator()(int i, int j, int k, int, VectorType *v)
  {
    double vx, vy, vz;
    _Input->Get(i, j, k, vx, vy, vz);
    Add(v, vx, vy, vz);
  }

protected:
  const BSplineFreeFormTransformationSV *_Input;  ///< Input transformation
  BaseImage                             *_Output; ///< Output image
};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Levels
// =============================================================================

// -----------------------------------------------------------------------------
inline BSplineFreeFormTransformationSV *
MultiLevelStationaryVelocityTransformation::SVFFD(int n)
{
  if (n < 0) n += _NumberOfLevels;
  return reinterpret_cast<BSplineFreeFormTransformationSV *>(_LocalTransformation[n]);
}

// -----------------------------------------------------------------------------
inline const BSplineFreeFormTransformationSV *
MultiLevelStationaryVelocityTransformation::SVFFD(int n) const
{
  if (n < 0) n += _NumberOfLevels;
  return reinterpret_cast<BSplineFreeFormTransformationSV *>(_LocalTransformation[n]);
}

// -----------------------------------------------------------------------------
inline BSplineFreeFormTransformationSV *
MultiLevelStationaryVelocityTransformation::ActiveSVFFD()
{
  for (int l = _NumberOfLevels - 1; l >= 0; --l) {
    if (LocalTransformationIsActive(l)) return SVFFD(l);
  }
  return _NumberOfLevels > 0 ? SVFFD(_NumberOfLevels - 1) : NULL;
}

// -----------------------------------------------------------------------------
inline const BSplineFreeFormTransformationSV *
MultiLevelStationaryVelocityTransformation::ActiveSVFFD() const
{
  for (int l = _NumberOfLevels - 1; l >= 0; --l) {
    if (LocalTransformationIsActive(l)) return SVFFD(l);
  }
  return _NumberOfLevels > 0 ? SVFFD(_NumberOfLevels - 1) : NULL;
}

// =============================================================================
// Point transformation
// =============================================================================

// -----------------------------------------------------------------------------
inline double MultiLevelStationaryVelocityTransformation
::UpperIntegrationLimit(double t, double t0) const
{
  const BSplineFreeFormTransformationSV *svffd = ActiveSVFFD();
  return (svffd ? svffd->UpperIntegrationLimit(t, t0) : (t - t0));
}

// -----------------------------------------------------------------------------
inline int MultiLevelStationaryVelocityTransformation
::NumberOfStepsForIntervalLength(double T) const
{
  const BSplineFreeFormTransformationSV *svffd = ActiveSVFFD();
  return (svffd ? svffd->NumberOfStepsForIntervalLength(T) : 32);
}

// -----------------------------------------------------------------------------
inline void MultiLevelStationaryVelocityTransformation
::NumberOfStepsForIntervalLength(double T, int n) const
{
  const BSplineFreeFormTransformationSV *svffd = ActiveSVFFD();
  if (svffd) svffd->NumberOfStepsForIntervalLength(T, n);
}

// -----------------------------------------------------------------------------
inline bool MultiLevelStationaryVelocityTransformation
::UseScalingAndSquaring() const
{
  const BSplineFreeFormTransformationSV *svffd = ActiveSVFFD();
  return (svffd ? svffd->IntegrationMethod() == FFDIM_SS ||
                  svffd->IntegrationMethod() == FFDIM_FastSS
                : false);
}

// -----------------------------------------------------------------------------
inline bool MultiLevelStationaryVelocityTransformation
::FastScalingAndSquaring() const
{
  const BSplineFreeFormTransformationSV *svffd = ActiveSVFFD();
  return (svffd ? svffd->IntegrationMethod() == FFDIM_FastSS : false);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
bool MultiLevelStationaryVelocityTransformation
::CanUseScalingAndSquaring(const GenericImage<VoxelType> &d) const
{
  const BSplineFreeFormTransformationSV *svffd = ActiveSVFFD();
  return (svffd || (svffd->Z() <= 1 && d.Z() <= 1) ||
                   (svffd->Z() >  1 && d.Z() >  1));
}

// -----------------------------------------------------------------------------
inline double MultiLevelStationaryVelocityTransformation
::MaxScaledVelocity() const
{
  const BSplineFreeFormTransformationSV *svffd = ActiveSVFFD();
  return (svffd ? svffd->MaxScaledVelocity() : .0);
}

// -----------------------------------------------------------------------------
inline bool MultiLevelStationaryVelocityTransformation
::RequiresCachingOfDisplacements() const
{
  return this->UseScalingAndSquaring();
}

// -----------------------------------------------------------------------------
template <class VectorType>
void MultiLevelStationaryVelocityTransformation
::Velocity(int m, int n, GenericImage<VectorType> &v, bool coeff) const
{
  MIRTK_START_TIMING();
  if (n < 0 || n > _NumberOfLevels) n = _NumberOfLevels;
  const int l1 = (m < 0 ? 0 : m);
  const ImageAttributes &grid = v.Attributes();
  // Determine local SV FFDs whose coefficients can be added without prior deconvolution
  Array<bool> add_coeff(n, false);
  bool convert_to_coeff = false;
  if (coeff) {
    for (int l = l1; l < n; ++l) {
      add_coeff[l] = grid.EqualInSpace(SVFFD(l)->Attributes());
    }
  }
  // Add global velocities
  if (m < 0 && !_LogA.IsIdentity()) {
    ParallelForEachVoxel(EvaluateGlobalSVFFD(_LogA, &v), grid, v);
    convert_to_coeff = coeff;
  } else {
    v = .0;
  }
  // Add local velocities (which require deconvolution)
  for (int l = l1; l < n; ++l) {
    if (!add_coeff[l]) {
      ParallelForEachVoxel(AddBSplineSVFFD(SVFFD(l), &v), grid, v);
      convert_to_coeff = coeff;
    }
  }
  // Convert to cubic B-spline coefficients if needed and requested
  if (convert_to_coeff) ConvertToCubicBSplineCoefficients(v);
  // Add local velocities (which do not require deconvolution)
  for (int l = l1; l < n; ++l) {
    if (add_coeff[l]) {
      ParallelForEachVoxel(AddDOFsOfBSplineSVFFD(SVFFD(l), &v), grid, v);
    }
  }
  MIRTK_DEBUG_TIMING(3, "MultiLevelStationaryVelocityTransformation::Velocity");
}

template <> void MultiLevelStationaryVelocityTransformation
::Velocity(int, int, GenericImage<float>  &, bool) const;

template <> void MultiLevelStationaryVelocityTransformation
::Velocity(int, int, GenericImage<double> &, bool) const;

// -----------------------------------------------------------------------------
template <class VoxelType>
void MultiLevelStationaryVelocityTransformation
::ScalingAndSquaring(int m, int n, GenericImage<VoxelType> &d, double T,
                     const WorldCoordsImage * /*unused*/) const
{
  if (n < 0 || n > _NumberOfLevels) n = _NumberOfLevels;
  if (m < 0 && _GlobalTransformation.IsIdentity()) m = 0;
  if (m >= 0 && n == 0) return;
  VelocityToDisplacementFieldSS<VoxelType> exp;
  GenericImage<VoxelType>                  v, tmp;
  // Use only vector fields defined at control points with B-spline interpolation.
  // This results in an approximate solution due to the error at each squaring.
  if (n > 0 && FastScalingAndSquaring()) {
    v.Initialize(SVFFD(n-1)->Attributes(), 3);
    Velocity(m, n, v, true);
    exp.Interpolation(Interpolation_BSpline);
    exp.ComputeInterpolationCoefficients(false);
  // Evaluate velocities at output voxels beforehand and use
  // linear interpolation of dense vector fields during squaring
  } else {
    v.Initialize(d.Attributes(), 3);
    Velocity(m, n, v, false);
    exp.Interpolation(Interpolation_Linear);
  }
  // Exponentiate velocity field
  exp.UpperIntegrationLimit(T);
  exp.NumberOfSteps(NumberOfStepsForIntervalLength(T));
  exp.MaxScaledVelocity(static_cast<VoxelType>(MaxScaledVelocity()));
  exp.Upsample(false); // better, but too expensive
  exp.Input(0, &v);    // velocity field to be exponentiated
  exp.Input(1, &d);    // input displacement field (may be zero)
  exp.Output(&tmp);    // result is exp(v) o d
  exp.Run();
  d.CopyFrom(tmp);
  // Update number of integration steps to use
  //NumberOfStepsForIntervalLength(T, pow(2.0, exp.NumberOfSquaringSteps()));
}

// -----------------------------------------------------------------------------
inline void MultiLevelStationaryVelocityTransformation
::RKE1(int m, int n, double &x, double &y, double &z, double T) const
{
  if (n < 0 || n > _NumberOfLevels) n = _NumberOfLevels;
  // According to
  //   Bossa, M., Zacur, E., & Olmos, S. (2008). Algorithms for
  //   computing the group exponential of diffeomorphisms: Performance evaluation.
  //   In 2008 IEEE Computer Society Conference on Computer Vision and Pattern
  //   Recognition Workshops (pp. 1–8). IEEE. doi:10.1109/CVPRW.2008.4563005
  // it is more accurate to sum the displacements instead of updating the
  // transformed point directly because the components of the sum have similar
  // order of magnitude. This promises a lower numerical error. Therefore,
  // we sum the displacements dx, dy, dz and add them to the initial point
  // at each iteration in order to evaluate the next velocity.
  const int N = NumberOfStepsForIntervalLength(T);
  if (N > 0) {
    const int    l1 = (m < 0 ? 0 : m);
    const double dt = T / static_cast<double>(N);
    double xi, yi, zi, vx, vy, vz, dx = .0, dy = .0, dz = .0;
    for (int i = 0; i < N; ++i) {
      // Add current displacement to initial point
      xi = x + dx, yi = y + dy, zi = z + dz;
      // Add displacement of global step
      if (m < 0) {
        dx += (_LogA(0, 0) * xi + _LogA(0, 1) * yi + _LogA(0, 2) * zi + _LogA(0, 3)) * dt;
        dy += (_LogA(1, 0) * xi + _LogA(1, 1) * yi + _LogA(1, 2) * zi + _LogA(1, 3)) * dt;
        dz += (_LogA(2, 0) * xi + _LogA(2, 1) * yi + _LogA(2, 2) * zi + _LogA(2, 3)) * dt;
      }
      // Add displacements of local steps
      for (int l = l1; l < n; ++l) {
        vx = xi, vy = yi, vz = zi;
        SVFFD(l)->WorldToLattice(vx, vy, vz);
        SVFFD(l)->Evaluate      (vx, vy, vz);
        dx += vx * dt;
        dy += vy * dt;
        dz += vz * dt;
      }
    }
    x += dx, y += dy, z += dz;
  }
}

// -----------------------------------------------------------------------------
inline void MultiLevelStationaryVelocityTransformation::Transform(int m, int n, double &x, double &y, double &z, double t, double t0) const
{
  RKE1(m, n, x, y, z, + UpperIntegrationLimit(t, t0));
}

// -----------------------------------------------------------------------------
inline bool MultiLevelStationaryVelocityTransformation::Inverse(int m, int n, double &x, double &y, double &z, double t, double t0) const
{
  RKE1(m, n, x, y, z, - UpperIntegrationLimit(t, t0));
  return true;
}

// -----------------------------------------------------------------------------
inline void MultiLevelStationaryVelocityTransformation
::Displacement(int m, int n, GenericImage<double> &d, double t, double t0, const WorldCoordsImage *wc) const
{
  double T;
  if ((T = + UpperIntegrationLimit(t, t0))) {
    // Use scaling and squaring method to efficiently compute displacements when possible
    if (UseScalingAndSquaring() && CanUseScalingAndSquaring(d)) {
      ScalingAndSquaring(m, n, d, T, wc);
    // Evaluate transformation at each voxel separately using explicit integration
    } else {
      MultiLevelTransformation::Displacement(m, n, d, t, t0, wc);
    }
  }
}

// -----------------------------------------------------------------------------
inline void MultiLevelStationaryVelocityTransformation
::Displacement(int m, int n, GenericImage<float> &d, double t, double t0, const WorldCoordsImage *wc) const
{
  double T;
  if ((T = + UpperIntegrationLimit(t, t0))) {
    // Use scaling and squaring method to efficiently compute displacements when possible
    if (UseScalingAndSquaring() && CanUseScalingAndSquaring(d)) {
      ScalingAndSquaring(m, n, d, T, wc);
    // Evaluate transformation at each voxel separately using explicit integration
    } else {
      MultiLevelTransformation::Displacement(m, n, d, t, t0, wc);
    }
  }
}

// -----------------------------------------------------------------------------
inline int MultiLevelStationaryVelocityTransformation
::InverseDisplacement(int m, int n, GenericImage<double> &d, double t, double t0, const WorldCoordsImage *wc) const
{
  double T;
  if ((T = - UpperIntegrationLimit(t, t0))) {
    // Use scaling and squaring method to efficiently compute displacements when possible
    if (UseScalingAndSquaring() && CanUseScalingAndSquaring(d)) {
      ScalingAndSquaring(m, n, d, T, wc);
    // Evaluate transformation at each voxel separately using explicit integration
    } else {
      MultiLevelTransformation::InverseDisplacement(m, n, d, t, t0, wc);
    }
  }
  return 0;
}

// -----------------------------------------------------------------------------
inline int MultiLevelStationaryVelocityTransformation
::InverseDisplacement(int m, int n, GenericImage<float> &d, double t, double t0, const WorldCoordsImage *wc) const
{
  double T;
  if ((T = - UpperIntegrationLimit(t, t0))) {
    // Use scaling and squaring method to efficiently compute displacements when possible
    if (UseScalingAndSquaring() && CanUseScalingAndSquaring(d)) {
      ScalingAndSquaring(m, n, d, T, wc);
    // Evaluate transformation at each voxel separately using explicit integration
    } else {
      MultiLevelTransformation::InverseDisplacement(m, n, d, t, t0, wc);
    }
  }
  return 0;
}


} // namespace mirtk

#endif // MIRTK_MultiLevelStationaryVelocityTransformation_H
