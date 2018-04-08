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

#ifndef MIRTK_BSplineFreeFormTransformationSV_H
#define MIRTK_BSplineFreeFormTransformationSV_H

#include "mirtk/BSplineFreeFormTransformation3D.h"

#include "mirtk/Math.h"
#include "mirtk/ImageFunction.h"
#include "mirtk/VoxelFunction.h"
#include "mirtk/FFDIntegrationMethod.h"


namespace mirtk {


/**
 * Free-form transformation parameterized by a stationary velocity field.
 *
 * This class implements a free-form transformation which is represented
 * by a stationary velocity field (2D/3D). The displacement field is obtained
 * by exponentiating the velocity field, resulting in a diffeomorphic
 * transformation. A stationary velocity field can be further approximated
 * from a given displacement field using the group logarithm of diffeomorphisms.
 */
class BSplineFreeFormTransformationSV : public BSplineFreeFormTransformation3D
{
  mirtkTransformationMacro(BSplineFreeFormTransformationSV);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Cross-sectional upper integration limit, i.e., time difference between
  /// target and source if both volumes have the same temporal origin.
  /// Otherwise, the difference between the temporal origins as given by the
  /// t and t0 parameters of LocalTransform et al. is used instead.
  ///
  /// \note Can be set to zero for the generation of a deformation movie using
  ///       the transformation tool based on ImageTransformation by varying
  ///       the -Tt parameter from zero to one and keeping -St fixed at one.
  ///       Use the doftool to set the --cross-sectional-time-interval to zero.
  mirtkPublicAttributeMacro(double, T);

  /// Temporal unit used to determine actual number of integration steps.
  /// If zero, the number of steps is independent of the temporal interval.
  mirtkPublicAttributeMacro(double, TimeUnit);

  /// Number of integration steps per _TimeUnit
  mirtkPublicAttributeMacro(int, NumberOfSteps);

  /// Maximum absolute value of scaled velocity when using the
  /// scaling-and-squaring for computing the displacement field.
  /// This introduces an adapative choice of the number of integration steps
  /// needed. At a minimum, NumberOfSteps are taken, but if needed, more
  /// steps added. If this attribute is set to zero, however, the number
  /// of integration steps always corresponds to the specified NumberOfSteps.
  mirtkPublicAttributeMacro(double, MaxScaledVelocity);

  /// Integration method used to obtain displacement from velocities
  mirtkPublicAttributeMacro(FFDIntegrationMethod, IntegrationMethod);

  /// Whether to evaluate BCH formulat for parametric gradient calculation
  /// on dense image lattice of non-parametric gradient using linear interpolation
  mirtkPublicAttributeMacro(bool, UseDenseBCHGrid);

  /// Use Lie derivative definition of Lie bracket based on vector field
  /// Jacobian matrices instead of difference of vector field composition
  mirtkPublicAttributeMacro(bool, LieDerivative);

  /// Number of BCH terms used by ParametricGradient. When set to an invalid
  /// value, i.e., no. of terms < 2, the gradient is computed using the specified
  /// integration method instead.
  mirtkPublicAttributeMacro(int, NumberOfBCHTerms);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

protected:

  /// Initialize extrapolation function
  void InitializeExtrapolator();

public:

  /// Default constructor
  BSplineFreeFormTransformationSV();

  /// Construct free-form transformation for given image domain and lattice spacing
  explicit BSplineFreeFormTransformationSV(const ImageAttributes &,
                                           double = -1, double = -1, double = -1);

  /// Construct free-form transformation for given target image and lattice spacing
  explicit BSplineFreeFormTransformationSV(const BaseImage &,
                                           double, double, double);

  /// Construct free-form transformation from existing 3D+t deformation field
  explicit BSplineFreeFormTransformationSV(const GenericImage<double> &, bool = false);

  /// Copy constructor
  BSplineFreeFormTransformationSV(const BSplineFreeFormTransformationSV &);

  /// Destructor
  virtual ~BSplineFreeFormTransformationSV();

  // Import other Initialize overloads
  using BSplineFreeFormTransformation3D::Initialize;

  /// Initialize free-form transformation constructed by default constructor
  virtual void Initialize(const ImageAttributes &);

  /// Initialize free-form transformation which approximates a given transformation
  virtual void Initialize(const ImageAttributes &, double, double, double,
                          const Transformation *);

  /// Subdivide FFD lattice
  virtual void Subdivide(bool = true, bool = true, bool = true, bool = true);

  // ---------------------------------------------------------------------------
  // Approximation/Interpolation

  using BSplineFreeFormTransformation3D::ApproximateAsNew;

  /// Approximate displacements: This function takes a set of points and a set
  /// of displacements and finds !new! parameters such that the resulting
  /// transformation approximates the displacements as good as possible.
  virtual void ApproximateDOFs(const double *, const double *, const double *, const double *,
                               const double *, const double *, const double *, int);

  /// Finds gradient of approximation error: This function takes a set of points
  /// and a set of errors. It finds a gradient w.r.t. the transformation parameters
  /// which minimizes the L2 norm of the approximation error and adds it to the
  /// input gradient with the given weight.
  virtual void ApproximateDOFsGradient(const double *, const double *, const double *, const double *,
                                       const double *, const double *, const double *, int,
                                       double *, double = 1.0) const;

  /// Get image domain on which this free-form transformation should be defined
  /// in order to reduce the error when the given transformation is approximated
  /// and applied to resample an image on the specified image grid.
  virtual ImageAttributes ApproximationDomain(const ImageAttributes &,
                                              const Transformation *);

  /// Approximate another transformation and return approximation error
  virtual double ApproximateAsNew(const ImageAttributes &, const Transformation *, int = 1, double = .0);

  /// Approximate displacements: This function takes a 3D displacement field
  /// and finds a stationary velocity field which, when exponentiated, approximates
  /// these displacements by computing the group logarithm of diffeomorphisms.
  /// The displacements are replaced by the residual displacements of the newly
  /// approximated transformation. Returns the approximation error of the resulting FFD.
  virtual double ApproximateAsNew(GenericImage<double> &, int = 1, double = .0);

  /// Approximate displacements: This function takes a 3D displacement field
  /// and finds a stationary velocity field which, when exponentiated, approximates
  /// these displacements by computing the group logarithm of diffeomorphisms.
  /// The displacements are replaced by the residual displacements of the newly
  /// approximated transformation. Returns the approximation error of the resulting FFD.
  virtual double ApproximateAsNew(GenericImage<double> &, bool, int = 3, int = 8);

  /// Interpolates displacements: This function takes a set of displacements defined at
  /// the control points and finds a stationary velocity field such that the
  /// resulting transformation interpolates these displacements.
  virtual void Interpolate(const double *, const double *, const double *);

  /// Interpolates velocities: This function takes a set of velocites defined at
  /// the control points and interpolates these velocities.
  virtual void InterpolateVelocities(const double *, const double *, const double *);

  /// Approximate velocities: This function takes a velocity field and initializes
  /// the cubic B-spline coefficients of this transformation such that the spline
  /// function approximates the given velocity field. Returns the approximation error
  /// of the velocity spline function.
  virtual double ApproximateVelocitiesAsNew(GenericImage<double> &);

  /// Combine transformations: This function takes a transformation and finds
  /// a stationary velocity field which, when exponentiated, approximates the
  /// composition of this and the given transformation.
  virtual void CombineWith(const Transformation *);

  /// Combine transformations: This function takes a transformation and finds
  /// a stationary velocity field which, when exponentiated, approximates the
  /// composition of this and the given transformation.
  virtual void CombineWith(const BSplineFreeFormTransformationSV *);

  /// Invert this transformation
  virtual void Invert();

  // ---------------------------------------------------------------------------
  // Parameters (non-DoFs)

  using BSplineFreeFormTransformation3D::Parameter;

  /// Set named (non-DoF) parameter from value as string
  virtual bool Set(const char *, const char *);

  /// Get (non-DoF) parameters as key/value as string map
  virtual ParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Evaluation (for use FreeFormTransformationRungeKutta integration)

  using BSplineFreeFormTransformation3D::Evaluate;
  using BSplineFreeFormTransformation3D::EvaluateJacobianWorld;
  using BSplineFreeFormTransformation3D::EvaluateJacobianDOFs;

  /// Evaluates the FFD at a point in lattice coordinates
  void Evaluate(double &, double &, double &, double) const;

  /// Calculates the Jacobian of the FFD at a point in lattice coordinates
  /// and converts the resulting Jacobian to derivatives w.r.t world coordinates
  void EvaluateJacobianWorld(Matrix &, double, double, double, double) const;

  /// Calculates the Jacobian of the 2D transformation w.r.t. the transformation parameters
  void EvaluateJacobianDOFs(TransformationJacobian &,
                            double, double) const;

  /// Calculates the Jacobian of the 3D transformation w.r.t. the transformation parameters
  void EvaluateJacobianDOFs(TransformationJacobian &,
                            double, double, double) const;

  /// Calculates the Jacobian of the transformation w.r.t. the transformation parameters
  void EvaluateJacobianDOFs(TransformationJacobian &,
                            double, double, double, double) const;

  // ---------------------------------------------------------------------------
  // Point transformation

public:

  using BSplineFreeFormTransformation3D::Displacement;

protected:

  /// Get number of integration steps for a given temporal integration interval
  int NumberOfStepsForIntervalLength(double) const;

  /// Set number of integration steps for a given temporal integration interval
  void NumberOfStepsForIntervalLength(double, int) const;

  /// Get length of steps for a given temporal integration interval
  double StepLengthForIntervalLength(double) const;

public:

  /// Scale velocities
  ///
  /// Alternatively, change integration inveral _T.
  void ScaleVelocities(double);

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

  /// Transforms a single point using a Runge-Kutta integration
  ///
  /// \param[in] T Upper integration limit, i.e., length of time interval.
  void IntegrateVelocities(double &, double &, double &, double T = 1.0) const;

  /// Compute displacement field using the scaling and squaring (SS) method
  ///
  /// \attention If an input/output displacement field is provided, the resulting
  ///            transformation is the composition exp(v) o d.
  ///
  /// \attention The scaling and squaring cannot be used to obtain the displacements
  ///            generated by a 3D velocity field on a 2D slice only. It always
  ///            must be applied to the entire domain of the velocity field.
  ///
  /// \param[in,out] d  Displacement field generated by stationary velocity field.
  /// \param[in]     T  Upper integration limit, i.e., length of time interval.
  /// \param[in]     wc Pre-computed world coordinates of image voxels.
  template <class VoxelType>
  void ScalingAndSquaring(GenericImage<VoxelType> *d,
                          double T = 1.0, const WorldCoordsImage *wc = NULL) const;

  /// Compute displacement field and derivatives using the scaling and squaring (SS) method
  ///
  /// \attention If an input/output displacement field is provided, the resulting
  ///            transformation is the composition exp(v) o d.
  ///
  /// \attention The scaling and squaring cannot be used to obtain the displacements
  ///            generated by a 3D velocity field on a 2D slice only. It always
  ///            must be applied to the entire domain of the velocity field.
  ///
  /// \param[in]     attr Attributes of output images.
  /// \param[in,out] d    Displacement field generated by stationary velocity field.
  /// \param[out]    dx   Partial derivatives of corresponding transformation w.r.t. x.
  /// \param[out]    dj   Determinant of Jacobian w.r.t. x, i.e., det(jac(\p dx)).
  /// \param[out]    lj   Log of determinant of Jacobian w.r.t. x, i.e., log(det(jac(\p dx))).
  /// \param[in]     T    Upper integration limit, i.e., length of time interval.
  /// \param[in]     wc   Pre-computed world coordinates of image voxels.
  template <class VoxelType>
  void ScalingAndSquaring(const ImageAttributes   &attr,
                          GenericImage<VoxelType> *d,
                          GenericImage<VoxelType> *dx,
                          GenericImage<VoxelType> *dj = NULL,
                          GenericImage<VoxelType> *lj = NULL,
                          double T = 1.0, const WorldCoordsImage *wc = NULL) const;

  /// Whether the caching of the transformation displacements is required
  /// (or preferred) by this transformation. For some transformations such as
  /// those parameterized by velocities, caching of the displacements for
  /// each target voxel results in better performance or is needed for example
  /// for the scaling and squaring method.
  virtual bool RequiresCachingOfDisplacements() const;

  /// Transforms a single point using the local transformation component only
  virtual void LocalTransform(double &, double &, double &, double = 0, double = NaN) const;

  /// Transforms a single point using the inverse of the local transformation only
  virtual bool LocalInverse(double &, double &, double &, double = 0, double = NaN) const;

  /// Get stationary velocity field
  void Velocity(GenericImage<float> &) const;

  /// Get stationary velocity field
  void Velocity(GenericImage<double> &) const;

  /// Calculates the displacement vectors for a whole image domain
  ///
  /// \attention The displacements are computed at the positions after applying the
  ///            current displacements at each voxel. These displacements are then
  ///            added to the current displacements. Therefore, set the input
  ///            displacements to zero if only interested in the displacements of
  ///            this transformation at the voxel positions.
  virtual void Displacement(GenericImage<double> &, double, double = NaN,
                            const WorldCoordsImage * = NULL) const;
  
  /// Calculates the displacement vectors for a whole image domain
  ///
  /// \attention The displacements are computed at the positions after applying the
  ///            current displacements at each voxel. These displacements are then
  ///            added to the current displacements. Therefore, set the input
  ///            displacements to zero if only interested in the displacements of
  ///            this transformation at the voxel positions.
  virtual void Displacement(GenericImage<float> &, double, double = NaN,
                            const WorldCoordsImage * = NULL) const;

  /// Calculates the inverse displacement vectors for a whole image domain
  ///
  /// \attention The displacements are computed at the positions after applying the
  ///            current displacements at each voxel. These displacements are then
  ///            added to the current displacements. Therefore, set the input
  ///            displacements to zero if only interested in the displacements of
  ///            this transformation at the voxel positions.
  ///
  /// \returns Always zero.
  virtual int InverseDisplacement(GenericImage<double> &, double, double = NaN,
                                  const WorldCoordsImage * = NULL) const;

  /// Calculates the inverse displacement vectors for a whole image domain
  ///
  /// \attention The displacements are computed at the positions after applying the
  ///            current displacements at each voxel. These displacements are then
  ///            added to the current displacements. Therefore, set the input
  ///            displacements to zero if only interested in the displacements of
  ///            this transformation at the voxel positions.
  ///
  /// \returns Always zero.
  virtual int InverseDisplacement(GenericImage<float> &, double, double = NaN,
                                  const WorldCoordsImage * = NULL) const;

  // ---------------------------------------------------------------------------
  // Derivatives

  using BSplineFreeFormTransformation3D::JacobianDOFs;
  using BSplineFreeFormTransformation3D::ParametricGradient;

  /// Calculates the Jacobian of the local transformation w.r.t world coordinates
  virtual void LocalJacobian(Matrix &, double, double, double, double = 0, double = NaN) const;

  /// Calculates the Hessian for each component of the local transformation w.r.t world coordinates
  virtual void LocalHessian(Matrix [3], double, double, double, double = 0, double = NaN) const;

  /// Calculates the Jacobian of the transformation w.r.t a control point
  virtual void JacobianDOFs(Matrix &, int, double, double, double, double = 0, double = NaN) const;

  /// Calculates the Jacobian of the transformation w.r.t a transformation parameters
  ///
  /// \attention This overriden base class function returns one column of the 3x3 Jacobian
  ///            matrix, i.e., only the partial derivatives of the transformation
  ///            w.r.t. the specified transformation parameter (not control point).
  ///            This is in accordance to the Transformation::JacobianDOFs
  ///            definition, which is violated, however, by most other FFDs.
  virtual void JacobianDOFs(double [3], int, double, double, double, double = 0, double = NaN) const;

  /// Calculates the Jacobian of the transformation w.r.t the transformation parameters
  virtual void JacobianDOFs(TransformationJacobian &, double, double, double, double = 0, double = NaN) const;

  /// Calculates derivatives of the Jacobian determinant at world point w.r.t. DoFs of a control point
  ///
  /// \param[out] dJ  Partial derivatives of Jacobian determinant at (x, y, z) w.r.t. DoFs of control point.
  /// \param[in]  cp  Index of control point w.r.t. whose DoFs the derivatives are computed.
  /// \param[in]  x   World coordinate along x axis at which to evaluate derivatives.
  /// \param[in]  y   World coordinate along y axis at which to evaluate derivatives.
  /// \param[in]  z   World coordinate along z axis at which to evaluate derivatives.
  /// \param[in]  adj Adjugate of Jacobian matrix evaluated at (x, y, z).
  virtual void JacobianDetDerivative(double dJ[3], const Matrix &adj,
                                     int cp, double x, double y, double z, double t = 0, double t0 = NaN,
                                     bool wrt_world = true, bool use_spacing = true) const;

protected:

  /// Evaluates the BCH formula s.t. u = log(exp(tau * v) o exp(eta * w))
  /// \sa Bossa, M., Hernandez, M., & Olmos, S. Contributions to 3D diffeomorphic
  ///     atlas estimation: application to brain images. MICCAI 2007, 10(Pt 1), 667–74.
  void EvaluateBCHFormula(int, CPImage &, double, const CPImage &, double, const CPImage &, bool = false) const;

  /// Applies the chain rule to convert spatial non-parametric gradient
  /// to a gradient w.r.t the parameters of this transformation given the
  /// temporal interval which the displacement gradient corresponds to
  void ParametricGradient(const GenericImage<double> *, double *,
                          const WorldCoordsImage *,
                          const WorldCoordsImage *,
                          double, double, double) const;

public:

  /// Applies the chain rule to convert spatial non-parametric gradient
  /// to a gradient w.r.t the parameters of this transformation
  virtual void ParametricGradient(const GenericImage<double> *, double *,
                                  const WorldCoordsImage *,
                                  const WorldCoordsImage *,
                                  double = -1, double = 1) const;

  /// Applies the chain rule to convert point-wise non-parametric gradient
  /// to a gradient w.r.t the parameters of this transformation.
  virtual void ParametricGradient(const PointSet &, const Vector3D<double> *,
                                  double *, double = 0, double = -1, double = 1) const;

  // ---------------------------------------------------------------------------
  // I/O

  // Do not hide methods of base class
  using BSplineFreeFormTransformation3D::Print;

  /// Prints the parameters of the transformation
  virtual void Print(ostream &, Indent = 0) const;

  /// Whether this transformation can read a file of specified type (i.e. format)
  virtual bool CanRead(TransformationType) const;

protected:

  /// Reads transformation parameters from a file stream
  virtual Cifstream &ReadDOFs(Cifstream &, TransformationType);

  /// Writes transformation parameters to a file stream
  virtual Cofstream &WriteDOFs(Cofstream &) const;

  // ---------------------------------------------------------------------------
  // Friends
  friend class EvaluateBSplineSVFFD2D;
  friend class EvaluateBSplineSVFFD3D;
  friend class PartialBSplineFreeFormTransformationSV;
  friend class MultiLevelStationaryVelocityTransformation;
};

////////////////////////////////////////////////////////////////////////////////
// Auxiliary voxel functions for subclasses
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
/// Evaluate global SV FFD at image voxels of vector field
class EvaluateGlobalSVFFD : public VoxelFunction
{
public:

  EvaluateGlobalSVFFD(const Matrix &logA, BaseImage *output)
  :
    _LogA(logA), _Output(output)
  {}

//  template <class T>
//  void Evaluate(Vector2D<T> &v, double x, double y, double z)
//  {
//    v._x = static_cast<T>(_LogA(0, 0) * x + _LogA(0, 1) * y + _LogA(0, 2) * z + _LogA(0, 3));
//    v._y = static_cast<T>(_LogA(1, 0) * x + _LogA(1, 1) * y + _LogA(1, 2) * z + _LogA(1, 3));
//  }

  template <class T>
  void Evaluate(Vector3D<T> &v, double x, double y, double z)
  {
    v._x = static_cast<T>(_LogA(0, 0) * x + _LogA(0, 1) * y + _LogA(0, 2) * z + _LogA(0, 3));
    v._y = static_cast<T>(_LogA(1, 0) * x + _LogA(1, 1) * y + _LogA(1, 2) * z + _LogA(1, 3));
    v._z = static_cast<T>(_LogA(2, 0) * x + _LogA(2, 1) * y + _LogA(2, 2) * z + _LogA(2, 3));
  }

  template <class VectorType>
  void operator()(int i, int j, int k, int, VectorType *v)
  {
    double x = i, y = j, z = k;
    _Output->ImageToWorld(x, y, z);
    Evaluate(*v, x, y, z);
  }

private:
  const Matrix _LogA;   ///< Input transformation
  BaseImage   *_Output; ///< Output image
};

// -----------------------------------------------------------------------------
/// Evaluate global SV FFD at image voxels of 3D vector field
class EvaluateGlobalSVFFD3D : public VoxelFunction
{
public:

  EvaluateGlobalSVFFD3D(const Matrix &logA, BaseImage *output)
  :
    _LogA  (logA),
    _Output(output),
    _y     (output->NumberOfSpatialVoxels()),
    _z     (2 * _y)
  {}

  template <class T>
  void operator()(int i, int j, int k, int, T *out)
  {
    double x = i, y = j, z = k;
    _Output->ImageToWorld(x, y, z);
    out[_x] = static_cast<T>(_LogA(0, 0) * x + _LogA(0, 1) * y + _LogA(0, 2) * z + _LogA(0, 3));
    out[_y] = static_cast<T>(_LogA(1, 0) * x + _LogA(1, 1) * y + _LogA(1, 2) * z + _LogA(1, 3));
    out[_z] = static_cast<T>(_LogA(2, 0) * x + _LogA(2, 1) * y + _LogA(2, 2) * z + _LogA(2, 3));
  }

private:
  const Matrix _LogA;   ///< Input transformation
  BaseImage   *_Output; ///< Output image

  static const int _x = 0; ///< Offset of x component
  int              _y;     ///< Offset of y component
  int              _z;     ///< Offset of z component
};

// -----------------------------------------------------------------------------
/// Evaluate B-spline SV FFD at image voxels
class EvaluateBSplineSVFFD : public VoxelFunction
{
public:

  EvaluateBSplineSVFFD(const BSplineFreeFormTransformationSV *input,
                       BaseImage                             *output)
  :
    _Input(input), _Output(output)
  {}

//  template <class T>
//  void Put(Vector2D<T> *v, double vx, double vy, double)
//  {
//    v->_x = vx, v->_y = vy;
//  }

  template <class T>
  void Put(Vector3D<T> *v, double vx, double vy, double vz)
  {
    v->_x = static_cast<T>(vx);
    v->_y = static_cast<T>(vy);
    v->_z = static_cast<T>(vz);
  }

  template <class VectorType>
  void operator()(int i, int j, int k, int, VectorType *v)
  {
    double vx = i, vy = j, vz = k;
    _Output->ImageToWorld  (vx, vy, vz);
    _Input ->WorldToLattice(vx, vy, vz);
    _Input ->Evaluate      (vx, vy, vz);
    Put(v, vx, vy, vz);
  }

protected:
  const BSplineFreeFormTransformationSV *_Input;  ///< Input transformation
  BaseImage                             *_Output; ///< Output image
};

// -----------------------------------------------------------------------------
/// Evaluate B-spline SV FFD at image voxels
class EvaluateBSplineSVFFD3D : public VoxelFunction
{
public:

  EvaluateBSplineSVFFD3D(const BSplineFreeFormTransformationSV *input, BaseImage *output)
  :
    _Input (input),
    _Output(output),
    _y     (output->NumberOfSpatialVoxels()),
    _z     (2 * _y)
  {}

  template <class T>
  void operator()(int i, int j, int k, int, T *out)
  {
    double x = i, y = j, z = k;
    _Output->ImageToWorld  (x, y, z);
    _Input ->WorldToLattice(x, y, z);
    _Input ->Evaluate      (x, y, z);
    out[_x] = static_cast<T>(x);
    out[_y] = static_cast<T>(y);
    out[_z] = static_cast<T>(z);
  }

private:
  const BSplineFreeFormTransformationSV *_Input;  ///< Input transformation
  BaseImage                             *_Output; ///< Output image
  
  static const int _x = 0; ///< Offset of x component
  int              _y;     ///< Offset of y component
  int              _z;     ///< Offset of z component
};

// -----------------------------------------------------------------------------
/// Evaluate and add B-spline SV FFD at image voxels
class AddBSplineSVFFD : public VoxelFunction
{
public:

  AddBSplineSVFFD(const BSplineFreeFormTransformationSV *input,
                  BaseImage                             *output)
  :
    _Input(input), _Output(output)
  {}

//  template <class T>
//  void Add(Vector2D<T> *v, double vx, double vy, double)
//  {
//    v->_x += static_cast<T>(vx);
//    v->_y += static_cast<T>(vy);
//  }

  template <class T>
  void Add(Vector3D<T> *v, double vx, double vy, double vz)
  {
    v->_x += static_cast<T>(vx);
    v->_y += static_cast<T>(vy);
    v->_z += static_cast<T>(vz);
  }

  template <class VectorType>
  void operator()(int i, int j, int k, int, VectorType *v)
  {
    double vx = i, vy = j, vz = k;
    _Output->ImageToWorld  (vx, vy, vz);
    _Input ->WorldToLattice(vx, vy, vz);
    _Input ->Evaluate      (vx, vy, vz);
    Add(v, vx, vy, vz);
  }

protected:
  const BSplineFreeFormTransformationSV *_Input;  ///< Input transformation
  BaseImage                             *_Output; ///< Output image
};

// -----------------------------------------------------------------------------
/// Evaluate and add 3D B-spline SV FFD at image voxels
class AddBSplineSVFFD3D : public VoxelFunction
{
public:

  AddBSplineSVFFD3D(const BSplineFreeFormTransformationSV *input, BaseImage *output)
  :
    _Input (input),
    _Output(output),
    _y     (output->NumberOfSpatialVoxels()),
    _z     (2 * _y)
  {}

  template <class T>
  void operator()(int i, int j, int k, int, T *out)
  {
    double x = i, y = j, z = k;
    _Output->ImageToWorld  (x, y, z);
    _Input ->WorldToLattice(x, y, z);
    _Input ->Evaluate      (x, y, z);
    out[_x] += static_cast<T>(x);
    out[_y] += static_cast<T>(y);
    out[_z] += static_cast<T>(z);
  }

private:
  const BSplineFreeFormTransformationSV *_Input;  ///< Input transformation
  BaseImage                             *_Output; ///< Output image

  static const int _x = 0; ///< Offset of x component
  int              _y;     ///< Offset of y component
  int              _z;     ///< Offset of z component
};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
inline void BSplineFreeFormTransformationSV::Evaluate(double &x, double &y, double &z, double) const
{
  BSplineFreeFormTransformation3D::Evaluate(x, y, z);
}

// -----------------------------------------------------------------------------
inline void BSplineFreeFormTransformationSV::EvaluateJacobianWorld(Matrix &jac, double x, double y, double z, double) const
{
  if (_z == 1) BSplineFreeFormTransformation3D::EvaluateJacobianWorld(jac, x, y);
  else         BSplineFreeFormTransformation3D::EvaluateJacobianWorld(jac, x, y, z);
}

// -----------------------------------------------------------------------------
inline void BSplineFreeFormTransformationSV
::EvaluateJacobianDOFs(TransformationJacobian &jac, double x, double y, double z, double) const
{
  if (_z == 1) EvaluateJacobianDOFs(jac, x, y);
  else         EvaluateJacobianDOFs(jac, x, y, z);
}

// -----------------------------------------------------------------------------
inline void BSplineFreeFormTransformationSV::Velocity(GenericImage<float> &v) const
{
  if (v.T() != 3) {
    Throw(ERR_InvalidArgument, __FUNCTION__, "Output image must have 3 channels!");
  }
  v.PutTSize(0.);
  ParallelForEachVoxel(EvaluateBSplineSVFFD3D(this, &v), v.Attributes(), v);
}

// -----------------------------------------------------------------------------
inline void BSplineFreeFormTransformationSV::Velocity(GenericImage<double> &v) const
{
  if (v.T() != 3) {
    Throw(ERR_InvalidArgument, __FUNCTION__, "Output image must have 3 channels!");
  }
  v.PutTSize(0.);
  ParallelForEachVoxel(EvaluateBSplineSVFFD3D(this, &v), v.Attributes(), v);
}

// =============================================================================
// Point transformation
// =============================================================================

// -----------------------------------------------------------------------------
inline bool BSplineFreeFormTransformationSV::RequiresCachingOfDisplacements() const
{
  return (_IntegrationMethod == FFDIM_SS || _IntegrationMethod == FFDIM_FastSS);
}

// -----------------------------------------------------------------------------
inline double BSplineFreeFormTransformationSV::UpperIntegrationLimit(double t, double t0) const
{
  double T = t - t0;
  return !IsNaN(T) && !AreEqual(T, 0.) ? T : _T; // if zero/NaN, return cross-sectional interval _T instead
}

// -----------------------------------------------------------------------------
inline int BSplineFreeFormTransformationSV::NumberOfStepsForIntervalLength(double T) const
{
  return (_TimeUnit > .0) ? static_cast<int>((abs(T) / _TimeUnit) * _NumberOfSteps + 0.5) : _NumberOfSteps;
}

// -----------------------------------------------------------------------------
inline void BSplineFreeFormTransformationSV::NumberOfStepsForIntervalLength(double T, int nsteps) const
{
  if (_TimeUnit > .0) nsteps = static_cast<int>((_TimeUnit / abs(T)) * nsteps + 0.5);
  const_cast<BSplineFreeFormTransformationSV *>(this)->_NumberOfSteps = nsteps;
}

// -----------------------------------------------------------------------------
inline double BSplineFreeFormTransformationSV::StepLengthForIntervalLength(double T) const
{
  if (_NumberOfSteps == 0) return .0;
  return ((_TimeUnit > .0) ? (_TimeUnit / static_cast<double>(_NumberOfSteps))
                           : ( abs(T)  / static_cast<double>(_NumberOfSteps)));
}

// =============================================================================
// Derivatives
// =============================================================================

// -----------------------------------------------------------------------------
inline void BSplineFreeFormTransformationSV
::ParametricGradient(const GenericImage<double> *in, double *out,
                     const WorldCoordsImage *i2w, const WorldCoordsImage *wc,
                     double t0, double w) const
{
  this->ParametricGradient(in, out, i2w, wc, in->GetTOrigin(), t0, w);
}


} // namespace mirtk

#endif // MIRTK_BSplineFreeFormTransformationSV_H
