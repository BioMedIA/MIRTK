/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2015 Imperial College London
 * Copyright 2008-2013 Daniel Rueckert, Julia Schnabel
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

#ifndef MIRTK_Transformation_H
#define MIRTK_Transformation_H

#include "mirtk/TransformationConfig.h" // for client code

#include "mirtk/Observable.h"

#include "mirtk/Status.h"
#include "mirtk/BSpline.h"
#include "mirtk/Cfstream.h"
#include "mirtk/Indent.h"
#include "mirtk/PointSet.h"
#include "mirtk/Matrix.h"
#include "mirtk/Vector3D.h"
#include "mirtk/GenericImage.h"

#include "mirtk/TransformationType.h"
#include "mirtk/TransformationJacobian.h"


namespace mirtk {



////////////////////////////////////////////////////////////////////////////////
// Abstract transformation class
////////////////////////////////////////////////////////////////////////////////

/**
 * Abstract base class for general transformations.
 *
 * This is the abstract base class which defines a common interface for all
 * transformations. Each derived class has to implement at least the abstract
 * methods and some of the virtual ones. Most other methods call these virtual
 * methods and should not be required to be overwritten in subclasses.
 *
 * The second time argument to the interface methods corresponds to the time
 * of the untransformed source image. It is only considered by some 3D+t
 * transformations, in particular those which parameterize the transformation
 * using a non-stationary velocity field.
 */
class Transformation : public Observable
{
  mirtkAbstractMacro(Transformation);

public:

  /// Type of transformation parameter value
  typedef double   DOFValue;

  /// Type of transforamtion parameter status
  typedef Status   DOFStatus;

protected:

  /// Number of transformation parameters
  int _NumberOfDOFs;

  /// Value of each transformation parameter
  DOFValue *_Param;

  /// Status of each transformation parameter (Active or Passive)
  DOFStatus *_Status;

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Default constructor
  Transformation(int = 0);

  /// Copy constructor
  Transformation(const Transformation &);

  /// Copy constructor
  Transformation(const Transformation &, int);

  /// Initialize transformation parameters
  void InitializeDOFs(int);

  /// Copy transformation parameters (DoFs) and their status
  void InitializeDOFs(const Transformation &, int = -1);

public:

  /// Static constructor. This function returns a pointer to a concrete
  /// new transformation of the specified type (e.g., TRANSFORMATION_RIGID).
  static Transformation *New(TransformationType);
  static Transformation *New(unsigned int type) {
    return New(static_cast<TransformationType>(type));
  }

  /// Static constructor. This function returns a pointer to a concrete
  /// transformation by copying the transformation passed to it.
  static Transformation *New(const Transformation *);

  /// Static constructor. This function returns a pointer to a concrete
  /// transformation by reading the transformation parameters from a file
  /// and creating the appropriate transformation.
  static Transformation *New(const char *);

  /// Default destructor.
  virtual ~Transformation();

  // ---------------------------------------------------------------------------
  // Transformation parameters (DoFs)

  /// Copy active transformation parameters (DoFs) from given
  /// transformation if possible and return \c false, otherwise
  virtual bool CopyFrom(const Transformation *);

  /// Get number of transformation parameters
  virtual int NumberOfDOFs() const;

  /// Get number of active transformation parameters
  int NumberOfActiveDOFs() const;

  /// Get number of passive transformation parameters
  int NumberOfPassiveDOFs() const;

  /// Get norm of the gradient vector
  virtual double DOFGradientNorm(const double *) const;

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

  /// Checks whether transformation depends on the same vector of parameters
  virtual bool HasSameDOFsAs(const Transformation *) const;

  /// Checks whether the transformation is an identity mapping
  virtual bool IsIdentity() const;

  /// Gets the spatial bounding box for a transformation parameter in image coordinates.
  /// The last parameter specifies what fraction of the bounding box to return.
  /// The default is 1 which equals 100% of the bounding box.
  virtual bool DOFBoundingBox(const Image *, int, int &, int &, int &,
                                                  int &, int &, int &, double = 1) const;

  /// Reset transformation
  virtual void Reset();

  // ---------------------------------------------------------------------------
  // Approximation

  /// Evaluates RMS error of transformation compared to another
  double EvaluateRMSError(const ImageAttributes &, const Transformation *) const;

  /// Evaluates RMS error of transformation compared to given displacement field
  ///
  /// This overloaded version of EvaluateRMSError is recommended for displacement
  /// fields defined on a regular lattice. It does not require memory for the
  /// explicit storage of the locations of each displacement vector. Moreover,
  /// if this transformation requires the caching of the displacements, the other
  /// overloads are either very slow or cannot be used.
  double EvaluateRMSError(const ImageAttributes &, double *, double *) const;

  /// Evaluates RMS error of transformation compared to given displacement field
  ///
  /// This overloaded version of EvaluateRMSError is recommended for displacement
  /// fields defined on a regular lattice. It does not require memory for the
  /// explicit storage of the locations of each displacement vector. Moreover,
  /// if this transformation requires the caching of the displacements, the other
  /// overloads are either very slow or cannot be used.
  double EvaluateRMSError(const ImageAttributes &, double *, double *, double *) const;

  /// Evaluates RMS error of transformation compared to displacement field
  double EvaluateRMSError(const double *, const double *, const double *, double,
                          double       *, double       *, double       *, int no) const;

  /// Evaluates RMS error of transformation compared to given displacement field
  double EvaluateRMSError(const double *,  const double *, const double *, const double *,
                          double       *, double        *, double       *, int no) const;

  /// Approximate another transformation and return approximation error
  virtual double Approximate(const ImageAttributes &, const Transformation *,
                             int = 1, double = .0);

  /// Approximate displacements: This function takes a set of points and a set
  /// of displacements and finds a transformation which approximates these
  /// displacements. After approximation, the displacements are replaced by the
  /// residual displacement errors at the points.
  virtual double Approximate(GenericImage<double> &, int = 1, double = .0);

  /// Approximate displacements: This function takes a set of points and a set
  /// of displacements and finds a transformation which approximates these
  /// displacements. After approximation, the displacements are replaced by the
  /// residual displacement errors at the points.
  virtual double Approximate(const ImageAttributes &, double *, double *, double *,
                             int = 1, double = .0);

  /// Approximate displacements: This function takes a set of points and a set
  /// of displacements and finds a transformation which approximates these
  /// displacements. After approximation, the displacements are replaced by the
  /// residual displacement errors at the points.
  virtual double Approximate(const double *, const double *, const double *,
                             double *,       double *,       double *, int,
                             int = 1, double = .0);

  /// Approximate displacements: This function takes a set of points and a set
  /// of displacements and finds a transformation which approximates these
  /// displacements. After approximation, the displacements are replaced by the
  /// residual displacement errors at the points.
  virtual double Approximate(const double *, const double *, const double *, const double *,
                             double *,       double *,       double *,       int,
                             int = 1, double = .0);

  /// Approximate another transformation and return approximation error
  virtual double ApproximateAsNew(const ImageAttributes &, const Transformation *,
                                  int = 1, double = .0);

  /// Approximate displacements: This function takes a set of points and a set
  /// of displacements and finds a !new! transformation which approximates these
  /// displacements. After approximation, the displacements are replaced by the
  /// residual displacement errors at the points.
  virtual double ApproximateAsNew(GenericImage<double> &, int = 1, double = .0);

  /// Approximate displacements: This function takes a set of points and a set
  /// of displacements and finds a !new! transformation which approximates these
  /// displacements. After approximation, the displacements are replaced by the
  /// residual displacement errors at the points.
  virtual double ApproximateAsNew(const ImageAttributes &, double *, double *, double *,
                                  int = 1, double = .0);

  /// Approximate displacements: This function takes a set of points and a set
  /// of displacements and finds a !new! transformation which approximates these
  /// displacements. After approximation, the displacements are replaced by the
  /// residual displacement errors at the points.
  virtual double ApproximateAsNew(const double *, const double *, const double *,
                                  double *,       double *,       double *, int,
                                  int = 1, double = .0);

  /// Approximate displacements: This function takes a set of points and a set
  /// of displacements and finds a !new! transformation which approximates these
  /// displacements. After approximation, the displacements are replaced by the
  /// residual displacement errors at the points.
  virtual double ApproximateAsNew(const double *, const double *, const double *, const double *,
                                  double *,       double *,       double *,       int,
                                  int = 1, double = .0);

  /// Finds gradient of approximation error: This function takes a set of points
  /// and a set of errors and finds a gradient to minimize the L2 norm of the error.
  virtual void ApproximateGradient(const ImageAttributes &,
                                   const double *, const double *, const double *,
                                   double *, double = 1.0) const;

  /// Finds gradient of approximation error: This function takes a set of points
  /// and a set of errors and finds a gradient to minimize the L2 norm of the error.
  virtual void ApproximateGradient(const double *, const double *, const double *,
                                   const double *, const double *, const double *, int,
                                   double *, double = 1.0) const;

  /// Finds gradient of approximation error: This function takes a set of points
  /// and a set of errors and finds a gradient to minimize the L2 norm of the error.
  virtual void ApproximateGradient(const double *, const double *, const double *, const double *,
                                   const double *, const double *, const double *, int,
                                   double *, double = 1.0) const;

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

  // ---------------------------------------------------------------------------
  // Parameters (non-DoFs)

  // Import other overloads
  using Observable::Parameter;

  /// Set named (non-DoF) parameter from value as string
  virtual bool Set(const char *, const char *);

  /// Get (non-DoF) parameters as key/value as string map
  virtual ParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Point transformation

  /// Whether the caching of the transformation displacements is required
  /// (or preferred) by this transformation. For some transformations such as
  /// those parameterized by velocities, caching of the displacements for
  /// each target voxel results in better performance or is needed for example
  /// for the scaling and squaring method.
  virtual bool RequiresCachingOfDisplacements() const;

  /// Transforms a single point using the global transformation component only
  virtual void GlobalTransform(double &, double &, double &, double = 0, double = NaN) const = 0;

  /// Transforms a single point using the local transformation component only
  virtual void LocalTransform(double &, double &, double &, double = 0, double = NaN) const = 0;

  /// Transforms a single point
  virtual void Transform(double &, double &, double &, double = 0, double = NaN) const = 0;

  /// Transforms a single point
  virtual void Transform(Point &, double = 0, double = NaN) const;

  /// Transforms a set of points
  virtual void Transform(PointSet &, double = 0, double = NaN) const;

  /// Transforms a set of points
  virtual void Transform(int, double *, double *, double *, double = 0, double = NaN) const;

  /// Transforms a set of points
  virtual void Transform(int, double *, double *, double *, const double *, double = NaN) const;

  /// Transforms world coordinates of image voxels
  virtual void Transform(WorldCoordsImage &, double = NaN) const;

  /// Calculates the displacement of a single point using the global transformation component only
  virtual void GlobalDisplacement(double &, double &, double &, double = 0, double = NaN) const;

  /// Calculates the displacement of a single point using the local transformation component only
  virtual void LocalDisplacement(double &, double &, double &, double = 0, double = NaN) const;

  /// Calculates the displacement of a single point
  virtual void Displacement(double &, double &, double &, double = 0, double = NaN) const;

  /// Calculates the displacement at specified lattice points
  ///
  /// \attention The displacements are computed at the positions after applying the
  ///            current displacements at each point. These displacements are then
  ///            added to the current displacements. Therefore, set the input
  ///            displacements to zero if only interested in the displacements of
  ///            this transformation at the lattice points.
  virtual void Displacement(const ImageAttributes &, double *, double *, double *) const;

  /// Calculates the displacement vectors for a whole image domain
  virtual void Displacement(GenericImage<double> &, double = NaN, const WorldCoordsImage * = NULL) const;

  /// Calculates the displacement vectors for a whole image domain
  virtual void Displacement(GenericImage<float> &, double = NaN, const WorldCoordsImage * = NULL) const;

  /// Calculates the displacement vectors for a whole image domain
  ///
  /// \attention The displacements are computed at the positions after applying the
  ///            current displacements at each voxel. These displacements are then
  ///            added to the current displacements. Therefore, set the input
  ///            displacements to zero if only interested in the displacements of
  ///            this transformation at the voxel positions.
  virtual void Displacement(GenericImage<double> &, double, double, const WorldCoordsImage * = NULL) const;

  /// Calculates the displacement vectors for a whole image domain
  ///
  /// \attention The displacements are computed at the positions after applying the
  ///            current displacements at each voxel. These displacements are then
  ///            added to the current displacements. Therefore, set the input
  ///            displacements to zero if only interested in the displacements of
  ///            this transformation at the voxel positions.
  virtual void Displacement(GenericImage<float> &, double, double, const WorldCoordsImage * = NULL) const;

  /// Whether this transformation implements a more efficient update of a given
  /// displacement field given the desired change of a transformation parameter
  virtual bool CanModifyDisplacement(int = -1) const;

  /// Updates the displacement vectors for a whole image domain
  ///
  /// \param[in]     dof Transformation parameter.
  /// \param[in]     dv  Change of transformation parameter value.
  /// \param[in,out] dx  Displacement field to be updated.
  /// \param[in]     t   Time point of start point.
  /// \param[in]     t0  Time point of end point.
  /// \param[in]     i2w Pre-computed world coordinates.
  virtual void DisplacementAfterDOFChange(int dof, double dv,
                                          GenericImage<double> &dx,
                                          double t, double t0 = NaN,
                                          const WorldCoordsImage *i2w = NULL) const;

  /// Transforms a single point using the inverse of the global transformation only
  virtual void GlobalInverse(double &, double &, double &, double = 0, double = NaN) const;

  /// Transforms a single point using the inverse of the local transformation only
  virtual bool LocalInverse(double &, double &, double &, double = 0, double = NaN) const;

  /// Transforms a single point using the inverse of the transformation
  virtual bool Inverse(double &, double &, double &, double = 0, double = NaN) const;

  /// Transforms a single point using the inverse of the transformation
  virtual bool Inverse(Point &, double = 0, double = NaN) const;

  /// Transforms a set of points using the inverse of the transformation
  ///
  /// \returns Number of points at which transformation is non-invertible.
  virtual int Inverse(PointSet &, double = 0, double = NaN) const;

  /// Calculates the displacement of a single point using the inverse of the global transformation only
  virtual void GlobalInverseDisplacement(double &, double &, double &, double = 0, double = NaN) const;

  /// Calculates the displacement of a single point using the inverse of the local transformation only
  virtual bool LocalInverseDisplacement(double &, double &, double &, double = 0, double = NaN) const;

  /// Calculates the displacement of a single point using the inverse of the transformation
  virtual bool InverseDisplacement(double &, double &, double &, double = 0, double = NaN) const;

  /// Calculates the inverse displacement at specified lattice points
  ///
  /// \attention The displacements are computed at the positions after applying the
  ///            current displacements at each point. These displacements are then
  ///            added to the current displacements. Therefore, set the input
  ///            displacements to zero if only interested in the displacements of
  ///            this transformation at the lattice points.
  ///
  /// \returns Number of points at which transformation is non-invertible.
  virtual int InverseDisplacement(const ImageAttributes &, double *, double *, double *) const;

  /// Calculates the inverse displacement vectors for a whole image domain
  ///
  /// \returns Number of points at which transformation is non-invertible.
  virtual int InverseDisplacement(GenericImage<double> &, double = NaN, const WorldCoordsImage * = NULL) const;

  /// Calculates the inverse displacement vectors for a whole image domain
  ///
  /// \returns Number of points at which transformation is non-invertible.
  virtual int InverseDisplacement(GenericImage<float> &, double = NaN, const WorldCoordsImage * = NULL) const;

  /// Calculates the inverse displacement vectors for a whole image domain
  ///
  /// \attention The displacements are computed at the positions after applying the
  ///            current displacements at each voxel. These displacements are then
  ///            added to the current displacements. Therefore, set the input
  ///            displacements to zero if only interested in the displacements of
  ///            this transformation at the voxel positions.
  ///
  /// \returns Number of points at which transformation is non-invertible.
  virtual int InverseDisplacement(GenericImage<double> &, double, double, const WorldCoordsImage * = NULL) const;

  /// Calculates the inverse displacement vectors for a whole image domain
  ///
  /// \attention The displacements are computed at the positions after applying the
  ///            current displacements at each voxel. These displacements are then
  ///            added to the current displacements. Therefore, set the input
  ///            displacements to zero if only interested in the displacements of
  ///            this transformation at the voxel positions.
  ///
  /// \returns Number of points at which transformation is non-invertible.
  virtual int InverseDisplacement(GenericImage<float> &, double, double, const WorldCoordsImage * = NULL) const;

  // ---------------------------------------------------------------------------
  // Derivatives

  /// Calculates the Jacobian of the global transformation w.r.t world coordinates
  virtual void GlobalJacobian(Matrix &, double, double, double, double = 0, double = NaN) const;

  /// Calculates the Jacobian of the local transformation w.r.t world coordinates
  virtual void LocalJacobian(Matrix &, double, double, double, double = 0, double = NaN) const;

  /// Calculates the Jacobian of the transformation w.r.t world coordinates
  virtual void Jacobian(Matrix &, double, double, double, double = 0, double = NaN) const;

  /// Calculates the determinant of the Jacobian of the global transformation w.r.t world coordinates
  virtual double GlobalJacobian(double, double, double, double = 0, double = NaN) const;

  /// Calculates the determinant of the Jacobian of the local transformation w.r.t world coordinates
  virtual double LocalJacobian(double, double, double, double = 0, double = NaN) const;

  /// Calculates the determinant of the Jacobian of the transformation w.r.t world coordinates
  virtual double Jacobian(double, double, double, double = 0, double = NaN) const;

  /// Calculates the Hessian for each component of the global transformation w.r.t world coordinates
  virtual void GlobalHessian(Matrix [3], double, double, double, double = 0, double = NaN) const;

  /// Calculates the Hessian for each component of the local transformation w.r.t world coordinates
  virtual void LocalHessian(Matrix [3], double, double, double, double = 0, double = NaN) const;

  /// Calculates the Hessian for each component of the transformation w.r.t world coordinates
  virtual void Hessian(Matrix [3], double, double, double, double = 0, double = NaN) const;

  /// Calculates the Jacobian of the transformation w.r.t a transformation parameter
  virtual void JacobianDOFs(double [3], int, double, double, double, double = 0, double = NaN) const;

  /// Calculates the derivative of the Jacobian of the transformation (w.r.t. world coordinates) w.r.t. a transformation parameter
  virtual void DeriveJacobianWrtDOF(Matrix &, int, double, double, double, double = 0, double = NaN) const;

  /// Applies the chain rule to convert spatial non-parametric gradient
  /// to a gradient w.r.t the parameters of this transformation.
  ///
  /// If the transformation itself is non-parametric, the gradient will be passed through
  /// unchanged. The default implementation uses the full Jacobian matrix computed for each
  /// DoF separately (i.e., calls JacobianDOFs for each DoF).
  ///
  /// For 4D transformations, the temporal coordinate t used for the computation of the Jacobian
  /// of the transformation w.r.t the transformation parameters for all spatial (x, y, z) voxel
  /// coordinates in world units, is assumed to correspond to the temporal origin of the given
  /// gradient image. For 4D transformations parameterized by velocities, a second time for the
  /// upper integration bound can be provided as last argument to this method. This last
  /// argument is ignored by transformations parameterized by displacements.
  ///
  /// \sa ImageSimilarityMetric::EvaluateGradient
  virtual void ParametricGradient(const GenericImage<double> *, double *,
                                  const WorldCoordsImage *,
                                  const WorldCoordsImage *,
                                  double = NaN, double = 1) const;

  /// Applies the chain rule to convert spatial non-parametric gradient
  /// to a gradient w.r.t the parameters of this transformation.
  void ParametricGradient(const GenericImage<double> *, double *,
                          const WorldCoordsImage *,
                          double = NaN, double = 1) const;

  /// Applies the chain rule to convert spatial non-parametric gradient
  /// to a gradient w.r.t the parameters of this transformation.
  void ParametricGradient(const GenericImage<double> *, double *,
                          double = NaN, double = 1) const;

  /// Applies the chain rule to convert spatial non-parametric gradient
  /// to a gradient w.r.t the parameters of this transformation.
  virtual void ParametricGradient(const GenericImage<double> **, int, double *,
                                  const WorldCoordsImage *,
                                  const WorldCoordsImage *,
                                  const double * = NULL, double = 1) const;

  /// Applies the chain rule to convert spatial non-parametric gradient
  /// to a gradient w.r.t the parameters of this transformation.
  void ParametricGradient(const GenericImage<double> **, int, double *,
                          const WorldCoordsImage *,
                          const double * = NULL, double = 1) const;

  /// Applies the chain rule to convert spatial non-parametric gradient
  /// to a gradient w.r.t the parameters of this transformation.
  void ParametricGradient(const GenericImage<double> **, int, double *,
                          const double * = NULL, double = 1) const;

  /// Applies the chain rule to convert point-wise non-parametric gradient
  /// to a gradient w.r.t the parameters of this transformation.
  virtual void ParametricGradient(const PointSet &, const Vector3D<double> *,
                                  double *, double = 0, double = NaN, double = 1) const;

  // ---------------------------------------------------------------------------
  // I/O

  /// Prints information about the transformation
  void Print(Indent = 0) const;

  /// Prints information about the transformation
  virtual void Print(ostream &os, Indent = 0) const = 0;

  /// Reads a transformation from a file
  virtual void Read(const char *);

  /// Writes a transformation to a file
  virtual void Write(const char *) const;

  /// Reads a transformation from a file stream
  virtual Cifstream &Read(Cifstream &);

  /// Writes a transformation to a file stream
  virtual Cofstream &Write(Cofstream &) const;

  /// Whether magic number in header of given file
  static bool CheckHeader(const char *);

  /// Whether this transformation can read a file of specified type (i.e. format)
  virtual bool CanRead(TransformationType) const;

protected:

  /// Reads transformation parameters from a file stream
  virtual Cifstream &ReadDOFs(Cifstream &, TransformationType);

  /// Writes transformation parameters to a file stream
  virtual Cofstream &WriteDOFs(Cofstream &) const;

public:

  // ---------------------------------------------------------------------------
  // Others

  /// Returns type ID corresponding to transformations of the named class
  static TransformationType TypeOfClass(const char *);

  /// Returns type ID of the instantiated transformation class
  virtual TransformationType TypeOfClass() const;

  /// Verifies that the transformation is well constructed
  /// according to class-specific rules
  virtual void Verify();

};

////////////////////////////////////////////////////////////////////////////////
// Auxiliary macros for transformation implementation
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
#define mirtkAbstractTransformationMacro(name)                                  \
  mirtkAbstractMacro(name)

// -----------------------------------------------------------------------------
#define mirtkTransformationMacro(name)                                          \
  mirtkObjectMacro(name)

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Transformation parameters
// =============================================================================

// -----------------------------------------------------------------------------
inline int Transformation::NumberOfDOFs() const
{
  return _NumberOfDOFs;
}

// -----------------------------------------------------------------------------
inline int Transformation::NumberOfActiveDOFs() const
{
  int nactive = 0;
  for (int dof = 0; dof < this->NumberOfDOFs(); ++dof) {
    if (this->GetStatus(dof) == Active) ++nactive;
  }
  return nactive;
}

// -----------------------------------------------------------------------------
inline int Transformation::NumberOfPassiveDOFs() const
{
  return this->NumberOfDOFs() - this->NumberOfActiveDOFs();
}

// -----------------------------------------------------------------------------
inline double Transformation::DOFGradientNorm(const double *gradient) const
{
  double norm, max = .0;
  for (int dof = 0; dof < _NumberOfDOFs; ++dof) {
    norm = abs(gradient[dof]);
    if (norm > max) max = norm;
  }
  return max;
}

// -----------------------------------------------------------------------------
inline void Transformation::Put(int idx, double x)
{
  if (_Param[idx] != static_cast<DOFValue>(x)) {
    _Param[idx] = static_cast<DOFValue>(x);
    this->Changed(true);
  }
}

// -----------------------------------------------------------------------------
inline double Transformation::Get(int idx) const
{
  return static_cast<double>(_Param[idx]);
}

// -----------------------------------------------------------------------------
inline void Transformation::Put(const DOFValue *x)
{
  for (int idx = 0; idx < _NumberOfDOFs; ++idx) {
    if (_Param[idx] != x[idx]) {
      _Param[idx] = x[idx];
      this->Changed(true);
    }
  }
}

// -----------------------------------------------------------------------------
inline void Transformation::Add(const DOFValue *dx)
{
  for (int idx = 0; idx < _NumberOfDOFs; ++idx) {
    if (dx[idx] != .0) {
      _Param[idx] += dx[idx];
      this->Changed(true);
    }
  }
}

// -----------------------------------------------------------------------------
inline double Transformation::Update(const DOFValue *dx)
{
  double delta, max_delta = .0;
  for (int idx = 0; idx < _NumberOfDOFs; ++idx) {
    _Param[idx] += dx[idx];
    delta = abs(dx[idx]);
    if (delta > max_delta) max_delta = delta;
  }
  if (max_delta > .0) this->Changed(true);
  return max_delta;
}

// -----------------------------------------------------------------------------
inline void Transformation::Get(DOFValue *x) const
{
  memcpy(x, _Param, _NumberOfDOFs * sizeof(DOFValue));
}

// -----------------------------------------------------------------------------
inline void Transformation::PutStatus(int idx, DOFStatus s)
{
  _Status[idx] = s;
}

// -----------------------------------------------------------------------------
inline Transformation::DOFStatus Transformation::GetStatus(int idx) const
{
  return _Status[idx];
}

// -----------------------------------------------------------------------------
inline bool Transformation::HasSameDOFsAs(const Transformation *t) const
{
  return (_Param == t->_Param);
}

// -----------------------------------------------------------------------------
inline bool HaveSameDOFs(const Transformation *t1, const Transformation *t2)
{
  // If the virtual HasSameDOFsAs function was not overriden by the specialized
  // type of one transformation, it is assumed that it uses the _Param memory to
  // store its transformation parameters. However, the other transformation
  // might be a specialized type which does not make use of it. In this case,
  // even if the two _Param pointers do not reference the same memory, the two
  // transformations may still use the same parameters because the other
  // transformation directly or indirectly wraps this transformation. Therefore,
  // the check with the transformations exchanged and the boolean OR (not AND).
  return t1->HasSameDOFsAs(t2) || t2->HasSameDOFsAs(t1);
}

// -----------------------------------------------------------------------------
inline bool Transformation::IsIdentity() const
{
  for (int i = 0; i < this->NumberOfDOFs(); ++i) {
    if (this->Get(i) != .0) return false;
  }
  return true;
}

// -----------------------------------------------------------------------------
inline bool Transformation
::DOFBoundingBox(const Image *image, int, int &i1, int &j1, int &k1,
                                          int &i2, int &j2, int &k2, double fraction) const
{
  i1 = j1 = k1 = 0;
  i2 = image->X() - 1, j2 = image->Y() -1, k2 = image->Z() - 1;
  return !image->IsEmpty();
}

// =============================================================================
// Point transformation
// =============================================================================

// -----------------------------------------------------------------------------
inline bool Transformation::RequiresCachingOfDisplacements() const
{
  return false;
}

// -----------------------------------------------------------------------------
inline void Transformation::Transform(Point &p, double t, double t0) const
{
  this->Transform(p._x, p._y, p._z, t, t0);
}

// -----------------------------------------------------------------------------
inline void Transformation::Transform(PointSet &pset, double t, double t0) const
{
  for (int i = 0; i < pset.Size(); i++) this->Transform(pset(i), t, t0);
}

// -----------------------------------------------------------------------------
inline bool Transformation::Inverse(Point &p, double t, double t0) const
{
  return this->Inverse(p._x, p._y, p._z, t, t0);
}

// -----------------------------------------------------------------------------
inline int Transformation::Inverse(PointSet &pset, double t, double t0) const
{
  int n = 0;
  for (int i = 0; i < pset.Size(); ++i) {
    if (!this->Inverse(pset(i), t, t0)) ++n;
  }
  return n;
}

// -----------------------------------------------------------------------------
inline void Transformation::GlobalDisplacement(double &x, double &y, double &z, double t, double t0) const
{
  const double u = x;
  const double v = y;
  const double w = z;

  this->GlobalTransform(x, y, z, t, t0);

  x -= u;
  y -= v;
  z -= w;
}

// -----------------------------------------------------------------------------
inline void Transformation::LocalDisplacement(double &x, double &y, double &z, double t, double t0) const
{
  const double u = x;
  const double v = y;
  const double w = z;

  this->LocalTransform(x, y, z, t, t0);

  x -= u;
  y -= v;
  z -= w;
}

// -----------------------------------------------------------------------------
inline void Transformation::Displacement(double &x, double &y, double &z, double t, double t0) const
{
  const double u = x;
  const double v = y;
  const double w = z;

  this->Transform(x, y, z, t, t0);

  x -= u;
  y -= v;
  z -= w;
}

// -----------------------------------------------------------------------------
inline void Transformation::GlobalInverseDisplacement(double &x, double &y, double &z, double t, double t0) const
{
  const double u = x;
  const double v = y;
  const double w = z;

  this->GlobalInverse(x, y, z, t, t0);

  x -= u;
  y -= v;
  z -= w;
}

// -----------------------------------------------------------------------------
inline bool Transformation::LocalInverseDisplacement(double &x, double &y, double &z, double t, double t0) const
{
  const double u = x;
  const double v = y;
  const double w = z;

  bool ok = this->LocalInverse(x, y, z, t, t0);

  x -= u;
  y -= v;
  z -= w;

  return ok;
}

// -----------------------------------------------------------------------------
inline bool Transformation::InverseDisplacement(double &x, double &y, double &z, double t, double t0) const
{
  const double u = x;
  const double v = y;
  const double w = z;

  bool ok = this->Inverse(x, y, z, t, t0);

  x -= u;
  y -= v;
  z -= w;

  return ok;
}

// =============================================================================
// Derivatives
// =============================================================================

// -----------------------------------------------------------------------------
inline void Transformation::GlobalJacobian(Matrix &, double, double, double, double, double) const
{
  cerr << this->NameOfClass() << "::GlobalJacobian: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
inline void Transformation::LocalJacobian(Matrix &, double, double, double, double, double) const
{
  cerr << this->NameOfClass() << "::LocalJacobian: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
inline void Transformation::Jacobian(Matrix &, double, double, double, double, double) const
{
  cerr << this->NameOfClass() << "::Jacobian: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
inline void Transformation::DeriveJacobianWrtDOF(Matrix &, int, double, double, double, double, double) const
{
  cerr << this->NameOfClass() << "::DeriveJacobianWrtDOF: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
inline double Transformation::GlobalJacobian(double x, double y, double z, double t, double t0) const
{
  Matrix jac(3, 3);
  this->GlobalJacobian(jac, x, y, z, t, t0);
  return jac.Det3x3();
}

// -----------------------------------------------------------------------------
inline double Transformation::LocalJacobian(double x, double y, double z, double t, double t0) const
{
  Matrix jac(3, 3);
  this->LocalJacobian(jac, x, y, z, t, t0);
  return jac.Det3x3();
}

// -----------------------------------------------------------------------------
inline double Transformation::Jacobian(double x, double y, double z, double t, double t0) const
{
  Matrix jac(3, 3);
  this->Jacobian(jac, x, y, z, t, t0);
  return jac.Det3x3();
}

// -----------------------------------------------------------------------------
inline void Transformation::GlobalHessian(Matrix [3], double, double, double, double, double) const
{
  cerr << this->NameOfClass() << "::GlobalHessian: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
inline void Transformation::LocalHessian(Matrix [3], double, double, double, double, double) const
{
  cerr << this->NameOfClass() << "::LocalHessian: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
inline void Transformation::Hessian(Matrix [3], double, double, double, double, double) const
{
  cerr << this->NameOfClass() << "::Hessian: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
inline void Transformation::JacobianDOFs(double [3], int, double, double, double, double, double) const
{
  cerr << this->NameOfClass() << "::JacobianDOFs: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
inline void Transformation
::ParametricGradient(const GenericImage<double> *in, double *out,
                     const WorldCoordsImage *i2w, double t0, double w) const
{
  this->ParametricGradient(in, out, i2w, NULL, t0, w);
}

// -----------------------------------------------------------------------------
inline void Transformation
::ParametricGradient(const GenericImage<double> *in, double *out, double t0, double w) const
{
  this->ParametricGradient(in, out, NULL, NULL, t0, w);
}

// -----------------------------------------------------------------------------
inline void Transformation
::ParametricGradient(const GenericImage<double> **in, int n, double *out,
                     const WorldCoordsImage *i2w, const WorldCoordsImage *wc,
                     const double *t0, double w) const
{
  for (int i = 0; i < n; ++i) {
    this->ParametricGradient(in[i], out, i2w, wc, t0 ? t0[i] : 1.0, w);
  }
}

// -----------------------------------------------------------------------------
inline void Transformation
::ParametricGradient(const GenericImage<double> **in, int n, double *out,
                     const WorldCoordsImage *i2w, const double *t0, double w) const
{
  this->ParametricGradient(in, n, out, i2w, NULL, t0, w);
}

// -----------------------------------------------------------------------------
inline void Transformation
::ParametricGradient(const GenericImage<double> **in, int n, double *out, const double *t0, double w) const
{
  this->ParametricGradient(in, n, out, NULL, NULL, t0, w);
}

// =============================================================================
// Others
// =============================================================================

// -----------------------------------------------------------------------------
inline void Transformation::Print(Indent indent) const
{
  this->Print(cout, indent);
}

// -----------------------------------------------------------------------------
inline TransformationType Transformation::TypeOfClass() const
{
  return Transformation::TypeOfClass(this->NameOfClass());
}

////////////////////////////////////////////////////////////////////////////////
// Auxiliary functions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
/// Check whether a named file is an MIRTK transformation file
///
/// \param[in] name File name.
///
/// \returns Whether the named file exists and stores an MIRTK transformation.
bool IsTransformation(const char *name);


} // namespace mirtk

#endif // MIRTK_Transformation_H
