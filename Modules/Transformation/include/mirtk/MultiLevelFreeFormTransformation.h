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

#ifndef MIRTK_MultiLevelFreeFormTransformation_H
#define MIRTK_MultiLevelFreeFormTransformation_H

#include "mirtk/MultiLevelTransformation.h"


namespace mirtk {


/**
 * Class for multi-level FFD where global and local transformations are summed up.
 *
 * T_mffd(x) = T_global(x) + sum_i T_local^i(x)
 */
class MultiLevelFreeFormTransformation : public MultiLevelTransformation
{
  mirtkTransformationMacro(MultiLevelFreeFormTransformation);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Default constructor
  MultiLevelFreeFormTransformation();

  /// Construct multi-level transformation given a rigid transformation
  MultiLevelFreeFormTransformation(const RigidTransformation &);

  /// Construct multi-level transformation given an affine transformation
  MultiLevelFreeFormTransformation(const AffineTransformation &);

  /// Copy constructor
  MultiLevelFreeFormTransformation(const MultiLevelFreeFormTransformation &);

  /// Destructor
  virtual ~MultiLevelFreeFormTransformation();

  // ---------------------------------------------------------------------------
  // Levels

  /// Combine local transformations on stack
  virtual void CombineLocalTransformation();

  /// Convert the global transformation from a matrix representation to a
  /// FFD and incorporate it with any existing local transformation
  virtual void MergeGlobalIntoLocalDisplacement();

  // ---------------------------------------------------------------------------
  // Bounding box

  /// Gets the spatial bounding box for a transformation parameter in image coordinates.
  /// The last parameter specifies what fraction of the bounding box to return.
  /// The default is 1 which equals 100% of the bounding box.
  bool DOFBoundingBox(const Image *, int, int &, int &, int &,
                                          int &, int &, int &, double = 1) const;

  // ---------------------------------------------------------------------------
  // Approximation

  /// Approximate displacements: This function takes a set of points and a set
  /// of displacements and finds a transformation which approximates these
  /// displacements. After approximation, the displacements are replaced by the
  /// residual displacement errors at the points.
  ///
  /// \note This function does not change the global transformation and passive levels.
  ///       Use ApproximateAsNew to also approximate a new global transformation.
  virtual double Approximate(const ImageAttributes &, double *, double *, double *,
                             int = 1, double = .0);

  /// Approximate displacements: This function takes a set of points and a set
  /// of displacements and finds a transformation which approximates these
  /// displacements. After approximation, the displacements are replaced by the
  /// residual displacement errors at the points.
  ///
  /// \note This function does not change the global transformation and passive levels.
  ///       Use ApproximateAsNew to also approximate a new global transformation.
  virtual double Approximate(const double *, const double *, const double *,
                             double *,       double *,       double *, int,
                             int = 1, double = .0);

  /// Approximate displacements: This function takes a set of points and a set
  /// of displacements and finds a transformation which approximates these
  /// displacements. After approximation, the displacements are replaced by the
  /// residual displacement errors at the points.
  ///
  /// \note This function does not change the global transformation and passive levels.
  ///       Use ApproximateAsNew to also approximate a new global transformation.
  virtual double Approximate(const double *, const double *, const double *, const double *,
                             double *,       double *,       double *,       int,
                             int = 1, double = .0);

  /// Approximate displacements: This function takes a set of points and a set
  /// of displacements and finds a !new! transformation which approximates these
  /// displacements. After approximation, the displacements are replaced by the
  /// residual displacement errors at the points.
  ///
  /// \note This function also modifies the global transformation.
  ///       Use Reset and Approximate instead if this is not desired.
  virtual double ApproximateAsNew(const ImageAttributes &, double *, double *, double *,
                                  int = 1, double = .0);

  /// Approximate displacements: This function takes a set of points and a set
  /// of displacements and finds a !new! transformation which approximates these
  /// displacements. After approximation, the displacements are replaced by the
  /// residual displacement errors at the points.
  ///
  /// \note This function also modifies the global transformation.
  ///       Use Reset and Approximate instead if this is not desired.
  virtual double ApproximateAsNew(const double *, const double *, const double *,
                                  double *,       double *,       double *, int,
                                  int = 1, double = .0);

  /// Approximate displacements: This function takes a set of points and a set
  /// of displacements and finds a !new! transformation which approximates these
  /// displacements. After approximation, the displacements are replaced by the
  /// residual displacement errors at the points.
  ///
  /// \note This function also modifies the global transformation.
  ///       Use Reset and Approximate instead if this is not desired.
  virtual double ApproximateAsNew(const double *, const double *, const double *, const double *,
                                  double *,       double *,       double *,       int,
                                  int = 1, double = .0);

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
  // Point transformation

  // Do not hide base class methods
  using MultiLevelTransformation::LocalTransform;
  using MultiLevelTransformation::Transform;
  using MultiLevelTransformation::Displacement;
  using MultiLevelTransformation::InverseDisplacement;

  /// Transforms a single point using the local transformation component only
  virtual void LocalTransform(int, int, double &, double &, double &, double = 0, double = NaN) const;

  /// Transforms a single point
  virtual void Transform(int, int, double &, double &, double &, double = 0, double = NaN) const;

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
                                          double t, double t0 = -1,
                                          const WorldCoordsImage *i2w = NULL) const;

  /// Calculates the inverse displacement vectors for a whole image domain
  ///
  /// \attention The displacements are computed at the positions after applying the
  ///            current displacements at each voxel. These displacements are then
  ///            added to the current displacements. Therefore, set the input
  ///            displacements to zero if only interested in the displacements of
  ///            this transformation at the voxel positions.
  ///
  /// \returns Number of points for which transformation is non-invertible.
  virtual int InverseDisplacement(int, int, GenericImage<double> &, double, double = NaN, const WorldCoordsImage * = NULL) const;

  /// Calculates the inverse displacement vectors for a whole image domain
  ///
  /// \attention The displacements are computed at the positions after applying the
  ///            current displacements at each voxel. These displacements are then
  ///            added to the current displacements. Therefore, set the input
  ///            displacements to zero if only interested in the displacements of
  ///            this transformation at the voxel positions.
  ///
  /// \returns Number of points for which transformation is non-invertible.
  virtual int InverseDisplacement(int, int, GenericImage<float> &, double, double = NaN, const WorldCoordsImage * = NULL) const;

  // ---------------------------------------------------------------------------
  // Derivatives

  // Do not hide base class methods
  using MultiLevelTransformation::LocalJacobian;
  using MultiLevelTransformation::Jacobian;
  using MultiLevelTransformation::LocalHessian;
  using MultiLevelTransformation::Hessian;
  using MultiLevelTransformation::DeriveJacobianWrtDOF;

  /// Calculates the Jacobian of the transformation w.r.t world coordinates
  virtual void Jacobian(int, int, Matrix &, double, double, double, double = 0, double = NaN) const;

  /// Calculates the Hessian for each component of the transformation w.r.t world coordinates
  virtual void Hessian(int, int, Matrix [3], double, double, double, double = 0, double = NaN) const;

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
                                  const WorldCoordsImage * = NULL,
                                  const WorldCoordsImage * = NULL,
                                  double = NaN, double = 1) const;

  /// Applies the chain rule to convert point-wise non-parametric gradient
  /// to a gradient w.r.t the parameters of this transformation.
  virtual void ParametricGradient(const PointSet &, const Vector3D<double> *,
                                  double *, double = 0, double = NaN, double = 1) const;

  // ---------------------------------------------------------------------------
  // I/O

  // Do not hide methods of base class
  using MultiLevelTransformation::Print;

  /// Prints the parameters of the transformation
  virtual void Print(ostream &, Indent = 0) const;

  // ---------------------------------------------------------------------------
  // Backwards compatibility

  /// Bending of multi-level free-form deformation
  double Bending(double, double, double) const;

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Bounding box
// =============================================================================

// -----------------------------------------------------------------------------
inline bool MultiLevelFreeFormTransformation
::DOFBoundingBox(const Image *image, int dof, int &i1, int &j1, int &k1,
                                                  int &i2, int &j2, int &k2, double fraction) const
{
  const FreeFormTransformation *ffd;
  DOFIndexToLocalTransformation(this, dof, ffd, dof);
  return ffd->DOFBoundingBox(image, dof, i1, j1, k1, i2, j2, k2, fraction);
}


} // namespace mirtk

#endif // MIRTK_MultiLevelFreeFormTransformation_H
