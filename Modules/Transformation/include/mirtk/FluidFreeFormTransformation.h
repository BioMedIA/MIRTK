/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2017 Imperial College London
 * Copyright 2008-2013 Daniel Rueckert, Julia Schnabel
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

#ifndef MIRTK_FluidFreeFormTransformation_H
#define MIRTK_FluidFreeFormTransformation_H

#include "mirtk/MultiLevelTransformation.h"

#include "mirtk/Memory.h" // swap


namespace mirtk {


/**
 * Class for multi-level FFD where global and local transformations are composed.
 *
 * T_mffd(x) = T_affine ° T_local^N ° T_local^N-1 ° ... ° T_local^1 ° T_global(x)
 *
 * where T_global is the initial global affine transformation, T_local are
 * the non-rigid free-form deformations (FFD), and T_affine is an optional
 * affine transformation which is applied after the deformations.
 */
class FluidFreeFormTransformation : public MultiLevelTransformation
{
  mirtkTransformationMacro(FluidFreeFormTransformation);

  // ---------------------------------------------------------------------------
  // Attributes

protected:

  /// Affine transformation applied after the local free-form deformations.
  AffineTransformation _AffineTransformation;

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Default constructor
  FluidFreeFormTransformation();

  /// Construct multi-level transformation given a rigid transformation
  FluidFreeFormTransformation(const RigidTransformation &);

  /// Construct multi-level transformation given an affine transformation
  FluidFreeFormTransformation(const AffineTransformation &);

  /// Copy constructor
  FluidFreeFormTransformation(const FluidFreeFormTransformation &);

  /// Destructor
  virtual ~FluidFreeFormTransformation();

  // ---------------------------------------------------------------------------
  // Levels

  /// Compose current transformation with another
  virtual void PushTransformation(const Transformation *, ImageAttributes *attr = nullptr);

  /// Combine local transformations on stack
  virtual void CombineLocalTransformation();

  /// Convert the global transformation from a matrix representation to a
  /// FFD and incorporate it with any existing local transformation
  virtual void MergeGlobalIntoLocalDisplacement();

  /// Get affine transformation applied after the local free-form deformations
  AffineTransformation *GetAffineTransformation();

  /// Get affine transformation applied after the local free-form deformations
  const AffineTransformation *GetAffineTransformation() const;

protected:

  /// Checks whether a given transformation is supported as local transformation
  virtual void CheckTransformation(FreeFormTransformation *) const;

public:

  // ---------------------------------------------------------------------------
  // Bounding box

  /// Gets the spatial bounding box for a transformation parameter in image coordinates.
  /// The last parameter specifies what fraction of the bounding box to return.
  /// The default is 1 which equals 100% of the bounding box.
  virtual bool DOFBoundingBox(const Image *, int, int &, int &, int &,
                                                  int &, int &, int &, double = 1) const;

  // ---------------------------------------------------------------------------
  // Approximation

  /// Approximate displacements: This function takes a set of points and a set
  /// of displacements and finds a transformation which approximates these
  /// displacements. After approximation, the displacements are replaced by the
  /// residual displacement errors at the points.
  ///
  /// \note This function does not change the global transformation and the
  ///       initial passive local levels of this multi-level FFD. Use
  ///       ApproximateAsNew to also approximate a new global transformation.
  virtual double Approximate(const ImageAttributes &, double *, double *, double *,
                             int = 1, double = .0);

  /// Approximate displacements: This function takes a set of points and a set
  /// of displacements and finds a transformation which approximates these
  /// displacements. After approximation, the displacements are replaced by the
  /// residual displacement errors at the points.
  ///
  /// \note This function does not change the global transformation and the
  ///       initial passive local levels of this multi-level FFD. Use
  ///       ApproximateAsNew to also approximate a new global transformation.
  virtual double Approximate(const double *, const double *, const double *,
                             double *,       double *,       double *, int,
                             int = 1, double = .0);

  /// Approximate displacements: This function takes a set of points and a set
  /// of displacements and finds a transformation which approximates these
  /// displacements. After approximation, the displacements are replaced by the
  /// residual displacement errors at the points.
  ///
  /// \note This function does not change the global transformation and the
  ///       initial passive local levels of this multi-level FFD. Use
  ///       ApproximateAsNew to also approximate a new global transformation.
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
  // Transformation parameters (DoFs)

  /// Copy active transformation parameters (DoFs) from given
  /// transformation if possible and return \c false, otherwise
  virtual bool CopyFrom(const Transformation *);

  /// Checks whether transformation is an identity mapping
  virtual bool IsIdentity() const;

  /// Reset transformation (does not remove local transformations)
  virtual void Reset();

  /// Reset transformation and remove all local transformations
  virtual void Clear();

  // ---------------------------------------------------------------------------
  // Point transformation

  using MultiLevelTransformation::Transform;
  using MultiLevelTransformation::Inverse;
  using MultiLevelTransformation::Displacement;
  using MultiLevelTransformation::InverseDisplacement;

  /// Transforms a single point
  virtual void Transform(int, int, double &, double &, double &, double = 0, double = -1) const;

  /// Transforms a single point using the inverse of the transformation
  virtual bool Inverse(int, int, double &, double &, double &, double = 0, double = -1) const;

  /// Calculates the displacement vectors for a whole image domain
  ///
  /// \attention The displacements are computed at the positions after applying the
  ///            current displacements at each voxel. These displacements are then
  ///            added to the current displacements. Therefore, set the input
  ///            displacements to zero if only interested in the displacements of
  ///            this transformation at the voxel positions.
  virtual void Displacement(int, int, GenericImage<double> &, double, double = -1, const WorldCoordsImage * = NULL) const;

  /// Calculates the displacement vectors for a whole image domain
  ///
  /// \attention The displacements are computed at the positions after applying the
  ///            current displacements at each voxel. These displacements are then
  ///            added to the current displacements. Therefore, set the input
  ///            displacements to zero if only interested in the displacements of
  ///            this transformation at the voxel positions.
  virtual void Displacement(int, int, GenericImage<float> &, double, double = -1, const WorldCoordsImage * = NULL) const;

  /// Calculates the inverse displacement vectors for a whole image domain
  ///
  /// \attention The displacements are computed at the positions after applying the
  ///            current displacements at each voxel. These displacements are then
  ///            added to the current displacements. Therefore, set the input
  ///            displacements to zero if only interested in the displacements of
  ///            this transformation at the voxel positions.
  virtual int InverseDisplacement(int, int, GenericImage<double> &, double, double = -1, const WorldCoordsImage * = NULL) const;

  /// Calculates the inverse displacement vectors for a whole image domain
  ///
  /// \attention The displacements are computed at the positions after applying the
  ///            current displacements at each voxel. These displacements are then
  ///            added to the current displacements. Therefore, set the input
  ///            displacements to zero if only interested in the displacements of
  ///            this transformation at the voxel positions.
  virtual int InverseDisplacement(int, int, GenericImage<float> &, double, double = -1, const WorldCoordsImage * = NULL) const;

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
                                  double = -1, double = 1) const;

  /// Applies the chain rule to convert point-wise non-parametric gradient
  /// to a gradient w.r.t the parameters of this transformation.
  virtual void ParametricGradient(const PointSet &, const Vector3D<double> *,
                                  double *, double = 0, double = -1, double = 1) const;

  // ---------------------------------------------------------------------------
  // Derivatives

  using MultiLevelTransformation::Jacobian;
  using MultiLevelTransformation::Hessian;

  /// Calculates the Jacobian of the transformation w.r.t world coordinates
  virtual void Jacobian(int, int, Matrix &, double, double, double, double = 0, double = -1) const;

  /// Calculates the Hessian of the transformation w.r.t world coordinates
  virtual void Hessian(int, int, Matrix [3], double, double, double, double = 0, double = -1) const;

  // ---------------------------------------------------------------------------
  // I/O

  // Do not hide methods of base class
  using MultiLevelTransformation::Print;

  /// Prints the parameters of the transformation
  virtual void Print(ostream &, Indent = 0) const;

  /// Whether this transformation can read a file of specified type (i.e. format)
  virtual bool CanRead(TransformationType) const;

protected:

  /// Reads transformation parameters from a file stream
  virtual Cifstream &ReadDOFs(Cifstream &, TransformationType);

  /// Writes transformation parameters to a file stream
  virtual Cofstream &WriteDOFs(Cofstream &) const;

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Levels
// =============================================================================

// -----------------------------------------------------------------------------
inline AffineTransformation *FluidFreeFormTransformation::GetAffineTransformation()
{
  return &_AffineTransformation;
}

// -----------------------------------------------------------------------------
inline  const AffineTransformation *FluidFreeFormTransformation::GetAffineTransformation() const
{
  return &_AffineTransformation;
}

// =============================================================================
// Bounding box
// =============================================================================

// -----------------------------------------------------------------------------
inline bool FluidFreeFormTransformation
::DOFBoundingBox(const Image *image, int dof, int &i1, int &j1, int &k1,
                                              int &i2, int &j2, int &k2, double fraction) const
{
  // Get local transformation and corresponding parameter index
  const FreeFormTransformation *ffd;
  const int n = DOFIndexToLocalTransformation(this, dof, ffd, dof);

  // Calculate bounding box in world coordinates
  double x1, y1, z1, x2, y2, z2;
  ffd->BoundingBox(ffd->DOFToIndex(dof), x1, y1, z1, x2, y2, z2, fraction);

  // Apply inverse transformations which preceed this local transformation
  this->Inverse(-1, n, x1, y1, z1);
  this->Inverse(-1, n, x2, y2, z2);

  // Convert to image coordinates
  image->WorldToImage(x1, y1, z1);
  image->WorldToImage(x2, y2, z2);

  if (x2 < x1) swap(x1, x2);
  if (y2 < y1) swap(y1, y2);
  if (z2 < z1) swap(z1, z2);

  // Round up/down to nearest voxel in image domain
  i1 = static_cast<int>(ceil (x1));
  i2 = static_cast<int>(floor(x2));
  j1 = static_cast<int>(ceil (y1));
  j2 = static_cast<int>(floor(y2));
  k1 = static_cast<int>(ceil (z1));
  k2 = static_cast<int>(floor(z2));
  // When both indices are outside in opposite directions,
  // use the full range [0, N[. If they are both outside in
  // the same direction, the condition i1 <= i2 is false which
  // indicates that the bounding box is empty in this case
  i1 = (i1 < 0 ?  0 : (i1 >= image->X() ? image->X()     : i1));
  i2 = (i2 < 0 ? -1 : (i2 >= image->X() ? image->X() - 1 : i2));
  j1 = (j1 < 0 ?  0 : (j1 >= image->Y() ? image->Y()     : j1));
  j2 = (j2 < 0 ? -1 : (j2 >= image->Y() ? image->Y() - 1 : j2));
  k1 = (k1 < 0 ?  0 : (k1 >= image->Z() ? image->Z()     : k1));
  k2 = (k2 < 0 ? -1 : (k2 >= image->Z() ? image->Z() - 1 : k2));
  return i1 <= i2 && j1 <= j2 && k1 <= k2;
}


} // namespace mirtk

#endif // MIRTK_FluidFreeFormTransformation_H
