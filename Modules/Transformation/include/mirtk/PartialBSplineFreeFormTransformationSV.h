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

#ifndef MIRTK_PartialBSplineFreeFormTransformationSV_H
#define MIRTK_PartialBSplineFreeFormTransformationSV_H

#include "mirtk/Transformation.h"
#include "mirtk/BSplineFreeFormTransformationSV.h"


namespace mirtk {


/**
 * Decorator for SV FFD transformation.
 *
 * This decorator wraps a SV FFD transformation but defines it's own attributes
 * regarding integration and whether or not to invert the transformation by
 * default. As the SV FFD is parameterized by stationary velocities, inverting
 * the transformation is simply achieved by a negative upper integration limit.
 * Two instances of this decorator are in particular used by
 * SymmetricBSplineFreeFormTransformation for the two half transformations
 * to deform both input images.
 */
class PartialBSplineFreeFormTransformationSV : public Transformation
{
  mirtkTransformationMacro(PartialBSplineFreeFormTransformationSV);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Pointer to decorated transformation
  mirtkPublicAggregateMacro(BSplineFreeFormTransformationSV, Transformation);

  /// Fraction of decorated transformation, negative value corresponds to inverse
  mirtkPublicAttributeMacro(double, Fraction);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Constructor
  PartialBSplineFreeFormTransformationSV(BSplineFreeFormTransformationSV * = NULL, double = 1.0);

  /// Copy constructor
  PartialBSplineFreeFormTransformationSV(const PartialBSplineFreeFormTransformationSV &);

  /// Destructor
  virtual ~PartialBSplineFreeFormTransformationSV();

  // ---------------------------------------------------------------------------
  // Transformation parameters (DoFs)

  /// Copy active transformation parameters (DoFs) from given
  /// transformation if possible and return \c false, otherwise
  virtual bool CopyFrom(const class Transformation *);

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

  /// Checks whether transformation depends on the same vector of parameters
  virtual bool HasSameDOFsAs(const class Transformation *) const;

  /// Checks whether the transformation is an identity mapping
  virtual bool IsIdentity() const;

  // ---------------------------------------------------------------------------
  // Parameters (non-DoFs)

  // Import other overloads
  using Transformation::Parameter;

  /// Set named (non-DoF) parameter from value as string
  virtual bool Set(const char *, const char *);

  /// Get (non-DoF) parameters as key/value as string map
  virtual ParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Point transformation

private:

  /// Partial upper integration limit used given the temporal origin of both target
  /// and source images. If both images have the same temporal origin, the sign
  /// of the fraction determines the direction of integration. Otherwise, if the
  /// temporal origin of the two images differs, the signed difference between
  /// these determines the direction of integration.
  double UpperIntegrationLimit(double t, double t0) const;

public:

  /// Whether the caching of the transformation displacements is required
  /// (or preferred) by this transformation. For some transformations such as
  /// those parameterized by velocities, caching of the displacements for
  /// each target voxel results in better performance or is needed for example
  /// for the scaling and squaring method.
  virtual bool RequiresCachingOfDisplacements() const;

  /// Transforms a single point using the global transformation component only
  virtual void GlobalTransform(double &, double &, double &, double = 0, double = NaN) const;

  /// Transforms a single point using the local transformation component only
  virtual void LocalTransform(double &, double &, double &, double = 0, double = NaN) const;

  /// Transforms a single point
  virtual void Transform(double &, double &, double &, double = 0, double = NaN) const;

  /// Transforms a single point using the inverse of the global transformation only
  virtual void GlobalInverse(double &, double &, double &, double = 0, double = NaN) const;

  /// Transforms a single point using the inverse of the local transformation only
  virtual bool LocalInverse(double &, double &, double &, double = 0, double = NaN) const;

  /// Transforms a single point using the inverse of the transformation
  virtual bool Inverse(double &, double &, double &, double = 0, double = NaN) const;

  /// Calculates the displacement vectors for a whole image domain
  ///
  /// \attention The displacements are computed at the positions after applying the
  ///            current displacements at each voxel. These displacements are then
  ///            added to the current displacements. Therefore, set the input
  ///            displacements to zero if only interested in the displacements of
  ///            this transformation at the voxel positions.
  virtual void Displacement(GenericImage<double> &, double, double = NaN, const WorldCoordsImage * = NULL) const;

  /// Calculates the displacement vectors for a whole image domain
  ///
  /// \attention The displacements are computed at the positions after applying the
  ///            current displacements at each voxel. These displacements are then
  ///            added to the current displacements. Therefore, set the input
  ///            displacements to zero if only interested in the displacements of
  ///            this transformation at the voxel positions.
  virtual void Displacement(GenericImage<float> &, double, double = NaN, const WorldCoordsImage * = NULL) const;

  /// Calculates the inverse displacement vectors for a whole image domain
  ///
  /// \attention The displacements are computed at the positions after applying the
  ///            current displacements at each voxel. These displacements are then
  ///            added to the current displacements. Therefore, set the input
  ///            displacements to zero if only interested in the displacements of
  ///            this transformation at the voxel positions.
  ///
  /// \returns Always zero.
  virtual int InverseDisplacement(GenericImage<double> &, double, double = NaN, const WorldCoordsImage * = NULL) const;

  /// Calculates the inverse displacement vectors for a whole image domain
  ///
  /// \attention The displacements are computed at the positions after applying the
  ///            current displacements at each voxel. These displacements are then
  ///            added to the current displacements. Therefore, set the input
  ///            displacements to zero if only interested in the displacements of
  ///            this transformation at the voxel positions.
  ///
  /// \returns Always zero.
  virtual int InverseDisplacement(GenericImage<float> &, double, double = NaN, const WorldCoordsImage * = NULL) const;

  // ---------------------------------------------------------------------------
  // Derivatives
  using Transformation::ParametricGradient;
  using Transformation::GlobalJacobian;
  using Transformation::LocalJacobian;
  using Transformation::Jacobian;

  /// Calculates the Jacobian of the global transformation w.r.t world coordinates
  virtual void GlobalJacobian(Matrix &, double, double, double, double = 0, double = NaN) const;

  /// Calculates the Jacobian of the local transformation w.r.t world coordinates
  virtual void LocalJacobian(Matrix &, double, double, double, double = 0, double = NaN) const;

  /// Calculates the Jacobian of the transformation w.r.t world coordinates
  virtual void Jacobian(Matrix &, double, double, double, double = 0, double = NaN) const;

  /// Calculates the Hessian for each component of the global transformation w.r.t world coordinates
  virtual void GlobalHessian(Matrix [3], double, double, double, double = 0, double = NaN) const;

  /// Calculates the Hessian for each component of the local transformation w.r.t world coordinates
  virtual void LocalHessian(Matrix [3], double, double, double, double = 0, double = NaN) const;

  /// Calculates the Hessian for each component of the transformation w.r.t world coordinates
  virtual void Hessian(Matrix [3], double, double, double, double = 0, double = NaN) const;

  /// Calculates the Jacobian of the transformation w.r.t the transformation parameters
  virtual void JacobianDOFs(double [3], int, double, double, double, double = 0, double = NaN) const;

  /// Applies the chain rule to convert spatial non-parametric gradient
  /// to a gradient w.r.t the parameters of this transformation.
  virtual void ParametricGradient(const GenericImage<double> *, double *,
                                  const WorldCoordsImage *,
                                  const WorldCoordsImage *,
                                  double = NaN, double = 1) const;


  /// Applies the chain rule to convert spatial non-parametric gradient
  /// to a gradient w.r.t the parameters of this transformation.
  virtual void ParametricGradient(const GenericImage<double> **, int, double *,
                                  const WorldCoordsImage *,
                                  const WorldCoordsImage *,
                                  const double * = NULL, double = 1) const;

  // ---------------------------------------------------------------------------
  // I/O

  // Do not hide methods of base class
  using Transformation::Print;
  using Transformation::Read;
  using Transformation::Write;

  /// Prints information about the transformation
  virtual void Print(ostream &, Indent = 0) const;

  /// Reads transformation from a file stream
  virtual Cifstream &Read(Cifstream &);

  /// Writes transformation to a file stream
  virtual Cofstream &Write(Cofstream &) const;

};


} // namespace mirtk

#endif // MIRTK_PartialBSplineFreeFormTransformationSV_H
