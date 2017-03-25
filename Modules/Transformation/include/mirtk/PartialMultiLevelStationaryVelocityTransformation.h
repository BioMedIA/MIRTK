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

#ifndef MIRTK_PartialMultiLevelStationaryVelocityTransformation_H
#define MIRTK_PartialMultiLevelStationaryVelocityTransformation_H

#include "mirtk/MultiLevelTransformation.h"
#include "mirtk/MultiLevelStationaryVelocityTransformation.h"


namespace mirtk {


class MultiLevelStationaryVelocityTransformation;


/**
 * Decorator for SV MFFD transformation.
 *
 * This decorator wraps a SV MFFD transformation but defines its own attributes
 * regarding integration and whether or not to invert the transformation by
 * default. As the SV MFFD is parameterized by a sum of stationary velocities,
 * inverting the transformation is simply achieved by a negative upper
 * integration limit.
 */
class PartialMultiLevelStationaryVelocityTransformation : public MultiLevelTransformation
{
  mirtkTransformationMacro(PartialMultiLevelStationaryVelocityTransformation);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Pointer to decorated transformation
  mirtkPublicAggregateMacro(MultiLevelStationaryVelocityTransformation, Transformation);

  /// Fraction of decorated transformation, negative value corresponds to inverse
  mirtkPublicAttributeMacro(double, Fraction);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Constructor
  PartialMultiLevelStationaryVelocityTransformation(MultiLevelStationaryVelocityTransformation * = NULL, double = 1.0);

  /// Copy constructor
  PartialMultiLevelStationaryVelocityTransformation(const PartialMultiLevelStationaryVelocityTransformation &);

  /// Destructor
  virtual ~PartialMultiLevelStationaryVelocityTransformation();

  // ---------------------------------------------------------------------------
  // Levels

  /// Returns the number of levels
  int NumberOfLevels() const;

  /// Gets global transformation
  AffineTransformation *GetGlobalTransformation();

  // Get global transformation
  const AffineTransformation *GetGlobalTransformation() const;

  /// Gets local transformation
  FreeFormTransformation *GetLocalTransformation(int);

  /// Gets local transformation
  const FreeFormTransformation *GetLocalTransformation(int) const;

  /// Put local transformation and return pointer to previous one (needs to be deleted if not used)
  FreeFormTransformation *PutLocalTransformation(FreeFormTransformation *, int, bool = true);

  /// Push local transformation on stack (append transformation)
  void PushLocalTransformation(FreeFormTransformation *, bool = true);

  /// Insert local transformation
  void InsertLocalTransformation(FreeFormTransformation *, int = 0, bool = true);

  /// Pop local transformation from stack (remove last transformation)
  FreeFormTransformation *PopLocalTransformation();

  /// Remove local transformation and return the pointer (need to be deleted if not used)
  FreeFormTransformation *RemoveLocalTransformation(int = 0);

  /// Combine local transformations on stack
  virtual void CombineLocalTransformation();

  /// Convert the global transformation from a matrix representation to a
  /// FFD and incorporate it with any existing local transformation
  virtual void MergeGlobalIntoLocalDisplacement();

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
  using MultiLevelTransformation::LocalTransform;
  using MultiLevelTransformation::Transform;
  using MultiLevelTransformation::LocalInverse;
  using MultiLevelTransformation::Inverse;
  using MultiLevelTransformation::Displacement;
  using MultiLevelTransformation::InverseDisplacement;

  /// Whether the caching of the transformation displacements is required
  /// (or preferred) by this transformation. For some transformations such as
  /// those parameterized by velocities, caching of the displacements for
  /// each target voxel results in better performance or is needed for example
  /// for the scaling and squaring method.
  virtual bool RequiresCachingOfDisplacements() const;

  /// Transforms a single point using the global transformation component only
  virtual void GlobalTransform(double &, double &, double &, double = 0, double = NaN) const;

  /// Transforms a single point using the local transformation component only
  virtual void LocalTransform(int, int, double &, double &, double &, double = 0, double = NaN) const;

  /// Transforms a single point
  virtual void Transform(int, int, double &, double &, double &, double = 0, double = NaN) const;

  /// Transforms a single point using the inverse of the global transformation only
  virtual void GlobalInverse(double &, double &, double &, double = 0, double = NaN) const;

  /// Transforms a single point using the inverse of the local transformation only
  virtual bool LocalInverse(int, int, double &, double &, double &, double = 0, double = NaN) const;

  /// Transforms a single point using the inverse of the transformation
  virtual bool Inverse(int, int, double &, double &, double &, double = 0, double = NaN) const;

  /// Calculates the displacement vectors for a whole image domain
  ///
  /// \attention The displacements are computed at the positions after applying the
  ///            current displacements at each voxel. These displacements are then
  ///            added to the current displacements. Therefore, set the input
  ///            displacements to zero if only interested in the displacements of
  ///            this transformation at the voxel positions.
  virtual void Displacement(int, int, GenericImage<double> &, double, double = 1, const WorldCoordsImage * = NULL) const;

  /// Calculates the displacement vectors for a whole image domain
  ///
  /// \attention The displacements are computed at the positions after applying the
  ///            current displacements at each voxel. These displacements are then
  ///            added to the current displacements. Therefore, set the input
  ///            displacements to zero if only interested in the displacements of
  ///            this transformation at the voxel positions.
  virtual void Displacement(int, int, GenericImage<float> &, double, double = 1, const WorldCoordsImage * = NULL) const;

  /// Calculates the inverse displacement vectors for a whole image domain
  ///
  /// \attention The displacements are computed at the positions after applying the
  ///            current displacements at each voxel. These displacements are then
  ///            added to the current displacements. Therefore, set the input
  ///            displacements to zero if only interested in the displacements of
  ///            this transformation at the voxel positions.
  ///
  /// \returns Always zero.
  virtual int InverseDisplacement(int, int, GenericImage<double> &, double, double = 1, const WorldCoordsImage * = NULL) const;

  /// Calculates the inverse displacement vectors for a whole image domain
  ///
  /// \attention The displacements are computed at the positions after applying the
  ///            current displacements at each voxel. These displacements are then
  ///            added to the current displacements. Therefore, set the input
  ///            displacements to zero if only interested in the displacements of
  ///            this transformation at the voxel positions.
  ///
  /// \returns Always zero.
  virtual int InverseDisplacement(int, int, GenericImage<float> &, double, double = 1, const WorldCoordsImage * = NULL) const;

  // ---------------------------------------------------------------------------
  // Derivatives
  using MultiLevelTransformation::ParametricGradient;

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
  using MultiLevelTransformation::Print;
  using MultiLevelTransformation::Read;
  using MultiLevelTransformation::Write;

  /// Prints information about the transformation
  virtual void Print(ostream &, Indent = 0) const;

  /// Reads a transformation from a file stream
  virtual Cifstream &Read(Cifstream &);

  /// Writes a transformation to a file stream
  virtual Cofstream &Write(Cofstream &) const;

};


} // namespace mirtk

#endif // MIRTK_PartialMultiLevelStationaryVelocityTransformation_H
