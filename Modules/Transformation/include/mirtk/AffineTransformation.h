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

#ifndef MIRTK_AffineTransformation_H
#define MIRTK_AffineTransformation_H

#include "mirtk/SimilarityTransformation.h"

#include "mirtk/Matrix.h"

#include <cstdlib>
#include <iostream>


namespace mirtk {


/**
 * Class for affine transformations.
 *
 * This class defines and implements affine transformations. In addition to
 * the rigid body transformation parameters, affine transformations are
 * parameterized by three scaling and and six skewing parameters. The three
 * scaling parameters define the scaling along the axis of the coordinate
 * transformations. The six skewing parameters define the skewing angles in
 * different planes. Note that the six skewing parameters are not independent.
 */
class AffineTransformation : public SimilarityTransformation
{
  mirtkTransformationMacro(AffineTransformation);

  // Wrapper class used for inverse consistent affine registration
  // which depends on the same parameters as the affine transformation
  friend class InverseAffineTransformation;

  // ---------------------------------------------------------------------------
  // Data members

protected:

  /// Tangent of shear angle sxy
  double _tansxy;

  /// Tangent of shear angle sxz
  double _tansxz;

  /// Tangent of shear angle syz
  double _tansyz;

  /// Update cached tangents of shearing angles
  void UpdateShearingTangent();

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Default constructor with given number of parameters
  AffineTransformation(int);

  /// Copy constructor with given number of parameters
  AffineTransformation(const RigidTransformation &, int);

  /// Copy constructor with given number of parameters
  AffineTransformation(const SimilarityTransformation &, int);

  /// Copy constructor with given number of parameters
  AffineTransformation(const AffineTransformation &, int);

public:

  /// Default constructor
  AffineTransformation();

  /// Copy constructor
  AffineTransformation(const RigidTransformation &);

  /// Copy constructor
  AffineTransformation(const SimilarityTransformation &);

  /// Copy constructor
  AffineTransformation(const AffineTransformation &);

  /// Destructor
  virtual ~AffineTransformation();

  // ---------------------------------------------------------------------------
  // Approximation

  /// Approximate displacements: This function takes a set of points and a set
  /// of displacements and finds !new! parameters such that the resulting
  /// transformation approximates the displacements as good as possible.
  virtual void ApproximateDOFs(const double *, const double *, const double *, const double *,
                               const double *, const double *, const double *, int);

  // ---------------------------------------------------------------------------
  // Parameters (non-DoFs)

  using Transformation::Parameter;

  /// Set named (non-DoF) parameter from value as string
  virtual bool Set(const char *, const char *);

  /// Get (non-DoF) parameters as key/value as string map
  virtual ParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Transformation parameters

  /// Construct transformation matrix based on parameters
  static Matrix DOFs2Matrix(const double *);

  /// Extract parameters from transformation matrix
  static void Matrix2DOFs(const Matrix &, double *);

  /// Update transformation matrix after change of parameter
  virtual void UpdateMatrix();
  
  /// Update transformation parameters after change of matrix
  virtual void UpdateDOFs();

  /// Copy active transformation parameters (DoFs) from given
  /// transformation if possible and return \c false, otherwise
  virtual bool CopyFrom(const Transformation *);

  /// Puts global scaling factor
  void PutScale(double);

  /// Not supported for affine transformation!
  double GetScale() const;

  /// Puts scaling factor along the x-axis
  void PutScaleX(double);

  /// Gets scaling factor along the x-axis
  double GetScaleX() const;

  /// Puts scaling factor along the y-axis
  void PutScaleY(double);

  /// Gets scaling factor along the y-axis
  double GetScaleY() const;

  /// Puts scaling factor along the z-axis
  void PutScaleZ(double);

  /// Gets scaling factor along the z-axis
  double GetScaleZ() const;

  /// Puts y-dependent skewing angle in the x direction (in degrees)
  void PutShearXY(double);

  /// Gets y-dependent skewing angle in the x direction (in degrees)
  double GetShearXY() const;

  /// Puts z-dependent skewing angle in the y direction (in degrees)
  void PutShearYZ(double);

  /// Gets z-dependent skewing angle in the y direction (in degrees)
  double GetShearYZ() const;

  /// Puts z-dependent skewing angle in the x direction (in degrees)
  void PutShearXZ(double);

  /// Gets z-dependent skewing angle in the x direction (in degrees)
  double GetShearXZ() const;

  /// Set status of translation parameters to active/passive
  void AllowTranslations(bool);

  /// Whether all translation parameters are active
  bool AllowTranslations() const;

  /// Set status of rotation parameters to active/passive
  void AllowRotations(bool);

  /// Whether all rotation parameters are active
  bool AllowRotations() const;

  /// Set status of scaling parameters to active/passive
  void AllowScaling(bool);

  /// Whether all scaling parameters are active
  bool AllowScaling() const;

  /// Set status of shearing parameters to active/passive
  void AllowShearing(bool);

  /// Whether all shearing parameters are active
  bool AllowShearing() const;

  // ---------------------------------------------------------------------------
  // Derivatives
  using Transformation::JacobianDOFs;

  /// Calculates the Jacobian of the transformation w.r.t the parameters
  virtual void JacobianDOFs(double [3], int, double, double, double, double = 0, double = -1) const;

  /// Calculates the derivative of the Jacobian of the transformation (w.r.t. world coordinates) w.r.t. a transformation parameter
  virtual void DeriveJacobianWrtDOF(Matrix &, int, double, double, double, double = 0, double = -1) const;

  // ---------------------------------------------------------------------------
  // I/O

  // Do not hide methods of base class
  using SimilarityTransformation::Print;
  using SimilarityTransformation::Write;

  /// Prints the parameters of the transformation
  virtual void Print(ostream &, Indent = 0) const;

  /// Whether this transformation can read a file of specified type (i.e. format)
  virtual bool CanRead(TransformationType) const;

  /// Writes a transformation to a file stream
  virtual Cofstream &Write(Cofstream &) const;

protected:

  /// Reads a transformation from a file stream
  virtual Cifstream &ReadDOFs(Cifstream &, TransformationType);

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Transformation parameters
// =============================================================================

// -----------------------------------------------------------------------------
inline void AffineTransformation::PutScaleX(double sx)
{
  Put(SX, sx);
}

// -----------------------------------------------------------------------------
inline double AffineTransformation::GetScaleX() const
{
  return Get(SX);
}

// -----------------------------------------------------------------------------
inline void AffineTransformation::PutScaleY(double sy)
{
  Put(SY, sy);
}

// -----------------------------------------------------------------------------
inline double AffineTransformation::GetScaleY() const
{
  return Get(SY);
}

// -----------------------------------------------------------------------------
inline void AffineTransformation::PutScaleZ(double sz)
{
  Put(SZ, sz);
}

// -----------------------------------------------------------------------------
inline double AffineTransformation::GetScaleZ() const
{
  return Get(SZ);
}

// -----------------------------------------------------------------------------
inline void AffineTransformation::PutScale(double s)
{
  _Param[SX] = _Param[SY] = _Param[SZ] = s;
  this->Update(MATRIX);
}

// -----------------------------------------------------------------------------
inline double AffineTransformation::GetScale() const
{
  cerr << this->NameOfClass() << "::GetScale: Use GetScaleX, GetScaleY, and/or GetScaleZ!" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
inline void AffineTransformation::PutShearXY(double sxy)
{
  Put(SXY, sxy);
}

// -----------------------------------------------------------------------------
inline double AffineTransformation::GetShearXY() const
{
  return Get(SXZ);
}

// -----------------------------------------------------------------------------
inline void AffineTransformation::PutShearYZ(double syz)
{
  Put(SYZ, syz);
}

// -----------------------------------------------------------------------------
inline double AffineTransformation::GetShearYZ() const
{
  return Get(SYZ);
}

// -----------------------------------------------------------------------------
inline void AffineTransformation::PutShearXZ(double sxz)
{
  Put(SXZ, sxz);
}

// -----------------------------------------------------------------------------
inline double AffineTransformation::GetShearXZ() const
{
  return Get(SXZ);
}

// -----------------------------------------------------------------------------
inline void AffineTransformation::AllowTranslations(bool b)
{
  _Status[TX] = _Status[TY] = _Status[TZ] = (b ? Active : Passive);
}

// -----------------------------------------------------------------------------
inline bool AffineTransformation::AllowTranslations() const
{
  return _Status[TX] == Active && _Status[TY] == Active && _Status[TZ] == Active;
}

// -----------------------------------------------------------------------------
inline void AffineTransformation::AllowRotations(bool b)
{
  _Status[RX] = _Status[RY] = _Status[RZ] = (b ? Active : Passive);
}

// -----------------------------------------------------------------------------
inline bool AffineTransformation::AllowRotations() const
{
  return _Status[RX] == Active && _Status[RY] == Active && _Status[RZ] == Active;
}

// -----------------------------------------------------------------------------
inline void AffineTransformation::AllowScaling(bool b)
{
  _Status[SX] = _Status[SY] = _Status[SZ] = (b ? Active : Passive);
}

// -----------------------------------------------------------------------------
inline bool AffineTransformation::AllowScaling() const
{
  return _Status[SX] == Active && _Status[SY] == Active && _Status[SZ] == Active;
}

// -----------------------------------------------------------------------------
inline void AffineTransformation::AllowShearing(bool b)
{
  _Status[SXY] = _Status[SYZ] = _Status[SXZ] = (b ? Active : Passive);
}

// -----------------------------------------------------------------------------
inline bool AffineTransformation::AllowShearing() const
{
  return _Status[SXY] == Active && _Status[SYZ] == Active && _Status[SXZ] == Active;
}


} // namespace mirtk

#endif // MIRTK_AffineTransformation_H
