/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2015 Imperial College London
 * Copyright 2015 Stefan Pszczolkowski Parraguez
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

#ifndef MIRTK_BSplineFreeFormTransformationStatistical_H
#define MIRTK_BSplineFreeFormTransformationStatistical_H

#include "mirtk/BSplineFreeFormTransformation3D.h"

#include "mirtk/Matrix.h"
#include "mirtk/Vector.h"
#include "mirtk/Vector3D.h"

#include <string>


namespace mirtk {


/**
 * Class for free-form transformations based on tensor product B-splines.
 *
 * This class implements 3D statistical free-form transformation using B-splines.
 */
class BSplineFreeFormTransformationStatistical : public BSplineFreeFormTransformation3D
{
  mirtkTransformationMacro(BSplineFreeFormTransformationStatistical);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Basis vectors (colums of the Matrix object)
  mirtkReadOnlyAttributeMacro(Matrix, BasisVectors);

  /// Mean vector
  mirtkReadOnlyAttributeMacro(mirtk::Vector, MeanVector);

  /// Name of file from which statistical deformation model was read
  mirtkAttributeMacro(string, ModelFile);

public:

  /// Default constructor
  BSplineFreeFormTransformationStatistical();

  /// Constructor based on a basis vectors matrix and a mean vector
  BSplineFreeFormTransformationStatistical(const ImageAttributes &,
                                           CPStatus ****,
                                           const Matrix &,
                                           const Vector &);

  /// Copy Constructor
  BSplineFreeFormTransformationStatistical(const BSplineFreeFormTransformationStatistical &);

  /// Destructor
  virtual ~BSplineFreeFormTransformationStatistical();

  // Import other Initialize overloads
  using BSplineFreeFormTransformation3D::Initialize;

  /// Initialize free-form transformation
  virtual void Initialize(const ImageAttributes &);

  // ---------------------------------------------------------------------------
  // Approximation/Interpolation

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

  /// Interpolates displacements: This function takes a set of displacements defined
  /// at the control points and finds a FFD which interpolates these displacements.
  virtual void Interpolate(const double *, const double *, const double * = NULL);

  // ---------------------------------------------------------------------------
  // Lattice

  using BSplineFreeFormTransformation3D::CropPadPassiveCPs;

  /// Crop/pad lattice to discard passive control points at the boundary,
  /// keeping only a layer of passive control points of given width.
  /// The DoF values of passive control points are optionally reset to zero.
  virtual bool CropPadPassiveCPs(int, int, int = 0, int = 0, bool = false);

  // ---------------------------------------------------------------------------
  // Transformation parameters (DOFs)

  /// Get norm of the gradient vector
  virtual double DOFGradientNorm(const double *) const;

  /// Puts a transformation parameter
  virtual void Put(int, DOFValue);

  /// Puts transformation parameters
  virtual void Put(const DOFValue *);

  /// Add change to transformation parameters
  virtual void Add(const DOFValue *);

  // ---------------------------------------------------------------------------
  // Update

  /// Update control point displacements after change of parameters
  void UpdateCPs();

  /// Update parameters after change of control point displacements
  void UpdateDOFs();

  // ---------------------------------------------------------------------------
  // Parameters (non-DoFs)

  using BSplineFreeFormTransformation3D::Parameter;

  /// Set named (non-DoF) parameter from value as string
  virtual bool Set(const char *, const char *);

  /// Get (non-DoF) parameters as key/value as string map
  virtual ParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Derivatives

  /// Applies the chain rule to convert spatial non-parametric gradient
  /// to a gradient w.r.t the parameters of this transformation
  virtual void ParametricGradient(const GenericImage<double> *, double *,
                                  const WorldCoordsImage *,
                                  const WorldCoordsImage *,
                                  double = 1, double = 1) const;

  // ---------------------------------------------------------------------------
  // Properties

  /// Calculates the gradient of the bending energy w.r.t the transformation parameters
  virtual void BendingEnergyGradient(double *, double = 1, bool = false, bool = true, bool = true) const;

  // ---------------------------------------------------------------------------
  // I/O

  // Do not hide methods of base class
  using BSplineFreeFormTransformation3D::Print;
  using BSplineFreeFormTransformation3D::Write;

  /// Prints the parameters of the transformation
  virtual void Print(ostream &, Indent = 0) const;

  /// Writes a transformation to a file stream
  virtual Cofstream &Write(Cofstream &) const;

  /// Reads statistical deformation model from a file
  virtual void ReadSDM(const char *);

  /// Writes statistical deformation model to a file
  virtual void WriteSDM(const char *);

  // ---------------------------------------------------------------------------
  // Other

  /// Verifies that the transformation is well constructed
  /// according to class-specific rules
  virtual void Verify();

protected:

  /// Reads transformation parameters from a file stream
  virtual Cifstream &ReadDOFs(Cifstream &, TransformationType);

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Transformation parameters
// =============================================================================

// -----------------------------------------------------------------------------
inline double BSplineFreeFormTransformationStatistical::DOFGradientNorm(const double *gradient) const
{
  return Transformation::DOFGradientNorm(gradient);
}

// -----------------------------------------------------------------------------
inline void BSplineFreeFormTransformationStatistical::Put(int i, DOFValue x)
{
  Transformation::Put(i, x);
  this->UpdateCPs();
}

// -----------------------------------------------------------------------------
inline void BSplineFreeFormTransformationStatistical::Put(const DOFValue *x)
{
  Transformation::Put(x);
  this->UpdateCPs();
}

// -----------------------------------------------------------------------------
inline void BSplineFreeFormTransformationStatistical::Add(const DOFValue *dx)
{
  Transformation::Add(dx);
  this->UpdateCPs();
}


} // namespace mirtk

#endif // MIRTK_BSplineFreeFormTransformationStatistical_H
