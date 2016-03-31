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

#ifndef MIRTK_PartialAffineTransformation_H
#define MIRTK_PartialAffineTransformation_H

#include "mirtk/AffineTransformation.h"
#include "mirtk/EventDelegate.h"


namespace mirtk {


/**
 * Class for partial linear transformation.
 *
 * An instance of this class decorates either a rigid, similarity, or an
 * affine transformation and represents only a fraction of this transformation,
 * i.e., T(x) = A^t x = exp(t * log(A)) x. A negative fraction corresponds to
 * the inverse transformation. The ParametricGradient function computes the
 * update of the parameters of the decorated transformation. Instances of this
 * class are used by the image registration filter for an inverse consistent
 * and symmetric linear registration.
 *
 * Note: For a symmetric inverse consistent linear registration, the use of
 *       AffineTransformation and InverseAffineTransformation is more
 *       efficient then two PartialAffineTransformation instances.
 *       For example, use the \c ireg energy function setting
 *       "-NMI(I1 o T^-1, I2 o T)" instead of "-NMI(I1 o T^-0.5, I2 o T^0.5)".
 *       The resulting transformation has to be squared, i.e., applied twice
 *       to obtain the full transformation between I1 and I2, however.
 *
 * \sa InverseAffineTransformation
 */
class PartialAffineTransformation : public AffineTransformation
{
  mirtkTransformationMacro(PartialAffineTransformation);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Pointer to decorated transformation
  mirtkReadOnlyAggregateMacro(HomogeneousTransformation, Transformation);

  /// Fraction of decorated transformation, negative value corresponds to inverse
  mirtkPublicAttributeMacro(double, Fraction);

  /// Observes changes of decorated transformation
  EventDelegate _TransformationObserver;

  /// Update parameters and matrix after change of decorated transformation
  void OnTransformationChanged();

public:

  /// Set decorated rigid, similarity, or affine transformation
  void Transformation(HomogeneousTransformation *);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Constructor
  PartialAffineTransformation(HomogeneousTransformation * = NULL, double = 1.0);

  /// Destructor
  virtual ~PartialAffineTransformation();

  // ---------------------------------------------------------------------------
  // Transformation parameters

  /// Checks whether transformation depends on the same vector of parameters
  virtual bool HasSameDOFsAs(const class Transformation *) const;

  // ---------------------------------------------------------------------------
  // Derivatives

  using AffineTransformation::ParametricGradient;

  /// Applies the chain rule to convert spatial non-parametric gradient
  /// to a gradient w.r.t the parameters of this transformation
  virtual void ParametricGradient(const GenericImage<double> *, double *,
                                  const WorldCoordsImage * = NULL,
                                  const WorldCoordsImage * = NULL,
                                  double = 0, double = 1) const;

};


} // namespace mirtk

#endif // MIRTK_PartialAffineTransformation_H
