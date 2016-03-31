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

#ifndef MIRTK_PointSetDistance_H
#define MIRTK_PointSetDistance_H

#include "mirtk/DataFidelity.h"

#include "mirtk/Array.h"
#include "mirtk/Vector3D.h"
#include "mirtk/PointSetDistanceMeasure.h"
#include "mirtk/RegisteredPointSet.h"


namespace mirtk {


/**
 * Base class for point set distance measures
 */
class PointSetDistance : public DataFidelity
{
  mirtkAbstractMacro(PointSetDistance);

  // ---------------------------------------------------------------------------
  // Types
public:

  /// Type of gradient w.r.t a single transformed data point
  typedef Vector3D<double> GradientType;

  // ---------------------------------------------------------------------------
  // Attributes
protected:

  /// First point set
  mirtkPublicAggregateMacro(RegisteredPointSet, Target);

  /// Second point set
  mirtkPublicAggregateMacro(RegisteredPointSet, Source);

  /// Memory for (non-parametric) gradient w.r.t points of target
  mirtkComponentMacro(GradientType, GradientWrtTarget);

  /// Memory for (non-parametric) gradient w.r.t points of source
  mirtkComponentMacro(GradientType, GradientWrtSource);

  /// Whether Update has not been called since initialization
  mirtkAttributeMacro(bool, InitialUpdate);

  /// Allocate memory for (non-parametric) gradient w.r.t points of target
  void AllocateGradientWrtTarget(int);

  /// Allocate memory for (non-parametric) gradient w.r.t points of source
  void AllocateGradientWrtSource(int);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
protected:

  /// Constructor
  PointSetDistance(const char * = "", double = 1.0);

  /// Copy constructor
  PointSetDistance(const PointSetDistance &, int = -1, int = -1);

  /// Assignment operator
  PointSetDistance &operator =(const PointSetDistance &);

  /// Copy attributes from other point set distance measure
  void CopyAttributes(const PointSetDistance &, int = -1, int = -1);

public:

  /// Instantiate specified similarity measure
  static PointSetDistance *New(PointSetDistanceMeasure,
                               const char * = "", double = 1.0);

  /// Destructor
  virtual ~PointSetDistance();

  // ---------------------------------------------------------------------------
  // Initialization

protected:

  /// Initialize distance measure once input and parameters have been set
  void Initialize(int, int);

  /// Reinitialize distance measure after change of input topology
  ///
  /// This function is called in particular when an input surface has been
  /// reparameterized, e.g., by a local remeshing filter.
  void Reinitialize(int, int);

public:

  /// Initialize distance measure once input and parameters have been set
  virtual void Initialize();

  /// Reinitialize distance measure after change of input topology
  ///
  /// This function is called in particular when an input surface has been
  /// reparameterized, e.g., by a local remeshing filter.
  virtual void Reinitialize();

  // ---------------------------------------------------------------------------
  // Evaluation

  /// Update moving input points and internal state of distance measure
  virtual void Update(bool = true);

protected:

  /// Compute non-parametric gradient w.r.t points of given data set
  ///
  /// \param[in]  target   Transformed point set.
  /// \param[out] gradient Non-parametric gradient of point set distance measure.
  virtual void NonParametricGradient(const RegisteredPointSet *target,
                                     GradientType             *gradient) = 0;

  /// Convert non-parametric gradient of point set distance measure into
  /// gradient w.r.t transformation parameters
  ///
  /// This function calls Transformation::ParametricGradient of the
  /// transformation to apply the chain rule in order to obtain the gradient
  /// of the distance measure w.r.t the transformation parameters.
  /// It adds the weighted gradient to the final registration energy gradient.
  ///
  /// \param[in]     target      Transformed point set.
  /// \param[in]     np_gradient Point-wise non-parametric gradient.
  /// \param[in,out] gradient    Gradient to which the computed parametric gradient
  ///                            is added, after multiplication by the given \p weight.
  /// \param[in]     weight      Weight of point set distance measure.
  virtual void ParametricGradient(const RegisteredPointSet *target,
                                  const GradientType       *np_gradient,
                                  double                   *gradient,
                                  double                    weight);

  /// Evaluate gradient of point set distance measure
  ///
  /// This function calls the virtual NonParametricGradient function to be
  /// implemented by subclasses for each transformed input point set to obtain
  /// the gradient of the point set distance measure. It then converts this gradient
  /// into a gradient w.r.t the transformation parameters using the ParametricGradient.
  ///
  /// If both target and source point sets are transformed by different transformations,
  /// the resulting gradient vector contains first the derivative values w.r.t
  /// the parameters of the target transformation followed by those computed
  /// w.r.t the parameters of the source transformation. If both point sets are
  /// transformed by the same transformation, the sum of the derivative values
  /// is added to the resulting gradient vector. This is in particular the case
  /// for a velocity based transformation model which is applied to deform both
  /// point sets "mid-way". Otherwise, only one input point set is transformed
  /// (usually the target) and the derivative values of only the respective
  /// transformation parameters added to the gradient vector.
  ///
  /// \sa NonParametricGradient, ParametricGradient
  ///
  /// \param[in,out] gradient Gradient to which the computed gradient of the
  ///                         point set distance measure is added after
  ///                         multiplying by the given similarity \p weight.
  /// \param[in]     step     Step length for finite differences (unused).
  /// \param[in]     weight   Weight of point set distance measure.
  virtual void EvaluateGradient(double *gradient, double step, double weight);

  // ---------------------------------------------------------------------------
  // Debugging
public:

  /// Write input of data fidelity term
  virtual void WriteDataSets(const char *, const char *, bool = true) const;

  /// Write gradient of data fidelity term w.r.t each transformed input
  virtual void WriteGradient(const char *, const char *) const;

protected:

  /// Write gradient of data fidelity term w.r.t each transformed input
  virtual void WriteGradient(const char *,
                             const RegisteredPointSet *,
                             const GradientType *,
                             const Array<int> * = NULL) const;

};


} // namespace mirtk

#endif // MIRTK_PointSetDistance_H
