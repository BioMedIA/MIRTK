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

#ifndef MIRTK_SurfaceDistance_H
#define MIRTK_SurfaceDistance_H

#include "mirtk/PointSetDistance.h"

#include "mirtk/Array.h"


namespace mirtk {


/**
 * Base class for point set surface distance measures
 *
 * Subclasses of this base class are evaluating the distance of the surfaces
 * of the input point sets. Note that also other, more generic, point set
 * distance types which are subclasses of PointSetDistance may be used
 * for surface meshes as well, given that the input point sets are surfaces.
 */
class SurfaceDistance : public PointSetDistance
{
  mirtkAbstractMacro(SurfaceDistance);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

protected:

  /// Constructor
  SurfaceDistance(const char * = "", double = 1.0);

  /// Copy constructor
  SurfaceDistance(const SurfaceDistance &, int = -1, int = -1);

  /// Assignment operator
  SurfaceDistance &operator =(const SurfaceDistance &);

  /// Copy attributes of surface distance term
  void CopyAttributes(const SurfaceDistance &, int = -1, int = -1);

public:

  /// Destructor
  virtual ~SurfaceDistance();

  // ---------------------------------------------------------------------------
  // Initialization

  /// Initialize distance measure once input and parameters have been set
  virtual void Initialize();

  /// Reinitialize distance measure after change of input topology
  ///
  /// This function is called in particular when an input surface has been
  /// reparameterized, e.g., by a local remeshing filter.
  virtual void Reinitialize();

  // ---------------------------------------------------------------------------
  // Evaluation

protected:

  /// Convert non-parametric gradient of polydata distance measure into
  /// gradient w.r.t transformation parameters
  ///
  /// This function calls Transformation::ParametricGradient of the
  /// transformation to apply the chain rule in order to obtain the gradient
  /// of the distance measure w.r.t the transformation parameters.
  /// It adds the weighted gradient to the final registration energy gradient.
  ///
  /// \param[in]     target      Transformed data set.
  /// \param[in]     np_gradient Point-wise non-parametric gradient.
  /// \param[in,out] gradient    Gradient to which the computed parametric gradient
  ///                            is added, after multiplication by the given \p weight.
  /// \param[in]     weight      Weight of polydata distance measure.
  virtual void ParametricGradient(const RegisteredPointSet *target,
                                  const GradientType       *np_gradient,
                                  double                   *gradient,
                                  double                    weight);

  // ---------------------------------------------------------------------------
  // Debugging

public:

  /// Write input of data fidelity term
  virtual void WriteDataSets(const char *, const char *, bool = true) const;

  // Import public overload
  using PointSetDistance::WriteGradient;

protected:

  /// Write gradient of data fidelity term w.r.t each transformed input
  virtual void WriteGradient(const char *,
                             const RegisteredPointSet *,
                             const GradientType *,
                             const Array<int> * = NULL) const;

};


} // namespace mirtk

#endif // MIRTK_SurfaceDistance_H
