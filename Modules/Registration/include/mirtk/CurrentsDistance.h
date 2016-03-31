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

#ifndef MIRTK_CurrentsDistance_H
#define MIRTK_CurrentsDistance_H

#include "mirtk/PointSetDistance.h"

#include "vtkSmartPointer.h"
#include "vtkPointSet.h"
#include "vtkPolyData.h"


namespace mirtk {


/**
 * Currents distance measure of
 * points (0-currents), curves (1-currents), or surfaces (2-currents)
 *
 * This implementation is based on the currents distance similarity
 * computations of the Deformetrica registration software package.
 *
 * \sa http://www.deformetrica.org/
 */
class CurrentsDistance : public PointSetDistance
{
  mirtkEnergyTermMacro(CurrentsDistance, EM_CurrentsDistance);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Current representation of target data set
  mirtkAttributeMacro(vtkSmartPointer<vtkPolyData>, TargetCurrent);

  /// Current representation of source data set
  mirtkAttributeMacro(vtkSmartPointer<vtkPolyData>, SourceCurrent);

  /// Sigma value of currents kernel
  mirtkPublicAttributeMacro(double, Sigma);

  /// Whether to ensure symmetry of currents dot product
  mirtkPublicAttributeMacro(bool, Symmetric);

  /// Sum of squared norm of fixed (i.e., untransformed) data set(s)
  mirtkAttributeMacro(double, TargetNormSquared);

  // ---------------------------------------------------------------------------
  // Currents representation
protected:

  /// Convert data set to current
  static vtkSmartPointer<vtkPolyData> ToCurrent(vtkPointSet *);

  /// Convert surface mesh to current
  static vtkSmartPointer<vtkPolyData> SurfaceToCurrent(vtkPolyData *);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
public:

  /// Constructor
  CurrentsDistance(const char * = "", double = 1.0);

  /// Copy constructor
  CurrentsDistance(const CurrentsDistance &);

  /// Assignment operator
  CurrentsDistance &operator =(const CurrentsDistance &);

  /// Destructor
  virtual ~CurrentsDistance();

  // ---------------------------------------------------------------------------
  // Initialization

protected:

  /// Common (re-)initialization code of this class (must be non-virtual function!)
  void Init();

public:

  /// Initialize distance measure after input and parameters were set
  virtual void Initialize();

  /// Reinitialize distance measure after input topology changed
  virtual void Reinitialize();

  // ---------------------------------------------------------------------------
  // Parameters

protected:

  /// Set parameter value from string
  virtual bool SetWithPrefix(const char *, const char *);

  /// Set parameter value from string
  virtual bool SetWithoutPrefix(const char *, const char *);

public:

  // Import other overloads
  using PointSetDistance::Parameter;

  /// Get parameter key/value as string map
  virtual ParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Evaluation

  /// Update moving input points and internal state of distance measure
  virtual void Update(bool);

protected:

  /// Evaluate unweighted energy term
  virtual double Evaluate();

  /// Compute non-parametric gradient w.r.t points of given data set
  ///
  /// \param[in]  target   Transformed data set.
  /// \param[out] gradient Non-parametric gradient of polydata distance measure.
  virtual void NonParametricGradient(const RegisteredPointSet *target,
                                     GradientType                 *gradient);

  // ---------------------------------------------------------------------------
  // Debugging
public:

  /// Write input of data fidelity term
  virtual void WriteDataSets(const char *, const char *, bool = true) const;

};


} // namespace mirtk

#endif // MIRTK_CurrentsDistance_H
