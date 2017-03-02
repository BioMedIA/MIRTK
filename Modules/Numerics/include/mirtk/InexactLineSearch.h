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

#ifndef MIRTK_InexactLineSearch_H
#define MIRTK_InexactLineSearch_H

#include "mirtk/LineSearch.h"


namespace mirtk {


/**
 * Searches sufficiently optimal step length along search direction
 *
 * This local optimizer implements an inexact line search with adaptive step
 * size control, increasing the step size while steps are accepted, and
 * decreasing it when a step did not yield a sufficient improvement.
 */
class InexactLineSearch : public LineSearch
{
  mirtkAbstractMacro(InexactLineSearch);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Maximum number consecutive rejected steps
  mirtkPublicAttributeMacro(int, MaxRejectedStreak);

  /// Whether to start new search using step length of previous search
  mirtkPublicAttributeMacro(bool, ReusePreviousStepLength);

  /// Whether (reused) step length is strictly limited by [min, max]
  ///
  /// 0: Step length is allowed to exceed max which may be the case when the
  ///    previously accumulated total step length is being reused.
  /// 1: Incremental steps are strictly limited to [min, max].
  /// 2: Accumulated total step length is strictly limited to [min, max].
  /// 3: Equivalent to 2, but indicates that both incremental and accumulated
  ///    step length ranges have been restricted.
  mirtkPublicAttributeMacro(int, StrictStepLengthRange);

  /// Whether to refuse any function parameter sign changes
  ///
  /// If \c false for a particular function parameter, the line search sets
  /// the function parameter to zero whenever the sign of the parameter would
  /// change when taking a full step along the scaled gradient direction.
  mirtkPublicAggregateMacro(bool, AllowSignChange);

protected:

  /// Previous function parameters
  double *_CurrentDoFValues;

  /// Line search direction scaled by current step length
  double *_ScaledDirection;

private:

  /// Copy attributes of this class from another instance
  void CopyAttributes(const InexactLineSearch &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
public:

  using LineSearch::Function;

  /// Constructor
  InexactLineSearch(ObjectiveFunction * = NULL);

  /// Copy constructor
  InexactLineSearch(const InexactLineSearch &);

  /// Assignment operator
  InexactLineSearch &operator =(const InexactLineSearch &);

  /// Destructor
  virtual ~InexactLineSearch();

  /// Set objective function
  virtual void Function(ObjectiveFunction *);

  // ---------------------------------------------------------------------------
  // Parameters
  using LineSearch::Parameter;

  /// Set parameter value from string
  virtual bool Set(const char *, const char *);

  /// Get parameters as key/value as string map
  virtual ParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Optimization
protected:

  /// Take step in search direction
  ///
  /// \returns Maximum change of DoF
  double Advance(double);

  /// Revert step in search direction
  void Retreat(double);

  /// Evaluate objective function for a given step length
  ///
  /// This function takes a step in the given direction with the specified
  /// step length and evaluates the objective function value at this point.
  /// It then restores the function parameters by reverting the step again.
  ///
  /// \returns Objective function value at new point.
  double Value(double, double * = NULL);

};


} // namespace mirtk

#endif // MIRTK_InexactLineSearch_H
