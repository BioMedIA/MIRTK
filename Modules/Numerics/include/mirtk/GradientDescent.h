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

#ifndef MIRTK_GradientDescent_H
#define MIRTK_GradientDescent_H

#include "mirtk/LocalOptimizer.h"
#include "mirtk/LineSearch.h"
#include "mirtk/EventDelegate.h"


namespace mirtk {



/**
 * Minimizes objective function using gradient descent
 */
class GradientDescent : public LocalOptimizer
{
  mirtkOptimizerMacro(GradientDescent, OM_GradientDescent);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Maximum number of restarts after upgrade of energy function
  mirtkPublicAttributeMacro(int, NumberOfRestarts);

  /// Maximum streak of unsuccessful restarts without improvement
  mirtkPublicAttributeMacro(int, NumberOfFailedRestarts);

  /// Line search strategy
  mirtkPublicAttributeMacro(enum LineSearchStrategy, LineSearchStrategy);

  /// (Possible) Line search parameter
  mirtkAttributeMacro(ParameterList, LineSearchParameter);

  /// Line search optimization method
  mirtkReadOnlyAggregateMacro(class LineSearch, LineSearch);

  /// Whether line search object is owned by this optimizer
  mirtkAttributeMacro(bool, LineSearchOwner);

  /// Forwards line search event messages to observers of optimization
  EventDelegate _EventDelegate;

  /// Copy attributes of this class from another instance
  void CopyAttributes(const GradientDescent &);

protected:

  /// Allocated memory for line search direction
  double *_Gradient;

  /// Whether to allow function parameter sign to change
  ///
  /// If \c false for a particular function parameter, the line search sets
  /// the function parameter to zero whenever the sign of the parameter would
  /// change when taking a full step along the scaled gradient direction.
  bool *_AllowSignChange;

  // ---------------------------------------------------------------------------
  // Construction/Destruction
public:

  /// Constructor
  GradientDescent(ObjectiveFunction * = NULL);

  /// Copy constructor
  GradientDescent(const GradientDescent &);

  /// Assignment operator
  GradientDescent &operator =(const GradientDescent &);

  /// Destructor
  virtual ~GradientDescent();

  // Import overloads from base class
  using LocalOptimizer::Function;

  /// Set objective function
  virtual void Function(ObjectiveFunction *);

  /// Set line search object
  virtual void LineSearch(class LineSearch *, bool = false);

  // ---------------------------------------------------------------------------
  // Parameters
  using LocalOptimizer::Parameter;

  /// Set parameter value from string
  virtual bool Set(const char *, const char *);

  /// Get parameters as key/value as string map
  virtual ParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Execution

  /// Initialize gradient descent
  ///
  /// This member funtion is implicitly called by Run. It can, however, be
  /// called prior to Run explicitly in order to be able to set up the line
  /// search instance. Otherwise, use the generic Set member function to
  /// change the line search parameters and simply have Run call Initialize.
  virtual void Initialize();

  /// Optimize objective function using gradient descent
  virtual double Run();

protected:

  /// Compute descent direction
  virtual void Gradient(double *, double = .0, bool * = NULL);

  /// Finalize gradient descent
  virtual void Finalize();

};


} // namespace mirtk

#endif // MIRTK_GradientDescent_H
