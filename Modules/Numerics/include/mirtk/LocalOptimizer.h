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

#ifndef MIRTK_LocalOptimizer_H
#define MIRTK_LocalOptimizer_H

#include "mirtk/Observable.h"

#include "mirtk/Math.h"
#include "mirtk/Array.h"
#include "mirtk/Queue.h"
#include "mirtk/OptimizationMethod.h"
#include "mirtk/ObjectiveFunction.h"
#include "mirtk/ObjectFactory.h"


namespace mirtk {


// Forward declaration of cyclic dependency
class StoppingCriterion;


/**
 * Minimizes objective function without guarantee of global solution
 */
class LocalOptimizer : public Observable
{
  mirtkAbstractMacro(LocalOptimizer);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Objective function
  mirtkPublicAggregateMacro(ObjectiveFunction, Function);

  /// Maximum number of iterative steps
  mirtkPublicAttributeMacro(int, NumberOfSteps);

  /// First convergence criterium: Required minimum change of objective function
  mirtkPublicAttributeMacro(double, Epsilon);

  /// Second convergence criterium: Required maximum change of DoFs
  mirtkPublicAttributeMacro(double, Delta);

  /// Last n objective function values
  mirtkReadOnlyAttributeMacro(Deque<double>, LastValues);

  /// Current slope of curve fit to last n function values
  mirtkReadOnlyAttributeMacro(double, LastValuesSlope);

  /// Number of last objective function values to store
  mirtkPublicAttributeMacro(int, NumberOfLastValues);

  /// Whether optimization converged within maximum number of steps
  mirtkReadOnlyAttributeMacro(bool, Converged);

  /// List of stopping criteria
  Array<class StoppingCriterion *> _StoppingCriteria;

  /// Copy attributes of this class from another instance
  void CopyAttributes(const LocalOptimizer &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
protected:

  /// Constructor
  LocalOptimizer(ObjectiveFunction * = NULL);

  /// Copy constructor
  LocalOptimizer(const LocalOptimizer &);

  /// Assignment operator
  LocalOptimizer &operator =(const LocalOptimizer &);

public:

  /// Type of optimizer factory
  typedef ObjectFactory<enum OptimizationMethod, LocalOptimizer> FactoryType;

  /// Get global optimizer factory instance
  static FactoryType &Factory();

  /// Construct optimizer
  static LocalOptimizer *New(enum OptimizationMethod, ObjectiveFunction * = NULL);

  /// Optimization method implemented by this optimizer
  virtual enum OptimizationMethod OptimizationMethod() const = 0;

  /// Destructor
  virtual ~LocalOptimizer();

  // ---------------------------------------------------------------------------
  // Parameters

  // Import other overloads
  using Observable::Parameter;

  /// Set parameter value from string
  virtual bool Set(const char *, const char *);

  /// Get parameters as key/value as string map
  virtual ParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Stopping criteria

  /// Get number of stopping criteria
  int NumberOfStoppingCriteria() const;

  /// Add stopping criterion and take over ownership of the object
  void AddStoppingCriterion(StoppingCriterion *);

  /// Remove stopping criterion and revoke ownership of the object
  void RemoveStoppingCriterion(StoppingCriterion *);

  /// Get the n-th stopping criterion
  class StoppingCriterion *StoppingCriterion(int);

  /// Get the n-th stopping criterion
  const class StoppingCriterion *StoppingCriterion(int) const;

  /// Delete all stopping criteria
  void ClearStoppingCriteria();

protected:

  /// Test whether a given function value is better than another
  ///
  /// \param[in] prev  Objective function value at previous iteration.
  /// \param[in] value Objective function value at current  iteration.
  bool IsImprovement(double prev, double value) const;

  /// Test stopping criteria
  ///
  /// The objective function must be up-to-date when this function is called.
  /// This is usually the case anyway because the current objective function
  /// value after a change of the parameters must be evaluated before this
  /// function is called to be able to provide this value as argument.
  ///
  /// \note The objective function value may be infinite in case of a non-parametric
  ///       deformable surface model. In this case, stopping criteria are based
  ///       only on the current surface geometry or last node displacements.
  ///       Stopping criteria based on the objective function value should
  ///       never be fulfilled in this case and always return \c false.
  ///
  /// \param[in] iter  Current number of iterations.
  /// \param[in] value Objective function value at current iteration.
  /// \param[in] delta Last change of objective function parameters.
  ///
  /// \returns True when at least one stopping criterion is fulfilled.
  virtual bool Converged(int iter, double value, const double *delta);

  // ---------------------------------------------------------------------------
  // Optimization
public:

  /// Initialize optimization
  virtual void Initialize();

  /// Run optimization
  /// \returns Value of local minimum (maximum) of objective function
  virtual double Run() = 0;

};

////////////////////////////////////////////////////////////////////////////////
// Auxiliary macros for optimizer implementation
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
#define mirtkOptimizerMacro(name, method)                                      \
  mirtkObjectMacro(name);                                                      \
public:                                                                        \
  /** Optimization method implemented by this optimizer */                     \
  static enum OptimizationMethod ID() { return method; }                       \
  /** Optimization method implemented by this optimizer */                     \
  virtual enum OptimizationMethod OptimizationMethod() const { return method; }\
private:

// -----------------------------------------------------------------------------
/// Register object type with factory singleton
#define mirtkRegisterOptimizerMacro(type)                                      \
  mirtk::LocalOptimizer::Factory()                                             \
      .Register(type::ID(), type::NameOfType(),                                \
                mirtk::New<mirtk::LocalOptimizer, type>)

// -----------------------------------------------------------------------------
#define mirtkAutoRegisterOptimizerMacro(type)                                  \
  mirtkAutoRegisterObjectTypeMacro(mirtk::LocalOptimizer::Factory(),           \
                                   mirtk::OptimizationMethod, type::ID(),      \
                                   mirtk::LocalOptimizer,     type)

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
inline bool LocalOptimizer::IsImprovement(double prev, double value) const
{
  if (IsInf(prev))  return !IsNaN(value);
  if (IsInf(value)) return false;
  return prev - value > _Epsilon;
}


} // namespace mirtk

#endif // MIRTK_LocalOptimizer_H
