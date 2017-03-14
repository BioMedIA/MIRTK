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

#ifndef MIRTK_StoppingCriterion_H
#define MIRTK_StoppingCriterion_H

#include "mirtk/Object.h"
#include "mirtk/ObjectiveFunction.h"


namespace mirtk {


// Forward declaration of cyclic dependency
class LocalOptimizer;


/**
 * Stopping criterion for iterative optimization
 */
class StoppingCriterion : public Object
{
  mirtkAbstractMacro(StoppingCriterion);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Optimizer to which this stopping criterion belongs to
  mirtkPublicAggregateMacro(const LocalOptimizer, Optimizer);

  /// Objective function
  mirtkPublicAggregateMacro(const ObjectiveFunction, Function);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const StoppingCriterion &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
protected:

  /// Constructor
  StoppingCriterion(const ObjectiveFunction * = NULL);

  /// Constructor
  StoppingCriterion(const LocalOptimizer *);

  /// Copy constructor
  StoppingCriterion(const StoppingCriterion &);

  /// Assignment operator
  StoppingCriterion &operator =(const StoppingCriterion &);

public:

  /// Create new copy of this instance
  virtual StoppingCriterion *New() const = 0;

  /// Destructor
  virtual ~StoppingCriterion();

  /// Initialize stopping criterion after input and parameters are set
  virtual void Initialize();

  // ---------------------------------------------------------------------------
  // Evaluation

  /// Test stopping criterion
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
  /// \returns Whether stopping criterion is fulfilled.
  virtual bool Fulfilled(int iter, double value, const double *delta) = 0;

  // ---------------------------------------------------------------------------
  // Logging

  /// Print current stopping criterion status / value
  ///
  /// This function must be called after Fulfilled, which should update any
  /// cached values that are needed by this function to avoid a costly
  /// reevaluation of the stopping criterion.
  ///
  /// \note The printed string is expected to be considerably short and must
  ///       not end with a newline or space characters.
  virtual void Print(ostream &) const;

};


} // namespace mirtk

#endif // MIRTK_StoppingCriterion_H
