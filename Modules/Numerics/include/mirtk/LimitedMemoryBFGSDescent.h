/*
 * Medical Image Registration ToolKit (MIRTK) LBFGS Library
 *
 * Copyright 2013-2015 Imperial College London
 * Copyright 2013-2015 Andreas Schuh
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

#ifndef MIRTK_LimitedMemoryBFGSDescent_H
#define MIRTK_LimitedMemoryBFGSDescent_H

#include "mirtk/LocalOptimizer.h"


namespace mirtk {


/**
 * Minimizes objective function using L-BFGS
 */
class LimitedMemoryBFGSDescent : public LocalOptimizer
{
  mirtkOptimizerMacro(LimitedMemoryBFGSDescent, OM_LBFGS);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Maximum number of iterations
  mirtkPublicAttributeMacro(int, NumberOfIterations);

  /// Minimum length of steps
  mirtkPublicAttributeMacro(double, MinStepLength);

  /// Maximum length of steps
  mirtkPublicAttributeMacro(double, MaxStepLength);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const LimitedMemoryBFGSDescent &);

public:

  /// Current line search progress
  LineSearchStep _CurrentStep;

  // ---------------------------------------------------------------------------
  // Construction/Destruction
public:

  /// Constructor
  LimitedMemoryBFGSDescent(ObjectiveFunction * = NULL);

  /// Copy constructor
  LimitedMemoryBFGSDescent(const LimitedMemoryBFGSDescent &);

  /// Assignment operator
  LimitedMemoryBFGSDescent &operator =(const LimitedMemoryBFGSDescent &);

  /// Destructor
  virtual ~LimitedMemoryBFGSDescent();

  // ---------------------------------------------------------------------------
  // Parameters

  // Import other overloads
  using LocalOptimizer::Parameter;

  /// Set parameter value from string
  virtual bool Set(const char *, const char *);

  /// Get parameters as key/value as string map
  virtual ParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Execution

  /// Optimize objective function using gradient descent
  virtual double Run();

};


} // namespace mirtk

#endif // MIRTK_LimitedMemoryBFGSDescent_H
