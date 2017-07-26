/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2016 Imperial College London
 * Copyright 2013-2016 Andreas Schuh
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

#ifndef MIRTK_GenericRegistrationLogger_H
#define MIRTK_GenericRegistrationLogger_H

#include "mirtk/Observer.h"


namespace mirtk {


/**
 * Prints progress of registration to output stream
 */
class GenericRegistrationLogger : public Observer
{
  mirtkObjectMacro(GenericRegistrationLogger);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Whether to always report raw, i.e., unweighted and unnormalized, energy term
  /// value in parentheses after the weighted (and normalized) energy term value.
  mirtkPublicAttributeMacro(bool, AlwaysReportRawValue);

  /// Verbosity level
  mirtkPublicAttributeMacro(int, Verbosity);

  /// Output stream for progress report
  mirtkPublicAggregateMacro(ostream, Stream);

  /// Whether to use SGR commands for colored terminal output
  mirtkPublicAttributeMacro(bool, Color);

  /// Whether to flush stream buffer after each printed message
  mirtkPublicAttributeMacro(bool, FlushBuffer);

  int _NumberOfIterations;    ///< Number of actual line search iterations
  int _NumberOfSteps;         ///< Number of iterative line search steps
  int _NumberOfGradientSteps; ///< Number of gradient descent steps

  // ---------------------------------------------------------------------------
  // Construction/Destruction
private:

  /// Copy construction
  /// \note Intentionally not implemented.
  GenericRegistrationLogger(const GenericRegistrationLogger &);

  /// Assignment operator
  /// \note Intentionally not implemented.
  GenericRegistrationLogger &operator =(const GenericRegistrationLogger &);

public:

  /// Constructor
  GenericRegistrationLogger(ostream * = &cout);

  /// Destructor
  ~GenericRegistrationLogger();

  /// Handle event and print message to output stream
  void HandleEvent(Observable *, Event, const void *);

};


} // namespace mirtk

#endif // MIRTK_GenericRegistrationLogger_H
