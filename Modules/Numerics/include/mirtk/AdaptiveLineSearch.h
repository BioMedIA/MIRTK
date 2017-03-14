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

#ifndef MIRTK_AdaptiveLineSearch_H
#define MIRTK_AdaptiveLineSearch_H

#include "mirtk/InexactLineSearch.h"


namespace mirtk {


/**
 * Searches sufficiently optimal step length along search direction
 *
 * This local optimizer implements an inexact line search with adaptive step
 * size control, increasing the step size while steps are accepted, and
 * decreasing it when a step did not yield a sufficient improvement.
 */
class AdaptiveLineSearch : public InexactLineSearch
{
  mirtkLineSearchMacro(AdaptiveLineSearch, LS_Adaptive);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Increase factor of step size if accepted
  mirtkPublicAttributeMacro(double, StepLengthRise);

  /// Decrease factor of step size if rejected
  mirtkPublicAttributeMacro(double, StepLengthDrop);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const AdaptiveLineSearch &other);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
public:

  /// Constructor
  AdaptiveLineSearch(ObjectiveFunction * = NULL);

  /// Copy constructor
  AdaptiveLineSearch(const AdaptiveLineSearch &);

  /// Assignment operator
  AdaptiveLineSearch &operator =(const AdaptiveLineSearch &);

  /// Destructor
  virtual ~AdaptiveLineSearch();

  // ---------------------------------------------------------------------------
  // Parameters
  using InexactLineSearch::Parameter;

  /// Set parameter value from string
  virtual bool Set(const char *, const char *);

  /// Get parameters as key/value as string map
  virtual ParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Optimization

  /// Initialize optimization
  virtual void Initialize();

  /// Make optimal step along search direction
  virtual double Run();

};


} // namespace mirtk

#endif // MIRTK_AdaptiveLineSearch_H
