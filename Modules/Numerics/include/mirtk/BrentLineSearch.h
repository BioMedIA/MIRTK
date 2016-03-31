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

#ifndef MIRTK_BrentLineSearch_H
#define MIRTK_BrentLineSearch_H

#include "mirtk/InexactLineSearch.h"


namespace mirtk {


/**
 * Performs a line search using Brent's method
 *
 * This line search strategy uses Brent's method to find an optimal step length
 * within a given search interval. It initially brackets the minimum when the
 * user specified step length range does not need to be strictly followed.
 */
class BrentLineSearch : public InexactLineSearch
{
  mirtkLineSearchMacro(BrentLineSearch, LS_Brent);

  // ---------------------------------------------------------------------------
  // Constants

  static const double _MinificationRatio;
  static const double _MagnificationRatio;
  static const double _MaxMagnificationRatio;

  // ---------------------------------------------------------------------------
  // Attributes

  /// Tolerance value for interval
  mirtkPublicAttributeMacro(double, Tolerance);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const BrentLineSearch &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
public:

  /// Constructor
  BrentLineSearch(ObjectiveFunction * = NULL);

  /// Copy constructor
  BrentLineSearch(const BrentLineSearch &);

  /// Assignment operator
  BrentLineSearch &operator =(const BrentLineSearch &);

  /// Destructor
  virtual ~BrentLineSearch();

  // ---------------------------------------------------------------------------
  // Parameters
  using InexactLineSearch::Parameter;

  /// Set parameter value from string
  virtual bool Set(const char *, const char *);

  /// Get parameters as key/value as string map
  virtual ParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Optimization

  /// Make optimal step along search direction
  virtual double Run();

protected:

  /// Bracket extremum of objective function
  ///
  /// This function finds an initial triplet (a, b, c) of step lengths which
  /// brackets a minimum (or maximum if \p _Descent is \c false) of the
  /// objective function.
  ///
  /// \param[out] a     Minimum step length.
  /// \param[out] b     Step length for which objective function achieves an
  ///                   extremum in the step length interval [a, b].
  /// \param[out] c     Maximum step length.
  /// \param[out] delta Maximum change of objective function parameter for
  ///                   step length equal to \p b. Not returned if \c NULL.
  ///
  /// \returns Objective function value for step length \p b.
  double BracketExtremum(double &a, double &b, double &c, double *delta = NULL);

};


} // namespace mirtk

#endif // MIRTK_BrentLineSearch_H
