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

#ifndef MIRTK_LineSearch_H
#define MIRTK_LineSearch_H

#include "mirtk/LocalOptimizer.h"

#include "mirtk/Math.h"


namespace mirtk {


/// Enumeration of available line search strategies
enum LineSearchStrategy
{
  LS_None,     ///< No line search
  // Add new enumeration values below
  LS_Adaptive, ///< Inexact line search with adaptive step length
  LS_Brent,    ///< Numerical recipes linmin function using Brent's method
  LS_LinMin = LS_Brent, ///< Alias for LS_Brent
  // Add new enumeration values above
  LS_Last
};


/**
 * Finds optimal step length along given search direction
 */
class LineSearch : public LocalOptimizer
{
  mirtkAbstractMacro(LineSearch);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Direction of line search
  mirtkPublicAggregateMacro(double, Direction);

  /// Whether to perform line search in opposite direction
  mirtkPublicAttributeMacro(bool, Revert);

  /// Initial objective function value
  mirtkPublicAttributeMacro(double, CurrentValue);

  /// Maximum number of search iterations
  mirtkPublicAttributeMacro(int, NumberOfIterations);

  /// Minimum length of step
  mirtkPublicAttributeMacro(double, MinStepLength);

  /// Maximum length of step
  mirtkPublicAttributeMacro(double, MaxStepLength);

  /// Unit of step length, e.g., maximum gradient norm
  mirtkPublicAttributeMacro(double, StepLengthUnit);

  /// Initial/final/optimal step length
  mirtkPublicAttributeMacro(double, StepLength);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const LineSearch &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
protected:

  /// Constructor
  LineSearch(ObjectiveFunction * = NULL);

  /// Copy constructor
  LineSearch(const LineSearch &);

  /// Assignment operator
  LineSearch &operator =(const LineSearch &);

public:

  /// Instantiate line search implementing specified strategy
  static LineSearch *New(LineSearchStrategy &, ObjectiveFunction * = NULL);

  /// Destructor
  virtual ~LineSearch();

  /// Line search strategy implemented by this line search
  virtual LineSearchStrategy Strategy() const = 0;

  // ---------------------------------------------------------------------------
  // Parameters

  // Import other overloads
  using LocalOptimizer::Parameter;

  /// Set parameter value from string
  virtual bool Set(const char *, const char *);

  /// Get parameters as key/value as string map
  virtual ParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Optimization

  /// Initialize optimization
  virtual void Initialize();

  /// Make optimal step along search direction
  /// \returns New value of objective function or previous if no step successful
  virtual double Run() = 0;
};

////////////////////////////////////////////////////////////////////////////////
// Auxiliary macros for line search implementation
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
#define mirtkLineSearchMacro(name, strategy)                                   \
  mirtkObjectMacro(name);                                                      \
  public:                                                                      \
    /** Optimization method implemented by this optimizer */                   \
    virtual mirtk::OptimizationMethod OptimizationMethod() const               \
    {                                                                          \
      return OM_LineSearch;                                                    \
    }                                                                          \
    /** Line search strategy implemented by this line search */                \
    virtual mirtk::LineSearchStrategy Strategy() const { return strategy; }    \
  private:

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
template <>
inline string ToString(const LineSearchStrategy &m, int w, char c, bool left)
{
  const char *str;
  switch (m) {
    case LS_None:     str = "None"; break;
    case LS_Adaptive: str = "Adaptive"; break;
    case LS_Brent:    str = "Brent"; break;
    default:          str = "Unknown"; break;
  }
  return ToString(str, w, c, left);
}

// -----------------------------------------------------------------------------
template <>
inline bool FromString(const char *str, LineSearchStrategy &m)
{
  m = LS_None;
  if (strcmp(str, "None") == 0) return true;
  if (strcmp(str, "linmin") == 0 || strcmp(str, "LinMin") == 0) m = LS_Brent;
  if (m == LS_None) {
    m = static_cast<LineSearchStrategy>(LS_Last - 1);
    while (m != LS_None) {
      if (ToString(m) == str) break;
      m = static_cast<LineSearchStrategy>(m - 1);
    }
  }
  return m != LS_None;
}

// -----------------------------------------------------------------------------
// For use by subclass Run implementation to get current function value
inline double LineSearch::Run()
{
  if (IsNaN(_CurrentValue)) {
    Function()->Update(false);
    _CurrentValue = Function()->Value();
  }
  return _CurrentValue;
}


} // namespace mirtk

#endif // MIRTK_LineSearch_H
