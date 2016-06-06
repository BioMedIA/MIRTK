/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2016 Imperial College London
 * Copyright 2016 Andreas Schuh
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

#ifndef MIRTK_PolynomialSolvers_H
#define MIRTK_PolynomialSolvers_H

#include "mirtk/Assert.h"

// Forward declaration of Boost polynomial argument such that this
// source file does not require itself the Boost library
namespace boost { namespace math { namespace tools {
template <class T> class polynomial;
} } } // namespace boost::math::tools


namespace mirtk {


// =============================================================================
// Root finding
// =============================================================================

//------------------------------------------------------------------------------
/// Compute roots of cubic polynomial
///
/// This function solves the equation
///
///   \f$x^3 + a x^2 + b x + c = 0/f$
///
/// \param[out] x Solutions to the cubic equation.
/// \param[in]  a Quadratic coefficient.
/// \param[in]  b Linear coefficient.
/// \param[in]  c Constant coefficient.
///
/// \returns Number of real solutions.
///
/// \retval 1 x[0] contains real solution, and
///           x1[1] +/- i x[2] is pair of complex conjugated roots.
/// \retval 2 First two x array entries are the real solutions.
/// \retval 3 All three x array entries contain the solutions.
///
/// \sa https://www.easycalculation.com/algebra/learn-cubic-equation.php
int SolveCubicEquation(double x[3], double a, double b, double c);

//------------------------------------------------------------------------------
/// Compute roots of cubic polynomial
///
/// This function solves the equation
///
///   \f$a x^3 + b x^2 + c x + d = 0/f$
///
/// \param[out] x Solutions to the cubic equation.
/// \param[in]  a Cubic coefficient.
/// \param[in]  b Quadratic coefficient.
/// \param[in]  c Linear coefficient.
/// \param[in]  d Constant coefficient.
///
/// \returns Number of real solutions.
///
/// \retval 1 x[0] contains real solution, and
///           x1[1] +/- i x[2] is pair of complex conjugated roots.
/// \retval 2 First two x array entries are the real solutions.
/// \retval 3 All three x array entries contain the solutions.
///
/// \sa https://www.easycalculation.com/algebra/learn-cubic-equation.php
inline int SolveCubicEquation(double x[3], double a, double b, double c, double d)
{
  mirtkAssert(a != 0., "cubic coefficient is non-zero");
  return SolveCubicEquation(x, b / a, c / a, d / a);
}

//------------------------------------------------------------------------------
/// Compute roots of cubic polynomial
///
/// This function solves the equation
///
///   \f$a x^3 + b x^2 + c x + d = 0/f$
///
/// \param[out] x Solutions to the cubic equation.
/// \param[in]  p Cubic polynomial.
///
/// \returns Number of real solutions.
///
/// \retval 1 x[0] contains real solution, and
///           x1[1] +/- i x[2] is pair of complex conjugated roots.
/// \retval 2 First two x array entries are the real solutions.
/// \retval 3 All three x array entries contain the solutions.
///
/// \sa https://www.easycalculation.com/algebra/learn-cubic-equation.php
template <class T>
inline int SolveCubicEquation(double x[3], const boost::math::tools::polynomial<T> &p)
{
  mirtkAssert(p.degree() == 3,    "polynomial degree must be 3");
  mirtkAssert(double(p[3]) != 0., "cubic coefficient is non-zero");
  return SolveCubicEquation(x, p[3], p[2], p[1], p[0]);
}

// =============================================================================
// Minimization
// =============================================================================

//------------------------------------------------------------------------------
/// Find minimum of 4th degree polynomial
///
/// This function minimizes the polynomial function
///
///   \f$x^4 + a x^3 + b x^2 + c x + d/f$
///
/// \param[in] a Cubic coefficient.
/// \param[in] b Quadratic coefficient.
/// \param[in] c Linear coefficient.
///
/// \returns Value for which function attains its minimum value.
double MinimumOf4thDegreePolynomial(double a, double b, double c);

//------------------------------------------------------------------------------
/// Find minimum of 4th degree polynomial
///
/// This function minimizes the polynomial function
///
///   \f$x^4 + a x^3 + b x^2 + c x + d/f$
///
/// \param[in] a Cubic coefficient.
/// \param[in] b Quadratic coefficient.
/// \param[in] c Linear coefficient.
/// \param[in] d Constant coefficient.
///
/// \returns Value for which function attains its minimum value.
inline double MinimumOf4thDegreePolynomial(double a, double b, double c, double)
{
  return MinimumOf4thDegreePolynomial(a, b, c);
}

//------------------------------------------------------------------------------
/// Find minimum of 4th degree polynomial
///
/// This function minimizes the polynomial function
///
///   \f$a x^4 + b x^3 + c x^2 + d x + e/f$
///
/// \param[in] a Quartic coefficient.
/// \param[in] b Cubic coefficient.
/// \param[in] c Quadratic coefficient.
/// \param[in] d Linear coefficient.
/// \param[in] e Constant coefficient.
///
/// \returns Value for which function attains its minimum value.
inline double MinimumOf4thDegreePolynomial(double a, double b, double c, double d, double)
{
  mirtkAssert(a != 0., "quartic coefficient is non-zero");
  return MinimumOf4thDegreePolynomial(b / a, c / a, d / a);
}

//------------------------------------------------------------------------------
/// Find minimum of 4th degree polynomial
///
/// This function minimizes the polynomial function
///
///   \f$a x^4 + b x^3 + c x^2 + d x + e/f$
///
/// \param[in] p Quartic polynomial.
///
/// \returns Value for which function attains its minimum value.
template <class T>
inline double MinimumOf4thDegreePolynomial(const boost::math::tools::polynomial<T> &p)
{
  mirtkAssert(p.degree() == 4, "polynomial degree must be 4");
  return MinimumOf4thDegreePolynomial(p[4], p[3], p[2], p[1], p[0]);
}


} // namespace mirtk

#endif // MIRTK_PolynomialSolvers_H
