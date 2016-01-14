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

#ifndef MIRTK_Polynomial_H
#define MIRTK_Polynomial_H

#include <mirtkObject.h>

#include <mirtkAssert.h>
#include <mirtkVector.h>
#include <mirtkMatrix.h>
#include <mirtkPointSet.h>


namespace mirtk {


/**
 * N-dimensional polynomial regression model.
 *
 * This implementation is based on the polyfitn and polyvaln MATLAB functions
 * written and distributed on MathCentral by John D'Errico.
 *
 * \sa http://uk.mathworks.com/matlabcentral/fileexchange/34765-polyfitn
 */
class Polynomial : public Object
{
  mirtkObjectMacro(Polynomial);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Order of the polynomial model, i.e., the highest exponent
  mirtkReadOnlyAttributeMacro(int, Order);

  /// Dimension of the independent variable space
  mirtkReadOnlyAttributeMacro(int, Dimension);

  /// Exponents of the model for each independent variable
  ///
  /// Each column corresponds to one variable, and each
  /// row to the exponents of each term for one variable.
  mirtkReadOnlyAttributeMacro(Matrix, ModelTerms);

  /// Coefficients of each model term
  mirtkReadOnlyAttributeMacro(Vector, Coefficients);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const Polynomial &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
public:

  /// Constructor
  ///
  /// \param[in] order Order of the polynomial, i.e., highest model term exponent.
  Polynomial(int order = 2);

  /// Construct complete polynomial order
  ///
  /// Use Coefficients(const Vector &) to set the coefficients of the model terms.
  ///
  /// \param[in] p     Dimension of independent variable space.
  /// \param[in] order Order of the polynomial, i.e., highest model term exponent.
  Polynomial(int p, int order);

  /// Copy constructor
  Polynomial(const Polynomial &);

  /// Assignment operator
  Polynomial &operator =(const Polynomial &);

  /// Destructor
  virtual ~Polynomial();

  /// Number of model terms/coefficients
  int NumberOfTerms() const;

  /// Get exponent of independent variable in specified model term
  ///
  /// \param[in] i Index of model term.
  /// \param[in] j Index of independent variable.
  ///
  /// \returns Exponent of j-th independent variable in i-th model term.
  int Exponent(int i, int j = 0) const;

  /// Whether i-th term is a constant term
  ///
  /// \param[in] i Index of model term.
  ///
  /// \returns Whether all exponents of i-th term are zero.
  bool IsConstant(int i) const;

  /// Set coefficients of polynomial model terms
  ///
  /// \param[in] coeff Coefficient vector of length NumberOfTerms().
  void Coefficients(const Vector &coeff);

  /// Set coefficient of i-th term
  ///
  /// \param[in] i Index of model term.
  /// \param[in] c Model term coefficient.
  void Coefficient(int i, double c);

  /// Get coefficient of i-th term
  ///
  /// \param[in] i Index of model term.
  ///
  /// \returns Coefficient of i-th term.
  double Coefficient(int i) const;

protected:

  /// Build complete polynomial regression model of given order
  ///
  /// \param[in] p     Dimension of independent variable space.
  /// \param[in] order Order of polynomial model, i.e., highest model term exponent.
  ///                  When non-positive, the current order of the model is used.
  ///
  /// \returns Number of model terms.
  int Initialize(int p, int order = 0);

  // ---------------------------------------------------------------------------
  // Regression
public:

  /// Fits polynomial regression model to one or more independent variables
  ///
  /// \param[in] x     Values of independent variables, a matrix of size (n x p), where
  ///                  n is the number of data points and
  ///                  p is the dimension of the independent variable space.
  /// \param[in] y     Values of dependent variable.
  /// \param[in] order Order of polynomial model, i.e., highest model term exponent.
  ///                  When non-positive, the current order of the model is used.
  ///
  /// \returns RMS error of the regression.
  double Fit(const Matrix &x, const Vector &y, int order = 0);

  /// Fits polynomial regression model to one or more independent variables
  ///
  /// \param[in] x     Values of independent variables.
  ///                  The dimension of the independent variable space is 2 or 3.
  /// \param[in] y     Values of dependent variable.
  /// \param[in] order Order of polynomial model, i.e., highest model term exponent.
  ///                  When non-positive, the current order of the model is used.
  /// \param[in] twoD  Whether to ignore z coordinate.
  ///
  /// \returns RMS error of the regression.
  double Fit(const PointSet &x, const Vector &y, int order = 0, bool twoD = false);

  /// Fits polynomial regression model to one independent variable
  ///
  /// \param[in] x     Values of independent variable.
  /// \param[in] y     Values of dependent variable.
  /// \param[in] order Order of polynomial model, i.e., highest model term exponent.
  ///                  When non-positive, the current order of the model is used.
  ///
  /// \returns RMS error of the regression.
  double Fit(const Vector &x, const Vector &y, int order = 0);

  // ---------------------------------------------------------------------------
  // Evaluation

  /// Evaluates polynomial for n p-dimensional data points
  ///
  /// \param[in] x Matrix of size (n x p), where n is the number of data points
  ///              to evaluate, and p is the dimension of the independent variable
  ///              space of the polynomial model.
  ///
  /// \returns Regressed value.
  Vector Evaluate(const Matrix &x) const;

  /// Evaluates polynomial for n 2- or 3-dimensional data points
  ///
  /// \param[in] x    Data points. The dimension of the model must be 2 or 3.
  /// \param[in] twoD Whether to ignore z coordinate.
  ///
  /// \returns Regressed value.
  Vector Evaluate(const PointSet &x, bool twoD = false) const;

  /// Evaluates polynomial for n 1-dimensional data points
  ///
  /// \param[in] x Values of independent variable.
  ///               The dimension of the model must be 1.
  ///
  /// \returns Regressed values.
  Vector Evaluate(const Vector &x) const;

  /// Evaluates polynomial for one 3D point
  ///
  /// \param[in] x Data point. The dimension of the model must be 3.
  ///
  /// \returns Regressed value.
  double Evaluate(const Point &x) const;

  /// Evaluates polynomial for one independent variable value
  ///
  /// \param[in] x Data point. The dimension of the model must be 1.
  ///
  /// \returns Regressed value.
  double Evaluate(double x) const;

  // ---------------------------------------------------------------------------
  // Debugging

  /// Print polynomial model equation
  ostream &Print(ostream &os, Indent = 0) const;

  /// Print polynomial model equation
  void Print(Indent = 0) const;

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
int Polynomial::NumberOfTerms() const
{
  return _ModelTerms.Rows();
}

// -----------------------------------------------------------------------------
int Polynomial::Exponent(int i, int j) const
{
  mirtkAssert(i < NumberOfTerms(), "model term index is within bounds");
  mirtkAssert(j < _Dimension,      "variable index is within bounds");
  return static_cast<int>(_ModelTerms(i, j));
}

// -----------------------------------------------------------------------------
bool Polynomial::IsConstant(int i) const
{
  mirtkAssert(i < NumberOfTerms(), "model term index is within bounds");
  return (_ModelTerms.RowSum(i) == .0);
}

// -----------------------------------------------------------------------------
void Polynomial::Coefficients(const Vector &coeff)
{
  mirtkAssert(coeff.Rows() == NumberOfTerms(), "coefficient vector has element for each term");
  _Coefficients = coeff;
}

// -----------------------------------------------------------------------------
void Polynomial::Coefficient(int i, double c)
{
  mirtkAssert(i < NumberOfTerms(), "model term index is within bounds");
  _Coefficients(i) = c;
}

// -----------------------------------------------------------------------------
double Polynomial::Coefficient(int i) const
{
  mirtkAssert(i < NumberOfTerms(), "model term index is within bounds");
  return _Coefficients(i);
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
Vector Polynomial::Evaluate(const PointSet &x, bool twoD) const
{
  mirtkAssert(twoD ? _Dimension == 2 : _Dimension == 3,
      "dimension of independent variable space must be " << (twoD ? 2 : 3));
  return Evaluate(Matrix(x, twoD));
}

// -----------------------------------------------------------------------------
Vector Polynomial::Evaluate(const Vector &x) const
{
  mirtkAssert(_Dimension == 1, "dimension of independent variable space must be 1");
  return Evaluate(Matrix(x));
}

// -----------------------------------------------------------------------------
double Polynomial::Evaluate(const Point &p) const
{
  mirtkAssert(_Dimension == 3, "dimension of independent variable space must be 3");
  Matrix x(1, 3);
  x(0, 0) = p._x;
  x(0, 1) = p._y;
  x(0, 2) = p._z;
  return Evaluate(x)(0);
}

// -----------------------------------------------------------------------------
double Polynomial::Evaluate(double x) const
{
  mirtkAssert(_Dimension == 1, "dimension of independent variable space must be 1");
  Matrix m(1, 1);
  m(0, 0) = x;
  return Evaluate(m)(0);
}

// =============================================================================
// Debugging
// =============================================================================

// -----------------------------------------------------------------------------
void Polynomial::Print(Indent indent) const
{
  Print(cout, indent);
}


} // namespace mirtk

#endif // MIRTK_Polynomial_H
