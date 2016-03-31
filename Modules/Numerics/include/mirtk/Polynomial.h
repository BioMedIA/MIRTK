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

#include "mirtk/Object.h"

#include "mirtk/Assert.h"
#include "mirtk/Array.h"
#include "mirtk/Status.h"
#include "mirtk/Vector.h"
#include "mirtk/Vector3.h"
#include "mirtk/Matrix.h"
#include "mirtk/Matrix3x3.h"
#include "mirtk/PointSet.h"
#include "mirtk/OrderedSet.h"


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

  /// Status of coefficients, only active coefficients are free
  mirtkAttributeMacro(Array<enum Status>, Status);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const Polynomial &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
public:

  /// Constructor
  ///
  /// \param[in] order Order of the polynomial, i.e., highest model term exponent.
  Polynomial(int order = 2);

  /// Construct complete polynomial
  ///
  /// \param[in] p     Dimension of independent variable space.
  /// \param[in] order Order of the polynomial, i.e., highest model term exponent.
  /// \param[in] coeff Coefficients of each model term. If not specified, the
  ///                  coefficients are initialized to zero.
  ///                  Use Coefficients(const Vector &) in this case to set the
  ///                  coefficients of the model terms after model construction.
  Polynomial(int p, int order, const Vector &coeff = Vector());

  /// Copy constructor
  Polynomial(const Polynomial &);

  /// Assignment operator
  Polynomial &operator =(const Polynomial &);

  /// Destructor
  virtual ~Polynomial();

  /// Build complete polynomial regression model of given order
  ///
  /// \param[in] p     Dimension of independent variable space.
  /// \param[in] order Order of polynomial model, i.e., highest model term exponent.
  ///                  When non-positive, the current order of the model is used.
  ///
  /// \returns Number of model terms.
  int Initialize(int p, int order = 0);

  /// Number of model terms/coefficients
  int NumberOfTerms() const;

  /// Number of passive model terms/coefficients
  int NumberOfPassiveTerms() const;

  /// Number of active model terms/coefficients
  int NumberOfActiveTerms() const;

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

  /// Set status of i-th coefficient
  ///
  /// \param[in] i Index of coefficient.
  /// \param[in] s Status of coefficient.
  void Status(int i, enum Status s);

  /// Get status of i-th coefficient
  ///
  /// \param[in] i Index of coefficient.
  ///
  /// \returns Status of i-th coefficient.
  enum Status Status(int i) const;

  /// Set all model coefficients to zero
  void SetCoefficientsToZero();

  /// Sets coefficient(s) of constant model terms
  ///
  /// \param[in] value  Coefficient value.
  /// \param[in] status Status of constant model terms. Set to Passive by
  ///                   default, i.e., constant terms are fixed for model fitting.
  void SetConstantCoefficient(double value, enum Status status = Passive);

  // ---------------------------------------------------------------------------
  // Regression
public:

  /// Fits polynomial regression model to one or more independent variables
  ///
  /// \param[in] x Values of independent variables, a matrix of size (n x p), where
  ///              n is the number of data points and
  ///              p is the dimension of the independent variable space.
  /// \param[in] y Values of dependent variable.
  ///
  /// \returns RMS error of the regression.
  double Fit(const Matrix &x, const Vector &y);

  /// Fits polynomial regression model to one or more independent variables
  ///
  /// \param[in] x      Values of independent variables, a matrix of size (n x p), where
  ///                   n is the number of data points and
  ///                   p is the dimension of the independent variable space.
  /// \param[in] subset Indices of points to use for fitting.
  /// \param[in] y      Values of dependent variable.
  ///
  /// \returns RMS error of the regression.
  double Fit(const Matrix &x, const Vector &y, const Array<int> &subset);

  /// Fits polynomial regression model to one or more independent variables
  ///
  /// \param[in] x      Values of independent variables, a matrix of size (n x p), where
  ///                   n is the number of data points and
  ///                   p is the dimension of the independent variable space.
  /// \param[in] subset Indices of points to use for fitting.
  /// \param[in] y      Values of dependent variable.
  ///
  /// \returns RMS error of the regression.
  double Fit(const Matrix &x, const Vector &y, const OrderedSet<int> &subset);

  /// Fits polynomial regression model to one or more independent variables
  ///
  /// \param[in] x    Values of independent variables.
  ///                 The dimension of the independent variable space is 2 or 3.
  /// \param[in] y    Values of dependent variable.
  /// \param[in] twoD Whether to ignore z coordinate.
  ///
  /// \returns RMS error of the regression.
  double Fit(const PointSet &x, const Vector &y, bool twoD = false);

  /// Fits polynomial regression model to one independent variable
  ///
  /// \param[in] x Values of independent variable.
  /// \param[in] y Values of dependent variable.
  ///
  /// \returns RMS error of the regression.
  double Fit(const Vector &x, const Vector &y);

  /// Fits polynomial surface model to 3D point cloud
  ///
  /// The dependent variable of the surface model is the algebraic distance of a
  /// point from the surface. This function sets the coefficients of constant model
  /// term(s) to 1 and excludes them from the fit to avoid the trivial solution
  /// of the homogeneous linear system of equations.
  ///
  /// \param[in] x Points on surface.
  ///
  /// \returns RMS error of the regression.
  double FitSurface(const PointSet &x);

  /// Fits polynomial surface model to subset of 3D point cloud
  ///
  /// The dependent variable of the surface model is the algebraic distance of a
  /// point from the surface. This function sets the coefficients of constant model
  /// term(s) to 1 and excludes them from the fit to avoid the trivial solution
  /// of the homogeneous linear system of equations.
  ///
  /// \param[in] x      Points on surface.
  /// \param[in] subset Indices of points to use for fitting.
  ///
  /// \returns RMS error of the regression.
  double FitSurface(const PointSet &x, const Array<int> &subset);

  /// Fits polynomial surface model to subset of 3D point cloud
  ///
  /// The dependent variable of the surface model is the algebraic distance of a
  /// point from the surface. This function sets the coefficients of constant model
  /// term(s) to 1 and excludes them from the fit to avoid the trivial solution
  /// of the homogeneous linear system of equations.
  ///
  /// \param[in] x      Points on surface.
  /// \param[in] subset Indices of points to use for fitting.
  ///
  /// \returns RMS error of the regression.
  double FitSurface(const PointSet &x, const OrderedSet<int> &subset);

  /// Fits polynomial surface model to 3D point cloud
  ///
  /// The dependent variable of the surface model is the algebraic distance of a
  /// point from the surface.
  ///
  /// \param[in] x Points on surface.
  /// \param[in] n Surface normals at input points.
  /// \param[in] c Offset along normal direction for points inside and outside
  ///              the surface to be modelled.
  ///
  /// \returns RMS error of the regression.
  double FitSurface(const PointSet &x, const PointSet &n, double c = 1.0);

  /// Fits polynomial surface model to subset of 3D point cloud
  ///
  /// The dependent variable of the surface model is the distance of a point from
  /// the surface.
  ///
  /// \param[in] x      Points on surface.
  /// \param[in] n      Surface normals at input points.
  /// \param[in] subset Indices of points to use for fitting.
  /// \param[in] c      Offset along normal direction for points inside and outside
  ///                   the surface to be modelled.
  ///
  /// \returns RMS error of the regression.
  double FitSurface(const PointSet &x, const PointSet &n,
                    const Array<int> &subset, double c = 1.0);

  /// Fits polynomial surface model to subset of 3D point cloud
  ///
  /// The dependent variable of the surface model is the distance of a point from
  /// the surface.
  ///
  /// \param[in] x      Points on surface.
  /// \param[in] n      Surface normals at input points.
  /// \param[in] subset Indices of points to use for fitting.
  /// \param[in] c      Offset along normal direction for points inside and outside
  ///                   the surface to be modelled.
  ///
  /// \returns RMS error of the regression.
  double FitSurface(const PointSet &x, const PointSet &n,
                    const OrderedSet<int> &subset, double c = 1.0);

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
  /// \param[in] x Values of independent variable. The dimension of the model must be 1.
  ///
  /// \returns Regressed values.
  Vector Evaluate(const Vector &x) const;

  /// Evaluates polynomial for one 3D point
  ///
  /// \param[in] x    Data point. The dimension of the model must be 2 or 3.
  /// \param[in] twoD Whether to ignore z coordinate.
  ///
  /// \returns Regressed value.
  double Evaluate(const Point &x, bool twoD = false) const;

  /// Evaluates polynomial for one 3D point
  ///
  /// \param[in] x First  point coordinate.
  /// \param[in] y Second point coordinate.
  /// \param[in] z Third  point coordinate.
  ///
  /// \returns Regressed value.
  double Evaluate(double x, double y, double z) const;

  /// Evaluates polynomial for one 2D point
  ///
  /// \param[in] x First  point coordinate.
  /// \param[in] y Second point coordinate.
  ///
  /// \returns Regressed value.
  double Evaluate(double x, double y) const;

  /// Evaluates polynomial for one independent variable value
  ///
  /// \param[in] x Data point. The dimension of the model must be 1.
  ///
  /// \returns Regressed value.
  double Evaluate(double x) const;

  // ---------------------------------------------------------------------------
  // 1st order derivatives

  /// Evaluate 1st order partial derivative of polynomial
  ///
  /// \param[in] j Index of independent variable w.r.t. which the derivative is taken.
  /// \param[in] x Matrix of size (n x p), where n is the number of data points
  ///              to evaluate, and p is the dimension of the independent variable
  ///              space of the polynomial model.
  ///
  /// \returns Values of 1st order derivative evaluated at given data points.
  Vector Evaluate1stOrderDerivative(int j, const Matrix &x) const;

  /// Evaluate 1st order partial derivative of polynomial
  ///
  /// \param[in] j    Index of independent variable w.r.t. which the derivative is taken.
  /// \param[in] x    Data points. The dimension of the model must be 2 or 3.
  /// \param[in] twoD Whether to ignore z coordinate.
  ///
  /// \returns Values of 1st order derivative evaluated at given data points.
  Vector Evaluate1stOrderDerivative(int j, const PointSet &x, bool twoD = false) const;

  /// Evaluate 1st order partial derivative of polynomial
  ///
  /// \param[in] j Index of independent variable w.r.t. which the derivative is taken.
  /// \param[in] x Values of independent variable. The dimension of the model must be 1.
  ///
  /// \returns Values of 1st order derivative evaluated at given data points.
  Vector Evaluate1stOrderDerivative(int j, const Vector &x) const;

  /// Evaluate 1st order partial derivative of polynomial
  ///
  /// \param[in] j    Index of independent variable w.r.t. which the derivative is taken.
  /// \param[in] x    Data point. The dimension of the model must be 2 or 3.
  /// \param[in] twoD Whether to ignore z coordinate.
  ///
  /// \returns Value of 1st order derivative evaluated at given data point.
  double Evaluate1stOrderDerivative(int j, const Point &x, bool twoD = false) const;

  /// Evaluate 1st order partial derivative of polynomial
  ///
  /// \param[in] j Index of independent variable w.r.t. which the derivative is taken.
  /// \param[in] x Data point. The dimension of the model must be 1.
  ///
  /// \returns Value of 1st order derivative evaluated at given data point.
  double Evaluate1stOrderDerivative(int j, double x) const;

  /// Evaluate gradient of 3D polynomial
  ///
  /// \param[in] x Data point. The dimension of the model must be 3.
  ///
  /// \returns Values of 1st order derivatives evaluated at given data point.
  Vector3 EvaluateGradient(const Point &x) const;

  // ---------------------------------------------------------------------------
  // 2nd order derivatives

  /// Evaluate 2nd order partial derivative of polynomial
  ///
  /// \param[in] x  Matrix of size (n x p), where n is the number of data points
  ///               to evaluate, and p is the dimension of the independent variable
  ///               space of the polynomial model.
  /// \param[in] j1 Index of independent variable w.r.t. which the 1st derivative is taken.
  /// \param[in] j2 Index of independent variable w.r.t. which the 2nd derivative is taken.
  ///
  /// \returns Values of 2nd order derivative evaluated at given data points.
  Vector Evaluate2ndOrderDerivative(int j1, int j2, const Matrix &x) const;

  /// Evaluate 2nd order partial derivative of polynomial
  ///
  /// \param[in] j1   Index of independent variable w.r.t. which the 1st derivative is taken.
  /// \param[in] j2   Index of independent variable w.r.t. which the 2nd derivative is taken.
  /// \param[in] x    Data points. The dimension of the model must be 2 or 3.
  /// \param[in] twoD Whether to ignore z coordinate.
  ///
  /// \returns Values of 2nd order derivative evaluated at given data points.
  Vector Evaluate2ndOrderDerivative(int j1, int j2, const PointSet &x, bool twoD = false) const;

  /// Evaluate 2nd order partial derivative of polynomial
  ///
  /// \param[in] j1 Index of independent variable w.r.t. which the 1st derivative is taken.
  /// \param[in] j2 Index of independent variable w.r.t. which the 2nd derivative is taken.
  /// \param[in] x  Values of independent variable. The dimension of the model must be 1.
  ///
  /// \returns Values of 2nd order derivative evaluated at given data points.
  Vector Evaluate2ndOrderDerivative(int j1, int j2, const Vector &x) const;

  /// Evaluate 2nd order partial derivative of polynomial
  ///
  /// \param[in] j1   Index of independent variable w.r.t. which the 1st derivative is taken.
  /// \param[in] j2   Index of independent variable w.r.t. which the 2nd derivative is taken.
  /// \param[in] x    Data point. The dimension of the model must be 2 or 3.
  /// \param[in] twoD Whether to ignore z coordinate.
  ///
  /// \returns Value of 2nd order derivative evaluated at given data point.
  double Evaluate2ndOrderDerivative(int j1, int j2, const Point &x, bool twoD = false) const;

  /// Evaluate 2nd order partial derivative of polynomial
  ///
  /// \param[in] j1 Index of independent variable w.r.t. which the 1st derivative is taken.
  /// \param[in] j2 Index of independent variable w.r.t. which the 2nd derivative is taken.
  /// \param[in] x  Data point. The dimension of the model must be 1.
  ///
  /// \returns Value of 2nd order derivative evaluated at given data point.
  double Evaluate2ndOrderDerivative(int j1, int j2, double x) const;

  /// Evaluate Hessian of 3D polynomial
  ///
  /// \param[in] x Data point. The dimension of the model must be 3.
  ///
  /// \returns Values of 2nd order derivatives evaluated at given data point.
  Matrix3x3 EvaluateHessian(const Point &x) const;

  // ---------------------------------------------------------------------------
  // Implicit surface

  /// Compute Taubin distance of a point to the implicit polynomial surface
  ///
  /// The Taubin distance is the agebraic distance, i.e., the value of the
  /// implicit polynomial function at \p x, divided the norm of its gradient.
  /// It is a first-order approximation to the Euclidean distance of the point
  /// to the closest point on the surface.
  ///
  /// \param[in] x Data point. The dimension of the model must be 3.
  ///
  /// \returns Taubin distance of the point to the surface.
  double EvaluateTaubinDistance(const Point &x) const;

  /// Evaluate Gaussian curvature of implicit polynomial surface
  ///
  /// Ron Goldman, Curvature formulas for implicit curves and surfaces,
  /// Computer Aided Geometric Design 22 (2005) 632–658.
  ///
  /// \sa http://www.cgeo.ulg.ac.be/CAO/Goldman_Curvature_formulas_implicit.pdf
  ///
  /// \param[in] x Data point. The dimension of the model must be 3.
  ///
  /// \returns Gaussian curvature of implicit surface at \p x.
  double EvaluateGaussianCurvature(const Point &x) const;

  /// Evaluate mean curvature of implicit polynomial surface
  ///
  /// Ron Goldman, Curvature formulas for implicit curves and surfaces,
  /// Computer Aided Geometric Design 22 (2005) 632–658.
  ///
  /// \sa http://www.cgeo.ulg.ac.be/CAO/Goldman_Curvature_formulas_implicit.pdf
  ///
  /// \param[in] x Data point. The dimension of the model must be 3.
  ///
  /// \returns Mean curvature of implicit surface at \p x.
  double EvaluateMeanCurvature(const Point &x) const;

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
inline int Polynomial::NumberOfTerms() const
{
  return _ModelTerms.Rows();
}

// -----------------------------------------------------------------------------
inline int Polynomial::NumberOfPassiveTerms() const
{
  int n = 0;
  for (int i = 0; i < NumberOfTerms(); ++i) {
    if (Status(i) == Passive) ++n;
  }
  return n;
}

// -----------------------------------------------------------------------------
inline int Polynomial::NumberOfActiveTerms() const
{
  return NumberOfTerms() - NumberOfPassiveTerms();
}

// -----------------------------------------------------------------------------
inline int Polynomial::Exponent(int i, int j) const
{
  mirtkAssert(i < NumberOfTerms(), "model term index is within bounds");
  mirtkAssert(j < _Dimension,      "variable index is within bounds");
  return static_cast<int>(_ModelTerms(i, j));
}

// -----------------------------------------------------------------------------
inline bool Polynomial::IsConstant(int i) const
{
  mirtkAssert(i < NumberOfTerms(), "model term index is within bounds");
  return (_ModelTerms.RowSum(i) == .0);
}

// -----------------------------------------------------------------------------
inline void Polynomial::Coefficients(const Vector &coeff)
{
  mirtkAssert(coeff.Rows() == NumberOfTerms(), "coefficient vector has element for each term");
  _Coefficients = coeff;
}

// -----------------------------------------------------------------------------
inline void Polynomial::Coefficient(int i, double c)
{
  mirtkAssert(i < NumberOfTerms(), "model term index is within bounds");
  _Coefficients(i) = c;
}

// -----------------------------------------------------------------------------
inline double Polynomial::Coefficient(int i) const
{
  mirtkAssert(i < NumberOfTerms(), "model term index is within bounds");
  return _Coefficients(i);
}

// -----------------------------------------------------------------------------
inline void Polynomial::Status(int i, enum Status s)
{
  mirtkAssert(i < NumberOfTerms(), "coefficient index is within bounds");
  _Status[i] = s;
}

// -----------------------------------------------------------------------------
inline enum Status Polynomial::Status(int i) const
{
  mirtkAssert(i < NumberOfTerms(), "coefficient index is within bounds");
  return _Status[i];
}

// -----------------------------------------------------------------------------
inline void Polynomial::SetCoefficientsToZero()
{
  _Coefficients = .0;
}

// =============================================================================
// Regression
// =============================================================================

// -----------------------------------------------------------------------------
inline double Polynomial::Fit(const PointSet &x, const Vector &y, bool twoD)
{
  return Fit(Matrix(x, twoD), y);
}

// -----------------------------------------------------------------------------
inline double Polynomial::Fit(const Vector &x, const Vector &y)
{
  return Fit(Matrix(x), y);
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
inline Vector Polynomial::Evaluate(const PointSet &x, bool twoD) const
{
  mirtkAssert(twoD ? _Dimension == 2 : _Dimension == 3,
      "dimension of independent variable space must be " << (twoD ? 2 : 3));
  return Evaluate(Matrix(x, twoD));
}

// -----------------------------------------------------------------------------
inline Vector Polynomial::Evaluate(const Vector &x) const
{
  mirtkAssert(_Dimension == 1, "dimension of independent variable space must be 1");
  return Evaluate(Matrix(x));
}

// -----------------------------------------------------------------------------
inline double Polynomial::Evaluate(const Point &p, bool twoD) const
{
  mirtkAssert(twoD ? _Dimension == 2 : _Dimension == 3,
      "dimension of independent variable space must be " << (twoD ? 2 : 3));
  Matrix x(1, twoD ? 2 : 3);
  x(0, 0) = p._x;
  x(0, 1) = p._y;
  if (!twoD) x(0, 2) = p._z;
  return Evaluate(x)(0);
}

// -----------------------------------------------------------------------------
inline double Polynomial::Evaluate(double x, double y, double z) const
{
  mirtkAssert(_Dimension == 3, "dimension of independent variable space must be 3");
  Matrix m(1, 3);
  m(0, 0) = x;
  m(0, 1) = y;
  m(0, 2) = z;
  return Evaluate(m)(0);
}

// -----------------------------------------------------------------------------
inline double Polynomial::Evaluate(double x, double y) const
{
  mirtkAssert(_Dimension == 2, "dimension of independent variable space must be 2");
  Matrix m(1, 2);
  m(0, 0) = x;
  m(0, 1) = y;
  return Evaluate(m)(0);
}

// -----------------------------------------------------------------------------
inline double Polynomial::Evaluate(double x) const
{
  mirtkAssert(_Dimension == 1, "dimension of independent variable space must be 1");
  Matrix m(1, 1);
  m(0, 0) = x;
  return Evaluate(m)(0);
}

// =============================================================================
// 1st order derivatives
// =============================================================================

// -----------------------------------------------------------------------------
inline Vector Polynomial::Evaluate1stOrderDerivative(int j1, const PointSet &x, bool twoD) const
{
  mirtkAssert(twoD ? _Dimension == 2 : _Dimension == 3,
      "dimension of independent variable space must be " << (twoD ? 2 : 3));
  return Evaluate1stOrderDerivative(j1, Matrix(x, twoD));
}

// -----------------------------------------------------------------------------
inline Vector Polynomial::Evaluate1stOrderDerivative(int j1, const Vector &x) const
{
  mirtkAssert(_Dimension == 1, "dimension of independent variable space must be 1");
  return Evaluate1stOrderDerivative(j1, Matrix(x));
}

// -----------------------------------------------------------------------------
inline double Polynomial::Evaluate1stOrderDerivative(int j1, const Point &p, bool twoD) const
{
  mirtkAssert(twoD ? _Dimension == 2 : _Dimension == 3,
      "dimension of independent variable space must be " << (twoD ? 2 : 3));
  Matrix x(1, twoD ? 2 : 3);
  x(0, 0) = p._x;
  x(0, 1) = p._y;
  if (!twoD) x(0, 2) = p._z;
  return Evaluate1stOrderDerivative(j1, x)(0);
}

// -----------------------------------------------------------------------------
inline double Polynomial::Evaluate1stOrderDerivative(int j1, double x) const
{
  mirtkAssert(_Dimension == 1, "dimension of independent variable space must be 1");
  Matrix m(1, 1);
  m(0, 0) = x;
  return Evaluate1stOrderDerivative(j1, m)(0);
}

// =============================================================================
// 2nd order derivatives
// =============================================================================

// -----------------------------------------------------------------------------
inline Vector Polynomial::Evaluate2ndOrderDerivative(int j1, int j2, const PointSet &x, bool twoD) const
{
  mirtkAssert(twoD ? _Dimension == 2 : _Dimension == 3,
      "dimension of independent variable space must be " << (twoD ? 2 : 3));
  return Evaluate2ndOrderDerivative(j1, j2, Matrix(x, twoD));
}

// -----------------------------------------------------------------------------
inline Vector Polynomial::Evaluate2ndOrderDerivative(int j1, int j2, const Vector &x) const
{
  mirtkAssert(_Dimension == 1, "dimension of independent variable space must be 1");
  return Evaluate2ndOrderDerivative(j1, j2, Matrix(x));
}

// -----------------------------------------------------------------------------
inline double Polynomial::Evaluate2ndOrderDerivative(int j1, int j2, const Point &p, bool twoD) const
{
  mirtkAssert(twoD ? _Dimension == 2 : _Dimension == 3,
      "dimension of independent variable space must be " << (twoD ? 2 : 3));
  Matrix x(1, twoD ? 2 : 3);
  x(0, 0) = p._x;
  x(0, 1) = p._y;
  if (!twoD) x(0, 2) = p._z;
  return Evaluate2ndOrderDerivative(j1, j2, x)(0);
}

// -----------------------------------------------------------------------------
inline double Polynomial::Evaluate2ndOrderDerivative(int j1, int j2, double x) const
{
  mirtkAssert(_Dimension == 1, "dimension of independent variable space must be 1");
  Matrix m(1, 1);
  m(0, 0) = x;
  return Evaluate2ndOrderDerivative(j1, j2, m)(0);
}

// =============================================================================
// Debugging
// =============================================================================

// -----------------------------------------------------------------------------
inline void Polynomial::Print(Indent indent) const
{
  Print(cout, indent);
}


} // namespace mirtk

#endif // MIRTK_Polynomial_H
