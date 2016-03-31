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

#include "mirtk/Polynomial.h"

#include "mirtk/Assert.h"
#include "mirtk/Math.h"
#include "mirtk/Eigen.h"

#include "Eigen/QR"


namespace mirtk {


// =============================================================================
// Auxiliaries
// =============================================================================

// -----------------------------------------------------------------------------
// build the exponent array recursively
static Matrix BuildCompleteModel(int p, int order)
{
  Matrix m;
  if (p > 0) {
    if (order <= 0) {
      m.Initialize(1, p);
    } else if (p == 1) {
      m.Initialize(order + 1, 1);
      for (int i = 0; i <= order; ++i) {
        m(i, 0) = order - i;
      }
    } else {
      Matrix m2;
      for (int i = order; i >= 0; --i) {
        const int r0 = m.Rows(); // before Resize
        m2 = BuildCompleteModel(p - 1, order - i);
        m.Resize(r0 + m2.Rows(), p);
        for (int r = 0; r < m2.Rows(); ++r) {
          m(r0 + r, 0) = i;
          for (int c = 0; c < m2.Cols(); ++c) {
            m(r0 + r, c + 1) = m2(r, c);
          }
        }
      }
    }
  }
  return m;
}

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
Polynomial::Polynomial(int order)
:
  _Order(order),
  _Dimension(0)
{
}

// -----------------------------------------------------------------------------
Polynomial::Polynomial(int p, int order, const Vector &coeff)
:
  _Order(0),
  _Dimension(0)
{
  Initialize(p, order);
  if (coeff.Rows() > 0) Coefficients(coeff);
}

// -----------------------------------------------------------------------------
void Polynomial::CopyAttributes(const Polynomial &other)
{
  _Order        = other._Order;
  _Dimension    = other._Dimension;
  _ModelTerms   = other._ModelTerms;
  _Coefficients = other._Coefficients;
  _Status       = other._Status;
}

// -----------------------------------------------------------------------------
Polynomial::Polynomial(const Polynomial &other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
Polynomial &Polynomial::operator =(const Polynomial &other)
{
  if (this != &other) {
    Object::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
Polynomial::~Polynomial()
{
}

// -----------------------------------------------------------------------------
int Polynomial::Initialize(int p, int order)
{
  if (p <= 0) {
    cerr << this->NameOfType() << "::Initialize: Dimension of independent variable space must be positive" << endl;
    exit(1);
  }
  if ((order > 0 && order != _Order) || p != _Dimension) {
    // Set order of polynomial model
    if (order > 0) {
      _Order = order;
    } else if (_Order <= 0) {
      cerr << this->NameOfType() << "::Initialize: Order of polynomial not set during construction" << endl;
      exit(1);
    }
    _Dimension  = p;
    _ModelTerms = BuildCompleteModel(_Dimension, _Order);
  }
  // Set coefficients to zero
  _Coefficients.Initialize(NumberOfTerms());
  // Set status of coefficients to active
  _Status.resize(NumberOfTerms());
  for (int i = 0; i < NumberOfTerms(); ++i) {
    _Status[i] = Active;
  }
  // Return number of terms
  return NumberOfTerms();
}

// -----------------------------------------------------------------------------
void Polynomial::SetConstantCoefficient(double value, enum Status status)
{
  for (int i = 0; i < NumberOfTerms(); ++i) {
    if (IsConstant(i)) {
      Coefficient(i, value);
      Status(i, status);
    }
  }
}

// =============================================================================
// Regression
// =============================================================================

// -----------------------------------------------------------------------------
double Polynomial::Fit(const Matrix &x, const Vector &y)
{
  const int n  = x.Rows();
  const int p  = x.Cols();
  const int nt = NumberOfTerms();

  if (n <= 0) {
    cerr << this->NameOfType() << "::Fit: No input points given!" << endl;
    exit(1);
  }
  if (y.Rows() != n) {
    cerr << this->NameOfType() << "::Fit: x and y have differing number of rows!" << endl;
    exit(1);
  }

  // Scale x to unit variance
  Vector sigma(p);
  Matrix xs(x);
  for (int i = 0; i < p; ++i) {
    sigma(i) = sqrt(x.ColVar(i));
    if (sigma(i) == .0) sigma(i) = 1.0;
    else xs.ScaleCol(i, 1.0 / sigma(i));
  }

  // Build design matrix
  Eigen::MatrixXd A(n, NumberOfActiveTerms());
  Eigen::VectorXd b = VectorToEigen(y);
  Vector scale(static_cast<int>(A.cols()));
  double t;

  A.setOnes(), scale = 1.0;
  for (int i = 0, c = 0; i < nt; ++i) {
    if (Status(i) == Active) {
      for (int j = 0; j < p;  ++j) {
        if (_ModelTerms(i, j) != 0) {
          for (int k = 0; k < n; ++k) {
            A(k, c) *= pow(xs(k, j), _ModelTerms(i, j));
          }
          scale(c) /= pow(sigma(j), _ModelTerms(i, j));
        }
      }
      ++c;
    } else {
      for (int k = 0; k < n; ++k) {
        t = 1.0;
        for (int j = 0; j < p;  ++j) {
          t *= pow(x(k, j), _ModelTerms(i, j));
        }
        b(k) -= t * _Coefficients(i);
      }
    }
  }

  // Estimate model using (column) pivoted QR for stability
  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(A);
  const Eigen::VectorXd coeff = qr.solve(b);

  // Recover the scaling
  for (int i = 0, c = 0; i < NumberOfTerms(); ++i) {
    if (Status(i) == Active) {
      _Coefficients(i) = scale(c) * coeff[c];
      ++c;
    }
  }

  // Return RMS error
  Eigen::VectorXd ypred = A * coeff;
  if (NumberOfPassiveTerms() > 0) {
    ypred += VectorToEigen(y) - b;
  }
  double sum2 = .0, delta;
  for (int i = 0; i < n; ++i) {
    delta = y(i) - ypred(i);
    sum2 += delta * delta;
  }
  return sqrt(sum2 / n);
}

// -----------------------------------------------------------------------------
double Polynomial::Fit(const Matrix &x, const Vector &y, const Array<int> &subset)
{
  Matrix _x(static_cast<int>(subset.size()), x.Cols());
  Vector _y(_x.Rows());
  for (int i = 0; i < _x.Rows(); ++i) {
    for (int j = 0; j < _x.Cols(); ++j) {
      _x(i, j) = x(subset[i], j);
    }
    _y(i) = y(subset[i]);
  }
  return Fit(_x, _y);
}

// -----------------------------------------------------------------------------
double Polynomial::Fit(const Matrix &x, const Vector &y, const OrderedSet<int> &subset)
{
  Matrix _x(static_cast<int>(subset.size()), x.Cols());
  Vector _y(_x.Rows());
  OrderedSet<int>::const_iterator it = subset.begin();
  for (int i = 0; i < _x.Rows(); ++i, ++it) {
    for (int j = 0; j < _x.Cols(); ++j) {
      _x(i, j) = x(*it, j);
    }
    _y(i) = y(*it);
  }
  return Fit(_x, _y);
}

// -----------------------------------------------------------------------------
double Polynomial::FitSurface(const PointSet &x, const PointSet &n, double c)
{
  mirtkAssert(x.Size() == n.Size(), "one normal for each point");
  Matrix m(3 * x.Size(), 3);
  Vector d(m.Rows());
  for (int i = 0, j = 0; j < x.Size(); ++j) {
    const Point &point  = x(j);
    const Point &normal = n(j);
    // point inside
    m(i, 0) = point._x - c * normal._x;
    m(i, 1) = point._y - c * normal._y;
    m(i, 2) = point._z - c * normal._z;
    d(i) = -c;
    ++i;
    // point on surface
    m(i, 0) = point._x;
    m(i, 1) = point._y;
    m(i, 2) = point._z;
    d(i) = .0;
    ++i;
    // point outside
    m(i, 0) = point._x + c * normal._x;
    m(i, 1) = point._y + c * normal._y;
    m(i, 2) = point._z + c * normal._z;
    d(i) = c;
    ++i;
  }
  return Fit(m, d);
}

// -----------------------------------------------------------------------------
double Polynomial::FitSurface(const PointSet &x, const PointSet &n,
                              const Array<int> &subset, double c)
{
  mirtkAssert(x.Size() == n.Size(), "one normal for each point");
  Matrix m(3 * static_cast<int>(subset.size()), 3);
  Vector d(m.Rows());
  int i = 0;
  for (Array<int>::const_iterator it = subset.begin(); it != subset.end(); ++it) {
    const Point &point  = x(*it);
    const Point &normal = n(*it);
    // point inside
    m(i, 0) = point._x - c * normal._x;
    m(i, 1) = point._y - c * normal._y;
    m(i, 2) = point._z - c * normal._z;
    d(i) = -c;
    ++i;
    // point on surface
    m(i, 0) = point._x;
    m(i, 1) = point._y;
    m(i, 2) = point._z;
    d(i) = .0;
    ++i;
    // point outside
    m(i, 0) = point._x + c * normal._x;
    m(i, 1) = point._y + c * normal._y;
    m(i, 2) = point._z + c * normal._z;
    d(i) = c;
    ++i;
  }
  return Fit(m, d);
}

// -----------------------------------------------------------------------------
double Polynomial::FitSurface(const PointSet &x, const PointSet &n,
                              const OrderedSet<int> &subset, double c)
{
  mirtkAssert(x.Size() == n.Size(), "one normal for each point");
  Matrix m(3 * static_cast<int>(subset.size()), 3);
  Vector d(m.Rows());
  int i = 0;
  for (OrderedSet<int>::const_iterator it = subset.begin(); it != subset.end(); ++it) {
    const Point &point  = x(*it);
    const Point &normal = n(*it);
    // point inside
    m(i, 0) = point._x - c * normal._x;
    m(i, 1) = point._y - c * normal._y;
    m(i, 2) = point._z - c * normal._z;
    d(i) = -c;
    ++i;
    // point on surface
    m(i, 0) = point._x;
    m(i, 1) = point._y;
    m(i, 2) = point._z;
    d(i) = .0;
    ++i;
    // point outside
    m(i, 0) = point._x + c * normal._x;
    m(i, 1) = point._y + c * normal._y;
    m(i, 2) = point._z + c * normal._z;
    d(i) = c;
    ++i;
  }
  return Fit(m, d);
}

// -----------------------------------------------------------------------------
double Polynomial::FitSurface(const PointSet &x, const Array<int> &subset)
{
  Matrix m(static_cast<int>(subset.size()), 3);
  Vector d(m.Rows());
  int i = 0;
  for (Array<int>::const_iterator it = subset.begin(); it != subset.end(); ++it) {
    const Point &point = x(*it);
    m(i, 0) = point._x;
    m(i, 1) = point._y;
    m(i, 2) = point._z;
    d(i)    = .0;
    ++i;
  }
  SetConstantCoefficient(1.0);
  return Fit(m, d);
}

// -----------------------------------------------------------------------------
double Polynomial::FitSurface(const PointSet &x, const OrderedSet<int> &subset)
{
  Matrix m(static_cast<int>(subset.size()), 3);
  Vector d(m.Rows());
  int i = 0;
  for (OrderedSet<int>::const_iterator it = subset.begin(); it != subset.end(); ++it) {
    const Point &point = x(*it);
    m(i, 0) = point._x;
    m(i, 1) = point._y;
    m(i, 2) = point._z;
    d(i)    = .0;
    ++i;
  }
  SetConstantCoefficient(1.0);
  return Fit(m, d);
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
Vector Polynomial::Evaluate(const Matrix &x) const
{
  mirtkAssert(_Dimension == x.Cols(), "independent variable dimensions must match");
  double t;
  Vector y(x.Rows());
  for (int k = 0; k < x.Rows(); ++k) {
    y(k) = .0;
    for (int i = 0; i < _ModelTerms.Rows(); ++i) {
      t = 1.0;
      for (int j = 0; j < _Dimension; ++j) {
        t *= pow(x(k, j), _ModelTerms(i, j));
      }
      y(k) += t * _Coefficients(i);
    }
  }
  return y;
}

// -----------------------------------------------------------------------------
Vector Polynomial::Evaluate1stOrderDerivative(int j1, const Matrix &x) const
{
  mirtkAssert(_Dimension == x.Cols(), "independent variable dimensions must match");
  double t;
  Vector y(x.Rows());
  for (int k = 0; k < x.Rows(); ++k) {
    y(k) = .0;
    for (int i = 0; i < _ModelTerms.Rows(); ++i) {
      t = 1.0;
      for (int j = 0; j < _Dimension; ++j) {
        if (_ModelTerms(i, j) != .0) {
          if (j == j1) {
            t *= _ModelTerms(i, j) * pow(x(k, j), _ModelTerms(i, j) - 1);
          } else {
            t *= pow(x(k, j), _ModelTerms(i, j));
          }
        }
      }
      y(k) += t * _Coefficients(i);
    }
  }
  return y;
}

// -----------------------------------------------------------------------------
Vector3 Polynomial::EvaluateGradient(const Point &p) const
{
  mirtkAssert(_Dimension == 3, "independent variable dimensions must match");
  Vector3 g(.0), t;
  for (int i = 0; i < _ModelTerms.Rows(); ++i) {
    t[0] = pow(p._x, _ModelTerms(i, 0) - 1) * _ModelTerms(i, 0)
         * pow(p._y, _ModelTerms(i, 1))
         * pow(p._z, _ModelTerms(i, 2));
    t[1] = pow(p._x, _ModelTerms(i, 0))
         * pow(p._y, _ModelTerms(i, 1) - 1) * _ModelTerms(i, 1)
         * pow(p._z, _ModelTerms(i, 2));
    t[2] = pow(p._x, _ModelTerms(i, 0))
         * pow(p._y, _ModelTerms(i, 1))
         * pow(p._z, _ModelTerms(i, 2) - 1) * _ModelTerms(i, 2);
    g += (t *= _Coefficients(i));
  }
  return g;
}

// -----------------------------------------------------------------------------
Vector Polynomial::Evaluate2ndOrderDerivative(int j1, int j2, const Matrix &x) const
{
  mirtkAssert(_Dimension == x.Cols(), "independent variable dimensions must match");
  double t;
  Vector y(x.Rows());
  for (int k = 0; k < x.Rows(); ++k) {
    y(k) = .0;
    for (int i = 0; i < _ModelTerms.Rows(); ++i) {
      t = 1.0;
      for (int j = 0; j < _Dimension; ++j) {
        if (_ModelTerms(i, j) != .0) {
          if (j == j1 || j == j2) {
            if (j1 == j2) {
              t *= (_ModelTerms(i, j) - 1) * _ModelTerms(i, j) * pow(x(k, j), _ModelTerms(i, j) - 2);
            } else {
              t *= _ModelTerms(i, j) * pow(x(k, j), _ModelTerms(i, j) - 1);
            }
          } else {
            t *= pow(x(k, j), _ModelTerms(i, j));
          }
        }
      }
      y(k) += t * _Coefficients(i);
    }
  }
  return y;
}

// -----------------------------------------------------------------------------
Matrix3x3 Polynomial::EvaluateHessian(const Point &p) const
{
  mirtkAssert(_Dimension == 3, "independent variable dimensions must match");
  Matrix3x3 h;
  for (int i = 0; i < _ModelTerms.Rows(); ++i) {
    h[0][0] += _Coefficients(i)
             * pow(p._x, _ModelTerms(i, 0) - 2) * _ModelTerms(i, 0) * (_ModelTerms(i, 0) - 1)
             * pow(p._y, _ModelTerms(i, 1))
             * pow(p._z, _ModelTerms(i, 2));
    h[0][1] += _Coefficients(i)
             * pow(p._x, _ModelTerms(i, 0) - 1) * _ModelTerms(i, 0)
             * pow(p._y, _ModelTerms(i, 1) - 1) * _ModelTerms(i, 1)
             * pow(p._z, _ModelTerms(i, 2));
    h[0][2] += _Coefficients(i)
             * pow(p._x, _ModelTerms(i, 0) - 1) * _ModelTerms(i, 0)
             * pow(p._y, _ModelTerms(i, 1))
             * pow(p._z, _ModelTerms(i, 2) - 1) * _ModelTerms(i, 2);
    h[1][1] += _Coefficients(i)
             * pow(p._x, _ModelTerms(i, 0))
             * pow(p._y, _ModelTerms(i, 1) - 2) * _ModelTerms(i, 1) * (_ModelTerms(i, 1) - 1)
             * pow(p._z, _ModelTerms(i, 2));
    h[1][2] += _Coefficients(i)
             * pow(p._x, _ModelTerms(i, 0))
             * pow(p._y, _ModelTerms(i, 1) - 1) * _ModelTerms(i, 1)
             * pow(p._z, _ModelTerms(i, 2) - 1) * _ModelTerms(i, 2);
    h[2][2] += _Coefficients(i)
             * pow(p._x, _ModelTerms(i, 0))
             * pow(p._y, _ModelTerms(i, 1))
             * pow(p._z, _ModelTerms(i, 2) - 2) * _ModelTerms(i, 2) * (_ModelTerms(i, 2) - 1);
  }
  h[1][0] = h[0][1];
  h[2][0] = h[0][2];
  h[2][1] = h[1][2];
  return h;
}

// =============================================================================
// Implicit surface
// =============================================================================

// -----------------------------------------------------------------------------
double Polynomial::EvaluateTaubinDistance(const Point &p) const
{
  mirtkAssert(_Dimension == 3, "independent variable dimensions must match");
  return Evaluate(p) / EvaluateGradient(p).Length();
}

// -----------------------------------------------------------------------------
double Polynomial::EvaluateGaussianCurvature(const Point &p) const
{
  mirtkAssert(_Dimension == 3, "independent variable dimensions must match");
  Vector3   g = EvaluateGradient(p);
  Matrix3x3 h = EvaluateHessian(p).Adjoint();
  const double l2 = g.SquaredLength();
  return g.Dot(h * g) / (l2 * l2);
}

// -----------------------------------------------------------------------------
double Polynomial::EvaluateMeanCurvature(const Point &p) const
{
  mirtkAssert(_Dimension == 3, "independent variable dimensions must match");
  Vector3   g = EvaluateGradient(p);
  Matrix3x3 h = EvaluateHessian(p);
  const double l2 = g.SquaredLength();
  return (g.Dot(h * g) - l2 * h.Trace()) / (2.0 * sqrt(l2) * l2);
}

// =============================================================================
// Debugging
// =============================================================================

// -----------------------------------------------------------------------------
ostream &Polynomial::Print(ostream &os, Indent indent) const
{
  ostringstream line;
  line << indent << "p(x) =";
  if (NumberOfTerms() == 0) {
    line << " undefined";
  } else {
    const int   varcount           = 3;
    const char *varnames[varcount] = {"x", "y", "z"};
    for (int i = 0; i < NumberOfTerms(); ++i) {
      if (i > 0) {
        if (line.str().size() > 70) {
          os << line.str() << "\n";
          line.str("");
          line << indent << "     +";
        } else {
          line << " +";
        }
      }
      line << " " << _Coefficients(i);
      for (int j = 0; j < _Dimension; ++j) {
        if (_ModelTerms(i, j) == .0) continue;
        line << " ";
        if (_Dimension <= varcount) {
          line << varnames[j];
        } else {
          line << "x" << ToString(j);
        }
        if (_ModelTerms(i, j) != 1.0) {
          line << "^" << ToString(_ModelTerms(i, j));
        }
      }
    }
  }
  os << line.str() << "\n";
  return os;
}


} // namespace mirtk
