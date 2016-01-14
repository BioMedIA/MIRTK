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

#include <mirtkPolynomial.h>

#include <mirtkAssert.h>
#include <mirtkMath.h>
#include <mirtkEigen.h>

#include <Eigen/QR>


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
Polynomial::Polynomial(int p, int order)
:
  _Order(0),
  _Dimension(0)
{
  Initialize(p, order);
}

// -----------------------------------------------------------------------------
void Polynomial::CopyAttributes(const Polynomial &other)
{
  _Order        = other._Order;
  _Dimension    = other._Dimension;
  _ModelTerms   = other._ModelTerms;
  _Coefficients = other._Coefficients;
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
  // Return number of terms
  return NumberOfTerms();
}

// =============================================================================
// Regression
// =============================================================================

// -----------------------------------------------------------------------------
double Polynomial::Fit(const Matrix &x, const Vector &y, int order)
{
  const int n  = x.Rows();
  const int p  = x.Cols();
  const int nt = Initialize(p, order);

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
  Eigen::MatrixXd A(n, nt);
  Vector scale(nt);

  A.setOnes(), scale = 1.0;
  for (int i = 0; i < nt; ++i)
  for (int j = 0; j < p;  ++j) {
    if (_ModelTerms(i, j) != 0) {
      for (int k = 0; k < n; ++k) {
        A(k, i) *= pow(xs(k, j), _ModelTerms(i, j));
      }
      scale(i) /= pow(sigma(j), _ModelTerms(i, j));
    }
  }

  // Estimate model using (column) pivoted QR for stability
  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(A);
  const Eigen::VectorXd coeff = qr.solve(VectorToEigen(y));
  const Eigen::VectorXd ypred = A * coeff;

  // Recover the scaling
  _Coefficients = EigenToVector(coeff);
  for (int i = 0; i < nt; ++i) {
    _Coefficients(i) *= scale(i);
  }

  // Return RMS error
  double sum2 = .0, delta;
  for (int i = 0; i < n; ++i) {
    delta = y(i) - ypred(i);
    sum2 += delta * delta;
  }
  return sqrt(sum2 / n);
}

// -----------------------------------------------------------------------------
double Polynomial::Fit(const PointSet &x, const Vector &y, int order, bool twoD)
{
  return Fit(Matrix(x, twoD), y, order);
}

// -----------------------------------------------------------------------------
double Polynomial::Fit(const Vector &x, const Vector &y, int order)
{
  return Fit(Matrix(x), y, order);
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
      line << " a" << ToString(i);
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
