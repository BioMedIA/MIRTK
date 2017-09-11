/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2015 Imperial College London
 * Copyright 2008-2015 Daniel Rueckert, Julia Schnabel
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

#include "mirtk/Matrix.h"

#include "mirtk/Assert.h"
#include "mirtk/Array.h"
#include "mirtk/Memory.h"
#include "mirtk/Math.h"
#include "mirtk/Vector.h"
#include "mirtk/Vector3D.h"
#include "mirtk/Indent.h"
#include "mirtk/PointSet.h"
#include "mirtk/String.h"
#include "mirtk/Stream.h"

#include "mirtk/NumericsConfig.h"
#if MIRTK_Numerics_WITH_MATLAB
#  include "mirtk/Matlab.h"
#endif

#include "mirtk/Eigen.h"

#ifdef _MSC_VER
#  pragma warning(push)
#  pragma warning(disable: 4244) // eigen\src/QR/FullPivHouseholderQR.h(125):
 // warning C4244: 'argument':
 // conversion from '__int64' to 'int',
 // possible loss of data
#endif
#include "Eigen/QR" // included by Eigen/SVD
#ifdef _MSC_VER
#  pragma warning(pop)
#endif

#include "Eigen/LU"
#include "Eigen/SVD"
#include "Eigen/Eigenvalues"


namespace mirtk {


// =============================================================================
// Auxiliary macros
// =============================================================================

// -----------------------------------------------------------------------------
#define MIRTK_CHECK_MATRIX_SIZE(op, m) \
  while (_rows != m._rows || _cols != m._cols) { \
    cerr << "Matrix::" op ": Right-hand matrix must have same size" << endl; \
    exit(1); \
  }

// -----------------------------------------------------------------------------
#define MIRTK_CHECK_MATRIX_IS_SQUARE(op) \
   while (_rows != _cols) { \
    cerr << "Matrix::" op ": Matrix must be square" << endl; \
    exit(1); \
  }

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
Matrix::Matrix()
:
  _rows  (0),
  _cols  (0),
  _matrix(NULL),
  _owner (false)
{
}

// -----------------------------------------------------------------------------
Matrix::Matrix(int rows, int cols, double *data)
:
  _rows  (rows),
  _cols  (cols < 0 ? rows : cols),
  _matrix(NULL),
  _owner (false)
{
  if (_rows * _cols > 0) {
    Allocate(_matrix, _rows, _cols, data);
    _owner = (this->RawPointer() != data);
    if (_owner) {
      memset(this->RawPointer(), 0, _rows * _cols * sizeof(double));
    }
  }
}

// -----------------------------------------------------------------------------
Matrix::Matrix(const Vector &v)
:
  _rows(v.Rows()),
  _cols(1),
  _matrix(NULL),
  _owner(true)
{
  if (_rows > 0) {
    Allocate(_matrix, _rows, _cols);
    memcpy(this->RawPointer(), v.RawPointer(), _rows * sizeof(double));
  } else {
    _owner = false;
  }
}

// -----------------------------------------------------------------------------
Matrix::Matrix(const PointSet &pset, bool twoD)
:
  _rows(pset.Size()),
  _cols(twoD ? 2 : 3),
  _matrix(NULL),
  _owner(true)
{
  if (_rows > 0) {
    Allocate(_matrix, _rows, _cols);
    double *x = _matrix[0];
    double *y = _matrix[1];
    if (twoD) {
      for (int i = 0; i < pset.Size(); ++i, ++x, ++y) {
        *x = pset(i)._x;
        *y = pset(i)._y;
      }
    } else {
      double *z = _matrix[2];
      for (int i = 0; i < pset.Size(); ++i, ++x, ++y, ++z) {
        *x = pset(i)._x;
        *y = pset(i)._y;
        *z = pset(i)._z;
      }
    }
  } else {
    _owner = false;
    _cols  = 0;
  }
}

// -----------------------------------------------------------------------------
Matrix::Matrix(const Matrix& m)
:
  Object(m),
  _rows  (m._rows),
  _cols  (m._cols),
  _matrix(NULL),
  _owner (false)
{
  if (_rows * _cols > 0) {
    Allocate(_matrix, _rows, _cols);
    _owner = true;
    memcpy(this->RawPointer(), m.RawPointer(), _rows * _cols * sizeof(double));
  }
}

// -----------------------------------------------------------------------------
Matrix::~Matrix()
{
  Clear();
}

// -----------------------------------------------------------------------------
Matrix& Matrix::operator =(const Matrix& m)
{
  if (_rows != m._rows || _cols != m._cols) {
    _rows = m._rows;
    _cols = m._cols;
    if (_owner) Deallocate(_matrix);
    Allocate(_matrix, _rows, _cols);
    _owner = true;
  }
  if (_matrix && this->RawPointer() != m.RawPointer()) {
    memcpy(this->RawPointer(), m.RawPointer(), _rows * _cols * sizeof(double));
  }
  return *this;
}

// -----------------------------------------------------------------------------
void Matrix::Initialize(int rows, int cols, double *data)
{
  if (cols < 0) cols = rows;
  if (_rows != rows || _cols != cols || (data && this->RawPointer() != data)) {
    if (_owner) Deallocate(_matrix);
    _rows = rows;
    _cols = cols;
    Allocate(_matrix, _rows, _cols, data);
    _owner = (this->RawPointer() != data);
  }
  if (_owner) {
    memset(this->RawPointer(), 0, _rows * _cols * sizeof(double));
  }
}

// -----------------------------------------------------------------------------
void Matrix::Resize(int rows, int cols)
{
  if (rows <= 0 || cols <= 0) {
    Clear();
  } else if (_rows != rows || _cols != cols) {
    double **matrix = Allocate<double>(rows, cols);
    const int m = min(rows, _rows);
    const int n = min(cols, _cols);
    for (int c = 0; c < n; ++c) {
      for (int r = 0; r < m; ++r) {
        matrix[c][r] = _matrix[c][r];
      }
      for (int r = m; r < rows; ++r) {
        matrix[c][r] = .0;
      }
    }
    for (int c = n; c < cols; ++c) {
      for (int r = 0; r < rows; ++r) {
        matrix[c][r] = .0;
      }
    }
    if (_owner) Deallocate(_matrix);
    _matrix = matrix;
    _owner  = true;
    _rows   = rows;
    _cols   = cols;
  }
}

// -----------------------------------------------------------------------------
void Matrix::Clear()
{
  if (_owner) Deallocate(_matrix);
  _rows = 0;
  _cols = 0;
}

// -----------------------------------------------------------------------------
void Matrix::Zero()
{
  memset(this->RawPointer(), 0, _rows * _cols * sizeof(double));
}

// =============================================================================
// Element access
// =============================================================================

// -----------------------------------------------------------------------------
void Matrix::operator()(Matrix &m, int rows, int cols)
{
  if ((m._rows + rows > _rows) || (m._cols + cols > _cols)) {
    cerr << "Matrix::operator(): Invalid range" << endl;
    exit(1);
  }
  for (int j = cols; j < m._cols + cols; ++j)
  for (int i = rows; i < m._rows + rows; ++i) {
    _matrix[j][i] = m._matrix[j - cols][i - rows];
  }
}

// -----------------------------------------------------------------------------
Matrix Matrix::operator()(int rows1, int cols1, int rows2, int cols2) const
{
  if ((rows1 < 0) || (rows2 > _rows) || (cols1 < 0) || (cols2 > _cols)) {
    cerr << "Matrix::operator(): Invalid range" << endl;
    exit(1);
  }

  // Create new matrix
  Matrix m(rows2-rows1, cols2-cols1);

  for (int j = cols1; j < cols2; ++j)
  for (int i = rows1; i < rows2; ++i) {
    m._matrix[j-cols1][i-rows1] = _matrix[j][i];
  }

  return m;
}

// =============================================================================
// Matrix-vector operations
// =============================================================================

// -----------------------------------------------------------------------------
Vector Matrix::operator *(const Vector &v) const
{
  if (_cols != v.Rows()) {
    cerr << "Matrix::operator*(Vector): Size mismatch" << endl;
    exit(1);
  }

  Vector p(_rows);
  double sum;
  for (int r = 0; r < _rows; ++r) {
    sum = .0;
    for (int c = 0; c < _cols; ++c) {
      sum += _matrix[c][r] * v(c);
    }
    p(r) = sum;
  }

  return p;
}

// =============================================================================
// Matrix-matrix operations
// =============================================================================

// -----------------------------------------------------------------------------
Matrix &Matrix::operator -=(const Matrix &m)
{
  MIRTK_CHECK_MATRIX_SIZE("operator -=", m);
  const int n = this->NumberOfElements();
  const double *ptr1 = m   . RawPointer();
  double       *ptr2 = this->RawPointer();
  for (int i = 0; i < n; ++i, ++ptr1, ++ptr2) (*ptr2) -= (*ptr1);
  return *this;
}

// -----------------------------------------------------------------------------
Matrix &Matrix::operator +=(const Matrix &m)
{
  MIRTK_CHECK_MATRIX_SIZE("operator +=", m);
  const int n = this->NumberOfElements();
  const double *ptr1 = m   . RawPointer();
  double       *ptr2 = this->RawPointer();
  for (int i = 0; i < n; ++i, ++ptr1, ++ptr2) (*ptr2) += (*ptr1);
  return *this;
}

// -----------------------------------------------------------------------------
Matrix& Matrix::operator *=(const Matrix& m)
{
  (*this) = (*this) * m;
  return *this;
}

// -----------------------------------------------------------------------------
Matrix Matrix::operator -(const Matrix& m) const
{
  MIRTK_CHECK_MATRIX_SIZE("operator -", m);
  Matrix m2(*this);
  m2 -= m;
  return m2;
}

// -----------------------------------------------------------------------------
Matrix Matrix::operator +(const Matrix& m) const
{
  MIRTK_CHECK_MATRIX_SIZE("operator +", m);
  Matrix m2(*this);
  m2 += m;
  return m2;
}

// -----------------------------------------------------------------------------
Matrix Matrix::operator *(const Matrix& m) const
{
  if (_cols != m.Rows()) {
    cerr << "Matrix::operator *: Matrix size mismatch" << endl;
    exit(1);
  }
  Matrix m2(_rows, m._cols);
  for (int i = 0; i < _rows; ++i) {
    for (int j = 0; j < m._cols; ++j) {
      m2._matrix[j][i] = 0.0;
      for (int k = 0; k < _cols; ++k) {
        m2._matrix[j][i] += _matrix[k][i] * m._matrix[j][k];
      }
    }
  }
  return m2;
}

// =============================================================================
// Matrix functions
// =============================================================================

// -----------------------------------------------------------------------------
double Matrix::RowSum(int r) const
{
  double sum = .0;
  for (int c = 0; c < _cols; ++c) {
    sum += _matrix[c][r];
  }
  return sum;
}

// -----------------------------------------------------------------------------
double Matrix::RowMean(int r) const
{
  return RowSum(r) / _cols;
}

// -----------------------------------------------------------------------------
double Matrix::RowVar(int r) const
{
  if (_cols == 0) return numeric_limits<double>::quiet_NaN();
  if (_cols <  2) return .0;
  double mean = .0, var = .0, delta;
  for (int c = 0; c < _cols; ++c) {
    delta = _matrix[c][r] - mean;
    mean += delta / (c + 1);
    var  += delta * (_matrix[c][r] - mean);
  }
  var /= (_cols - 1);
  return var;
}

// -----------------------------------------------------------------------------
double Matrix::RowStd(int r) const
{
  return sqrt(RowVar(r));
}

// -----------------------------------------------------------------------------
double Matrix::ColSum(int c) const
{
  double sum = .0, *d = _matrix[c];
  for (int r = 0; r < _rows; ++r, ++d) {
    sum += *d;
  }
  return sum;
}

// -----------------------------------------------------------------------------
double Matrix::ColMean(int c) const
{
  return ColSum(c) / _rows;
}

// -----------------------------------------------------------------------------
double Matrix::ColVar(int c) const
{
  if (_rows == 0) return numeric_limits<double>::quiet_NaN();
  if (_rows <  2) return .0;
  double mean = .0, var = .0, delta, *d = _matrix[c];
  for (int r = 0; r < _rows; ++r, ++d) {
    delta = *d - mean;
    mean += delta / (r + 1);
    var  += delta * (*d - mean);
  }
  var /= (_rows - 1);
  return var;
}

// -----------------------------------------------------------------------------
double Matrix::ColStd(int c) const
{
  return sqrt(ColVar(c));
}

// -----------------------------------------------------------------------------
Matrix Matrix::Exp() const
{
  // Matrix exponential via Pade approximation.
  // See Golub and Van Loan, Matrix Computations, Algorithm 11.3-1.

  // Number of iterations (determines accuracy).
  const int q = 6;

  // First this matrix is scaled by a power of 2 so its norm is < 1/2.
  // j is the index of the power required.
  double norm = InfinityNorm();
  int e = (int) ceil(log(norm) / log(2.0));
  int j = max(0, 1 + e);

  Matrix A = (*this) / (pow(2.0, j));
  Matrix D(_rows, _cols);
  Matrix N(_rows, _cols);
  Matrix X(_rows, _cols);

  D.Ident();
  N.Ident();
  X.Ident();

  double c             = 1.0;
  int    minusOnePower = 1;
  for (int k = 1; k <= q; ++k) {
    c = c * (q - k + 1) / ((double) k * (2*q - k + 1));
    X = A * X;
    N = N +  X * c;
    minusOnePower *= -1;
    D = D + X * minusOnePower * c;
  }

  D.Invert();
  X = D * N;

  // Squaring steps
  for (int k = 1; k <= j; ++k) X = X * X;
  return X;
}

// -----------------------------------------------------------------------------
Matrix Matrix::Log() const
{
  MIRTK_CHECK_MATRIX_IS_SQUARE("Log");

  Vector eigval;
  Matrix eigvec, tmp;
  Eigenvalues(eigvec, eigval, tmp);
  for (int i = 0; i < eigval.Rows(); ++i) {
    if (eigval(i) <= 0) {
      cerr << "Matrix::Log: Found non-positive eigenvalue: e(" << (i+1) << ") = " << eigval(i) << endl;
      this->Print(2);
      exit(1);
    }
  }

  const int    maxit = 100;
  const double tol   = 0.00001;
  int          i, k, n;

  Matrix A(*this);
  Matrix I(_rows, _cols);
  Matrix Z(_rows, _cols);
  Matrix X(_rows, _cols);
  Matrix D(_rows, _cols);

  I.Ident();

  D = A - I;
  k = 0, n = 0;
  while (D.InfinityNorm() > 0.5 && n < maxit) {
    A = A.Sqrt(), ++k;
    D = A - I,    ++n;
  }

  A = X = Z = I - A;
  i = 1, n = 0;
  while (Z.InfinityNorm() > tol && n < maxit) {
    Z = Z * A,       ++i;
    X = X + (Z / i), ++n;
  }

  X = X * -pow(2.0, k);
  return X;
}

// -----------------------------------------------------------------------------
Matrix Matrix::Sqrt() const
{
  const int    maxit = 100;
  const double tol   = 0.0001;

  Matrix X   (*this);
  Matrix Y   (_rows, _cols);
  Matrix D   (_rows, _cols);
  Matrix invX(_rows, _cols);
  Matrix invY(_rows, _cols);

  Y.Ident();
  D = X * X - (*this);

  for (int i = 0; i < maxit && D.InfinityNorm() > tol; ++i) {
    invX = X.Inverse();
    invY = Y.Inverse();
    X = (X + invY) * 0.5;
    Y = (Y + invX) * 0.5;
    D = X * X - (*this);
  }

  return X;
}

// -----------------------------------------------------------------------------
double Matrix::Det() const
{
  MIRTK_CHECK_MATRIX_IS_SQUARE("Det");

  double d, sign;
  int i;
  Matrix lu, perm;

  this->LU(lu, perm, sign);

  d = sign;
  for (i = 0; i < _rows; i++) {
    d *= lu(i, i);
  }

  return d;
}

// -----------------------------------------------------------------------------
Matrix &Matrix::Invert(bool use_svd_if_singular)
{
  MIRTK_CHECK_MATRIX_IS_SQUARE("Invert");

  double sign;
  int i, j, p, q;
  Matrix lu, perm, A, inv;
  Vector b, x, y;

  if (this->Det() == 0) {
    if (use_svd_if_singular) return this->SVDInvert();
    cerr << "Matrix::Invert: Zero determinant, matrix is singular" << endl;
    exit(1);
  }

  this->LU(lu, perm, sign);

  A.Initialize(_rows, _rows);
  b.Initialize(_rows);
  x.Initialize(_rows);
  y.Initialize(_rows);

  for (j = 0; j < _rows; j++) {
    for (i = 0; i < _rows; i++) {
      b(i) = 0.0;
    }
    b(j) = 1.0;

    // Forward substitution
    for (p = 0; p < _rows; p++)
    {
      y(p) = b(p);
      for (q = 0; q < p; q++)
      {
        y(p) -= lu(p, q) * y(q);
      }
    }

    // Back substitution
    for (p = _rows - 1; p >= 0; p--)
    {
      x(p) = y(p);
      for (q = p + 1; q < _rows; q++)
      {
        x(p) -= lu(p, q) * x(q);
      }
      x(p) /= lu(p, p);
    }

    // Write column j
    for (i = 0; i < _rows; i++) {
      A(i, j) = x(i);
    }
  }

  // Multiply with permutation matrix to get actual inverse
  inv = A * perm;

  // Copy into *this and return
  return (*this = inv);
}

// -----------------------------------------------------------------------------
Matrix &Matrix::SVDInvert()
{
  MIRTK_CHECK_MATRIX_IS_SQUARE("SVDInvert");

  Matrix U, V;
  Vector w;

  this->SVD(U, w, V);

  double wmax = .0;
  for (int r = 0; r < w.Rows(); ++r) {
    wmax = max(wmax, abs(w(r)));
  }
  double wmin = 1e-4 * wmax;
  for (int r = 0; r < w.Rows(); ++r) {
    if (abs(w(r)) <= wmin) w(r) = .0;
    else                   w(r) = 1.0 / w(r);
  }

  return ((*this) = Matrix(V * w) * U.Transpose());
}

// -----------------------------------------------------------------------------
Matrix &Matrix::PseudoInvert()
{
  if (_rows < _cols) {
    return Transpose().PseudoInvert().Transpose();
  }
  Matrix       &m = *this;
  const Matrix mT = m.Transposed();
  return ((*this) = (mT * m).Invert() * mT);
}

// -----------------------------------------------------------------------------
Matrix &Matrix::Adjugate(double &d)
{
  MIRTK_CHECK_MATRIX_IS_SQUARE("Adjugate");

  double sign;
  int i, j, p, q;
  Matrix lu, perm, A, inv;
  Vector b, x, y;

  d = this->Det();

  if (d == 0) {
    cerr << "Matrix::Invert: Zero determinant" << endl;
    exit(1);
  }

  this->LU(lu, perm, sign);

  A.Initialize(_rows, _rows);
  b.Initialize(_rows);
  x.Initialize(_rows);
  y.Initialize(_rows);

  for (j = 0; j < _rows; j++) {
    for (i = 0; i < _rows; i++) {
      b(i) = 0.0;
    }
    b(j) = 1.0;

    // Forward substitution
    for (p = 0; p < _rows; p++)
    {
      y(p) = b(p);
      for (q = 0; q < p; q++)
      {
        y(p) -= lu(p, q) * y(q);
      }
    }

    // Back substitution
    for (p = _rows - 1; p >= 0; p--)
    {
      x(p) = y(p);
      for (q = p + 1; q < _rows; q++)
      {
        x(p) -= lu(p, q) * x(q);
      }
      x(p) /= lu(p, p);
    }

    // Write column j
    for (i = 0; i < _rows; i++) {
      A(i, j) = x(i) * d;
    }
  }

  // Multiply with permutation matrix to get actual inverse
  inv = A * perm;

  // Copy into *this and return
  return (*this = inv);
}

// -----------------------------------------------------------------------------
Matrix &Matrix::Transpose()
{
  if (_rows == _cols) {
    for (int i = 0; i < _rows; ++i)
    for (int j = 0; j < i;     ++j) {
      swap(_matrix[i][j], _matrix[j][i]);
    }
  } else {
    Matrix tmp(_rows, _cols);
    for (int j = 0; j < _cols; ++j)
    for (int i = 0; i < _rows; ++i) {
      tmp._matrix[j][i] = _matrix[j][i];
    }
    Deallocate(_matrix);
    Allocate(_matrix, tmp._cols, tmp._rows);
    _rows = tmp._cols;
    _cols = tmp._rows;
    for (int j = 0; j < _cols; ++j)
    for (int i = 0; i < _rows; ++i) {
      _matrix[j][i] = tmp._matrix[i][j];
    }
  }

  return *this;
}

// -----------------------------------------------------------------------------
Matrix &Matrix::PermuteRows(Array<int> idx)
{
  mirtkAssert(idx.size() <= static_cast<size_t>(_cols), "valid permutation");
  for (int r1 = 0; r1 < static_cast<int>(idx.size()); ++r1) {
    int r2 = idx[r1];
    if (r2 == r1) continue;
    for (int c = 0; c < _cols; ++c) swap(_matrix[c][r1], _matrix[c][r2]);
    for (int r = r1 + 1; r < _rows; ++r) {
      if (idx[r] == r1) {
        swap(idx[r], idx[r1]);
        break;
      }
    }
  }
  return *this;
}

// -----------------------------------------------------------------------------
Matrix &Matrix::PermuteCols(Array<int> idx)
{
  mirtkAssert(idx.size() <= static_cast<size_t>(_cols), "valid permutation");
  for (int c1 = 0; c1 < static_cast<int>(idx.size()); ++c1) {
    int c2 = idx[c1];
    if (c2 == c1) continue;
    for (int r = 0; r < _rows; ++r) swap(_matrix[c1][r], _matrix[c2][r]);
    for (int c = c1 + 1; c < _cols; ++c) {
      if (idx[c] == c1) {
        swap(idx[c], idx[c1]);
        break;
      }
    }
  }
  return *this;
}

// -----------------------------------------------------------------------------
void Matrix::LU(Matrix &lu, Matrix &perm, double &sign) const
{
  MIRTK_CHECK_MATRIX_IS_SQUARE("LU");
  Eigen::PartialPivLU<Eigen::MatrixXd> ppiv_lu(MatrixToEigen(*this));
  lu   = EigenToMatrix(ppiv_lu.matrixLU());
  perm = EigenToMatrix(ppiv_lu.permutationP());
  // Trick to determine the sign (either 1 or -1) of the LU decomposition
  // since the PartialPivLU class interface does not provide a way to get it
  sign = ppiv_lu.matrixLU().diagonal().prod();
  if (AreEqual(sign, 0.)) {
    sign = 1.;
  } else {
    sign = ppiv_lu.determinant() / sign;
  }
}

// -----------------------------------------------------------------------------
void Matrix::SVD(Matrix &u, Vector &w, Matrix &v) const
{
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(MatrixToEigen(*this),
                                        Eigen::ComputeThinU | Eigen::ComputeThinV);
  u = EigenToMatrix(svd.matrixU());
  v = EigenToMatrix(svd.matrixV());
  w = EigenToVector(svd.singularValues());
}

// -----------------------------------------------------------------------------
bool Matrix::Eigenvalues(Matrix &E1, Vector &e, Matrix &E2) const
{
  bool ok = true;

  if (this->IsSymmetric()) {

    this->SymmetricEigen(E1, e);
    E2 = E1.Inverse();

  } else {

    // TODO: Consider use of Eigen::EigenSolver instead.
    //       http://eigen.tuxfamily.org/dox/classEigen_1_1EigenSolver.html

    this->SVD(E1, e, E2);
    for (int i = 0; i < e.Rows(); ++i) e(i) *= e(i);

    // Singular value decomposition differs from eigenvalue decomposition
    // if matrix is not square and diagonalizable
    ok = this->IsDiagonalizable();

  }

  return ok;
}

// -----------------------------------------------------------------------------
void Matrix::SymmetricEigen(Matrix &E, Vector &e) const
{
  MIRTK_CHECK_MATRIX_IS_SQUARE("SymmetricEigen");
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig_solver(MatrixToEigen(*this));
  E = EigenToMatrix(eig_solver.eigenvectors());
  e = EigenToVector(eig_solver.eigenvalues());
}

// -----------------------------------------------------------------------------
void Matrix::LeastSquaresFit(const Vector &y, Vector &x) const
{
  Matrix u, v, m_inv;
  Vector w;

  // nmatrix should be rows
  // ma should be cols

  if (y.Rows() != _rows) {
    cerr << "Matrix::LeastSquaresFit: Y has wrong dimensions" << endl;
    exit(1);
  }

  if (x.Rows() != _cols) {
    cerr << "Matrix::LeastSquaresFit: X has wrong dimensions" << endl;
    exit(1);
  }

  // Calculate least squares fit via SVD
  this->SVD(u, w, v);

  double wmax = .0;
  for (int j = 0; j < _cols; ++j) {
    if (w(j) > wmax) wmax = w(j);
  }
  const double thresh = EIGEN_TOL * wmax;
  for (int j = 0; j < _cols; ++j) {
    if (w(j) < thresh) w(j) = .0;
  }

  // Compute u_new = {w}^-1 * u^T
  u.Transpose();
  for (int i = 0; i < _cols; ++i)
  for (int j = 0; j < _rows; ++j) {
    if (w(i) != 0.0) {
      u(i, j) /= w(i);
    } else {
      u(i, j) = 0.0;
    }
  }

  // Compute x = V * u_new * y
  x = (v * u) * y;
}

// -----------------------------------------------------------------------------
Matrix &Matrix::Ident()
{
  for (int c = 0; c < _cols; ++c)
  for (int r = 0; r < _rows; ++r) {
    _matrix[c][r] = static_cast<double>(r == c);
  }
  return *this;
}

// -----------------------------------------------------------------------------
bool Matrix::IsIdentity() const
{
  if (_rows == 0 || _cols == 0 || _rows != _cols) return false;

  for (int c = 0; c < _cols; ++c)
  for (int r = 0; r < _rows; ++r) {
    if (!AreEqual(_matrix[c][r], r == c ? 1. : 0., 1e-6)) return false;
  }

  return true;
}

// -----------------------------------------------------------------------------
bool Matrix::IsSquare() const
{
  return _rows == _cols;
}

// -----------------------------------------------------------------------------
bool Matrix::IsSymmetric() const
{
  if (_rows == 0 || _cols == 0 || _rows != _cols) return false;
  for (int c = 0; c < _cols; ++c)
  for (int r = 0; r < _rows; ++r) {
    if (!AreEqual(_matrix[c][r], _matrix[r][c], 1e-6)) return false;
  }
  return true;
}

// -----------------------------------------------------------------------------
bool Matrix::IsDiagonalizable() const
{
  if (_rows == 0 || _cols == 0 || _rows != _cols) return false;
  Matrix AT  = (*this); AT.Transpose();
  Matrix ATA = AT * (*this);
  Matrix AAT = (*this) * AT;
  for (int c = 0; c < _cols; ++c)
  for (int r = 0; r < _rows; ++r) {
    if (!AreEqual(AAT(r, c), ATA(r, c), 1e-6)) return false;
  }
  return true;
}

// -----------------------------------------------------------------------------
void Matrix::MakeSymmetric()
{
  if (_rows != _cols) {
    cerr << "Matrix::MakeSymmetric: Matrix must be square" << endl;
    exit(1);
  }
  for (int c = 0; c != _cols; ++c)
  for (int r = 0; r != _rows; ++r) {
    _matrix[c][r] = _matrix[r][c] = 0.5 * (_matrix[c][r] + _matrix[r][c]);
  }
}

// =============================================================================
// I/O
// =============================================================================

// -----------------------------------------------------------------------------
ostream& operator <<(ostream& os, const Matrix &m)
{
  const int n = m._rows * m._cols;
  os << "irtkMatrix " << m._rows << " x " << m._cols << endl;
  double *data  = new double[n];
  int index = 0;
  for (int j = 0; j < m._cols; ++j)
  for (int i = 0; i < m._rows; ++i, ++index) {
    data[index] = m._matrix[j][i];
  }
  if (GetByteOrder() == LittleEndian) {
    swap64((char *)data, (char *)data, n);
  }
  os.write((char *)data, n * sizeof(double));
  delete[] data;
  return os;
}

// -----------------------------------------------------------------------------
istream& operator >>(istream& is, Matrix &m)
{
  // Read keyword
  char keyword[11];
  is.read(keyword, 11);
  if (strncmp(keyword, "irtkMatrix ", 11) != 0) {
    cerr << "Matrix: File does not appear to be a binary (M)IRTK Matrix file" << endl;
    exit(1);
  }

  // Read matrix size
  int cols = 0, rows = 0;
  if (!(is >> rows) || is.get() != ' ' || is.get() != 'x' || is.get() != ' ' || !(is >> cols) || cols <= 0 || rows <= 0) {
    cerr << "Matrix: Could not read matrix size from binary file" << endl;
    exit(1);
  }
  const int n = rows * cols;

  // Skip remaining characters (e.g., comment) on first line
  char buffer[256];
  while (is.get() != '\n')  {
    is.get(buffer, 256);
  }

  // Read binary data
  double *data  = new double[n];
  is.read((char *)data, n * sizeof(double));
  if (is.fail()) {
    delete[] data;
    cerr << "Matrix: File contains fewer values than expected: #rows = " << rows << ", #cols = " << cols << endl;
    exit(1);
  }
  if (GetByteOrder() == LittleEndian) {
    swap64((char *)data, (char *)data, n);
  }

  // Convert data
  m.Initialize(rows, cols);
  int index = 0;
  for (int j = 0; j < m._cols; ++j)
  for (int i = 0; i < m._rows; ++i, ++index) {
    m._matrix[j][i] = data[index];
  }
  delete[] data;

  return is;
}

// -----------------------------------------------------------------------------
Cofstream& operator <<(Cofstream &to, const Matrix &m)
{
  to.WriteAsChar("irtkMatrix", 11);
  to.WriteAsInt(&m._rows, 1);
  to.WriteAsInt(&m._cols, 1);
  to.WriteAsDouble(m.RawPointer(), m._rows * m._cols);
  return to;
}

// -----------------------------------------------------------------------------
Cifstream& operator >>(Cifstream &from, Matrix &m)
{
  char keyword[11];
  if (!from.ReadAsChar(keyword, 11) || strncmp(keyword, "irtkMatrix", 11) != 0) {
    cerr << "Matrix: File does not appear to be a binary (M)IRTK Matrix file" << endl;
    exit(1);
  }

  int rows = 0, cols = 0;
  if (!from.ReadAsInt(&rows, 1) || !from.ReadAsInt(&cols, 1) || rows <= 0 || cols <= 0) {
    cerr << "Matrix: Could not read matrix size from binary file" << endl;
    exit(1);
  }

  m.Initialize(rows, cols);
  if (!from.ReadAsDouble(m.RawPointer(), rows * cols)) {
    cerr << "Matrix: File contains fewer values than expected: #rows = " << rows << ", #cols = " << cols << endl;
    exit(1);
  }

  return from;
}

// -----------------------------------------------------------------------------
void Matrix::Print(Indent indent) const
{
  Print(cout, indent);
}

// -----------------------------------------------------------------------------
void Matrix::Print(ostream &os, Indent indent) const
{
  os << indent << "Matrix " << _rows << " x " << _cols << endl;
  ++indent;
  os.setf(ios::right);
  os.setf(ios::fixed);
  os.precision(4);
  for (int i = 0; i < _rows; i++) {
    os << indent;
    for (int j = 0; j < _cols; j++) {
      os << setw(15) << _matrix[j][i] << " ";
    }
    os << endl;
  }
  os.precision(6);
  os.unsetf(ios::right);
  os.unsetf(ios::fixed);
}

// -----------------------------------------------------------------------------
void Matrix::Read(const char *filename)
{
  // Open file stream
  ifstream from(filename, ios::in | ios::binary);

  // Check whether file opened ok
  if (!from) {
    cerr << "Matrix::Read: Can't open file " << filename << endl;
    exit(1);
  }

  // Read matrix
  from >> *this;
}

// -----------------------------------------------------------------------------
void Matrix::Write(const char *filename) const
{
  // Open file stream
  ofstream to(filename, ios::out | ios::binary);

  // Check whether file opened ok
  if (!to) {
    cerr << "Matrix::Write: Can't open file " << filename << endl;
    exit(1);
  }

  // Write matrix
  to << *this;
}

// -----------------------------------------------------------------------------
void Matrix::Import(const char *filename, int rows, int cols)
{
  // Open file stream
  ifstream from(filename);

  // Check whether file opened ok
  if (!from) {
    cerr << "Matrix::Read: Can't open file " << filename << endl;
    exit(1);
  }

  // Initialize matrix
  this->Initialize(rows, cols);

  // Read matrix
  for (int i = 0; i < _rows; ++i)
  for (int j = 0; j < _cols; ++j) {
    from >> _matrix[j][i];
  }
}

#if MIRTK_Numerics_WITH_MATLAB

// -----------------------------------------------------------------------------
inline mxArray *MatrixToMxArray(const Matrix &m)
{
  mxArray *mx = mxCreateDoubleMatrix(m.Rows(), m.Cols(), mxREAL);
  memcpy(mxGetPr(mx), m.RawPointer(), m.NumberOfElements() * sizeof(double));
  return mx;
}

// -----------------------------------------------------------------------------
bool Matrix::WriteMAT(const char *fname, const char *varname) const
{
  Matlab::Initialize();
  MATFile *fp = matOpen(fname, "w");
  if (fp == NULL) return false;
  mxArray *m = MatrixToMxArray(*this);
  if (matPutVariable(fp, varname, m) != 0) {
    mxDestroyArray(m);
    matClose(fp);
    return false;
  }
  mxDestroyArray(m);
  return (matClose(fp) == 0);
}

#endif // MIRTK_Numerics_WITH_MATLAB

////////////////////////////////////////////////////////////////////////////////
// Means
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
Matrix LogEuclideanMean(int n, const Matrix *matrices, const double *weights)
{
  mirtkAssert(n > 0, "at least one matrix must be given");
  if (n == 1) return matrices[0];

  double totalWeight;
  if (weights) {
    totalWeight = .0;
    for (int i = 0; i < n; ++i) {
      totalWeight += weights[i];
    }
  } else {
    totalWeight = static_cast<double>(n);
  }
  if (totalWeight <= .0) {
    cerr << "itkMatrix::LogEuclideanMean: Sum of weights must be positive" << endl;
    exit(1);
  }

  Matrix sumLogs(matrices[0].Rows(), matrices[0].Cols());
  if (weights) {
    for (int i = 0; i < n; ++i) {
      sumLogs += matrices[i].Log() * weights[i];
    }
  } else {
    for (int i = 0; i < n; ++i) {
      sumLogs += matrices[i].Log();
    }
  }
  sumLogs /= totalWeight;

  return sumLogs.Exp();
}

// -----------------------------------------------------------------------------
Matrix BiInvariantMean(int           n,
                       const Matrix *matrices,
                       const double *weights,
                       int           niter,
                       double        tolerance,
                       const Matrix *mu0)
{
  mirtkAssert(n > 0, "at least one matrix must be given");
  if (n == 1) return matrices[0];

  int    i, iter = 0;
  double totalWeight, normLogDeltaMu;
  Matrix mu, muInv, deltaMu, deltaM, sumLogs;

  // Normalize weights
  if (weights) {
    totalWeight = .0;
    for (i = 0; i < n; ++i) {
      totalWeight += weights[i];
    }
  } else {
    totalWeight = static_cast<double>(n);
  }
  if (totalWeight <= .0) {
    cerr << "itkMatrix::BiInvariantMean: Sum of weights must be positive" << endl;
    exit(1);
  }

  // Check that determinants of all matrices have the same sign
  double sign = sgn(matrices[0].Det());
  for (i = 1; i < n; ++i) {
    if (sgn(matrices[i].Det()) != sign) {
      cerr << "Matrix::BiInvariantMean: Sign of determinant is not the same for all matrices" << endl;
      exit(1);
    }
  }

  // Barycentric fixed point iteration
  mu = (mu0 ? *mu0 : matrices[0]);

  do {

    muInv = mu.Inverse();

    sumLogs.Initialize(mu.Rows(), mu.Cols()); // Reset sum to zero
    if (weights) {
      for (i = 0; i < n; ++i) {
        if (weights[i] == .0) continue;
        deltaM = muInv * matrices[i];
        sumLogs += deltaM.Log() * weights[i];
      }
    } else {
      for (i = 0; i < n; ++i) {
        deltaM = muInv * matrices[i];
        sumLogs += deltaM.Log();
      }
    }
    sumLogs /= totalWeight;
    deltaMu = sumLogs.Exp();

    mu = mu * deltaMu;

    normLogDeltaMu = deltaMu.Log().InfinityNorm();

  } while (normLogDeltaMu > tolerance && ++iter < niter);

  return mu;
}

////////////////////////////////////////////////////////////////////////////////
// Affine transformation matrices
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
// R = (Rz Ry Rx)^T
void RigidParametersToMatrix(double tx,    double ty,    double tz,
                             double cosrx, double cosry, double cosrz,
                             double sinrx, double sinry, double sinrz, Matrix &m)
{
  m.Initialize(4, 4);

  m(0, 0) = cosry * cosrz;
  m(0, 1) = cosry * sinrz;
  m(0, 2) = -sinry;
  m(0, 3) = tx;

  m(1, 0) = (sinrx * sinry * cosrz - cosrx * sinrz);
  m(1, 1) = (sinrx * sinry * sinrz + cosrx * cosrz);
  m(1, 2) = sinrx * cosry;
  m(1, 3) = ty;

  m(2, 0) = (cosrx * sinry * cosrz + sinrx * sinrz);
  m(2, 1) = (cosrx * sinry * sinrz - sinrx * cosrz);
  m(2, 2) = cosrx * cosry;
  m(2, 3) = tz;

  m(3, 0) = 0.0;
  m(3, 1) = 0.0;
  m(3, 2) = 0.0;
  m(3, 3) = 1.0;
}

// -----------------------------------------------------------------------------
void RigidParametersToMatrix(double tx,  double ty,  double tz,
                             double rx,  double ry,  double rz, Matrix &m)
{
  const double cosrx = cos(rx);
  const double cosry = cos(ry);
  const double cosrz = cos(rz);
  const double sinrx = sin(rx);
  const double sinry = sin(ry);
  const double sinrz = sin(rz);

  RigidParametersToMatrix(tx, ty, tz, cosrx, cosry, cosrz, sinrx, sinry, sinrz, m);
}

// -----------------------------------------------------------------------------
void MatrixToEulerAngles(const Matrix &m, double &rx,  double &ry,  double &rz)
{
  const double TOL = 0.000001;
  double tmp;

  tmp = asin(-1 * m(0, 2));

  // asin returns values for tmp in range -pi/2 to +pi/2, i.e. cos(tmp) >=
  // 0 so the division by cos(tmp) in the first part of the if clause was
  // not needed.
  if (fabs(cos(tmp)) > TOL) {
    rx = atan2(m(1,2), m(2,2));
    ry = tmp;
    rz = atan2(m(0,1), m(0,0));
  } else {
    //m(0,2) is close to +1 or -1
    rx = atan2(-1.0*m(0,2)*m(1,0), -1.0*m(0,2)*m(2,0));
    ry = tmp;
    rz = 0;
  }
}

// -----------------------------------------------------------------------------
void MatrixToRigidParameters(const Matrix &m,
                             double &tx,  double &ty,  double &tz,
                             double &rx,  double &ry,  double &rz)
{
  const double TOL = 0.000001;
  double tmp;

  tx = m(0, 3);
  ty = m(1, 3);
  tz = m(2, 3);

  tmp = asin(-1 * m(0, 2));

  // asin returns values for tmp in range -pi/2 to +pi/2, i.e. cos(tmp) >=
  // 0 so the division by cos(tmp) in the first part of the if clause was
  // not needed.
  if (fabs(cos(tmp)) > TOL) {
    rx = atan2(m(1,2), m(2,2));
    ry = tmp;
    rz = atan2(m(0,1), m(0,0));
  } else {
    //m(0,2) is close to +1 or -1
    rx = atan2(-1.0*m(0,2)*m(1,0), -1.0*m(0,2)*m(2,0));
    ry = tmp;
    rz = 0;
  }
}

// -----------------------------------------------------------------------------
void AffineParametersToMatrix(double tx,  double ty,  double tz,
                              double rx,  double ry,  double rz,
                              double sx,  double sy,  double sz,
                              double sxy, double sxz, double syz, Matrix &m)
{
  Matrix tmp(4, 4);

  // Construct rigid transformation matrix
  RigidParametersToMatrix(tx, ty, tz, rx, ry, rz, m);

  // Pre-multiply with shearing transformation
  tmp.Ident();
  tmp(0, 1) = tan(sxy);
  tmp(0, 2) = tan(sxz);
  tmp(1, 2) = tan(syz);
  m *= tmp;

  // Pre-multiply with scaling transformation
  tmp.Ident();
  tmp(0, 0) = sx;
  tmp(1, 1) = sy;
  tmp(2, 2) = sz;
  m *= tmp;
}

// -----------------------------------------------------------------------------
void MatrixToAffineParameters(const Matrix &m,
                              double &tx,  double &ty,  double &tz,
                              double &rx,  double &ry,  double &rz,
                              double &sx,  double &sy,  double &sz,
                              double &sxy, double &sxz, double &syz)
{
  const double TOL = 0.000001;
  double tansxy, tansxz, tansyz;

  if (abs(m(3, 3) - 1.0) > TOL) {
    cerr << "MatrixToAffineParameters: Value at m(3, 3) must equal 1." << endl;
    exit(1);
  }
  if (abs(m.Det()) < TOL) {
    cerr << "MatrixToAffineParameters: Matrix is singular or very close to singular!." << endl;
    exit(1);
  }

  // First Part Of Graphics Gems Code Ignored Because It Relates To
  // Perspective Transformation.
  if (abs(m(3, 0)) > TOL ||
      abs(m(3, 1)) > TOL ||
      abs(m(3, 2)) > TOL ) {
    cerr << "MatrixToAffineParameters: Matrix contains perspective distortion." << endl;
    exit(1);
  }

  Matrix copy(m);

  // Get scale and shear by manipulating the columns of the upper left 3x3
  // sub-matrix.
  Vector col_0, col_1, col_2;
  col_0.Initialize(3);
  col_1.Initialize(3);
  col_2.Initialize(3);
  for (int i = 0; i < 3; ++i) {
    col_0(i) = copy(i, 0);
    col_1(i) = copy(i, 1);
    col_2(i) = copy(i, 2);
  }

  // Compute X scale factor and normalize first col.
  sx = col_0.Norm();
  col_0 /= sx;

  // Compute XY shear factor and make 2nd col orthogonal to 1st.
  tansxy = col_0.ScalarProduct(col_1);
  col_1 = col_1 - col_0 * tansxy;

  // Actually, tansxy and col_1 are still to large by a factor of sy.
  // Now, compute Y scale and normalize 2nd col and rescale tansxy.
  sy = col_1.Norm();
  col_1  /= sy;
  tansxy /= sy;

  // Compute XZ and YZ shears, orthogonalize 3rd col
  tansxz = col_0.ScalarProduct(col_2);
  col_2 = col_2 - col_0 * tansxz;

  tansyz = col_1.ScalarProduct(col_2);
  col_2 = col_2 - col_1 * tansyz;

  // Actually, tansxz, tansyz and col_2 are still to large by a factor of
  // sz.  Next, get Z scale, normalize 3rd col and scale tansxz and tansyz.
  sz = col_2.Norm();
  col_2  /= sz;
  tansxz /= sz;
  tansyz /= sz;

  // At this point, the columns are orthonormal.  Check for a coordinate
  // system flip.  If the determinant is -1, then negate the matrix and the
  // scaling factors.
  Vector col_1_x_col_2;
  col_1_x_col_2.Initialize(3);
  col_1_x_col_2 = col_1.CrossProduct(col_2);

  if (col_0.ScalarProduct(col_1_x_col_2) < 0) {
    sx *= -1;
    sy *= -1;
    sz *= -1;
    col_0 *= -1;
    col_1 *= -1;
    col_2 *= -1;
  }

  // Retrieve the shear angles in degrees.
  sxy = atan(tansxy);
  sxz = atan(tansxz);
  syz = atan(tansyz);
  if (AreEqual(sxy, 0.)) sxy = 0.;
  if (AreEqual(sxz, 0.)) sxz = 0.;
  if (AreEqual(syz, 0.)) syz = 0.;

  // Now get the rigid transformation parameters.
  // Put the rotation matrix components into the upper left 3x3 submatrix.
  for (int i = 0; i < 3; ++i) {
    copy(i, 0) = col_0(i);
    copy(i, 1) = col_1(i);
    copy(i, 2) = col_2(i);
  }

  MatrixToRigidParameters(copy, tx, ty, tz, rx, ry, rz);
}

// -----------------------------------------------------------------------------
// Based on Insight Journal code: http://hdl.handle.net/10380/3299
//
// Solve Qa = C
//
//          | sum( q(1,i)*q(1,i) ) sum( q(1,i)*q(2,i) ) sum( q(1,i)*q(3,i) ) sum( q(1,i) ) |
// Q[4x4] = | sum( q(2,i)*q(1,i) ) sum( q(2,i)*q(2,i) ) sum( q(2,i)*q(3,i) ) sum( q(2,i) ) |
//          | sum( q(3,i)*q(1,i) ) sum( q(3,i)*q(2,i) ) sum( q(3,i)*q(3,i) ) sum( q(3,i) ) |
//          | sum( q(1,i) )        sum( q(2,i) )        sum( q(3,i) )        sum( q(4,i) ) |
//
//          | a(1,1) a(2,1) a(3,1) |
// a[4x3] = | a(1,2) a(2,2) a(3,2) |
//          | a(1,3) a(2,3) a(3,3) |
//          | a(1,4) a(2,4) a(3,4) |
//
//          | sum( q(1,i)*p(1,i) )  sum( q(1,i)*p(2,i) ) sum( q(1,i)*p(3,i) )  |
// C[4x3] = | sum( q(2,i)*p(1,i) )  sum( q(2,i)*p(2,i) ) sum( q(2,i)*p(3,i) )  |
//          | sum( q(3,i)*p(1,i) )  sum( q(3,i)*p(2,i) ) sum( q(3,i)*p(3,i) )  |
//          | sum( p(1,i) )         sum( p(2,i) )        sum( p(3,i) )         |
//
Matrix ApproximateAffineMatrix(const PointSet &target, const PointSet &source, const Vector &weight)
{
  // Check arguments
  const int no = source.Size();
  if (target.Size() != no) {
    cerr << "ApproximateAffineMatrix: Size mismatch between number of landark pairs!" << endl;
    exit(1);
  }
  if (no < 4) {
    cerr << "ApproximateAffineMatrix: Must have at least four points" << endl;
    exit(1);
  }

  double wnorm = .0;
  if (weight.Rows() != 0) {
    if (weight.Rows() != no) {
      cerr << "ApproximateAffineMatrix: Size mismatch between number of landmark pairs and weights" << endl;
      exit(1);
    }
    for (int i = 0; i < no; ++i) {
      wnorm += weight(i) * weight(i);
    }
    wnorm = sqrt(wnorm);
  }

  // Output matrix
  Matrix A(4, 4);
  A.Ident();

  typedef Eigen::Matrix<double, 4, Eigen::Dynamic> Matrix4xN;
  typedef Eigen::Matrix<double, 3, Eigen::Dynamic> Matrix3xN;
  typedef Eigen::Matrix<double, 4, 4>              Matrix4x4;
  typedef Eigen::Matrix<double, 4, 3>              Matrix4x3;
  typedef Eigen::Matrix<double, 4, 1>              Vector4;
  typedef Eigen::FullPivHouseholderQR<Matrix4x4>   DenseSolver;
 
  // Convert point sets to matrices and apply landmark weights
  Matrix4xN q(4, no);
  Matrix3xN p(3, no);

  double w;
  for (int i = 0; i < no; ++i) {
    // Normalized landmark weight
    w = (weight.Rows() != 0 ? weight(i) / wnorm : 1.0);
    // Target point
    q(0, i) = w * target(i)._x;
    q(1, i) = w * target(i)._y;
    q(2, i) = w * target(i)._z;
    q(3, i) = w;
    // Source point
    p(0, i) = w * source(i)._x;
    p(1, i) = w * source(i)._y;
    p(2, i) = w * source(i)._z;
  }

  // Solve Qa = C
  Matrix4x4 Q;
  Matrix4x3 C;
  Vector4   a;

  Q = Matrix4x4::Zero();
  C = Matrix4x3::Zero();
  for (int i = 0; i < no; ++i) {
    Q += q.col(i) * q.col(i).transpose();
    C += q.col(i) * p.col(i).transpose();
  }

  DenseSolver solver(Q);
  for (int i = 0; i < 3; ++i) {
    a = solver.solve(C.col(i));
    A(i, 0) = a(0), A(i, 1) = a(1), A(i, 2) = a(2), A(i, 3) = a(3);
  }

  return A;
}

// -----------------------------------------------------------------------------
Matrix ApproximateAffineMatrix(const PointSet &target, const PointSet &source)
{
  return ApproximateAffineMatrix(target, source, Vector());
}

// -----------------------------------------------------------------------------
void OrthoNormalize3x3(Matrix &mat)
{
  for (int i = 0; i < 3; i++) {
    Vector3D<double> vi(mat(0, i), mat(1, i), mat(2, i));
    // normalize column i
    vi.Normalize();
    // make other vectors linearly independent of column i
    for (int j = i + 1; j < 3; j++) {
      Vector3D<double> vj(mat(0, j), mat(1, j), mat(2, j));
      const double s = Vector3D<double>::DotProduct(vj, vi);
      mat(0, j) = vj._x - s * vi._x;
      mat(1, j) = vj._y - s * vi._y;
      mat(2, j) = vj._z - s * vi._z;
    }
  }
}


} // namespace mirtk
