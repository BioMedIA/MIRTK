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

#ifndef MIRTK_Matrix_H
#define MIRTK_Matrix_H

#include "mirtk/Object.h"

#include "mirtk/Math.h"
#include "mirtk/Array.h"
#include "mirtk/Pair.h"
#include "mirtk/Indent.h"
#include "mirtk/Cfstream.h"
#include "mirtk/Matrix3x3.h"


namespace mirtk {


// Forward declarations to reduce cyclic header dependencies
class PointSet;
class Vector;


/**
 * Dense matrix
 *
 * \sa SparseMatrix
 */

class Matrix : public Object
{
  mirtkObjectMacro(Matrix);

  // ---------------------------------------------------------------------------
  // Types
public:

  typedef double ElementType;

  // ---------------------------------------------------------------------------
  // Attributes
protected:

  /// Number of rows
  int _rows;

  /// Number of colums
  int _cols;

  /// Matrix values
  double **_matrix;

  /// Whether this instance owns the memory of the matrix elements
  bool _owner;

  // ---------------------------------------------------------------------------
  // Construction/destruction
public:

  /// Default constructor
  Matrix();

  /// Constructor for given number of rows and columns
  Matrix(int, int = -1, double * = NULL);

  /// Convert column vector to matrix
  Matrix(const Vector &);

  /// Convert point set to (n x 3) or (n x 2) matrix
  ///
  /// \param[in] twoD Discard z coordinate.
  Matrix(const PointSet &, bool twoD = false);

  /// Copy constructor
  Matrix(const Matrix &);

  /// Construct from 3x3 matrix
  explicit Matrix(const Matrix3x3 &);

  /// Destructor
  ~Matrix();

  /// Assign scalar value to all elements
  Matrix& operator =(double);

  /// Assignment operator
  Matrix& operator =(const Matrix &);

  /// Assignment operator
  Matrix& operator =(const Matrix3x3 &);

  /// Initialize matrix with number of rows and columns
  void Initialize(int, int = -1, double * = NULL);

  /// Resize matrix preserving existing entries
  void Resize(int, int = -1);

  /// Free memory
  void Clear();

  /// Set elements to zero
  void Zero();

  // ---------------------------------------------------------------------------
  // Indexing

  /// Get number of elements
  int NumberOfElements() const;

  /// Get matrix size in each dimension
  Pair<int, int> Size() const;

  /// Get number of rows
  int Rows() const;

  /// Get number of columns
  int Cols() const;

  /// Get linear index of element given its row and column indices
  int Index(int, int) const;

  /// Get row index of element given its linear index
  int RowIndex(int) const;

  /// Get column index of element given its linear index
  int ColIndex(int) const;

  /// Get row and column index of element given its linear index
  void SubIndex(int, int &, int &) const;

  /// Get row and column index of element given its linear index
  Pair<int, int> SubIndex(int) const;

  // ---------------------------------------------------------------------------
  // Element access

  /// Get pointer to linear memory which stores matrix elements in column-major order
  double *RawPointer(int r = 0, int c = 0);

  /// Get pointer to linear memory which stores matrix elements in column-major order
  const double *RawPointer(int r = 0, int c = 0) const;

  /// Get pointer to linear memory which stores matrix elements in column-major order
  /// \deprecated Use RawPointer instead.
  double *GetPointerToElements(int r = 0, int c = 0);

  /// Get pointer to linear memory which stores matrix elements in column-major order
  /// \deprecated Use RawPointer instead.
  const double *GetPointerToElements(int r = 0, int c = 0) const;

  /// Get pointer to matrix entries in specified column
  double *Col(int c);

  /// Get pointer to matrix entries in specified column
  const double *Col(int c) const;

  /// Get reference to element with specified linear index
  double &operator ()(int);

  /// Get const reference to element with specified linear index
  const double &operator ()(int) const;

  /// Get reference to element in specified row and column
  double &operator ()(int, int);

  /// Get const reference to element in specified row and column
  const double &operator ()(int, int) const;

  /// Set value of element with specified linear index
  void Put(int, double);

  /// Get value of element with specified linear index
  double Get(int) const;

  /// Set value of element in specified row and column
  void Put(int, int, double);

  /// Get value of element in specified row and column
  double Get(int, int) const;

  /// Get submatrix
  Matrix operator ()(int, int, int, int) const;

  /// Set submatrix
  void operator ()(Matrix &, int, int);

  // ---------------------------------------------------------------------------
  // Scalar matrix operations

  /// Add scalar to values in specified row
  Matrix &AddToRow(int, double);

  /// Add scalar to values in specified column
  Matrix &AddToCol(int, double);

  /// Scale row by a scalar
  Matrix &ScaleRow(int, double);

  /// Scale column by a scalar
  Matrix &ScaleCol(int, double);

  /// Subtraction of a double
  Matrix &operator -=(double);

  /// Addition of a double
  Matrix &operator +=(double);

  /// Multiplication with a double
  Matrix &operator *=(double);

  /// Division by a double
  Matrix &operator /=(double);

  /// Return result of subtraction of a double
  Matrix operator -(double) const;

  /// Return result of addition of a double
  Matrix operator +(double) const;

  /// Return result of multiplication with a double
  Matrix operator *(double) const;

  /// Return result of division by a double
  Matrix operator /(double) const;

  // ---------------------------------------------------------------------------
  // Matrix-vector operations

  /// Right-multiply matrix with vector
  Vector operator* (const Vector&) const;

  // ---------------------------------------------------------------------------
  // Matrix-matrix operations

  /// Matrix subtraction operator
  Matrix& operator -=(const Matrix&);

  /// Matrix addition operator
  Matrix& operator +=(const Matrix&);

  /// Matrix multiplication operator
  Matrix& operator *=(const Matrix&);

  /// Return result of matrix subtraction
  Matrix operator -(const Matrix&) const;

  /// Return result of matrix addition
  Matrix operator +(const Matrix&) const;

  /// Return result of matrix multiplication
  Matrix operator *(const Matrix&) const;

  // ---------------------------------------------------------------------------
  // Comparison operations

  /// Matrix equality
  bool operator ==(const Matrix &) const;

  /// Matrix inequality
  bool operator !=(const Matrix &) const;

  /// Find elements with specific value
  /// \returns Sorted linear indices of found elements
  Array<int> operator ==(double) const;

  /// Find elements with value not equal to specified value
  /// \returns Sorted linear indices of found elements
  Array<int> operator !=(double) const;

  /// Find elements with value less than specified value
  /// \returns Sorted linear indices of found elements
  Array<int> operator <(double) const;

  /// Find elements with value less or equal than specified value
  /// \returns Sorted linear indices of found elements
  Array<int> operator <=(double) const;

  /// Find elements with value greater than specified value
  /// \returns Sorted linear indices of found elements
  Array<int> operator >(double) const;

  /// Find elements with value greater or equal than specified value
  /// \returns Sorted linear indices of found elements
  Array<int> operator >=(double) const;

  // ---------------------------------------------------------------------------
  // Matrix functions

  /// Minimum value in specified row
  double RowMin(int) const;

  /// Maximum value in specified row
  double RowMax(int) const;

  /// Minimum and maximum value in specified row
  void RowRange(int, double &, double &) const;

  /// Sum of column values in specified row
  double RowSum(int) const;

  /// Mean of column values in specified row
  double RowMean(int) const;

  /// Variance of column values in specified row
  double RowVar(int) const;

  /// Standard deviation of column values in specified row
  double RowStd(int) const;

  /// Minimum value in specified column
  double ColMin(int) const;

  /// Maximum value in specified column
  double ColMax(int) const;

  /// Minimum and maximum value in specified column
  void ColRange(int, double &, double &) const;

  /// Sum of row values in specified column
  double ColSum(int) const;

  /// Mean of row values in specified column
  double ColMean(int) const;

  /// Variance of row values in specified column
  double ColVar(int) const;

  /// Standard deviation of row values in specified column
  double ColStd(int) const;

  /// Matrix exponential via Pade approximation
  /// (cf. Golub and Van Loan, Matrix Computations, Algorithm 11.3-1)
  Matrix Exp() const;

  /// Matrix logarithm
  Matrix Log() const;

  /// Matrix square root
  Matrix Sqrt() const;

  /// Calculate norm of matrix
  double Norm() const;

  /// Calculate trace of matrix
  double Trace() const;

  // The infinity norm is the maximum of the absolute value row sums.
  double InfinityNorm() const;

  /// Calculate determinant of matrix
  double Det() const;

  /// Calculate determinant of a 3x3 matrix
  double Det3x3() const;

  /// Get upper left 3x3 sub-matrix
  Matrix3x3 To3x3() const;

  /// Set to identity matrix
  Matrix &Ident();

  /// Returns true if the matrix is an identity matrix.
  bool IsIdentity() const;

  /// Whether matrix is square
  bool IsSquare() const;

  /// Whether matrix is symmetric
  bool IsSymmetric() const;

  /// Whether matrix is diagonalizable
  bool IsDiagonalizable() const;

  /// Make square matrix symmetric by adding its transpose and divide by 2
  void MakeSymmetric();

  /// Invert matrix
  Matrix &Invert(bool use_svd_if_singular = true);

  /// Get inverse matrix
  Matrix Inverse(bool use_svd_if_singular = true) const;

  /// Invert (singular) matrix using SVD
  Matrix &SVDInvert();

  /// Get inverse of (singular) matrix using SVD
  Matrix SVDInverse() const;

  /// Invert matrix using pseudo inverse
  Matrix &PseudoInvert();

  /// Get pseudo inverse of matrix
  Matrix PseudoInverse() const;

  /// Matrix inversion operator
  Matrix operator !();

  /// Adjugate matrix and return determinant
  Matrix &Adjugate(double &);

  /// Transpose matrix
  Matrix &Transpose();

  /// Get transposed matrix
  Matrix Transposed() const;

  /// Matrix transpose operator
  Matrix operator ~();

  /// Permute rows
  Matrix &PermuteRows(Array<int>);

  /// Permute columns
  Matrix &PermuteCols(Array<int>);

  /// Calculate LU decomposition of square matrix
  void LU(Matrix &, Matrix &, double &) const;

  /// Calculate singular value decomposition
  void SVD(Matrix &, Vector &, Matrix &) const;

  /// Calculate eigenvalues and eigenvectors of matrix
  ///
  /// This function chooses the appropriate algorithm
  /// according to matrix properties.
  ///
  /// \returns Whether eigenvalue decomposition exists.
  ///          If \c false, squared singular values are returned.
  bool Eigenvalues(Matrix &, Vector &, Matrix &) const;

  /// Calculate eigendecomposition of symmetric matrix
  void SymmetricEigen(Matrix &, Vector &) const;

  /// Calculate least square fit via SVD
  void LeastSquaresFit(const Vector &, Vector &) const;

  // ---------------------------------------------------------------------------
  // I/O

  /// Interface to output stream
  friend ostream& operator<< (ostream&, const Matrix&);

  /// Interface to input stream
  friend istream& operator>> (istream&, Matrix&);

  /// Interface to output stream
  friend Cofstream& operator<< (Cofstream&, const Matrix&);

  /// Interface to input stream
  friend Cifstream& operator>> (Cifstream&, Matrix&);

  /// Print matrix
  void Print(Indent = 0) const;

  /// Print matrix
  void Print(ostream &, Indent = 0) const;

  /// Read matrix from file
  void Read(const char *);

  /// Write matrix to file
  void Write(const char *) const;

  /// Import matrix from text file (requires no. of expected rows and cols)
  void Import(const char *, int, int);

  /// Write dense matrix to MAT-file
  ///
  /// \note Use only when MIRTK_Numerics_WITH_MATLAB is 1.
  bool WriteMAT(const char *, const char * = "A") const;

};


} // namespace mirtk

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

#include "mirtk/Point.h"


namespace mirtk {

// =============================================================================
// Construction
// =============================================================================

// -----------------------------------------------------------------------------
inline Matrix::Matrix(const Matrix3x3 &m)
:
  _rows(3),
  _cols(3),
  _owner(true)
{
  Allocate(_matrix, 3, 3);
  for (int c = 0; c < 3; ++c)
  for (int r = 0; r < 3; ++r) {
    _matrix[c][r] = m[r][c];
  }
}

// -----------------------------------------------------------------------------
inline Matrix &Matrix::operator =(const Matrix3x3 &m)
{
  if (_rows != 3 || _cols != 3) {
    Clear();
    Allocate(_matrix, 3, 3);
    _rows = _cols = 3;
    _owner = true;
  }
  for (int c = 0; c < 3; ++c)
  for (int r = 0; r < 3; ++r) {
    _matrix[c][r] = m[r][c];
  }
  return *this;
}

// =============================================================================
// Indexing
// =============================================================================

// -----------------------------------------------------------------------------
inline int Matrix::NumberOfElements() const
{
  return _rows * _cols;
}

// -----------------------------------------------------------------------------
inline Pair<int, int> Matrix::Size() const
{
  return MakePair(_rows, _cols);
}

// -----------------------------------------------------------------------------
inline int Matrix::Rows() const
{
  return _rows;
}

// -----------------------------------------------------------------------------
inline int Matrix::Cols() const
{
  return _cols;
}

// -----------------------------------------------------------------------------
inline int Matrix::Index(int r, int c) const
{
  return c * _cols + r;
}

// -----------------------------------------------------------------------------
inline int Matrix::RowIndex(int i) const
{
  return i % _cols;
}

// -----------------------------------------------------------------------------
inline int Matrix::ColIndex(int i) const
{
  return i / _cols;
}

// -----------------------------------------------------------------------------
inline void Matrix::SubIndex(int i, int &r, int &c) const
{
  r = RowIndex(i);
  c = ColIndex(i);
}

// -----------------------------------------------------------------------------
inline Pair<int, int> Matrix::SubIndex(int i) const
{
  return MakePair(RowIndex(i), ColIndex(i));
}

// =============================================================================
// Element access
// =============================================================================

// -----------------------------------------------------------------------------
inline double *Matrix::RawPointer(int r, int c)
{
  return &_matrix[c][r];
}

// -----------------------------------------------------------------------------
inline const double *Matrix::RawPointer(int r, int c) const
{
  return &_matrix[c][r];
}

// -----------------------------------------------------------------------------
inline double *Matrix::GetPointerToElements(int r, int c)
{
  return &_matrix[c][r];
}

// -----------------------------------------------------------------------------
inline const double *Matrix::GetPointerToElements(int r, int c) const
{
  return &_matrix[c][r];
}

// -----------------------------------------------------------------------------
inline double *Matrix::Col(int c)
{
  return _matrix[c];
}

// -----------------------------------------------------------------------------
inline const double *Matrix::Col(int c) const
{
  return _matrix[c];
}

// -----------------------------------------------------------------------------
inline double &Matrix::operator()(int i)
{
  return _matrix[0][i];
}

// -----------------------------------------------------------------------------
inline const double &Matrix::operator()(int i) const
{
  return const_cast<const double &>(const_cast<Matrix *>(this)->operator()(i));
}

// -----------------------------------------------------------------------------
inline double &Matrix::operator()(int r, int c)
{
  return _matrix[c][r];
}

// -----------------------------------------------------------------------------
inline const double &Matrix::operator()(int r, int c) const
{
  return const_cast<const double &>(const_cast<Matrix *>(this)->operator()(r, c));
}

// -----------------------------------------------------------------------------
inline void Matrix::Put(int i, double v)
{
  this->operator()(i) = v;
}

// -----------------------------------------------------------------------------
inline double Matrix::Get(int i) const
{
  return this->operator()(i);
}

// -----------------------------------------------------------------------------
inline void Matrix::Put(int r, int c, double v)
{
  this->operator()(r, c) = v;
}

// -----------------------------------------------------------------------------
inline double Matrix::Get(int r, int c) const
{
  return this->operator()(r, c);
}

// -----------------------------------------------------------------------------
inline Matrix3x3 Matrix::To3x3() const
{
  Matrix3x3 m;
  for (int r = 0; r < 3; ++r) {
    if (r < _rows) {
      for (int c = 0; c < 3; ++c) {
        if (c < _cols) {
          m[r][c] = _matrix[c][r];
        } else {
          m[r][c] = 0.;
        }
      }
    } else {
      memset(m[r], 0, 3 * sizeof(double));
    }
  }
  return m;
}

// =============================================================================
// Scalar matrix operations
// =============================================================================

// -----------------------------------------------------------------------------
inline Matrix &Matrix::AddToRow(int r, double s)
{
  for (int c = 0; c < _cols; ++c) _matrix[c][r] += s;
  return *this;
}

// -----------------------------------------------------------------------------
inline Matrix &Matrix::AddToCol(int c, double s)
{
  for (int r = 0; r < _rows; ++r) _matrix[c][r] += s;
  return *this;
}

// -----------------------------------------------------------------------------
inline Matrix &Matrix::ScaleRow(int r, double s)
{
  for (int c = 0; c < _cols; ++c) _matrix[c][r] *= s;
  return *this;
}

// -----------------------------------------------------------------------------
inline Matrix &Matrix::ScaleCol(int c, double s)
{
  for (int r = 0; r < _rows; ++r) _matrix[c][r] *= s;
  return *this;
}

// -----------------------------------------------------------------------------
inline Matrix &Matrix::operator =(double s)
{
  const int n = this->NumberOfElements();
  double   *p = this->RawPointer();
  for (int i = 0; i < n; ++i, ++p) (*p) = s;
  return *this;
}

// -----------------------------------------------------------------------------
inline Matrix &Matrix::operator -=(double x)
{
  const int n = this->NumberOfElements();
  double   *p = this->RawPointer();
  for (int i = 0; i < n; ++i, ++p) (*p) -= x;
  return *this;
}

// -----------------------------------------------------------------------------
inline Matrix &Matrix::operator +=(double x)
{
  const int n = this->NumberOfElements();
  double   *p = this->RawPointer();
  for (int i = 0; i < n; ++i, ++p) (*p) += x;
  return *this;
}

// -----------------------------------------------------------------------------
inline Matrix &Matrix::operator *=(double x)
{
  const int n = this->NumberOfElements();
  double   *p = this->RawPointer();
  for (int i = 0; i < n; ++i, ++p) (*p) *= x;
  return *this;
}

// -----------------------------------------------------------------------------
inline Matrix &Matrix::operator /=(double x)
{
  const int n = this->NumberOfElements();
  double   *p = this->RawPointer();
  for (int i = 0; i < n; ++i, ++p) (*p) /= x;
  return *this;
}

// -----------------------------------------------------------------------------
inline Matrix Matrix::operator -(double x) const
{
  return Matrix(*this) -= x;
}

// -----------------------------------------------------------------------------
inline Matrix Matrix::operator +(double x) const
{
  return Matrix(*this) += x;
}

// -----------------------------------------------------------------------------
inline Matrix Matrix::operator *(double x) const
{
  return Matrix(*this) *= x;
}

// -----------------------------------------------------------------------------
inline Matrix Matrix::operator /(double x) const
{
  return Matrix(*this) /= x;
}

// =============================================================================
// Comparison operations
// =============================================================================

// -----------------------------------------------------------------------------
inline bool Matrix::operator ==(const Matrix &m) const
{
  if ((m._rows != _rows) || (m._cols != _cols)) return false;
  const int n = this->NumberOfElements();
  const double *ptr1 = m   . RawPointer();
  const double *ptr2 = this->RawPointer();
  for (int i = 0; i < n; ++i, ++ptr1, ++ptr2) {
    if (!fequal(*ptr2, *ptr1)) return false;
  }
  return true;
}

// -----------------------------------------------------------------------------
inline bool Matrix::operator !=(const Matrix &m) const
{
  return !(*this == m);
}

// -----------------------------------------------------------------------------
inline Array<int> Matrix::operator ==(double x) const
{
  Array<int> idx;
  const int     n = this->NumberOfElements();
  const double *p = this->RawPointer();
  for (int i = 0; i < n; ++i, ++p) {
    if (*p == x) idx.push_back(i);
  }
  return idx;
}

// -----------------------------------------------------------------------------
inline Array<int> Matrix::operator !=(double x) const
{
  Array<int> idx;
  const int     n = this->NumberOfElements();
  const double *p = this->RawPointer();
  for (int i = 0; i < n; ++i, ++p) {
    if (*p != x) idx.push_back(i);
  }
  return idx;
}

// -----------------------------------------------------------------------------
inline Array<int> Matrix::operator <(double x) const
{
  Array<int> idx;
  const int     n = this->NumberOfElements();
  const double *p = this->RawPointer();
  for (int i = 0; i < n; ++i, ++p) {
    if (*p < x) idx.push_back(i);
  }
  return idx;
}

// -----------------------------------------------------------------------------
inline Array<int> Matrix::operator <=(double x) const
{
  Array<int> idx;
  const int     n = this->NumberOfElements();
  const double *p = this->RawPointer();
  for (int i = 0; i < n; ++i, ++p) {
    if (*p <= x) idx.push_back(i);
  }
  return idx;
}

// -----------------------------------------------------------------------------
inline Array<int> Matrix::operator >(double x) const
{
  Array<int> idx;
  const int     n = this->NumberOfElements();
  const double *p = this->RawPointer();
  for (int i = 0; i < n; ++i, ++p) {
    if (*p > x) idx.push_back(i);
  }
  return idx;
}

// -----------------------------------------------------------------------------
inline Array<int> Matrix::operator >=(double x) const
{
  Array<int> idx;
  const int     n = this->NumberOfElements();
  const double *p = this->RawPointer();
  for (int i = 0; i < n; ++i, ++p) {
    if (*p >= x) idx.push_back(i);
  }
  return idx;
}

// =============================================================================
// Matrix functions
// =============================================================================

// -----------------------------------------------------------------------------
inline double Matrix::RowMin(int r) const
{
  if (_rows == 0) return .0;
  double v = _matrix[0][r];
  for (int c = 1; c < _cols; ++c) v = min(v, _matrix[c][r]);
  return v;
}

// -----------------------------------------------------------------------------
inline double Matrix::RowMax(int r) const
{
  if (_rows == 0) return .0;
  double v = _matrix[0][r];
  for (int c = 1; c < _cols; ++c) v = max(v, _matrix[c][r]);
    return v;
}

// -----------------------------------------------------------------------------
inline void Matrix::RowRange(int r, double &v1, double &v2) const
{
  if (_rows > 0) {
    v1 = v2 = _matrix[0][r];
    for (int c = 1; c < _cols; ++c) {
      const double &v = _matrix[c][r];
      v1 = min(v1, v), v2 = max(v2, v);
    }
  } else {
    v1 = v2 = .0;
  }
}

// -----------------------------------------------------------------------------
inline double Matrix::ColMin(int c) const
{
  if (_cols == 0) return .0;
  double v = _matrix[c][0];
  for (int r = 1; r < _rows; ++r) v = min(v, _matrix[c][r]);
  return v;
}

// -----------------------------------------------------------------------------
inline double Matrix::ColMax(int c) const
{
  if (_cols == 0) return .0;
  double v = _matrix[c][0];
  for (int r = 1; r < _rows; ++r) v = min(v, _matrix[c][r]);
  return v;
}

// -----------------------------------------------------------------------------
inline void Matrix::ColRange(int c, double &v1, double &v2) const
{
  if (_cols > 0) {
    v1 = v2 = _matrix[c][0];
    for (int r = 1; r < _rows; ++r) {
      const double &v = _matrix[c][r];
      v1 = min(v1, v), v2 = max(v2, v);
    }
  } else {
    v1 = v2 = .0;
  }
}

// -----------------------------------------------------------------------------
inline double Matrix::Trace() const
{
  if (_rows == _cols) {
    double trace = 0;
    for (int j = 0; j < _cols; ++j)
    for (int i = 0; i < _rows; ++i) {
      trace += _matrix[j][i];
    }
    return trace;
  } else {
    cerr << "Matrix::Trace() matrix number of col != row" << endl;
    return 0;
  }
}

// -----------------------------------------------------------------------------
inline double Matrix::Norm() const
{
  int i, j;
  double norm = 0;

  // The norm of a matrix M is defined as trace(M * M~)
  for (j = 0; j < _cols; j++) {
    for (i = 0; i < _rows; i++) {
      norm += _matrix[j][i]*_matrix[j][i];
    }
  }
  return sqrt(norm);
}

// -----------------------------------------------------------------------------
inline double Matrix::InfinityNorm() const
{
  int i, j;
  double normInf = -1.0 * DBL_MAX;
  double sum;

  for (i = 0; i < _rows; ++i) {
    sum = 0;
    for (j = 0; j < _cols; ++j) {
      sum += abs(_matrix[j][i]);
    }
    if (sum > normInf)
      normInf = sum;
  }
  return normInf;
}

// -----------------------------------------------------------------------------
inline double Matrix::Det3x3() const
{
  return _matrix[0][0] * (_matrix[1][1] * _matrix[2][2] - _matrix[1][2] * _matrix[2][1])
       - _matrix[0][1] * (_matrix[1][0] * _matrix[2][2] - _matrix[1][2] * _matrix[2][0])
       + _matrix[0][2] * (_matrix[1][0] * _matrix[2][1] - _matrix[1][1] * _matrix[2][0]);
}

// -----------------------------------------------------------------------------
inline Matrix Matrix::Inverse(bool use_svd_if_singular) const
{
  return Matrix(*this).Invert(use_svd_if_singular);
}

// -----------------------------------------------------------------------------
inline Matrix Matrix::SVDInverse() const
{
  return Matrix(*this).SVDInvert();
}

// -----------------------------------------------------------------------------
inline Matrix Matrix::PseudoInverse() const
{
  return Matrix(*this).PseudoInvert();
}

// -----------------------------------------------------------------------------
inline Matrix Matrix::operator !()
{
  return Inverse();
}

// -----------------------------------------------------------------------------
inline Matrix Matrix::Transposed() const
{
  return Matrix(*this).Transpose();
}

// -----------------------------------------------------------------------------
inline Matrix Matrix::operator ~()
{
  return Matrix(*this).Transpose();
}

// =============================================================================
// Backwards compatibility
// =============================================================================

// -----------------------------------------------------------------------------
inline Matrix expm(const Matrix &m)
{
  return m.Exp();
}

// -----------------------------------------------------------------------------
inline Matrix logm(const Matrix &m)
{
  return m.Log();
}

// -----------------------------------------------------------------------------
inline Matrix sqrtm(const Matrix &m)
{
  return m.Sqrt();
}

////////////////////////////////////////////////////////////////////////////////
// Conversion to CUDA vector types
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
// Construct float3x3 matrix from upper left 3x3 sub-matrix of Matrix
MIRTKCU_HOST_API inline float3x3 make_float3x3(const Matrix &m)
{
  float3x3 f;
  f.a = make_float3(m(0, 0), m(0, 1), m(0, 2));
  f.b = make_float3(m(1, 0), m(1, 1), m(1, 2));
  f.c = make_float3(m(2, 0), m(2, 1), m(2, 2));
  return f;
}

// -----------------------------------------------------------------------------
// Construct float3x4 matrix from upper left 3x4 sub-matrix of Matrix
MIRTKCU_HOST_API inline float3x4 make_float3x4(const Matrix &m)
{
  float3x4 f;
  f.a = make_float4(m(0, 0), m(0, 1), m(0, 2), m(0, 3));
  f.b = make_float4(m(1, 0), m(1, 1), m(1, 2), m(1, 3));
  f.c = make_float4(m(2, 0), m(2, 1), m(2, 2), m(2, 3));
  return f;
}

// -----------------------------------------------------------------------------
// Construct float4x4 matrix from upper left 4x4 sub-matrix of Matrix
MIRTKCU_HOST_API inline float4x4 make_float4x4(const Matrix &m)
{
  float4x4 d;
  d.a = make_float4(m(0, 0), m(0, 1), m(0, 2), m(0, 3));
  d.b = make_float4(m(1, 0), m(1, 1), m(1, 2), m(1, 3));
  d.c = make_float4(m(2, 0), m(2, 1), m(2, 2), m(2, 3));
  d.d = make_float4(m(3, 0), m(3, 1), m(3, 2), m(3, 3));
  return d;
}

// -----------------------------------------------------------------------------
// Construct double3x3 matrix from upper left 3x3 sub-matrix of Matrix
MIRTKCU_HOST_API inline double3x3 make_double3x3(const Matrix &m)
{
  double3x3 d;
  d.a = make_double3(m(0, 0), m(0, 1), m(0, 2));
  d.b = make_double3(m(1, 0), m(1, 1), m(1, 2));
  d.c = make_double3(m(2, 0), m(2, 1), m(2, 2));
  return d;
}

// -----------------------------------------------------------------------------
// Construct double3x4 matrix from upper left 3x4 sub-matrix of Matrix
MIRTKCU_HOST_API inline double3x4 make_double3x4(const Matrix &m)
{
  double3x4 d;
  d.a = make_double4(m(0, 0), m(0, 1), m(0, 2), m(0, 3));
  d.b = make_double4(m(1, 0), m(1, 1), m(1, 2), m(1, 3));
  d.c = make_double4(m(2, 0), m(2, 1), m(2, 2), m(2, 3));
  return d;
}

// -----------------------------------------------------------------------------
// Construct double4x4 matrix from upper left 4x4 sub-matrix of Matrix
MIRTKCU_HOST_API inline double4x4 make_double4x4(const Matrix &m)
{
  double4x4 d;
  d.a = make_double4(m(0, 0), m(0, 1), m(0, 2), m(0, 3));
  d.b = make_double4(m(1, 0), m(1, 1), m(1, 2), m(1, 3));
  d.c = make_double4(m(2, 0), m(2, 1), m(2, 2), m(2, 3));
  d.d = make_double4(m(3, 0), m(3, 1), m(3, 2), m(3, 3));
  return d;
}

////////////////////////////////////////////////////////////////////////////////
// Means on Lie groups of transformation matrices
////////////////////////////////////////////////////////////////////////////////

/// Compute log-Euclidean mean of transformation matrices
///
/// \f[
///   mu = \exp\left( \sum w_i \log\left( matrices_i \right) \right)
/// \f]
///
/// \param[in] n        Number of matrices.
/// \param[in] matrices Transformation matrices.
/// \param[in] weights  Weights of transformation matrices.
///                     Uniform weighting if \c NULL.
///
/// \return Mean transformation matrix.
Matrix LogEuclideanMean(int n, const Matrix *matrices, const double *weights = NULL);

/// Compute exponential barycenter / bi-invariant mean of transformation matrices
///
/// This function finds the exponential barycenter of a set of rigid/affine
/// transformation matrices using a fixed point iteration (Gauss-Newton;
/// Barycentric fixed point iteration on Lie groups).
///
/// \f[
///   mu_{t+1} = mu_t \exp\left( \sum w_i \log\left( mu_t^{-1} matrices_i \right) \right)
/// \f]
///
/// \see Arsigny and Pennec, Exponential Barycenters of the Canonical Cartan
///      Connection and Invariant Means on Lie Groups,
///      Matrix Information Geometry (2012)
///
/// \param[in] n        Number of matrices.
/// \param[in] matrices Transformation matrices.
/// \param[in] weights  Weights of transformation matrices.
///                     Uniform weighting if \c NULL.
/// \param[in] niter    Maximum number of fixed point iterations.
/// \param[in] tol      Tolerance of residual infinity norm.
/// \param[in] mu0      Start value of fixed point iteration.
///                     First matrix if \c NULL.
///
/// \return Mean transformation matrix.
Matrix BiInvariantMean(int n, const Matrix *matrices,
                       const double *weights = NULL,
                       int niter = 20, double tol = 1e-12,
                       const Matrix *mu0 = NULL);

// -----------------------------------------------------------------------------
/// \deprecated Use BiInvariantMean instead
inline Matrix FrechetMean(const Matrix *matrices, const double *weights, int n,
                          int niter = 20, double tol = 1e-12,
                          const Matrix *mu0 = NULL)
{
  return BiInvariantMean(n, matrices, weights, niter, tol, mu0);
}

// -----------------------------------------------------------------------------
/// \deprecated Use BiInvariantMean instead
inline Matrix FrechetMean(const Matrix *matrices, int n,
                          int niter = 20, double tol = 1e-12,
                          const Matrix *mu0 = NULL)
{
  return BiInvariantMean(n, matrices, NULL, niter, tol, mu0);
}

////////////////////////////////////////////////////////////////////////////////
// Affine transformation matrices
////////////////////////////////////////////////////////////////////////////////

/// Muliply point by homogeneous transformation matrix
inline void Transform(const Matrix &m, double  x,  double  y,  double  z,
                                       double &mx, double &my, double &mz)
{
  mx = m(0, 0) * x + m(0, 1) * y + m(0, 2) * z + m(0, 3);
  my = m(1, 0) * x + m(1, 1) * y + m(1, 2) * z + m(1, 3);
  mz = m(2, 0) * x + m(2, 1) * y + m(2, 2) * z + m(2, 3);
}

/// Muliply point by homogeneous transformation matrix
inline void Transform(const Matrix &m, double &x, double &y, double &z)
{
  Transform(m, x, y, z, x, y, z);
}

/// Muliply point by homogeneous transformation matrix
inline Point Transform(const Matrix &m, const Point &p)
{
  Point p2;
  Transform(m, p._x, p._y, p._z, p2._x, p2._y, p2._z);
  return p2;
}

/// Muliply vector by 3x3 transformation matrix
inline void TransformVector(const Matrix &m, double  x,  double  y,  double  z,
                                             double &mx, double &my, double &mz)
{
  mx = m(0, 0) * x + m(0, 1) * y + m(0, 2) * z;
  my = m(1, 0) * x + m(1, 1) * y + m(1, 2) * z;
  mz = m(2, 0) * x + m(2, 1) * y + m(2, 2) * z;
}

/// Muliply vector by 3x3 transformation matrix
inline void TransformVector(const Matrix &m, double &x, double &y, double &z)
{
  TransformVector(m, x, y, z, x, y, z);
}

/// Muliply vector by 3x3 transformation matrix
inline void TransformVector(const Matrix &m, Vector &v)
{
  TransformVector(m, v(0), v(1), v(2));
}

/// Muliply vector by 3x3 transformation matrix
inline Vector TransformVector(const Matrix &m, const Vector &v)
{
  Vector u(v);
  TransformVector(m, u);
  return u;
}

/// Construct 4x4 homogeneous coordinate transformation matrix from rigid
/// transformation parameters. The output transformation matrix is the
/// composition of a rotation followed by a translation, i.e.,
///
///   T = Translate * Rotate, where Rotate = (Rz Ry Rx)^T
///
/// \note This overloaded function is mainly used by the linear transformation
///       classes such as in particular irtkRigidTransformation as these
///       store the cosine and sine values of the rotation angles as members.
///
/// \param[in]  tx    Translation along x axis.
/// \param[in]  ty    Translation along y axis.
/// \param[in]  tz    Translation along z axis.
/// \param[in]  cosrx Cosine of rotation around x axis in radians.
/// \param[in]  cosry Cosine of rotation around y axis in radians.
/// \param[in]  cosrz Cosine of rotation around z axis in radians.
/// \param[in]  sinrx Sine   of rotation around x axis in radians.
/// \param[in]  sinry Sine   of rotation around y axis in radians.
/// \param[in]  sinrz Sine   of rotation around z axis in radians.
/// \param[out] m     Homogeneous transformation matrix.
void RigidParametersToMatrix(double tx,     double ty,     double tz,
                             double cosrx,  double cosry,  double cosrz,
                             double sinrx,  double sinry,  double sinrz,
                             Matrix &m);

/// Construct 4x4 homogeneous coordinate transformation matrix from rigid
/// transformation parameters. The output transformation matrix is the
/// composition of a rotation followed by a translation, i.e.,
///
///   T = Translate * Rotate
///
/// \param[in]  tx Translation along x axis.
/// \param[in]  ty Translation along y axis.
/// \param[in]  tz Translation along z axis.
/// \param[in]  rx Rotation around x axis in radians.
/// \param[in]  ry Rotation around y axis in radians.
/// \param[in]  rz Rotation around z axis in radians.
/// \param[out] m  Homogeneous transformation matrix.
void RigidParametersToMatrix(double tx,  double ty,  double tz,
                             double rx,  double ry,  double rz, Matrix &m);

/// Construct 4x4 homogeneous coordinate transformation matrix from rigid
/// transformation parameters. The output transformation matrix is the
/// composition of a rotation followed by a translation, i.e.,
///
///   T = Translate * Rotate, where Rotate = (Rz Ry Rx)^T
///
/// \param[in]  tx Translation along x axis.
/// \param[in]  ty Translation along y axis.
/// \param[in]  tz Translation along z axis.
/// \param[in]  rx Rotation around x axis in radians.
/// \param[in]  ry Rotation around y axis in radians.
/// \param[in]  rz Rotation around z axis in radians.
///
/// \returns Homogeneous transformation matrix.
inline Matrix RigidParametersToMatrix(double tx,  double ty,  double tz,
                                      double rx,  double ry,  double rz)
{
  Matrix m;
  RigidParametersToMatrix(tx, ty, tz, rx, ry, rz, m);
  return m;
}

/// Extract Euler angles from 4x4 homogeneous coordinate transformation matrix.
/// The input transformation matrix is assumed to be the composition of a
/// rotation followed by a translation, i.e.,
///
///   T = Translate * Rotate, where Rotate = (Rz Ry Rx)^T
///
/// \param[in]  m  Homogeneous coordinate transformation matrix.
/// \param[out] rx Rotation around x axis in radians.
/// \param[out] ry Rotation around y axis in radians.
/// \param[out] rz Rotation around z axis in radians.
void MatrixToEulerAngles(const Matrix &m,
                         double &rx,  double &ry,  double &rz);

/// Extract rigid transformation parameters from 4x4 homogeneous coordinate
/// transformation matrix. The input transformation matrix is assumed to be
/// the composition of a rotation followed by a translation, i.e.,
///
///   T = Translate * Rotate, where Rotate = (Rz Ry Rx)^T
///
/// \param[in]  m  Homogeneous coordinate transformation matrix.
/// \param[out] tx Translation along x axis.
/// \param[out] ty Translation along y axis.
/// \param[out] tz Translation along z axis.
/// \param[out] rx Rotation around x axis in radians.
/// \param[out] ry Rotation around y axis in radians.
/// \param[out] rz Rotation around z axis in radians.
void MatrixToRigidParameters(const Matrix &m,
                             double &tx,  double &ty,  double &tz,
                             double &rx,  double &ry,  double &rz);

/// Construct 4x4 homogeneous coordinate transformation matrix from affine
/// transformation parameters. The output transformation matrix is the
/// composition of a shearing, followed by a scaling, followed by a
/// rotation, followed by a translation, i.e.,
///
///   T = Translate * Rotate * Scale * Shear, where Rotate = (Rz Ry Rx)^T
///
/// \param[in]  tx  Translation along x axis.
/// \param[in]  ty  Translation along y axis.
/// \param[in]  tz  Translation along z axis.
/// \param[in]  rx  Rotation around x axis in radians.
/// \param[in]  ry  Rotation around y axis in radians.
/// \param[in]  rz  Rotation around z axis in radians.
/// \param[in]  sx  Scaling of x axis (factor, not percentage).
/// \param[in]  sy  Scaling of y axis (factor, not percentage).
/// \param[in]  sz  Scaling of z axis (factor, not percentage).
/// \param[in]  sxy Skew between x and y axes in radians.
/// \param[in]  sxz Skew between x and z axes in radians.
/// \param[in]  syz Skew between y and z axes in radians.
/// \param[out] m  Homogeneous transformation matrix.
void AffineParametersToMatrix(double tx,  double ty,  double tz,
                              double rx,  double ry,  double rz,
                              double sx,  double sy,  double sz,
                              double sxy, double sxz, double syz, Matrix &m);

/// Construct 4x4 homogeneous coordinate transformation matrix from affine
/// transformation parameters. The output transformation matrix is the
/// composition a scaling, followed by a rotation, followed by a translation, i.e.,
///
///   T = Translate * Rotate * Scale, where Rotate = (Rz Ry Rx)^T
///
/// \param[in]  tx  Translation along x axis.
/// \param[in]  ty  Translation along y axis.
/// \param[in]  tz  Translation along z axis.
/// \param[in]  rx  Rotation around x axis in radians.
/// \param[in]  ry  Rotation around y axis in radians.
/// \param[in]  rz  Rotation around z axis in radians.
/// \param[in]  sx  Scaling of x axis (factor, not percentage).
/// \param[in]  sy  Scaling of y axis (factor, not percentage).
/// \param[in]  sz  Scaling of z axis (factor, not percentage).
/// \param[in]  sxy Skew between x and y axes in radians.
/// \param[in]  sxz Skew between x and z axes in radians.
/// \param[in]  syz Skew between y and z axes in radians.
///
/// \returns Homogeneous transformation matrix.
inline Matrix AffineParametersToMatrix(double tx,  double ty,  double tz,
                                       double rx,  double ry,  double rz,
                                       double sx,  double sy,  double sz,
                                       double sxy, double sxz, double syz)
{
  Matrix m(4, 4);
  AffineParametersToMatrix(tx, ty, tz, rx, ry, rz, sx, sy, sz, sxy, sxz, syz, m);
  return m;
}

/// Construct 4x4 homogeneous coordinate transformation matrix from affine
/// transformation parameters. The output transformation matrix is the
/// composition a scaling, followed by a rotation, followed by a translation, i.e.,
///
///   T = Translate * Rotate * Scale, where Rotate = (Rz Ry Rx)^T
///
/// \param[in]  tx  Translation along x axis.
/// \param[in]  ty  Translation along y axis.
/// \param[in]  tz  Translation along z axis.
/// \param[in]  rx  Rotation around x axis in radians.
/// \param[in]  ry  Rotation around y axis in radians.
/// \param[in]  rz  Rotation around z axis in radians.
/// \param[in]  sx  Scaling of x axis (factor, not percentage).
/// \param[in]  sy  Scaling of y axis (factor, not percentage).
/// \param[in]  sz  Scaling of z axis (factor, not percentage).
/// \param[out] m  Homogeneous transformation matrix.
inline void AffineParametersToMatrix(double tx,  double ty,  double tz,
                                     double rx,  double ry,  double rz,
                                     double sx,  double sy,  double sz, Matrix &m)
{
  AffineParametersToMatrix(tx, ty, tz, rx, ry, rz, sx, sy, sz, .0, .0, .0, m);
}

/// Construct 4x4 homogeneous coordinate transformation matrix from affine
/// transformation parameters. The output transformation matrix is the
/// composition a scaling, followed by a rotation, followed by a translation, i.e.,
///
///   T = Translate * Rotate * Scale, where Rotate = (Rz Ry Rx)^T
///
/// \param[in]  tx  Translation along x axis.
/// \param[in]  ty  Translation along y axis.
/// \param[in]  tz  Translation along z axis.
/// \param[in]  rx  Rotation around x axis in radians.
/// \param[in]  ry  Rotation around y axis in radians.
/// \param[in]  rz  Rotation around z axis in radians.
/// \param[in]  sx  Scaling of x axis (factor, not percentage).
/// \param[in]  sy  Scaling of y axis (factor, not percentage).
/// \param[in]  sz  Scaling of z axis (factor, not percentage).
///
/// \returns Homogeneous transformation matrix.
inline Matrix AffineParametersToMatrix(double tx,  double ty,  double tz,
                                       double rx,  double ry,  double rz,
                                       double sx,  double sy,  double sz)
{
  Matrix m(4, 4);
  AffineParametersToMatrix(tx, ty, tz, rx, ry, rz, sx, sy, sz, .0, .0, .0, m);
  return m;
}

/// Extract affine transformation parameters from 4x4 homogeneous coordinate
/// transformation matrix. The input transformation matrix is assumed to be
/// the composition of a shearing, followed by a scaling, followed by a
/// rotation, followed by a translation, i.e.,
///
///   T = Translate * Rotate * Scale * Shear, where Rotate = (Rz Ry Rx)^T
///
///  @sa Graphicx Gems II, with a 12 DoF model without perspective distortion
///      https://github.com/erich666/GraphicsGems/blob/master/gemsii/unmatrix.c
///
/// \param[in]  m   Homogeneous coordinate transformation matrix.
/// \param[out] tx  Translation along x axis.
/// \param[out] ty  Translation along y axis.
/// \param[out] tz  Translation along z axis.
/// \param[out] rx  Rotation around x axis in radians.
/// \param[out] ry  Rotation around y axis in radians.
/// \param[out] rz  Rotation around z axis in radians.
/// \param[out] sx  Scaling of x axis (factor, not percentage).
/// \param[out] sy  Scaling of y axis (factor, not percentage).
/// \param[out] sz  Scaling of z axis (factor, not percentage).
/// \param[out] sxy Skew between x and y axes in radians.
/// \param[out] sxz Skew between x and z axes in radians.
/// \param[out] syz Skew between y and z axes in radians.
void MatrixToAffineParameters(const Matrix &m,
                              double &tx,  double &ty,  double &tz,
                              double &rx,  double &ry,  double &rz,
                              double &sx,  double &sy,  double &sz,
                              double &sxy, double &sxz, double &syz);

/// Find affine transformation matrix which minimizes the mean squared distance
/// between two given sets of corresponding points (e.g., landmarks, fiducial markers)
///
/// \param[in] target (Transformed) Target point set.
/// \param[in] source Fixed source point set.
/// \param[in] weight Weight of corresponding point pair. All correspondences
///                   are weighted equally if vector is empty.
///
/// \returns Homogeneous transformation matrix of the target points.
Matrix ApproximateAffineMatrix(const PointSet &target,
                               const PointSet &source,
                               const Vector   &weight);

/// Find affine transformation matrix which minimizes the mean squared distance
/// between two given sets of corresponding points (e.g., landmarks, fiducial markers)
///
/// \param[in] target (Transformed) Target point set.
/// \param[in] source Fixed source point set.
///
/// \returns Homogeneous transformation matrix of the target points.
Matrix ApproximateAffineMatrix(const PointSet &target,
                               const PointSet &source);

/// Orthonormalize upper 3x3 matrix using stabilized Gram-Schmidt
void OrthoNormalize3x3(Matrix &);


} // namespace mirtk

#endif // MIRTK_Matrix_H
