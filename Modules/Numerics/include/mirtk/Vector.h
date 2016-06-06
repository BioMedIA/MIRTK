/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2015 Imperial College London
 * Copyright 2008-2015 Daniel Rueckert, Julia Schnabel
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

#ifndef MIRTK_Vector_H
#define MIRTK_Vector_H

#include "mirtk/Object.h"

#include "mirtk/Math.h"
#include "mirtk/Indent.h"
#include "mirtk/Cfstream.h"
#include "mirtk/Memory.h"
#include "mirtk/Array.h"


namespace mirtk {


// Forward declaration of specialized vector types included after
// declaration of Vector for definition of inline functions
template <class T> struct Vector3D;
template <class T> struct Vector4D;

/**
 * Vector class.
 */
class Vector : public Object
{
  mirtkObjectMacro(Vector);

  // ---------------------------------------------------------------------------
  // Data members

protected:

  /// Number of rows
  int _rows;

  /// Vector elements
  double *_vector;

public:

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Default constructor
  Vector();

  /// Constructor for given row dimensions
  explicit Vector(int);

  /// Constructor for given row dimensions
  Vector(int, double);

  /// Constructor for given row dimensions
  Vector(int, double *);

  /// Copy constructor
  Vector(const Vector &);

  /// Construct vector from 3D vector
  template <class T> Vector(const Vector3D<T> &);

  /// Construct vector from 4D vector
  template <class T> Vector(const Vector4D<T> &);

  /// Destructor
  ~Vector();

  /// Intialize matrix with number of rows
  void Initialize(int);

  /// Intialize matrix with number of rows
  void Initialize(int, double);

  /// Intialize matrix with number of rows
  void Initialize(int, double *);

  /// Change size of vector, preserving existing rows
  void Resize(int, double = .0);

  /// Free vector
  void Clear();

  /// Initialize from 3D vector
  template <class T> Vector &Put(const Vector3D<T> &);

  /// Initialize from 4D vector
  template <class T> Vector &Put(const Vector4D<T> &);

  /// Whether vector is non-empty, i.e., initialized and number of rows greater zero
  operator bool() const;

  // ---------------------------------------------------------------------------
  // Vector access functions

  /// Returns number of rows
  int Rows() const;

  /// Puts vector value
  void Put(int, double);

  /// Gets vector value
  const double &Get(int) const;

  // ---------------------------------------------------------------------------
  // Element access

  /// Get pointer to linear memory which stores vector elements
  double *RawPointer(int r = 0);

  /// Get pointer to linear memory which stores vector elements
  const double *RawPointer(int r = 0) const;

  /// Puts vector value
  double &operator()(int);

  /// Gets vector value
  const double &operator()(int) const;

  // ---------------------------------------------------------------------------
  // Vector/scalar operations

  /// Assignment of double
  Vector &operator =(double);

  /// Subtraction of a double
  Vector &operator-=(double);

  /// Addition of a double
  Vector &operator+=(double);

  /// Multiplication with a double
  Vector &operator*=(double);

  /// Division by a double
  Vector &operator/=(double);

  /// Return result of subtraction of a double
  Vector operator- (double) const;

  /// Return result of addition of a double
  Vector operator+ (double) const;

  /// Return result of multiplication with a double
  Vector operator* (double) const;

  /// Return result of division by a double
  Vector operator/ (double) const;

  // ---------------------------------------------------------------------------
  // Element-wise vector/vector operations

  /// Unary negation operator
  Vector operator -() const;

  /// Vector copy operator
  Vector &operator =(const Vector &);

  /// Vector subtraction operator
  Vector &operator-=(const Vector &);

  /// Vector addition operator
  Vector &operator+=(const Vector &);

  /// Vector componentwise multiplication operator (no scalar nor cross product)
  Vector &operator*=(const Vector &);

  /// Vector componentwise division operator
  Vector &operator/=(const Vector &);

  /// Return result for vector subtraction
  Vector operator- (const Vector &) const;

  /// Return result for vector addition
  Vector operator+ (const Vector &) const;

  /// Return result for componentwise vector multiplication (no scalar nor cross product)
  Vector operator* (const Vector &) const;

  /// Return result for componentwise vector division
  Vector operator/ (const Vector &) const;

  // ---------------------------------------------------------------------------
  // Comparison

  /// Comparison operator ==
  bool operator==(const Vector &) const;

  /// Comparison operator !=
  bool operator!=(const Vector &) const;

  /// Comparison operator <
  bool operator<(const Vector &) const;

  // ---------------------------------------------------------------------------
  // Vector products

  /// Scalar/dot product
  double ScalarProduct(const Vector &) const;

  /// Scalar/dot product
  double DotProduct(const Vector &) const;

  /// Vector/cross product
  Vector CrossProduct(const Vector &) const;

  // ---------------------------------------------------------------------------
  // Other vector functions

  /// Returns sum of a vector components
  double Sum() const;

  /// Compute mean value of vector components
  double Mean() const;

  /// Returns norm of a vector
  double Norm() const;

  /// Normalize vector
  Vector &Normalize();

  /// Replace each vector element by its inverse
  Vector &Inverse();

  /// Permute vector elements
  void PermuteRows(Array<int>);

  /// Permute vector elements (cf. PermuteRows)
  void Permute(const Array<int> &);

  // ---------------------------------------------------------------------------
  // I/O

  /// Print vector
  void Print(Indent = 0) const;

  /// Read vector from file
  void Read(const char *);

  /// Write vector to file
  void Write(const char *) const;

  /// Interface to output stream
  friend ostream &operator<< (ostream &, const Vector &);

  /// Interface to input stream
  friend istream &operator>> (istream &, Vector &);

  /// Interface to output stream
  friend Cofstream &operator<< (Cofstream &, const Vector &);

  /// Interface to input stream
  friend Cifstream &operator>> (Cifstream &, Vector &);

  /// Write vector to MAT-file
  ///
  /// \note Use only when MIRTK_Numerics_WITH_MATLAB is 1.
  bool WriteMAT(const char *, const char * = "A") const;

};


} // namespace mirtk

#include "mirtk/Vector3D.h"
#include "mirtk/Vector4D.h"

namespace mirtk {

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
inline Vector::Vector()
:
  _rows  (0),
  _vector(NULL)
{
}

// -----------------------------------------------------------------------------
inline Vector::Vector(int rows)
:
  _rows  (0),
  _vector(NULL)
{
  Initialize(rows);
}

// -----------------------------------------------------------------------------
inline Vector::Vector(int rows, double s)
:
  _rows  (0),
  _vector(NULL)
{
  Initialize(rows, s);
}

// -----------------------------------------------------------------------------
inline Vector::Vector(int rows, double *v)
:
  _rows  (0),
  _vector(NULL)
{
  Initialize(rows, v);
}

// -----------------------------------------------------------------------------
inline Vector::Vector(const Vector& v)
:
  Object(v),
  _rows  (0),
  _vector(NULL)
{
  Initialize(v._rows, v._vector);
}

// -----------------------------------------------------------------------------
template <class T>
inline Vector::Vector(const Vector3D<T> &v)
{
  _rows      = 3;
  _vector    = new double[3];
  _vector[0] = static_cast<double>(v._x);
  _vector[1] = static_cast<double>(v._y);
  _vector[2] = static_cast<double>(v._z);
}

// -----------------------------------------------------------------------------
template <class T>
inline Vector::Vector(const Vector4D<T> &v)
{
  _rows      = 4;
  _vector    = new double[4];
  _vector[0] = static_cast<double>(v._x);
  _vector[1] = static_cast<double>(v._y);
  _vector[2] = static_cast<double>(v._z);
  _vector[3] = static_cast<double>(v._t);
}

// -----------------------------------------------------------------------------
inline Vector::~Vector()
{
  delete[] _vector;
}

// -----------------------------------------------------------------------------
inline void Vector::Initialize(int rows)
{
  if (_rows != rows) {
    delete[] _vector;
    _rows   = rows;
    _vector = (_rows > 0) ? new double[_rows] : NULL;
  }
  memset(_vector, 0, _rows * sizeof(double));
}

// -----------------------------------------------------------------------------
inline void Vector::Initialize(int rows, double s)
{
  if (_rows != rows) {
    delete[] _vector;
    _rows   = rows;
    _vector = (_rows > 0) ? new double[_rows] : NULL;
  }
  for (int i = 0; i < _rows; i++) _vector[i] = s;
}

// -----------------------------------------------------------------------------
inline void Vector::Initialize(int rows, double *v)
{
  if (_rows != rows) {
    delete[] _vector;
    _rows   = rows;
    _vector = (_rows > 0) ? new double[_rows] : NULL;
  }
  if (v) memcpy(_vector, v, _rows * sizeof(double));
}

// -----------------------------------------------------------------------------
inline void Vector::Clear()
{
  delete[] _vector;
  _vector = NULL;
  _rows   = 0;
}

// -----------------------------------------------------------------------------
inline void Vector::Resize(int n, double value)
{
  if (n <= 0) {
    Clear();
  } else if (_rows != n) {
    double *vector = new double[n];
    const int m = min(n, _rows);
    for (int i = 0; i < m; ++i) vector[i] = _vector[i];
    for (int i = m; i < n; ++i) vector[i] = value;
    delete[] _vector;
    _vector = vector;
    _rows   = n;
  }
}

// -----------------------------------------------------------------------------
inline Vector::operator bool() const
{
  return _rows != 0;
}

// =============================================================================
// Access operators
// =============================================================================

// -----------------------------------------------------------------------------
inline int Vector::Rows() const
{
  return _rows;
}

// -----------------------------------------------------------------------------
inline double *Vector::RawPointer(int r)
{
  return &_vector[r];
}

// -----------------------------------------------------------------------------
inline const double *Vector::RawPointer(int r) const
{
  return &_vector[r];
}

// -----------------------------------------------------------------------------
inline void Vector::Put(int rows, double vector)
{
  _vector[rows] = vector;
}

// -----------------------------------------------------------------------------
template <class T>
inline Vector &Vector::Put(const Vector3D<T> &v)
{
  if (_rows != 3) {
    delete[] _vector;
    _rows      = 3;
    _vector    = new double[3];
  }
  _vector[0] = static_cast<double>(v._x);
  _vector[1] = static_cast<double>(v._y);
  _vector[2] = static_cast<double>(v._z);
  return *this;
}

// -----------------------------------------------------------------------------
template <class T>
inline Vector &Vector::Put(const Vector4D<T> &v)
{
  if (_rows != 4) {
    delete[] _vector;
    _rows      = 4;
    _vector    = new double[4];
  }
  _vector[0] = static_cast<double>(v._x);
  _vector[1] = static_cast<double>(v._y);
  _vector[2] = static_cast<double>(v._z);
  _vector[3] = static_cast<double>(v._t);
  return *this;
}

// -----------------------------------------------------------------------------
inline const double &Vector::Get(int rows) const
{
  return _vector[rows];
}

// -----------------------------------------------------------------------------
inline double &Vector::operator()(int rows)
{
  return _vector[rows];
}

// -----------------------------------------------------------------------------
inline const double &Vector::operator()(int rows) const
{
  return _vector[rows];
}

// =============================================================================
// Vector/scalar operations
// =============================================================================

// -----------------------------------------------------------------------------
inline Vector &Vector::operator =(double x)
{
  for (int i = 0; i < _rows; ++i) _vector[i] = x;
  return *this;
}

// -----------------------------------------------------------------------------
inline Vector &Vector::operator-=(double x)
{
  for (int i = 0; i < _rows; ++i) _vector[i] -= x;
  return *this;
}

// -----------------------------------------------------------------------------
inline Vector &Vector::operator+=(double x)
{
  for (int i = 0; i < _rows; ++i) _vector[i] += x;
  return *this;
}

// -----------------------------------------------------------------------------
inline Vector &Vector::operator*=(double x)
{
  for (int i = 0; i < _rows; ++i) _vector[i] *= x;
  return *this;
}

// -----------------------------------------------------------------------------
inline Vector &Vector::operator/=(double x)
{
  for (int i = 0; i < _rows; ++i) _vector[i] /= x;
  return *this;
}

// -----------------------------------------------------------------------------
inline Vector Vector::operator- (double x) const
{
  return (Vector(*this) -= x);
}

// -----------------------------------------------------------------------------
inline Vector Vector::operator+ (double x) const
{
  return (Vector(*this) += x);
}

// -----------------------------------------------------------------------------
inline Vector Vector::operator* (double x) const
{
  return (Vector(*this) *= x);
}

// -----------------------------------------------------------------------------
inline Vector Vector::operator/ (double x) const
{
  return (Vector(*this) /= x);
}

// -----------------------------------------------------------------------------
inline Vector operator- (double x, const Vector &v)
{
  return v - x;
}

// -----------------------------------------------------------------------------
inline Vector operator+ (double x, const Vector &v)
{
  return v + x;
}

// -----------------------------------------------------------------------------
inline Vector operator* (double x, const Vector &v)
{
  return v * x;
}

// =============================================================================
// Element-wise vector/vector operations
// =============================================================================

// -----------------------------------------------------------------------------
inline Vector Vector::operator -() const
{
  Vector negative(_rows);
  for (int i = 0; i < _rows; ++i) negative._vector[i] = -_vector[i];
  return negative;
}

// -----------------------------------------------------------------------------
inline Vector &Vector::operator =(const Vector &v)
{
  Initialize(v._rows, v._vector);
  return *this;
}

// -----------------------------------------------------------------------------
inline Vector &Vector::operator-=(const Vector &v)
{
  if (_rows != v._rows) {
    cerr << "Vector::operator-=: Size mismatch" << endl;
    exit(1);
  }
  for (int i = 0; i < _rows; ++i) _vector[i] -= v._vector[i];
  return *this;
}

// -----------------------------------------------------------------------------
inline Vector &Vector::operator+=(const Vector &v)
{
  if (_rows != v._rows) {
    cerr << "Vector::operator+=: Size mismatch" << endl;
    exit(1);
  }
  for (int i = 0; i < _rows; ++i) _vector[i] += v._vector[i];
  return *this;
}

// -----------------------------------------------------------------------------
inline Vector &Vector::operator*=(const Vector &v)
{
  if (_rows != v._rows) {
    cerr << "Vector::operator*=: Size mismatch" << endl;
    exit(1);
  }
  for (int i = 0; i < _rows; ++i) _vector[i] *= v._vector[i];
  return *this;
}

// -----------------------------------------------------------------------------
inline Vector &Vector::operator/=(const Vector &v)
{
  if (_rows != v._rows) {
    cerr << "Vector::operator/=: Size mismatch" << endl;
    exit(1);
  }
  for (int i = 0; i < _rows; ++i) _vector[i] /= v._vector[i];
  return *this;
}

// -----------------------------------------------------------------------------
inline Vector Vector::operator- (const Vector &v) const
{
  return (Vector(*this) -= v);
}

// -----------------------------------------------------------------------------
inline Vector Vector::operator+ (const Vector &v) const
{
  return (Vector(*this) += v);
}

// -----------------------------------------------------------------------------
inline Vector Vector::operator* (const Vector &v) const
{
  return (Vector(*this) *= v);
}

// -----------------------------------------------------------------------------
inline Vector Vector::operator/ (const Vector& v) const
{
  return (Vector(*this) /= v);
}

// =============================================================================
// Comparison
// =============================================================================

// -----------------------------------------------------------------------------
inline bool Vector::operator==(const Vector &v) const
{
  if (_rows != v._rows) return false;
  for (int i = 0; i < _rows; ++i) {
    if (!fequal(_vector[i], v._vector[i])) return false;
  }
  return true;
}

// -----------------------------------------------------------------------------
inline bool Vector::operator!=(const Vector &v) const
{
  return !(*this == v);
}

// -----------------------------------------------------------------------------
inline bool Vector::operator<(const Vector &v) const
{
  if (_rows > v._rows) return false;
  for (int i = 0; i < _rows; ++i) {
    if (_vector[i] >= v._vector[i]) return false;
  }
  return true;
}

// =============================================================================
// Vector products
// =============================================================================

// -----------------------------------------------------------------------------
inline double Vector::ScalarProduct(const Vector &v) const
{
  if (_rows != v._rows) {
    cerr << "Vector::ScalarProduct: Size mismatch" << endl;
    exit(1);
  }
  double s = .0;
  for (int i = 0; i < _rows; ++i) {
    s += _vector[i] * v._vector[i];
  }
  return s;
}

// -----------------------------------------------------------------------------
inline double ScalarProduct(const Vector &a, const Vector &b)
{
  return a.ScalarProduct(b);
}

// -----------------------------------------------------------------------------
inline double Vector::DotProduct(const Vector &v) const
{
  return ScalarProduct(v);
}

// -----------------------------------------------------------------------------
inline double DotProduct(const Vector &a, const Vector &b)
{
  return a.DotProduct(b);
}

// -----------------------------------------------------------------------------
inline Vector Vector::CrossProduct(const Vector &v) const
{
  if (_rows != v._rows) {
    cerr << "Vector::CrossProduct: Size mismatch" << endl;
    exit(1);
  }
  int    a, b;
  Vector c(_rows, (double*)NULL); // allocate without initialization
  for (int i = 0; i < _rows; ++i) {
    a = (i+1) % _rows;
    b = (i+2) % _rows;
    c._vector[i] = _vector[a] * v._vector[b] - _vector[b] * v._vector[a];
  }
  return c;
}

// -----------------------------------------------------------------------------
inline Vector CrossProduct(const Vector &a, const Vector &b)
{
  return a.CrossProduct(b);
}

// =============================================================================
// Functions
// =============================================================================

// -----------------------------------------------------------------------------
inline double Vector::Sum() const
{
  double sum = .0;
  for (int i = 0; i < _rows; i++) sum += _vector[i];
  return sum;
}

// -----------------------------------------------------------------------------
inline double Vector::Mean() const
{
  return Sum() / _rows;
}

// -----------------------------------------------------------------------------
inline double Vector::Norm() const
{
  double norm = .0;
  for (int i = 0; i < _rows; i++) norm += _vector[i] * _vector[i];
  return sqrt(norm);
}

// -----------------------------------------------------------------------------
inline Vector &Vector::Normalize()
{
  double norm = Norm();
  if (norm != .0) (*this) /= norm;
  return *this;
}

// -----------------------------------------------------------------------------
inline Vector &Vector::Inverse()
{
  for (int i = 0; i < _rows; i++) {
    if (_vector[i] != .0) _vector[i] = 1.0 / _vector[i];
  }
  return *this;
}

// -----------------------------------------------------------------------------
inline void Vector::Permute(const Array<int> &idx)
{
  PermuteRows(idx);
}


} // namespace mirtk

#endif // MIRTK_Vector_H
