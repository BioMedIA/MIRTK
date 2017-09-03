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

#ifndef MIRTK_Vector3D_H
#define MIRTK_Vector3D_H

#include "mirtk/Math.h"
#include "mirtk/Point.h"
#include "mirtk/Vector3.h"


namespace mirtk {


/// Represents a 3D vector
///
/// Must be a primitive type which can be treated as an array of three
/// values of type T such that sizeof(Vector<T>) == 3 * sizeof(T).
/// Thus, this primitive vector type may not have any other data members
/// besides the three vector components. This is required especially
/// when Vector3D is used as voxel type of an image and further an
/// externally allocated continuous block of 3 * sizeof(T) bytes used
/// internally by the image instance which only reinterprets the memory
/// as consecutive Vector3D<T> instances, i.e.,
///
/// \code
/// const int X   = 256;
/// const int Y   = 256;
/// const int Z   = 128;
/// const int num = X * Y * Z;
/// double *data = new double[3 * num];
/// GenericImage<Vector3D<double> > image(X, Y, Z, data);
/// \endcode
template <typename T>
struct Vector3D
{
  typedef T ComponentType;

  // ---------------------------------------------------------------------------
  // Attributes

  T _x; ///< The x component
  T _y; ///< The y component
  T _z; ///< The z component

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Default constructor
  Vector3D();

  /// Construct from scalar
  Vector3D(T);

  /// Construct from vector components
  Vector3D(T, T, T);

  /// Construct from C array
  explicit Vector3D(const T [3]);

  /// Construct from non-template 3D vector type
  explicit Vector3D(const Vector3 &);

  /// Construct from 3D point
  explicit Vector3D(const Point &);

  /// Assignment operator
  Vector3D &operator =(const Vector3D &);

  /// Assignment operator
  Vector3D &operator =(const Point &);

  /// Copy constructor
  template <typename T2> Vector3D(const Vector3D<T2> &);

  // ---------------------------------------------------------------------------
  // Accessors

  /// Number of vector components
  static int Rows() { return 3; }

  /// Set/get vector component at index 0: _x, 1: _y, or 2: _z
  T &operator ()(int);

  /// Get vector component at index 0: _x, 1: _y, or 2: _z
  T operator ()(int) const;

  /// Cast to C array pointer
  operator T *();

  /// Cast to C array pointer
  operator const T *() const;

  // ---------------------------------------------------------------------------
  // Vector/integral-valued scalar operators

  /// Assign integral valued scalar
  Vector3D &operator =(int);

  /// Add integral valued scalar
  Vector3D &operator +=(int);

  /// Subtract integral valued scalar
  Vector3D &operator -=(int);

  /// Multiply by integral valued scalar
  Vector3D &operator *=(int);

  /// Divide by integral valued scalar
  Vector3D &operator /=(int);

  /// Add integral valued scalar to vector
  Vector3D operator +(int) const;

  /// Subtract integral valued scalar to vector
  Vector3D operator -(int) const;

  /// Multiply vector by integral valued scalar
  Vector3D operator *(int) const;

  /// Divide vector by integral valued scalar
  Vector3D operator /(int) const;

  // ---------------------------------------------------------------------------
  // Vector/real-valued scalar operators

  /// Assign real valued scalar
  Vector3D &operator =(double);

  /// Add real valued scalar
  Vector3D &operator +=(double);

  /// Subtract real valued scalar
  Vector3D &operator -=(double);

  /// Multiply by real valued scalar
  Vector3D &operator *=(double);

  /// Divide by real valued scalar
  Vector3D &operator /=(double);

  /// Add real valued scalar to vector
  Vector3D operator +(double) const;

  /// Subtract real valued scalar to vector
  Vector3D operator -(double) const;

  /// Multiply vector by real valued scalar
  Vector3D operator *(double) const;

  /// Divide vector by real valued scalar
  Vector3D operator /(double) const;

  // ---------------------------------------------------------------------------
  // Vector/vector operators

  /// Unary negation operator
  Vector3D operator -() const;

  /// Assignment from other vector
  template <typename T2> Vector3D &operator =(const Vector3D<T2> &);

  /// Addition of other vector
  template <typename T2> Vector3D &operator +=(const Vector3D<T2> &);

  /// Subtraction of other vector
  template <typename T2> Vector3D &operator -=(const Vector3D<T2> &);

  /// Element-wise multiplication with other vector
  template <typename T2> Vector3D &operator *=(const Vector3D<T2> &);

  /// Element-wise division by other vector
  template <typename T2> Vector3D &operator /=(const Vector3D<T2> &);

  /// Addition of two vectors
  template <typename T2> Vector3D operator +(const Vector3D<T2> &) const;

  /// Subtraction of two vectors
  template <typename T2> Vector3D operator -(const Vector3D<T2> &) const;

  /// Element-wise multiplication of two vectors
  template <typename T2> Vector3D operator *(const Vector3D<T2> &) const;

  /// Element-wise division of two vectors
  template <typename T2> Vector3D operator /(const Vector3D<T2> &) const;

  // ---------------------------------------------------------------------------
  // Vector/integral-valued scalar comparison

  /// Element-wise equality comparison with inegral-valued scalar
  bool operator ==(int) const;

  /// Element-wise inequality comparison with inegral-valued scalar
  bool operator !=(int) const;

  /// Element-wise less than comparison to inegral-valued scalar
  bool operator <(int) const;

  /// Element-wise greater than comparison to inegral-valued scalar
  bool operator >(int) const;

  /// Element-wise less or equal than comparison to inegral-valued scalar
  bool operator <=(int) const;

  /// Element-wise greater or equal than comparison to inegral-valued scalar
  bool operator >=(int) const;

  // ---------------------------------------------------------------------------
  // Vector/real-valued scalar comparison

  /// Element-wise equality comparison with real-valued scalar
  bool operator ==(double) const;

  /// Element-wise inequality comparison with real-valued scalar
  bool operator !=(double) const;

  /// Element-wise less than comparison to real-valued scalar
  bool operator <(double) const;

  /// Element-wise greater than comparison to real-valued scalar
  bool operator >(double) const;

  /// Element-wise less or equal than comparison to real-valued scalar
  bool operator <=(double) const;

  /// Element-wise greater or equal than comparison to real-valued scalar
  bool operator >=(double) const;

  // ---------------------------------------------------------------------------
  // Vector/vector comparison

  /** Operator for testing equality of two vectors. */
  template <typename T2> bool operator ==(const Vector3D<T2> &) const;

  /** Operator for testing non-equality of two vector. */
  template <typename T2> bool operator !=(const Vector3D<T2> &) const;

  /** Operator for comparing sizes of vectors. */
  template <typename T2> bool operator <(const Vector3D<T2> &) const;

  /** Operator for comparing sizes of vectors. */
  template <typename T2> bool operator >(const Vector3D<T2> &) const;

  /** Operator for comparing sizes of vectors. */
  template <typename T2> bool operator <=(const Vector3D<T2> &) const;

  /** Operator for comparing sizes of vectors. */
  template <typename T2> bool operator >=(const Vector3D<T2> &) const;

  // ---------------------------------------------------------------------------
  // Other vector functions

  /// Compute squared length of vector
  double SquaredLength() const;

  /// Compute length of vector
  double Length() const;

  /// Normalize vector to length one
  void Normalize();

  /// Cross-product with other vector
  Vector3D CrossProduct(const Vector3D &) const;

  /// Dot-product with other vector
  double DotProduct(const Vector3D &) const;

  /// Cross-product of two vectors
  ///
  /// \deprecated Use v1.CrossProduct(v2) instead which is less to type because
  ///             it does not require the Vector<T>:: type prefix and further
  ///             puts the name of the operation in between the arguments.
  static Vector3D CrossProduct(const Vector3D &, const Vector3D &);

  /// Dot-product of two vectors
  ///
  /// \deprecated Use v1.DotProduct(v2) instead which is less to type because
  ///             it does not require the Vector<T>:: type prefix and further
  ///             puts the name of the operation in between the arguments.
  static double DotProduct(const Vector3D &, const Vector3D &);

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
template <typename T>
inline Vector3D<T>::Vector3D()
{
  _x = _y = _z = T();
}

// -----------------------------------------------------------------------------
template <typename T>
inline Vector3D<T>::Vector3D(T s)
{
  _x = _y = _z = s;
}

// -----------------------------------------------------------------------------
template <typename T>
inline Vector3D<T>::Vector3D(T x, T y, T z)
{
  _x = x;
  _y = y;
  _z = z;
}

// -----------------------------------------------------------------------------
template <typename T>
inline Vector3D<T>::Vector3D(const T v[3])
{
  _x = v[0];
  _y = v[1];
  _z = v[2];
}

// -----------------------------------------------------------------------------
template <typename T>
inline Vector3D<T>::Vector3D(const Vector3 &v)
{
  _x = static_cast<T>(v[0]);
  _y = static_cast<T>(v[1]);
  _z = static_cast<T>(v[2]);
}

// -----------------------------------------------------------------------------
template <typename T>
inline Vector3D<T>::Vector3D(const Point &p)
{
  _x = static_cast<T>(p._x);
  _y = static_cast<T>(p._y);
  _z = static_cast<T>(p._z);
}

// -----------------------------------------------------------------------------
template <typename T1> template <typename T2>
inline Vector3D<T1>::Vector3D(const Vector3D<T2> &v)
{
  _x = static_cast<T1>(v._x);
  _y = static_cast<T1>(v._y);
  _z = static_cast<T1>(v._z);
}

// -----------------------------------------------------------------------------
template <typename T>
inline Vector3D<T> &Vector3D<T>::operator =(const Vector3D &v)
{
  _x = static_cast<T>(v._x);
  _y = static_cast<T>(v._y);
  _z = static_cast<T>(v._z);
  return *this;
}

// -----------------------------------------------------------------------------
template <typename T>
inline Vector3D<T> &Vector3D<T>::operator =(const Point &p)
{
  _x = static_cast<T>(p._x);
  _y = static_cast<T>(p._y);
  _z = static_cast<T>(p._z);
  return *this;
}

// =============================================================================
// Accessors
// =============================================================================

// -----------------------------------------------------------------------------
template <typename T>
inline T& Vector3D<T>::operator ()(int i)
{
  switch (i) {
    case 0: return _x;
    case 1: return _y;
    case 2: return _z;
    default:
      cerr << "Vector3D::operator(): Invalid index " << i << endl;
      exit(1);
  }
}

// -----------------------------------------------------------------------------
template <typename T>
inline T Vector3D<T>::operator ()(int i) const
{
  return const_cast<Vector3D<T> *>(this)->operator()(i);
}

// -----------------------------------------------------------------------------
template <typename T>
inline Vector3D<T>::operator T *()
{
  return &_x;
}

// -----------------------------------------------------------------------------
template <typename T>
inline Vector3D<T>::operator const T *() const
{
  return &_x;
}

// =============================================================================
// 3D vector/integral-valued scalar operators
// =============================================================================

// -----------------------------------------------------------------------------
template <typename T>
inline Vector3D<T>& Vector3D<T>::operator =(int s)
{
  _x = static_cast<T>(s);
  _y = static_cast<T>(s);
  _z = static_cast<T>(s);
  return *this;
}

// -----------------------------------------------------------------------------
template <typename T>
inline Vector3D<T>& Vector3D<T>::operator +=(int s)
{
  _x = static_cast<T>(static_cast<double>(_x) + static_cast<double>(s));
  _y = static_cast<T>(static_cast<double>(_y) + static_cast<double>(s));
  _z = static_cast<T>(static_cast<double>(_z) + static_cast<double>(s));
  return *this;
}

// -----------------------------------------------------------------------------
template <typename T>
inline Vector3D<T>& Vector3D<T>::operator -=(int s)
{
  _x = static_cast<T>(static_cast<double>(_x) - static_cast<double>(s));
  _y = static_cast<T>(static_cast<double>(_y) - static_cast<double>(s));
  _z = static_cast<T>(static_cast<double>(_z) - static_cast<double>(s));
  return *this;
}

// -----------------------------------------------------------------------------
template <typename T>
inline Vector3D<T>& Vector3D<T>::operator *=(int s)
{
  _x = static_cast<T>(static_cast<double>(_x) * static_cast<double>(s));
  _y = static_cast<T>(static_cast<double>(_y) * static_cast<double>(s));
  _z = static_cast<T>(static_cast<double>(_z) * static_cast<double>(s));
  return *this;
}

// -----------------------------------------------------------------------------
template <typename T>
inline Vector3D<T>& Vector3D<T>::operator /=(int s)
{
  _x = static_cast<T>(static_cast<double>(_x) / static_cast<double>(s));
  _y = static_cast<T>(static_cast<double>(_y) / static_cast<double>(s));
  _z = static_cast<T>(static_cast<double>(_z) / static_cast<double>(s));
  return *this;
}

// -----------------------------------------------------------------------------
template <typename T>
inline Vector3D<T> Vector3D<T>::operator +(int s) const
{
  Vector3D<T> r(*this);
  r += s;
  return r;
}

// -----------------------------------------------------------------------------
template <typename T>
inline Vector3D<T> Vector3D<T>::operator -(int s) const
{
  Vector3D<T> r(*this);
  r -= s;
  return r;
}

// -----------------------------------------------------------------------------
template <typename T>
inline Vector3D<T> Vector3D<T>::operator *(int s) const
{
  Vector3D<T> r(*this);
  r *= s;
  return r;
}

// -----------------------------------------------------------------------------
template <typename T>
inline Vector3D<T> Vector3D<T>::operator /(int s) const
{
  Vector3D<T> r(*this);
  r /= s;
  return r;
}

// -----------------------------------------------------------------------------
template <typename T>
inline Vector3D<T> operator +(int s, const Vector3D<T> &v)
{
  return v + s;
}

// -----------------------------------------------------------------------------
template <typename T>
inline Vector3D<T> operator -(int s, const Vector3D<T> &v)
{
  return v - s;
}

// -----------------------------------------------------------------------------
template <typename T>
inline Vector3D<T> operator *(int s, const Vector3D<T> &v)
{
  return v * s;
}

// =============================================================================
// 3D vector/real-valued scalar operators
// =============================================================================

// -----------------------------------------------------------------------------
template <typename T>
inline Vector3D<T>& Vector3D<T>::operator =(double s)
{
  _x = static_cast<T>(s);
  _y = static_cast<T>(s);
  _z = static_cast<T>(s);
  return *this;
}

// -----------------------------------------------------------------------------
template <typename T>
inline Vector3D<T>& Vector3D<T>::operator +=(double s)
{
  _x = static_cast<T>(static_cast<double>(_x) + s);
  _y = static_cast<T>(static_cast<double>(_y) + s);
  _z = static_cast<T>(static_cast<double>(_z) + s);
  return *this;
}

// -----------------------------------------------------------------------------
template <typename T>
inline Vector3D<T>& Vector3D<T>::operator -=(double s)
{
  _x = static_cast<T>(static_cast<double>(_x) - s);
  _y = static_cast<T>(static_cast<double>(_y) - s);
  _z = static_cast<T>(static_cast<double>(_z) - s);
  return *this;
}

// -----------------------------------------------------------------------------
template <typename T>
inline Vector3D<T>& Vector3D<T>::operator *=(double s)
{
  _x = static_cast<T>(static_cast<double>(_x) * s);
  _y = static_cast<T>(static_cast<double>(_y) * s);
  _z = static_cast<T>(static_cast<double>(_z) * s);
  return *this;
}

// -----------------------------------------------------------------------------
template <typename T>
inline Vector3D<T>& Vector3D<T>::operator /=(double s)
{
  _x = static_cast<T>(static_cast<double>(_x) / s);
  _y = static_cast<T>(static_cast<double>(_y) / s);
  _z = static_cast<T>(static_cast<double>(_z) / s);
  return *this;
}

// -----------------------------------------------------------------------------
template <typename T>
inline Vector3D<T> Vector3D<T>::operator +(double s) const
{
  Vector3D<T> r(*this);
  r += s;
  return r;
}

// -----------------------------------------------------------------------------
template <typename T>
inline Vector3D<T> Vector3D<T>::operator -(double s) const
{
  Vector3D<T> r(*this);
  r -= s;
  return r;
}

// -----------------------------------------------------------------------------
template <typename T>
inline Vector3D<T> Vector3D<T>::operator *(double s) const
{
  Vector3D<T> r(*this);
  r *= s;
  return r;
}

// -----------------------------------------------------------------------------
template <typename T>
inline Vector3D<T> Vector3D<T>::operator /(double s) const
{
  Vector3D<T> r(*this);
  r /= s;
  return r;
}

// -----------------------------------------------------------------------------
template <typename T>
inline Vector3D<T> operator +(double s, const Vector3D<T> &v)
{
  return v + s;
}

// -----------------------------------------------------------------------------
template <typename T>
inline Vector3D<T> operator -(double s, const Vector3D<T> &v)
{
  return v - s;
}

// -----------------------------------------------------------------------------
template <typename T>
inline Vector3D<T> operator *(double s, const Vector3D<T> &v)
{
  return v * s;
}

// =============================================================================
// 3D vector/vector operators
// =============================================================================

// -----------------------------------------------------------------------------
template <typename T>
inline Vector3D<T> Vector3D<T>::operator -() const
{
  return Vector3D<T>(-_x, -_y, -_z);
}

// -----------------------------------------------------------------------------
template <typename T1> template <typename T2>
inline Vector3D<T1> &Vector3D<T1>::operator =(const Vector3D<T2> &v)
{
  _x = static_cast<T1>(v._x);
  _y = static_cast<T1>(v._y);
  _z = static_cast<T1>(v._z);
  return *this;
}

// -----------------------------------------------------------------------------
template <typename T1> template <typename T2>
inline Vector3D<T1> &Vector3D<T1>::operator +=(const Vector3D<T2>& v)
{
  _x = static_cast<T1>(static_cast<double>(_x) + static_cast<double>(v._x));
  _y = static_cast<T1>(static_cast<double>(_y) + static_cast<double>(v._y));
  _z = static_cast<T1>(static_cast<double>(_z) + static_cast<double>(v._z));
  return *this;
}

// -----------------------------------------------------------------------------
template <typename T1> template <typename T2>
inline Vector3D<T1> &Vector3D<T1>::operator -=(const Vector3D<T2> &v)
{
  _x = static_cast<T1>(static_cast<double>(_x) - static_cast<double>(v._x));
  _y = static_cast<T1>(static_cast<double>(_y) - static_cast<double>(v._y));
  _z = static_cast<T1>(static_cast<double>(_z) - static_cast<double>(v._z));
  return *this;
}

// -----------------------------------------------------------------------------
template <typename T1> template <typename T2>
inline Vector3D<T1> &Vector3D<T1>::operator *=(const Vector3D<T2> &v)
{
  _x = static_cast<T1>(static_cast<double>(_x) * static_cast<double>(v._x));
  _y = static_cast<T1>(static_cast<double>(_y) * static_cast<double>(v._y));
  _z = static_cast<T1>(static_cast<double>(_z) * static_cast<double>(v._z));
  return *this;
}

// -----------------------------------------------------------------------------
template <typename T1> template <typename T2>
inline Vector3D<T1> &Vector3D<T1>::operator /=(const Vector3D<T2> &v)
{
  _x = static_cast<T1>(static_cast<double>(_x) / static_cast<double>(v._x));
  _y = static_cast<T1>(static_cast<double>(_y) / static_cast<double>(v._y));
  _z = static_cast<T1>(static_cast<double>(_z) / static_cast<double>(v._z));
  return *this;
}

// -----------------------------------------------------------------------------
template <typename T1> template <typename T2>
inline Vector3D<T1> Vector3D<T1>::operator +(const Vector3D<T2> &v) const
{
  Vector3D<T1> r(*this);
  r += v;
  return r;
}

// -----------------------------------------------------------------------------
template <typename T1> template <typename T2>
inline Vector3D<T1> Vector3D<T1>::operator -(const Vector3D<T2> &v) const
{
  Vector3D<T1> r(*this);
  r -= v;
  return r;
}

// -----------------------------------------------------------------------------
template <typename T1> template <typename T2>
inline Vector3D<T1> Vector3D<T1>::operator *(const Vector3D<T2> &v) const
{
  Vector3D<T1> r(*this);
  r *= v;
  return r;
}

// -----------------------------------------------------------------------------
template <typename T1> template <typename T2>
inline Vector3D<T1> Vector3D<T1>::operator /(const Vector3D<T2> &v) const
{
  Vector3D<T1> r(*this);
  r /= v;
  return r;
}

// =============================================================================
// 3D Vector/integral-valued scalar comparison
// =============================================================================

// -----------------------------------------------------------------------------
template <typename T>
inline bool Vector3D<T>::operator ==(int s) const
{
  return (_x == s && _y == s && _z == s);
}

// -----------------------------------------------------------------------------
template <typename T>
inline bool Vector3D<T>::operator !=(int s) const
{
  return !(*this == s);
}

// -----------------------------------------------------------------------------
template <typename T>
inline bool Vector3D<T>::operator <(int s) const
{
  return (_x < s && _y < s && _z < s);
}

// -----------------------------------------------------------------------------
template <typename T>
inline bool Vector3D<T>::operator >(int s) const
{
  return (_x > s && _y > s && _z > s);
}

// -----------------------------------------------------------------------------
template <typename T>
inline bool Vector3D<T>::operator <=(int s) const
{
  return (_x <= s && _y <= s && _z <= s);
}

// -----------------------------------------------------------------------------
template <typename T>
inline bool Vector3D<T>::operator >=(int s) const
{
  return (_x >= s && _y >= s && _z >= s);
}

// =============================================================================
// 3D Vector/real-valued scalar comparison
// =============================================================================

// -----------------------------------------------------------------------------
template <typename T>
inline bool Vector3D<T>::operator ==(double s) const
{
  return (_x == s && _y == s && _z == s);
}

// -----------------------------------------------------------------------------
template <typename T>
inline bool Vector3D<T>::operator !=(double s) const
{
  return !(*this == s);
}

// -----------------------------------------------------------------------------
template <typename T>
inline bool Vector3D<T>::operator <(double s) const
{
  return (_x < s && _y < s && _z < s);
}

// -----------------------------------------------------------------------------
template <typename T>
inline bool Vector3D<T>::operator >(double s) const
{
  return (_x > s && _y > s && _z > s);
}

// -----------------------------------------------------------------------------
template <typename T>
inline bool Vector3D<T>::operator <=(double s) const
{
  return (_x <= s && _y <= s && _z <= s);
}

// -----------------------------------------------------------------------------
template <typename T>
inline bool Vector3D<T>::operator >=(double s) const
{
  return (_x >= s && _y >= s && _z >= s);
}

// =============================================================================
// 3D Vector/vector comparison
// =============================================================================

// -----------------------------------------------------------------------------
template <typename T1> template <typename T2>
inline bool Vector3D<T1>::operator ==(const Vector3D<T2> &v) const
{
  return ((_z == v._z) && (_y == v._y) && (_x == v._x));
}

// -----------------------------------------------------------------------------
template <typename T1> template <typename T2>
inline bool Vector3D<T1>::operator !=(const Vector3D<T2> &v) const
{
  return ((_z != v._z) || (_y != v._y) || (_x != v._x));
}

// -----------------------------------------------------------------------------
template <typename T1> template <typename T2>
inline bool Vector3D<T1>::operator <(const Vector3D<T2> &v) const
{
  return ((_z < v._z) ||
          ((_z == v._z) && (_y < v._y)) ||
          ((_z == v._z) && (_y == v._y) && (_x < v._x)));
}

// -----------------------------------------------------------------------------
template <typename T1> template <typename T2>
inline bool Vector3D<T1>::operator >(const Vector3D<T2> &v) const
{
  return ((_z > v._z) ||
          ((_z == v._z) && (_y > v._y)) ||
          ((_z == v._z) && (_y == v._y) && (_x > v._x)));
}

// -----------------------------------------------------------------------------
template <typename T1> template <typename T2>
inline bool Vector3D<T1>::operator <=(const Vector3D<T2> &v) const
{
  return ((*this < v) || (*this == v));
}

// -----------------------------------------------------------------------------
template <typename T1> template <typename T2>
inline bool Vector3D<T1>::operator >=(const Vector3D<T2> &v) const
{
  return ((*this > v) || (*this == v));
}

// =============================================================================
// 3D vector functions
// =============================================================================

// -----------------------------------------------------------------------------
template <typename T>
inline double Vector3D<T>::SquaredLength() const
{
  return static_cast<double>(_x*_x + _y*_y + _z*_z);
}

// -----------------------------------------------------------------------------
template <typename T>
inline double Vector3D<T>::Length() const
{
  return sqrt(SquaredLength());
}

// -----------------------------------------------------------------------------
template <typename T>
inline void Vector3D<T>::Normalize()
{
  const double length = Length();
  if (length != .0) (*this) /= length;
}

// -----------------------------------------------------------------------------
template<typename T>
inline Vector3D<T>
Vector3D<T>::CrossProduct(const Vector3D<T> &v) const
{
  return Vector3D<T>((_y * v._z - _z * v._y),
                     (_z * v._x - _x * v._z),
                     (_x * v._y - _y * v._x));
}

// -----------------------------------------------------------------------------
template<typename T>
inline Vector3D<T> Vector3D<T>::CrossProduct(const Vector3D<T> &v1,
                                             const Vector3D<T> &v2)
{
  return v1.CrossProduct(v2);
}

// -----------------------------------------------------------------------------
template <typename T>
inline double Vector3D<T>::DotProduct(const Vector3D<T> &v) const
{
  return (_x * v._x + _y * v._y + _z * v._z);
}

// -----------------------------------------------------------------------------
template<typename T>
inline double Vector3D<T>::DotProduct(const Vector3D<T> &v1,
                                      const Vector3D<T> &v2)
{
  return v1.DotProduct(v2);
}

// =============================================================================
// Indexed element access
// =============================================================================

// -----------------------------------------------------------------------------
template <class T>
inline T get(const mirtk::Vector3D<T> &v, int n)
{
  switch (n) {
    case 0: return v._x;
    case 1: return v._y;
    case 2: return v._z;
    default:
      cerr << "Invalid 3D vector element index: " << n << endl;
      cerr << "Set breakpoint in " << __FILE__ << ":" << __LINE__ << " to debug." << endl;
      exit(1);
  }
}

// -----------------------------------------------------------------------------
template <class T>
inline void put(mirtk::Vector3D<T> &v, int n, const T &value)
{
  switch (n) {
    case 0: v._x = value;
    case 1: v._y = value;
    case 2: v._z = value;
    default:
      cerr << "Invalid 3D vector element index: " << n << endl;
      cerr << "Set breakpoint in " << __FILE__ << ":" << __LINE__ << " to debug." << endl;
      exit(1);
  }
}

// =============================================================================
// Element-wise <cmath> functions
// =============================================================================

// -----------------------------------------------------------------------------
template <typename T>
Vector3D<T> pow(const Vector3D<T> &v, int e)
{
  return Vector3D<T>(pow(v._x, e), pow(v._y, e), pow(v._z, e));
}

// -----------------------------------------------------------------------------
template <typename T>
Vector3D<T> pow(const Vector3D<T> &v, double e)
{
  return Vector3D<T>(pow(v._x, e), pow(v._y, e), pow(v._z, e));
}

// -----------------------------------------------------------------------------
template <typename T>
Vector3D<T> sqrt(const Vector3D<T> &v)
{
  return Vector3D<T>(sqrt(v._x), sqrt(v._y), sqrt(v._z));
}


} // namespace mirtk

#endif // MIRTK_Vector3D_H
