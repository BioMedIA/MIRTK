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

#ifndef MIRTK_Vector4D_H
#define MIRTK_Vector4D_H

#include "mirtk/Math.h"
#include "mirtk/Stream.h"


namespace mirtk {


/// Represents a 4D vector
///
/// \note Must be a primitive type which can be treated as an array of
///       four values of type T without virtual function table et al.
///       Thus, keep this type as simple as possible!
template <typename T>
struct Vector4D
{
  typedef T ComponentType;

  // ---------------------------------------------------------------------------
  // Attributes

  T _x; ///< The x component
  T _y; ///< The y component
  T _z; ///< The z component
  T _t; ///< The t component

  // ---------------------------------------------------------------------------

  /** Constructor. */
  Vector4D();

  /** Constructor. */
  Vector4D(T x);

  /** Constructor. */
  Vector4D(T x, T y, T z, T t);

  /** Copy Constructor. */
  Vector4D(const Vector4D &);

  /** Assignment operator */
  Vector4D& operator=(T s);

  /** Assignment operator */
  Vector4D& operator=(const Vector4D& v);

  /** Unary negation operator. */
  Vector4D operator-() const;

  /** Operator for multiplying by a scalar. */
  template <typename S>
  Vector4D operator*(S s) const;

  /** Operator for adding a scalar to a vector. */
  template <typename S>
  Vector4D operator+(S s) const;

  /** Operator for adding two vectors. */
  Vector4D operator+(const Vector4D& v) const;

  /** Operator for subtracting a scalar to a vector. */
  template <typename S>
  Vector4D operator-(S s) const;

  /** Operator for subtraction. */
  Vector4D operator-(const Vector4D& v) const;

  /** Operator for multiplying by a scalar. */
  template <typename S>
  Vector4D& operator*=(S s);

  /** Operator for multiplying by a vector. */
  Vector4D& operator*=(const Vector4D& v);

  /** Operator for adding a scalar. */
  template <typename S>
  Vector4D& operator+=(S s);

  /** Operator for adding a vector. */
  Vector4D& operator+=(const Vector4D& v);

  /** Operator for subtracting a scalar. */
  template <typename S>
  Vector4D& operator-=(S s);

  /** Operator for subtracting a vector. */
  Vector4D& operator-=(const Vector4D& v);

  /** Operator for testing equality of two vectors. */
  bool operator==(const Vector4D& v) const;

  /** Operator for testing non-equality of two vector. */
  bool operator!=(const Vector4D& v) const;

  /** Operator for comparing sizes of vectors. */
  bool operator<(const Vector4D& v) const;

  /** Operator for comparing sizes of vectors. */
  bool operator>(const Vector4D& v) const;

  /** Operator for comparing sizes of vectors. */
  bool operator<=(const Vector4D& v) const;

  /** Operator for comparing sizes of vectors. */
  bool operator>=(const Vector4D& v) const;

  /** Operator for dividing one vector by another. */
  Vector4D& operator/=(const Vector4D& v);

  /** Operator for dividing one vector by another. */
  Vector4D operator/(const Vector4D& v) const;

  /** Operator for dividing a vector by a scalar. */
  template <typename S>
  Vector4D& operator/=(S s);

  /** Operator for dividing a vector by a scalar. */
  template <typename S>
  Vector4D operator/(S s) const;

  /** Normalizes the vector. */
  void Normalize();

  /** Returns the length of the vector. */
  double Length() const;

  /** Takes the dot-product of two vectors. */
  static double DotProduct(const Vector4D& v1, const Vector4D& v2);

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
template <typename T> inline Vector4D<T>::Vector4D()
{
  _x = 0;
  _y = 0;
  _z = 0;
  _t = 0;
}

// -----------------------------------------------------------------------------
template <typename T> inline Vector4D<T>::Vector4D(T x)
{
  _x = x;
  _y = x;
  _z = x;
  _t = x;
}

// -----------------------------------------------------------------------------
template <typename T> inline Vector4D<T>::Vector4D(T x, T y, T z, T t)
{
  _x = x;
  _y = y;
  _z = z;
  _t = t;
}

// -----------------------------------------------------------------------------
template <typename T> inline Vector4D<T>::Vector4D(const Vector4D &v)
{
  _x = v._x;
  _y = v._y;
  _z = v._z;
  _t = v._t;
}

// -----------------------------------------------------------------------------
template <typename T> inline Vector4D<T>& Vector4D<T>::operator=(T s)
{
  _x = s;
  _y = s;
  _z = s;
  _t = s;
  return *this;
}

// -----------------------------------------------------------------------------
template <typename T> inline Vector4D<T>& Vector4D<T>::operator=(const Vector4D<T>& v)
{
  _x = v._x;
  _y = v._y;
  _z = v._z;
  _t = v._t;
  return *this;
}

// -----------------------------------------------------------------------------
template <typename T> template <typename S> inline Vector4D<T> Vector4D<T>::operator*(S s) const
{
  Vector4D<T> r;

  r._x = _x*s;
  r._y = _y*s;
  r._z = _z*s;
  r._t = _t*s;

  return r;
}

// -----------------------------------------------------------------------------
template <typename T> template <typename S> inline Vector4D<T> Vector4D<T>::operator+(S s) const
{
  Vector4D<T> r;

  r._x = _x + s;
  r._y = _y + s;
  r._z = _z + s;
  r._t = _t + s;

  return r;
}

// -----------------------------------------------------------------------------
template <typename T> inline Vector4D<T> Vector4D<T>::operator+(const Vector4D<T>& v) const
{
  Vector4D<T> r;

  r._x = _x + v._x;
  r._y = _y + v._y;
  r._z = _z + v._z;
  r._t = _t + v._t;

  return r;
}

// -----------------------------------------------------------------------------
template <typename T> template <typename S> inline Vector4D<T> Vector4D<T>::operator-(S s) const
{
  Vector4D<T> r;

  r._x = _x - s;
  r._y = _y - s;
  r._z = _z - s;
  r._t = _t - s;

  return r;
}

// -----------------------------------------------------------------------------
template <typename T> inline Vector4D<T> Vector4D<T>::operator-(const Vector4D<T>& v) const
{
  Vector4D<T> r;

  r._x = _x - v._x;
  r._y = _y - v._y;
  r._z = _z - v._z;
  r._t = _t - v._t;

  return r;
}

// -----------------------------------------------------------------------------
template <typename T> template <typename S> inline Vector4D<T>& Vector4D<T>::operator*=(S s)
{
  _x *= s;
  _y *= s;
  _z *= s;
  _t *= s;

  return *this;
}

// -----------------------------------------------------------------------------
template <typename T> inline Vector4D<T>& Vector4D<T>::operator*=(const Vector4D<T>& v)
{
  _x *= v._x;
  _y *= v._y;
  _z *= v._z;
  _t *= v._t;

  return *this;
}

// -----------------------------------------------------------------------------
template <typename T> template <typename S> inline Vector4D<T>& Vector4D<T>::operator+=(S s)
{
  _x += s;
  _y += s;
  _z += s;
  _t += s;

  return *this;
}

// -----------------------------------------------------------------------------
template <typename T>
inline Vector4D<T> Vector4D<T>::operator -() const
{
  return Vector4D<T>(-_x, -_y, -_z, -_t);
}

// -----------------------------------------------------------------------------
template <typename T> inline Vector4D<T>& Vector4D<T>::operator+=(const Vector4D<T>& v)
{
  _x += v._x;
  _y += v._y;
  _z += v._z;
  _t += v._t; 

  return *this;
}

// -----------------------------------------------------------------------------
template <typename T> template <typename S> inline Vector4D<T>& Vector4D<T>::operator-=(S s)
{
  _x -= s;
  _y -= s;
  _z -= s;
  _t -= s;

  return *this;
}

// -----------------------------------------------------------------------------
template <typename T> inline Vector4D<T>& Vector4D<T>::operator-=(const Vector4D<T>& v)
{
  _x -= v._x;
  _y -= v._y;
  _z -= v._z;
  _t -= v._t;

  return *this;
}

// -----------------------------------------------------------------------------
template <typename T> inline bool Vector4D<T>::operator==(const Vector4D<T>& v) const
{
  return ((_z == v._z) && (_y == v._y) && (_x == v._x) && (_t == v._t));
}

// -----------------------------------------------------------------------------
template <typename T> inline bool Vector4D<T>::operator!=(const Vector4D<T>& v) const
{
  return ((_z != v._z) || (_y != v._y) || (_x != v._x) || (_t != v._t));
}

// -----------------------------------------------------------------------------
template <typename T> inline bool Vector4D<T>::operator<(const Vector4D<T>& v) const
{
  return ( (_t <  v._t) ||
          ((_t == v._t) && (_z <  v._z)) ||
          ((_t == v._t) && (_z == v._z) && (_y <  v._y)) ||
          ((_t == v._t) && (_z == v._z) && (_y == v._y) && (_x <  v._x)));
}

// -----------------------------------------------------------------------------
template <typename T> inline bool Vector4D<T>::operator>(const Vector4D<T>& v) const
{
  return ( (_t >  v._t) ||
          ((_t == v._t) && (_z >  v._z)) ||
          ((_t == v._t) && (_z == v._z) && (_y >  v._y)) ||
          ((_t == v._t) && (_z == v._z) && (_y == v._y) && (_x >  v._x)));
}

// -----------------------------------------------------------------------------
template <typename T> inline bool Vector4D<T>::operator<=(const Vector4D<T>& v) const
{
  return ((*this < v) || (*this == v));
}

// -----------------------------------------------------------------------------
template <typename T> inline bool Vector4D<T>::operator>=(const Vector4D<T>& v) const
{
  return ((*this > v) || (*this == v));
}

// -----------------------------------------------------------------------------
template <typename T> template <typename S> inline Vector4D<T>& Vector4D<T>::operator/=(S s)
{
  _x /= s;
  _y /= s;
  _z /= s;
  _t /= s;

  return *this;
}

// -----------------------------------------------------------------------------
template <typename T> template <typename S> inline Vector4D<T> Vector4D<T>::operator/(S s) const
{
  return Vector4D<T>(_x/s, _y/s, _z/s, _t/s);
}

// -----------------------------------------------------------------------------
template <typename T> inline void Vector4D<T>::Normalize()
{
  double length = sqrt(static_cast<double>(_x*_x + _y*_y + _z*_z + _t*_t));
  if (length != 0) {
    _x = static_cast<T>(_x/length);
    _y = static_cast<T>(_y/length);
    _z = static_cast<T>(_z/length);
    _t = static_cast<T>(_t/length);
  }
}

// -----------------------------------------------------------------------------
template <typename T> inline double Vector4D<T>::Length() const
{
  return sqrt(static_cast<double>(_x*_x + _y*_y + _z*_z + _t*_t));
}

// -----------------------------------------------------------------------------
template <typename T> inline double Vector4D<T>::DotProduct(const Vector4D<T>& v1, const Vector4D<T>& v2)
{
  return v1._x*v2._x + v1._y*v2._y + v1._z*v2._z + v1._t*v2._t;
}

// =============================================================================
// Indexed element access
// =============================================================================

// -----------------------------------------------------------------------------
template <class T>
inline T get(const Vector4D<T> &v, int n)
{
  switch (n) {
    case 0: return v._x;
    case 1: return v._y;
    case 2: return v._z;
    case 3: return v._t;
    default:
      cerr << "Invalid 4D vector element index: " << n << endl;
      cerr << "Set breakpoint in " << __FILE__ << ":" << __LINE__ << " to debug." << endl;
      exit(1);
  }
}

// -----------------------------------------------------------------------------
template <class T>
inline T put(Vector4D<T> &v, int n, const T &value)
{
  switch (n) {
    case 0: v._x = value;
    case 1: v._y = value;
    case 2: v._z = value;
    case 3: v._t = value;
    default:
      cerr << "Invalid 4D vector element index: " << n << endl;
      cerr << "Set breakpoint in " << __FILE__ << ":" << __LINE__ << " to debug." << endl;
      exit(1);
  }
}

// =============================================================================
// Element-wise <cmath> functions
// =============================================================================

// -----------------------------------------------------------------------------
template <typename T>
Vector4D<T> pow(const Vector4D<T> &v, int e)
{
  return Vector4D<T>(pow(v._x, e), pow(v._y, e), pow(v._z, e), pow(v._t, e));
}

// -----------------------------------------------------------------------------
template <typename T>
Vector4D<T> pow(const Vector4D<T> &v, double e)
{
  return Vector4D<T>(pow(v._x, e), pow(v._y, e), pow(v._z, e), pow(v._t, e));
}

// -----------------------------------------------------------------------------
template <typename T>
Vector4D<T> sqrt(const Vector4D<T> &v)
{
  return Vector4D<T>(sqrt(v._x), sqrt(v._y), sqrt(v._z), sqrt(v._t));
}


} // namespace mirtk

#endif // __IRTKVECTOR4D_H
