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

#ifndef MIRTK_VectorND_H
#define MIRTK_VectorND_H

#include "mirtk/Math.h"
#include "mirtk/Point.h"
#include "mirtk/Vector3.h"


namespace mirtk {


/// Represents a n-D vector of fixed dimensions
///
/// Must be a primitive type which can be treated as an array of n
/// values of type T such that sizeof(VectorND<n, T>) == n * sizeof(T).
/// Thus, this primitive vector type may not have any other data members
/// besides the n vector components. This is required especially
/// when VectorND is used as voxel type of an image and further an
/// externally allocated continuous block of n * sizeof(T) bytes used
/// internally by the image instance which only reinterprets the memory
/// as consecutive VectorND<n, T> instances, i.e.,
///
/// \code
/// const int X   = 256;
/// const int Y   = 256;
/// const int Z   = 128;
/// const int num = X * Y * Z;
/// const int dim = 9;
/// double *data = new double[dim * num];
/// GenericImage<VectorND<dim, double> > image(X, Y, Z, data);
/// \endcode
template <int N, typename T>
struct VectorND
{
  typedef T ComponentType;

  // ---------------------------------------------------------------------------
  // Attributes

  /// Vector components
  T _v[N];

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Construct from scalar
  VectorND(T = T(0));

  /// Construct from C array
  explicit VectorND(const T[N]);

  /// Construct from variable size vector type
  explicit VectorND(const Vector &);

  /// Copy constructor
  template <typename T2>
  VectorND(const VectorND<N, T2> &);

  // ---------------------------------------------------------------------------
  // Accessors

  /// Number of vector components
  static int Rows() { return N; }

  /// Set/get vector component at index
  T &operator ()(int);

  /// Get vector component at index
  T operator ()(int) const;

  /// Cast to C array pointer
  operator T *();

  /// Cast to C array pointer
  operator const T *() const;

  // ---------------------------------------------------------------------------
  // Vector/integral-valued scalar operators

  /// Assign integral valued scalar
  VectorND &operator =(int);

  /// Add integral valued scalar
  VectorND &operator +=(int);

  /// Subtract integral valued scalar
  VectorND &operator -=(int);

  /// Multiply by integral valued scalar
  VectorND &operator *=(int);

  /// Divide by integral valued scalar
  VectorND &operator /=(int);

  /// Add integral valued scalar to vector
  VectorND operator +(int) const;

  /// Subtract integral valued scalar to vector
  VectorND operator -(int) const;

  /// Multiply vector by integral valued scalar
  VectorND operator *(int) const;

  /// Divide vector by integral valued scalar
  VectorND operator /(int) const;

  // ---------------------------------------------------------------------------
  // Vector/real-valued scalar operators

  /// Assign real valued scalar
  VectorND &operator =(double);

  /// Add real valued scalar
  VectorND &operator +=(double);

  /// Subtract real valued scalar
  VectorND &operator -=(double);

  /// Multiply by real valued scalar
  VectorND &operator *=(double);

  /// Divide by real valued scalar
  VectorND &operator /=(double);

  /// Add real valued scalar to vector
  VectorND operator +(double) const;

  /// Subtract real valued scalar to vector
  VectorND operator -(double) const;

  /// Multiply vector by real valued scalar
  VectorND operator *(double) const;

  /// Divide vector by real valued scalar
  VectorND operator /(double) const;

  // ---------------------------------------------------------------------------
  // Vector/vector operators

  /// Unary negation operator
  VectorND operator -() const;

  /// Assignment from other vector
  template <typename T2> VectorND &operator =(const VectorND<N, T2> &);

  /// Addition of other vector
  template <typename T2> VectorND &operator +=(const VectorND<N, T2> &);

  /// Subtraction of other vector
  template <typename T2> VectorND &operator -=(const VectorND<N, T2> &);

  /// Element-wise multiplication with other vector
  template <typename T2> VectorND &operator *=(const VectorND<N, T2> &);

  /// Element-wise division by other vector
  template <typename T2> VectorND &operator /=(const VectorND<N, T2> &);

  /// Addition of two vectors
  template <typename T2> VectorND operator +(const VectorND<N, T2> &) const;

  /// Subtraction of two vectors
  template <typename T2> VectorND operator -(const VectorND<N, T2> &) const;

  /// Element-wise multiplication of two vectors
  template <typename T2> VectorND operator *(const VectorND<N, T2> &) const;

  /// Element-wise division of two vectors
  template <typename T2> VectorND operator /(const VectorND<N, T2> &) const;

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

  /// Operator for testing equality of two vectors
  template <typename T2> bool operator ==(const VectorND<N, T2> &) const;

  /// Operator for testing non-equality of two vector
  template <typename T2> bool operator !=(const VectorND<N, T2> &) const;

  /// Operator for ordering vectors
  template <typename T2> bool operator <(const VectorND<N, T2> &) const;

  /// Operator for ordering vectors
  template <typename T2> bool operator >(const VectorND<N, T2> &) const;

  /// Operator for ordering vectors
  template <typename T2> bool operator <=(const VectorND<N, T2> &) const;

  /// Operator for ordering vectors
  template <typename T2> bool operator >=(const VectorND<N, T2> &) const;

  // ---------------------------------------------------------------------------
  // Other vector functions

  /// Compute squared length of vector
  double SquaredLength() const;

  /// Compute length of vector
  double Length() const;

  /// Normalize vector to length one
  void Normalize();

  /// Dot-product with other vector
  double DotProduct(const VectorND &) const;

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
template <int N, typename T>
inline VectorND<N, T>::VectorND(T s)
{
  for (int i = 0; i < N; ++i) {
    _v[i] = s;
  }
}

// -----------------------------------------------------------------------------
template <int N, typename T>
inline VectorND<N, T>::VectorND(const T v[N])
{
  for (int i = 0; i < N; ++i) {
    _v[i] = v[i];
  }
}

// -----------------------------------------------------------------------------
template <int N, typename T>
inline VectorND<N, T>::VectorND(const Vector &v)
{
  if (v.Rows() != N) {
    Throw(ERR_InvalidArgument, __FUNCTION__, "Vector must have ", N, " components");
  }
  for (int i = 0; i < N; ++i) {
    _v[i] = static_cast<T>(v(i));
  }
}

// -----------------------------------------------------------------------------
template <int N, typename T1> template <typename T2>
inline VectorND<N, T1>::VectorND(const VectorND<N, T2> &v)
{
  for (int i = 0; i < N; ++i) {
    _v[i] = static_cast<T1>(v._v[i]);
  }
}

// =============================================================================
// Accessors
// =============================================================================

// -----------------------------------------------------------------------------
template <int N, typename T>
inline T& VectorND<N, T>::operator ()(int i)
{
  return _v[i];
}

// -----------------------------------------------------------------------------
template <int N, typename T>
inline T VectorND<N, T>::operator ()(int i) const
{
  return _v[i];
}

// -----------------------------------------------------------------------------
template <int N, typename T>
inline VectorND<N, T>::operator T *()
{
  return _v;
}

// -----------------------------------------------------------------------------
template <int N, typename T>
inline VectorND<N, T>::operator const T *() const
{
  return _v;
}

// =============================================================================
// Vector/integral-valued scalar operators
// =============================================================================

// -----------------------------------------------------------------------------
template <int N, typename T>
inline VectorND<N, T>& VectorND<N, T>::operator =(int s)
{
  for (int i = 0; i < N; ++i) {
    _v[i] = static_cast<T>(s);
  }
  return *this;
}

// -----------------------------------------------------------------------------
template <int N, typename T>
inline VectorND<N, T>& VectorND<N, T>::operator +=(int s)
{
  for (int i = 0; i < N; ++i) {
    _v[i] = static_cast<T>(static_cast<double>(_v[i]) + static_cast<double>(s));
  }
  return *this;
}

// -----------------------------------------------------------------------------
template <int N, typename T>
inline VectorND<N, T>& VectorND<N, T>::operator -=(int s)
{
  for (int i = 0; i < N; ++i) {
    _v[i] = static_cast<T>(static_cast<double>(_v[i]) - static_cast<double>(s));
  }
  return *this;
}

// -----------------------------------------------------------------------------
template <int N, typename T>
inline VectorND<N, T>& VectorND<N, T>::operator *=(int s)
{
  for (int i = 0; i < N; ++i) {
    _v[i] = static_cast<T>(static_cast<double>(_v[i]) * static_cast<double>(s));
  }
  return *this;
}

// -----------------------------------------------------------------------------
template <int N, typename T>
inline VectorND<N, T>& VectorND<N, T>::operator /=(int s)
{
  for (int i = 0; i < N; ++i) {
    _v[i] = static_cast<T>(static_cast<double>(_v[i]) / static_cast<double>(s));
  }
  return *this;
}

// -----------------------------------------------------------------------------
template <int N, typename T>
inline VectorND<N, T> VectorND<N, T>::operator +(int s) const
{
  VectorND<N, T> r(*this);
  r += s;
  return r;
}

// -----------------------------------------------------------------------------
template <int N, typename T>
inline VectorND<N, T> VectorND<N, T>::operator -(int s) const
{
  VectorND<N, T> r(*this);
  r -= s;
  return r;
}

// -----------------------------------------------------------------------------
template <int N, typename T>
inline VectorND<N, T> VectorND<N, T>::operator *(int s) const
{
  VectorND<N, T> r(*this);
  r *= s;
  return r;
}

// -----------------------------------------------------------------------------
template <int N, typename T>
inline VectorND<N, T> VectorND<N, T>::operator /(int s) const
{
  VectorND<N, T> r(*this);
  r /= s;
  return r;
}

// -----------------------------------------------------------------------------
template <int N, typename T>
inline VectorND<N, T> operator +(int s, const VectorND<N, T> &v)
{
  return v + s;
}

// -----------------------------------------------------------------------------
template <int N, typename T>
inline VectorND<N, T> operator -(int s, const VectorND<N, T> &v)
{
  return v - s;
}

// -----------------------------------------------------------------------------
template <int N, typename T>
inline VectorND<N, T> operator *(int s, const VectorND<N, T> &v)
{
  return v * s;
}

// =============================================================================
// Vector/real-valued scalar operators
// =============================================================================

// -----------------------------------------------------------------------------
template <int N, typename T>
inline VectorND<N, T>& VectorND<N, T>::operator =(double s)
{
  for (int i = 0; i < N; ++i) {
    _v[i] = static_cast<T>(s);
  }
  return *this;
}

// -----------------------------------------------------------------------------
template <int N, typename T>
inline VectorND<N, T>& VectorND<N, T>::operator +=(double s)
{
  for (int i = 0; i < N; ++i) {
    _v[i] = static_cast<T>(static_cast<double>(_v[i]) + s);
  }
  return *this;
}

// -----------------------------------------------------------------------------
template <int N, typename T>
inline VectorND<N, T>& VectorND<N, T>::operator -=(double s)
{
  for (int i = 0; i < N; ++i) {
    _v[i] = static_cast<T>(static_cast<double>(_v[i]) - s);
  }
  return *this;
}

// -----------------------------------------------------------------------------
template <int N, typename T>
inline VectorND<N, T>& VectorND<N, T>::operator *=(double s)
{
  for (int i = 0; i < N; ++i) {
    _v[i] = static_cast<T>(static_cast<double>(_v[i]) * s);
  }
  return *this;
}

// -----------------------------------------------------------------------------
template <int N, typename T>
inline VectorND<N, T>& VectorND<N, T>::operator /=(double s)
{
  for (int i = 0; i < N; ++i) {
    _v[i] = static_cast<T>(static_cast<double>(_v[i]) / s);
  }
  return *this;
}

// -----------------------------------------------------------------------------
template <int N, typename T>
inline VectorND<N, T> VectorND<N, T>::operator +(double s) const
{
  VectorND<N, T> r(*this);
  r += s;
  return r;
}

// -----------------------------------------------------------------------------
template <int N, typename T>
inline VectorND<N, T> VectorND<N, T>::operator -(double s) const
{
  VectorND<N, T> r(*this);
  r -= s;
  return r;
}

// -----------------------------------------------------------------------------
template <int N, typename T>
inline VectorND<N, T> VectorND<N, T>::operator *(double s) const
{
  VectorND<N, T> r(*this);
  r *= s;
  return r;
}

// -----------------------------------------------------------------------------
template <int N, typename T>
inline VectorND<N, T> VectorND<N, T>::operator /(double s) const
{
  VectorND<N, T> r(*this);
  r /= s;
  return r;
}

// -----------------------------------------------------------------------------
template <int N, typename T>
inline VectorND<N, T> operator +(double s, const VectorND<N, T> &v)
{
  return v + s;
}

// -----------------------------------------------------------------------------
template <int N, typename T>
inline VectorND<N, T> operator -(double s, const VectorND<N, T> &v)
{
  return v - s;
}

// -----------------------------------------------------------------------------
template <int N, typename T>
inline VectorND<N, T> operator *(double s, const VectorND<N, T> &v)
{
  return v * s;
}

// =============================================================================
// Vector/vector operators
// =============================================================================

// -----------------------------------------------------------------------------
template <int N, typename T>
inline VectorND<N, T> VectorND<N, T>::operator -() const
{
  VectorND<N, T> v(*this);
  for (int i = 0; i < N; ++i) {
    v._v[i] = -v._v[i];
  }
  return v;
}

// -----------------------------------------------------------------------------
template <int N, typename T1> template <typename T2>
inline VectorND<N, T1> &VectorND<N, T1>::operator =(const VectorND<N, T2> &v)
{
  for (int i = 0; i < N; ++i) {
    v._v[i] = static_cast<T1>(v._v[i]);
  }
  return *this;
}

// -----------------------------------------------------------------------------
template <int N, typename T1> template <typename T2>
inline VectorND<N, T1> &VectorND<N, T1>::operator +=(const VectorND<N, T2>& v)
{
  for (int i = 0; i < N; ++i) {
   _v[i] = static_cast<T1>(static_cast<double>(_v[i]) + static_cast<double>(v._v[i]));
  }
  return *this;
}

// -----------------------------------------------------------------------------
template <int N, typename T1> template <typename T2>
inline VectorND<N, T1> &VectorND<N, T1>::operator -=(const VectorND<N, T2> &v)
{
  for (int i = 0; i < N; ++i) {
   _v[i] = static_cast<T1>(static_cast<double>(_v[i]) - static_cast<double>(v._v[i]));
  }
  return *this;
}

// -----------------------------------------------------------------------------
template <int N, typename T1> template <typename T2>
inline VectorND<N, T1> &VectorND<N, T1>::operator *=(const VectorND<N, T2> &v)
{
  for (int i = 0; i < N; ++i) {
   _v[i] = static_cast<T1>(static_cast<double>(_v[i]) * static_cast<double>(v._v[i]));
  }
  return *this;
}

// -----------------------------------------------------------------------------
template <int N, typename T1> template <typename T2>
inline VectorND<N, T1> &VectorND<N, T1>::operator /=(const VectorND<N, T2> &v)
{
  for (int i = 0; i < N; ++i) {
   _v[i] = static_cast<T1>(static_cast<double>(_v[i]) / static_cast<double>(v._v[i]));
  }
  return *this;
}

// -----------------------------------------------------------------------------
template <int N, typename T1> template <typename T2>
inline VectorND<N, T1> VectorND<N, T1>::operator +(const VectorND<N, T2> &v) const
{
  VectorND<N, T1> r(*this);
  r += v;
  return r;
}

// -----------------------------------------------------------------------------
template <int N, typename T1> template <typename T2>
inline VectorND<N, T1> VectorND<N, T1>::operator -(const VectorND<N, T2> &v) const
{
  VectorND<N, T1> r(*this);
  r -= v;
  return r;
}

// -----------------------------------------------------------------------------
template <int N, typename T1> template <typename T2>
inline VectorND<N, T1> VectorND<N, T1>::operator *(const VectorND<N, T2> &v) const
{
  VectorND<N, T1> r(*this);
  r *= v;
  return r;
}

// -----------------------------------------------------------------------------
template <int N, typename T1> template <typename T2>
inline VectorND<N, T1> VectorND<N, T1>::operator /(const VectorND<N, T2> &v) const
{
  VectorND<N, T1> r(*this);
  r /= v;
  return r;
}

// =============================================================================
// Vector/real-valued scalar comparison
// =============================================================================

// -----------------------------------------------------------------------------
template <int N, typename T>
inline bool VectorND<N, T>::operator ==(double s) const
{
  for (int i = 0; i < N; ++i) {
    if (!AreEqual(static_cast<double>(_v[i]), s)) return false;
  }
}

// -----------------------------------------------------------------------------
template <int N, typename T>
inline bool VectorND<N, T>::operator !=(double s) const
{
  return !(*this == s);
}

// -----------------------------------------------------------------------------
template <int N, typename T>
inline bool VectorND<N, T>::operator <(double s) const
{
  for (int i = 0; i < N; ++i) {
    if (!(static_cast<double>(_v[i]) < s)) return false;
  }
  return true;
}

// -----------------------------------------------------------------------------
template <int N, typename T>
inline bool VectorND<N, T>::operator >(double s) const
{
  for (int i = 0; i < N; ++i) {
    if (!(static_cast<double>(_v[i]) > s)) return false;
  }
  return true;
}

// -----------------------------------------------------------------------------
template <int N, typename T>
inline bool VectorND<N, T>::operator <=(double s) const
{
  return !(*this > s);
}

// -----------------------------------------------------------------------------
template <int N, typename T>
inline bool VectorND<N, T>::operator >=(double s) const
{
  return !(*this < s);
}

// =============================================================================
// Vector/integral-valued scalar comparison
// =============================================================================

// -----------------------------------------------------------------------------
template <int N, typename T>
inline bool VectorND<N, T>::operator ==(int s) const
{
  return (*this == static_cast<double>(s));
}

// -----------------------------------------------------------------------------
template <int N, typename T>
inline bool VectorND<N, T>::operator !=(int s) const
{
  return !(*this == s);
}

// -----------------------------------------------------------------------------
template <int N, typename T>
inline bool VectorND<N, T>::operator <(int s) const
{
  return (*this < static_cast<double>(s));
}

// -----------------------------------------------------------------------------
template <int N, typename T>
inline bool VectorND<N, T>::operator >(int s) const
{
  return (*this > static_cast<double>(s));
}

// -----------------------------------------------------------------------------
template <int N, typename T>
inline bool VectorND<N, T>::operator <=(int s) const
{
  return !(*this > s);
}

// -----------------------------------------------------------------------------
template <int N, typename T>
inline bool VectorND<N, T>::operator >=(int s) const
{
  return !(*this < s);
}

// =============================================================================
// Vector/vector comparison
// =============================================================================

// -----------------------------------------------------------------------------
template <int N, typename T1> template <typename T2>
inline bool VectorND<N, T1>::operator ==(const VectorND<N, T2> &v) const
{
  for (int i = 0; i < N; ++i) {
    if (!AreEqual(static_cast<double>(_v[i]), static_cast<double>(v._v[i]))) return false;
  }
  return true;
}

// -----------------------------------------------------------------------------
template <int N, typename T1> template <typename T2>
inline bool VectorND<N, T1>::operator !=(const VectorND<N, T2> &v) const
{
  return !(*this == v);
}

// -----------------------------------------------------------------------------
template <int N, typename T1> template <typename T2>
inline bool VectorND<N, T1>::operator <(const VectorND<N, T2> &v) const
{
  for (int i = N - 1; i >= 0; --i) {
    if (static_cast<double>(_v[i]) < static_cast<double>(v._v[i])) return true;
    if (!AreEqual(static_cast<double>(_v[i]), static_cast<double>(v._v[i]))) break;
  }
  return false;
}

// -----------------------------------------------------------------------------
template <int N, typename T1> template <typename T2>
inline bool VectorND<N, T1>::operator >(const VectorND<N, T2> &v) const
{
  return !(*this == v) && !(*this < v);
}

// -----------------------------------------------------------------------------
template <int N, typename T1> template <typename T2>
inline bool VectorND<N, T1>::operator <=(const VectorND<N, T2> &v) const
{
  return (*this < v) || (*this == v);
}

// -----------------------------------------------------------------------------
template <int N, typename T1> template <typename T2>
inline bool VectorND<N, T1>::operator >=(const VectorND<N, T2> &v) const
{
  return (*this > v) || (*this == v);
}

// =============================================================================
// Vector functions
// =============================================================================

// -----------------------------------------------------------------------------
template <int N, typename T>
inline double VectorND<N, T>::DotProduct(const VectorND<N, T> &v) const
{
  double dp = 0.;
  for (int i = 0; i < N; ++i) {
    dp += static_cast<double>(_v[i]) * static_cast<double>(v._v[i]);
  }
  return dp;
}

// -----------------------------------------------------------------------------
template <int N, typename T>
inline double VectorND<N, T>::SquaredLength() const
{
  return DotProduct(*this);
}

// -----------------------------------------------------------------------------
template <int N, typename T>
inline double VectorND<N, T>::Length() const
{
  return sqrt(SquaredLength());
}

// -----------------------------------------------------------------------------
template <int N, typename T>
inline void VectorND<N, T>::Normalize()
{
  const double length = Length();
  if (!IsZero(length)) (*this) /= length;
}

// =============================================================================
// Indexed element access
// =============================================================================

// -----------------------------------------------------------------------------
template <int N, class T>
inline T get(const VectorND<N, T> &v, int n)
{
  return v(n);
}

// -----------------------------------------------------------------------------
template <int N, class T>
inline void put(VectorND<N, T> &v, int n, const T &value)
{
  v(n) = value;
}

// =============================================================================
// Element-wise <cmath> functions
// =============================================================================

// -----------------------------------------------------------------------------
template <int N, typename T>
VectorND<N, T> pow(const VectorND<N, T> &v, int e)
{
  VectorND<N, T> w;
  for (int i = 0; i < N; ++i) {
    w._v[i] = pow(v._v[i], e);
  }
  return w;
}

// -----------------------------------------------------------------------------
template <int N, typename T>
VectorND<N, T> pow(const VectorND<N, T> &v, double e)
{
  VectorND<N, T> w;
  for (int i = 0; i < N; ++i) {
    w._v[i] = pow(v._v[i], e);
  }
  return w;
}

// -----------------------------------------------------------------------------
template <int N, typename T>
VectorND<N, T> sqrt(const VectorND<N, T> &v)
{
  VectorND<N, T> w;
  for (int i = 0; i < N; ++i) {
    w._v[i] = sqrt(v._v[i]);
  }
  return w;
}


} // namespace mirtk

#endif // MIRTK_VectorND_H
