/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2016 Imperial College London
 * Copyright 2008-2015 Daniel Rueckert, Julia Schnabel
 * Copyright 2016      Andreas Schuh
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

#ifndef MIRTK_Point_H
#define MIRTK_Point_H

#include "mirtk/Object.h"
#include "mirtk/Math.h"


namespace mirtk {


// Forward declaration of Vector in particular, as its header
// includes mirtk/Vector3D.h which in turn includes this file again.
// Therefore, only include mirtk/Vector.h after Point is declared.
class Vector;
class Vector3;
class Matrix;


/**
 * 3D point.
 */
class Point
{
public:

  double _x; ///< x coordinate of Point
  double _y; ///< y coordinate of Point
  double _z; ///< z coordinate of Point

  //
  // Constructors and destructor
  //

  /// Constructor
  Point();

  /// Constructor with three coordinates
  Point(double, double, double);

  /// Constructor with three coordinates
  Point(const double [3]);

  /// Constructor with Point
  Point(const Point &);

  /// Constructor with Vector
  explicit Point(const Vector3 &);

  /// Constructor with Vector
  explicit Point(const Vector &);

  /// Default destructor
  virtual ~Point();

  //
  // Operators for indexed element access
  //

  /// Get reference to i-th point coordinate
  double &operator [](int i);

  /// Get const reference to i-th point coordinate
  const double &operator [](int i) const;

  /// Get reference to i-th point coordinate
  double &operator ()(int i);

  /// Get const reference to i-th point coordinate
  const double &operator ()(int i) const;

  /// Cast to C array pointer
  operator double *();

  /// Cast to C array pointer
  operator const double *() const;

  //
  // Operators for Point
  //

  /// Copy operator for point
  Point& operator =(const Point&);

  /// Substraction operator for point
  Point& operator-=(const Point&);

  /// Addition operator for point
  Point& operator+=(const Point&);

  /// Multiplication operator for point
  Point& operator*=(const Point&);

  /// Division operator for point
  Point& operator/=(const Point&);

  /// Return result of point substraction
  Point operator- (const Point&) const;

  /// Return result of point addition
  Point operator+ (const Point&) const;

  /// Return result of point multiplication
  Point operator* (const Point&) const;

  /// Return result of point division
  Point operator/ (const Point&) const;

  //
  // Operators for comparison
  //

  /// Comparison operator ==
  int operator==(const Point&) const;

  /// Comparison operator != (if USE_STL is defined, negate == operator)
  int operator!=(const Point&) const;

  /// Comparison operator <
  int operator<(const Point&) const;

  /// Comparison operator >
  int operator>(const Point&) const;

  //
  // Operators for double
  //

  /// Assign scalar value to all coordinates
  Point& operator =(double);

  /// Substraction of double
  Point& operator-=(double);

  /// Addition of double
  Point& operator+=(double);

  /// Multiplication with double
  Point& operator*=(double);

  /// Division by double
  Point& operator/=(double);

  // Return result of substraction of double
  Point  operator- (double) const;

  // Return result of addition of double
  Point  operator+ (double) const;

  // Return result of multiplication with double
  Point  operator* (double) const;

  // Return result of division by double
  Point  operator/ (double) const;

  //
  // Operators for Vector3
  //

  /// Copy operator for vectors
  Point& operator =(const Vector3&);

  /// Substraction operator for vectors
  Point& operator-=(const Vector3&);

  /// Addition operator for vectors
  Point& operator+=(const Vector3&);

  /// Multiplication operator for vectors (componentwise)
  Point& operator*=(const Vector3&);

  /// Division operator for vectors (componentwise)
  Point& operator/=(const Vector3&);

  // Return result of vector substraction
  Point operator- (const Vector3&) const;

  // Return result of vector addition
  Point operator+ (const Vector3&) const;

  // Return result of vector multiplication
  Point operator* (const Vector3&) const;

  // Return result of vector division
  Point operator/ (const Vector3&) const;

  //
  // Operators for Vector
  //

  /// Copy operator for vectors
  Point& operator =(const Vector&);

  /// Substraction operator for vectors
  Point& operator-=(const Vector&);

  /// Addition operator for vectors
  Point& operator+=(const Vector&);

  /// Multiplication operator for vectors (componentwise)
  Point& operator*=(const Vector&);

  /// Division operator for vectors (componentwise)
  Point& operator/=(const Vector&);

  // Return result of vector substraction
  Point operator- (const Vector&) const;

  // Return result of vector addition
  Point operator+ (const Vector&) const;

  // Return result of vector multiplication
  Point operator* (const Vector&) const;

  // Return result of vector division
  Point operator/ (const Vector&) const;

  //
  // Operators for Matrix
  //

  /// Point multiplication operator for matrices
  Point& operator*=(const Matrix&);

  /// Return result from Matrix multiplication
  Point operator* (const Matrix&) const;

  //
  // Distance methods
  //

  /// Squared distance from origin
  double SquaredDistance() const;

  /// Squared distance from point
  double SquaredDistance(const Point&) const;

  /// Distance from origin
  double Distance() const;

  /// Distance from point
  double Distance(const Point&) const;

};

} // namespace mirtk

#include "mirtk/Vector.h"
#include "mirtk/Vector3.h"
#include "mirtk/Matrix.h"

namespace mirtk {

////////////////////////////////////////////////////////////////////////////////
// Streaming operators
////////////////////////////////////////////////////////////////////////////////

/// Write point coordinates to output stream
ostream& operator <<(ostream&, const Point&);

/// Read point coordiantes from input stream
istream& operator >>(istream&, Point&);

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
inline Point::Point()
{
  _x = 0;
  _y = 0;
  _z = 0;
}

// -----------------------------------------------------------------------------
inline Point::Point(double x, double y, double z)
{
  _x = x;
  _y = y;
  _z = z;
}

// -----------------------------------------------------------------------------
inline Point::Point(const double p[3])
{
  _x = p[0];
  _y = p[1];
  _z = p[2];
}

// -----------------------------------------------------------------------------
inline Point::Point(const Point& p)
{
  _x = p._x;
  _y = p._y;
  _z = p._z;
}

// -----------------------------------------------------------------------------
inline Point::Point(const Vector3& v)
{
  _x = v._x;
  _y = v._y;
  _z = v._z;
}

// -----------------------------------------------------------------------------
inline Point::Point(const Vector& v)
{
  if ((v.Rows() < 0) || (v.Rows() > 3)) {
    cerr << "Point::Point(const Vector&) Illegal dimension: " << v.Rows() << endl;
    exit(1);
  } else {
    if (v.Rows() == 1) {
      _x = v(0);
      _y = 0;
      _z = 0;
    }
    if (v.Rows() == 2) {
      _x = v(0);
      _y = v(1);
      _z = 0;
    }
    if (v.Rows() == 3) {
      _x = v(0);
      _y = v(1);
      _z = v(2);
    }
  }
}

// -----------------------------------------------------------------------------
inline Point::~Point()
{
}

// -----------------------------------------------------------------------------
inline double &Point::operator [](int i)
{
  switch (i) {
    case 0: return _x;
    case 1: return _y;
    case 2: return _z;
    default:
      cerr << "Point::operator []: Invalid coorindate index: " << i << endl;
      exit(1);
  }
}

// -----------------------------------------------------------------------------
inline const double &Point::operator [](int i) const
{
  switch (i) {
    case 0: return _x;
    case 1: return _y;
    case 2: return _z;
    default:
      cerr << "Point::operator []: Invalid coorindate index: " << i << endl;
      exit(1);
  }
}

// -----------------------------------------------------------------------------
inline double &Point::operator ()(int i)
{
  return operator [](i);
}

// -----------------------------------------------------------------------------
inline const double &Point::operator ()(int i) const
{
  return operator [](i);
}

// -----------------------------------------------------------------------------
inline Point::operator double *()
{
  return &_x;
}

// -----------------------------------------------------------------------------
inline Point::operator const double *() const
{
  return &_x;
}

// -----------------------------------------------------------------------------
inline Point& Point::operator =(const Point& p)
{
  _x = p._x;
  _y = p._y;
  _z = p._z;
  return *this;
}

// -----------------------------------------------------------------------------
inline Point& Point::operator+=(const Point& p)
{
  _x += p._x;
  _y += p._y;
  _z += p._z;
  return *this;
}

// -----------------------------------------------------------------------------
inline Point& Point::operator-=(const Point& p)
{
  _x -= p._x;
  _y -= p._y;
  _z -= p._z;
  return *this;
}

// -----------------------------------------------------------------------------
inline Point& Point::operator*=(const Point& p)
{
  _x *= p._x;
  _y *= p._y;
  _z *= p._z;
  return *this;
}

// -----------------------------------------------------------------------------
inline Point& Point::operator/=(const Point& p)
{
  _x /= p._x;
  _y /= p._y;
  _z /= p._z;
  return *this;
}

// -----------------------------------------------------------------------------
inline Point Point::operator+ (const Point& p) const
{
  Point tmp;
  tmp._x = _x + p._x;
  tmp._y = _y + p._y;
  tmp._z = _z + p._z;
  return tmp;
}

// -----------------------------------------------------------------------------
inline Point Point::operator- (const Point& p) const
{
  Point tmp;
  tmp._x = _x - p._x;
  tmp._y = _y - p._y;
  tmp._z = _z - p._z;
  return tmp;
}

// -----------------------------------------------------------------------------
inline Point Point::operator* (const Point& p) const
{
  Point tmp;
  tmp._x = _x * p._x;
  tmp._y = _y * p._y;
  tmp._z = _z * p._z;
  return tmp;
}

// -----------------------------------------------------------------------------
inline Point Point::operator/ (const Point& p) const
{
  Point tmp;
  tmp._x = _x / p._x;
  tmp._y = _y / p._y;
  tmp._z = _z / p._z;
  return tmp;
}

// -----------------------------------------------------------------------------
inline int Point::operator==(Point const &p) const
{
  return ((_x == p._x) && (_y == p._y) && (_z == p._z));
}

// -----------------------------------------------------------------------------
inline int Point::operator!=(Point const &p) const
{
  return ((_x != p._x) || (_y != p._y) || (_z != p._z));
}

// -----------------------------------------------------------------------------
inline int Point::operator<(Point const& p) const
{
  return ((_x < p._x) && (_y < p._y) && (_z < p._z));
}

// -----------------------------------------------------------------------------
inline int Point::operator>(Point const& p) const
{
  return ((_x > p._x) && (_y > p._y) && (_z > p._z));
}

// -----------------------------------------------------------------------------
inline Point& Point::operator =(double x)
{
  _x = x;
  _y = x;
  _z = x;
  return *this;
}

// -----------------------------------------------------------------------------
inline Point& Point::operator+=(double x)
{
  _x += x;
  _y += x;
  _z += x;
  return *this;
}

// -----------------------------------------------------------------------------
inline Point& Point::operator-=(double x)
{
  _x -= x;
  _y -= x;
  _z -= x;
  return *this;
}

// -----------------------------------------------------------------------------
inline Point& Point::operator/=(double x)
{
  _x /= x;
  _y /= x;
  _z /= x;
  return *this;
}

// -----------------------------------------------------------------------------
inline Point& Point::operator*=(double x)
{
  _x *= x;
  _y *= x;
  _z *= x;
  return *this;
}

// -----------------------------------------------------------------------------
inline Point Point::operator+ (double x) const
{
  Point p;
  p._x = _x + x;
  p._y = _y + x;
  p._z = _z + x;
  return p;
}

// -----------------------------------------------------------------------------
inline Point Point::operator- (double x) const
{
  Point p;
  p._x = _x - x;
  p._y = _y - x;
  p._z = _z - x;
  return p;
}

// -----------------------------------------------------------------------------
inline Point Point::operator* (double x) const
{
  Point p;
  p._x = _x * x;
  p._y = _y * x;
  p._z = _z * x;
  return p;
}

// -----------------------------------------------------------------------------
inline Point Point::operator/ (double x) const
{
  Point p;
  p._x = _x / x;
  p._y = _y / x;
  p._z = _z / x;
  return p;
}

// -----------------------------------------------------------------------------
inline Point& Point::operator =(const Vector3& v)
{
  _x = v._x;
  _y = v._y;
  _z = v._z;
  return *this;
}

// -----------------------------------------------------------------------------
inline Point& Point::operator +=(const Vector3 &v)
{
  _x += v._x;
  _y += v._y;
  _z += v._z;
  return *this;
}

// -----------------------------------------------------------------------------
inline Point& Point::operator -=(const Vector3 &v)
{
  _x -= v._x;
  _y -= v._y;
  _z -= v._z;
  return *this;
}

// -----------------------------------------------------------------------------
inline Point& Point::operator *=(const Vector3 &v)
{
  _x *= v._x;
  _y *= v._y;
  _z *= v._z;
  return *this;
}

// -----------------------------------------------------------------------------
inline Point& Point::operator /=(const Vector3 &v)
{
  _x /= v._x;
  _y /= v._y;
  _z /= v._z;
  return *this;
}

// -----------------------------------------------------------------------------
inline Point Point::operator +(const Vector3 &v) const
{
  return (Point(*this) += v);
}

// -----------------------------------------------------------------------------
inline Point Point::operator -(const Vector3 &v) const
{
  return (Point(*this) -= v);
}

// -----------------------------------------------------------------------------
inline Point Point::operator *(const Vector3 &v) const
{
  return (Point(*this) *= v);
}

// -----------------------------------------------------------------------------
inline Point Point::operator /(const Vector3 &v) const
{
  return (Point(*this) /= v);
}

// -----------------------------------------------------------------------------
inline Point& Point::operator+=(const Vector& v)
{
  if (v.Rows() != 3) {
    cerr << "Point::operator+=(const Vector& v): Size mismatch" << endl;
    exit(1);
  }
  _x += v(0);
  _y += v(1);
  _z += v(2);
  return *this;
}

// -----------------------------------------------------------------------------
inline Point& Point::operator-=(const Vector& v)
{
  if (v.Rows() != 3) {
    cerr << "Point::operator-=(const Vector& v): Size mismatch" << endl;
    exit(1);
  }
  _x -= v(0);
  _y -= v(1);
  _z -= v(2);
  return *this;
}

// -----------------------------------------------------------------------------
inline Point& Point::operator*=(const Vector& v)
{
  if (v.Rows() != 3) {
    cerr << "Point::operator*=(const Vector& v): Size mismatch" << endl;
    exit(1);
  }
  _x *= v(0);
  _y *= v(1);
  _z *= v(2);
  return *this;
}

// -----------------------------------------------------------------------------
inline Point& Point::operator/=(const Vector& v)
{
  if (v.Rows() != 3) {
    cerr << "Point::operator/=(const Vector& v): Size mismatch" << endl;
    exit(1);
  }
  _x /= v(0);
  _y /= v(1);
  _z /= v(2);
  return *this;
}

// -----------------------------------------------------------------------------
inline Point Point::operator+ (const Vector& v) const
{
  Point tmp;
  if (v.Rows() != 3) {
    cerr << "Point::operator+(const Vector& v): Size mismatch" << endl;
    exit(1);
  }
  tmp._x = _x + v(0);
  tmp._y = _y + v(1);
  tmp._z = _z + v(2);
  return tmp;
}

// -----------------------------------------------------------------------------
inline Point Point::operator- (const Vector& v) const
{
  Point tmp;
  if (v.Rows() != 3) {
    cerr << "Point::operator-(const Vector& v): Size mismatch" << endl;
    exit(1);
  }
  tmp._x = _x - v(0);
  tmp._y = _y - v(1);
  tmp._z = _z - v(2);
  return tmp;
}

// -----------------------------------------------------------------------------
inline Point Point::operator* (const Vector& v) const
{
  Point tmp;
  if (v.Rows() != 3) {
    cerr << "Point::operator*(const Vector& v): Size mismatch" << endl;
    exit(1);
  }
  tmp._x = _x * v(0);
  tmp._y = _y * v(1);
  tmp._z = _z * v(2);
  return tmp;
}

// -----------------------------------------------------------------------------
inline Point Point::operator/ (const Vector& v) const
{
  Point tmp;
  if (v.Rows() != 3) {
    cerr << "Point::operator/(const Vector& v): Size mismatch" << endl;
    exit(1);
  }
  tmp._x = _x / v(0);
  tmp._y = _y / v(1);
  tmp._z = _z / v(2);
  return tmp;
}

// -----------------------------------------------------------------------------
inline Point Point::operator* (const Matrix& m) const
{
  Point tmp;
  if ((m.Rows() != 4) && (m.Cols() != 4)) {
    cerr << "Point::operator*(const Matrix& m): Size mismatch" << endl;
    exit(1);
  }
  tmp._x = _x * m(0, 0) + _y * m(0, 1) + _z * m(0, 2) + m(0, 3);
  tmp._y = _x * m(1, 0) + _y * m(1, 1) + _z * m(1, 2) + m(1, 3);
  tmp._z = _x * m(2, 0) + _y * m(2, 1) + _z * m(2, 2) + m(2, 3);
  return tmp;
}

// -----------------------------------------------------------------------------
inline Point& Point::operator*=(const Matrix& m)
{
  Point tmp;
  if ((m.Rows() != 4) && (m.Cols() != 4)) {
    cerr << "Point::operator*(const Matrix& m): Size mismatch" << endl;
    exit(1);
  }
  tmp._x = _x * m(0, 0) + _y * m(0, 1) + _z * m(0, 2) + m(0, 3);
  tmp._y = _x * m(1, 0) + _y * m(1, 1) + _z * m(1, 2) + m(1, 3);
  tmp._z = _x * m(2, 0) + _y * m(2, 1) + _z * m(2, 2) + m(2, 3);
  *this  = tmp;
  return *this;
}

// -----------------------------------------------------------------------------
inline double Point::SquaredDistance() const
{
  return _x*_x + _y*_y  + _z *_z;
}

// -----------------------------------------------------------------------------
inline double Point::SquaredDistance(const Point &p) const
{
  const double dx = _x - p._x;
  const double dy = _y - p._y;
  const double dz = _z - p._z;
  return dx*dx + dy*dy + dz*dz;
}

// -----------------------------------------------------------------------------
inline double Point::Distance() const
{
  return sqrt(SquaredDistance());
}

// -----------------------------------------------------------------------------
inline double Point::Distance(const Point &p) const
{
  return sqrt(SquaredDistance(p));
}


} // namespace mirtk

#endif // MIRTK_Point_H
