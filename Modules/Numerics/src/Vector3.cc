/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2015 Imperial College London
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

#include "mirtk/Vector3.h"

#include "mirtk/Assert.h"
#include "mirtk/Math.h"
#include "mirtk/Stream.h"
#include "mirtk/Point.h"


namespace mirtk {


//----------------------------------------------------------------------------
const Vector3 Vector3::ZERO(0,0,0);
const Vector3 Vector3::UNIT_X(1,0,0);
const Vector3 Vector3::UNIT_Y(0,1,0);
const Vector3 Vector3::UNIT_Z(0,0,1);

//----------------------------------------------------------------------------
Vector3::Vector3 ()
{
  // For efficiency in construction of large arrays of vectors, the
  // default constructor does not initialize the vector.
}

//----------------------------------------------------------------------------
Vector3::Vector3 (double fScalar)
{
  _x = _y = _z = fScalar;
}

//----------------------------------------------------------------------------
Vector3::Vector3 (double fX, double fY, double fZ)
{
  _x = fX;
  _y = fY;
  _z = fZ;
}

//----------------------------------------------------------------------------
Vector3::Vector3 (const double afCoordinate[3])
{
  _x = afCoordinate[0];
  _y = afCoordinate[1];
  _z = afCoordinate[2];
}

//----------------------------------------------------------------------------
Vector3::Vector3 (const Vector3& rkVector)
{
  _x = rkVector._x;
  _y = rkVector._y;
  _z = rkVector._z;
}

//----------------------------------------------------------------------------
Vector3::Vector3 (const Point& p)
{
  _x = p._x;
  _y = p._y;
  _z = p._z;
}

//----------------------------------------------------------------------------
double& Vector3::operator[] (int i)
{
  // assert:  0 <= i < 2; x, y, and z are packed into 3*sizeof(double)
  //          bytes
  mirtkAssert(&_y == (&_x + 1), "x, y, and z are packed into 3*sizeof(double) bytes");
  mirtkAssert(&_z == (&_x + 2), "x, y, and z are packed into 3*sizeof(double) bytes");
  return *(&_x + i);
}

//----------------------------------------------------------------------------
const double& Vector3::operator[] (int i) const
{
  // assert:  0 <= i < 2; x, y, and z are packed into 3*sizeof(double)
  //          bytes
  mirtkAssert(&_y == (&_x + 1), "x, y, and z are packed into 3*sizeof(double) bytes");
  mirtkAssert(&_z == (&_x + 2), "x, y, and z are packed into 3*sizeof(double) bytes");
  return *(&_x + i);
}

//----------------------------------------------------------------------------
Vector3::operator double* ()
{
  mirtkAssert(&_y == (&_x + 1), "x, y, and z are packed into 3*sizeof(double) bytes");
  mirtkAssert(&_z == (&_x + 2), "x, y, and z are packed into 3*sizeof(double) bytes");
  return &_x;
}

//----------------------------------------------------------------------------
Vector3::operator const double* () const
{
  mirtkAssert(&_y == (&_x + 1), "x, y, and z are packed into 3*sizeof(double) bytes");
  mirtkAssert(&_z == (&_x + 2), "x, y, and z are packed into 3*sizeof(double) bytes");
  return &_x;
}

//----------------------------------------------------------------------------
Vector3& Vector3::operator= (const Vector3& rkVector)
{
  _x = rkVector._x;
  _y = rkVector._y;
  _z = rkVector._z;
  return *this;
}

//----------------------------------------------------------------------------
Vector3& Vector3::operator= (double fScalar)
{
  _x = _y = _z = fScalar;
  return *this;
}

//----------------------------------------------------------------------------
bool Vector3::operator== (const Vector3& rkVector) const
{
  return (_x == rkVector._x && _y == rkVector._y && _z == rkVector._z);
}

//----------------------------------------------------------------------------
bool Vector3::operator!= (const Vector3& rkVector) const
{
  return (_x != rkVector._x || _y != rkVector._y || _z != rkVector._z);
}

//----------------------------------------------------------------------------
Vector3 Vector3::operator+ (const Vector3& rkVector) const
{
  Vector3 kSum;
  kSum._x = _x + rkVector._x;
  kSum._y = _y + rkVector._y;
  kSum._z = _z + rkVector._z;
  return kSum;
}

//----------------------------------------------------------------------------
Vector3 Vector3::operator- (const Vector3& rkVector) const
{
  Vector3 kDiff;
  kDiff._x = _x - rkVector._x;
  kDiff._y = _y - rkVector._y;
  kDiff._z = _z - rkVector._z;
  return kDiff;
}

//----------------------------------------------------------------------------
Vector3 Vector3::operator* (double fScalar) const
{
  Vector3 kProd;
  kProd._x = fScalar*_x;
  kProd._y = fScalar*_y;
  kProd._z = fScalar*_z;
  return kProd;
}

//----------------------------------------------------------------------------
Vector3 Vector3::operator- () const
{
  Vector3 kNeg;
  kNeg._x = -_x;
  kNeg._y = -_y;
  kNeg._z = -_z;
  return kNeg;
}

//----------------------------------------------------------------------------
Vector3 operator* (double fScalar, const Vector3& rkVector)
{
  Vector3 kProd;
  kProd._x = fScalar*rkVector._x;
  kProd._y = fScalar*rkVector._y;
  kProd._z = fScalar*rkVector._z;
  return kProd;
}

//----------------------------------------------------------------------------
Vector3& Vector3::operator+= (const Vector3& rkVector)
{
  _x += rkVector._x;
  _y += rkVector._y;
  _z += rkVector._z;
  return *this;
}

//----------------------------------------------------------------------------
Vector3& Vector3::operator-= (const Vector3& rkVector)
{
  _x -= rkVector._x;
  _y -= rkVector._y;
  _z -= rkVector._z;
  return *this;
}

//----------------------------------------------------------------------------
Vector3& Vector3::operator*= (double fScalar)
{
  _x *= fScalar;
  _y *= fScalar;
  _z *= fScalar;
  return *this;
}

//----------------------------------------------------------------------------
double Vector3::SquaredLength () const
{
  return _x*_x + _y*_y + _z*_z;
}

//----------------------------------------------------------------------------
double Vector3::Dot (const Vector3& rkVector) const
{
  return _x*rkVector._x + _y*rkVector._y + _z*rkVector._z;
}

//----------------------------------------------------------------------------
Vector3 Vector3::operator/ (double fScalar) const
{
  Vector3 kQuot;

  if ( fScalar != 0.0 ) {
    double fInvScalar = 1.0/fScalar;
    kQuot._x = fInvScalar*_x;
    kQuot._y = fInvScalar*_y;
    kQuot._z = fInvScalar*_z;
    return kQuot;
  } else {
    cerr << "Vector3::operator/: division by 0" << endl;
    exit(1);
  }
}

//----------------------------------------------------------------------------
Vector3& Vector3::operator/= (double fScalar)
{
  if ( fScalar != 0.0 ) {
    double fInvScalar = 1.0/fScalar;
    _x *= fInvScalar;
    _y *= fInvScalar;
    _z *= fInvScalar;
  } else {
    cerr << "Vector3::operator/=: division by 0" << endl;
    exit(1);
  }

  return *this;
}

//----------------------------------------------------------------------------
double Vector3::Length () const
{
  return sqrt(_x*_x + _y*_y + _z*_z);
}

//----------------------------------------------------------------------------
double Vector3::Normalize(double fTolerance)
{
  double fLength = Length();

  if ( fLength > fTolerance ) {
    double fInvLength = 1.0/fLength;
    _x *= fInvLength;
    _y *= fInvLength;
    _z *= fInvLength;
  } else {
    fLength = 0.0;
  }

  return fLength;
}

//----------------------------------------------------------------------------
double Vector3::Unitize(double fTolerance)
{
  return Normalize(fTolerance);
}

//----------------------------------------------------------------------------
Vector3 Vector3::Cross (const Vector3& rkVector) const
{
  Vector3 kCross;

  kCross._x = _y*rkVector._z - _z*rkVector._y;
  kCross._y = _z*rkVector._x - _x*rkVector._z;
  kCross._z = _x*rkVector._y - _y*rkVector._x;

  return kCross;
}

//----------------------------------------------------------------------------
Vector3 Vector3::UnitCross (const Vector3& rkVector) const
{
  Vector3 kCross;

  kCross._x = _y*rkVector._z - _z*rkVector._y;
  kCross._y = _z*rkVector._x - _x*rkVector._z;
  kCross._z = _x*rkVector._y - _y*rkVector._x;
  kCross.Unitize();

  return kCross;
}

//----------------------------------------------------------------------------
void Vector3::Orthonormalize (Vector3 akVector[3])
{
  // If the input vectors are v0, v1, and v2, then the Gram-Schmidt
  // orthonormalization produces vectors u0, u1, and u2 as follows,
  //
  //   u0 = v0/|v0|
  //   u1 = (v1-(u0*v1)u0)/|v1-(u0*v1)u0|
  //   u2 = (v2-(u0*v2)u0-(u1*v2)u1)/|v2-(u0*v2)u0-(u1*v2)u1|
  //
  // where |A| indicates length of vector A and A*B indicates dot
  // product of vectors A and B.

  // compute u0
  akVector[0].Unitize();

  // compute u1
  double fDot0 = akVector[0].Dot(akVector[1]);
  akVector[1] -= fDot0*akVector[0];
  akVector[1].Unitize();

  // compute u2
  double fDot1 = akVector[1].Dot(akVector[2]);
  fDot0 = akVector[0].Dot(akVector[2]);
  akVector[2] -= fDot0*akVector[0] + fDot1*akVector[1];
  akVector[2].Unitize();
}

//----------------------------------------------------------------------------
void Vector3::GenerateOrthonormalBasis(Vector3& rkU, Vector3& rkV,
                                       Vector3& rkW, bool bUnitLengthW)
{
  if ( !bUnitLengthW )
    rkW.Unitize();

  if (abs(rkW._x) >= abs(rkW._y) &&
      abs(rkW._x) >= abs(rkW._z)) {
    rkU._x = -rkW._y;
    rkU._y = +rkW._x;
    rkU._z = 0.0;
  } else {
    rkU._x = 0.0;
    rkU._y = +rkW._z;
    rkU._z = -rkW._y;
  }

  rkU.Unitize();
  rkV = rkW.Cross(rkU);
}


} // namespace mirtk
