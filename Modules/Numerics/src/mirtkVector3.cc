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

#include <mirtkVector3.h>

#include <mirtkMath.h>
#include <mirtkStream.h>


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
Vector3::Vector3 (double fX, double fY, double fZ)
{
  x = fX;
  y = fY;
  z = fZ;
}

//----------------------------------------------------------------------------
Vector3::Vector3 (double afCoordinate[3])
{
  x = afCoordinate[0];
  y = afCoordinate[1];
  z = afCoordinate[2];
}

//----------------------------------------------------------------------------
Vector3::Vector3 (const Vector3& rkVector)
{
  x = rkVector.x;
  y = rkVector.y;
  z = rkVector.z;
}

//----------------------------------------------------------------------------
double& Vector3::operator[] (int i) const
{
  // assert:  0 <= i < 2; x, y, and z are packed into 3*sizeof(double)
  //          bytes
  return (double&) *(&x + i);
}

//----------------------------------------------------------------------------
Vector3::operator double* ()
{
  return &x;
}

//----------------------------------------------------------------------------
Vector3& Vector3::operator= (const Vector3& rkVector)
{
  x = rkVector.x;
  y = rkVector.y;
  z = rkVector.z;
  return *this;
}

//----------------------------------------------------------------------------
bool Vector3::operator== (const Vector3& rkVector) const
{
  return ( x == rkVector.x && y == rkVector.y && z == rkVector.z );
}

//----------------------------------------------------------------------------
bool Vector3::operator!= (const Vector3& rkVector) const
{
  return ( x != rkVector.x || y != rkVector.y || z != rkVector.z );
}

//----------------------------------------------------------------------------
Vector3 Vector3::operator+ (const Vector3& rkVector) const
{
  Vector3 kSum;
  kSum.x = x + rkVector.x;
  kSum.y = y + rkVector.y;
  kSum.z = z + rkVector.z;
  return kSum;
}

//----------------------------------------------------------------------------
Vector3 Vector3::operator- (const Vector3& rkVector) const
{
  Vector3 kDiff;
  kDiff.x = x - rkVector.x;
  kDiff.y = y - rkVector.y;
  kDiff.z = z - rkVector.z;
  return kDiff;
}

//----------------------------------------------------------------------------
Vector3 Vector3::operator* (double fScalar) const
{
  Vector3 kProd;
  kProd.x = fScalar*x;
  kProd.y = fScalar*y;
  kProd.z = fScalar*z;
  return kProd;
}

//----------------------------------------------------------------------------
Vector3 Vector3::operator- () const
{
  Vector3 kNeg;
  kNeg.x = -x;
  kNeg.y = -y;
  kNeg.z = -z;
  return kNeg;
}

//----------------------------------------------------------------------------
Vector3 operator* (double fScalar, const Vector3& rkVector)
{
  Vector3 kProd;
  kProd.x = fScalar*rkVector.x;
  kProd.y = fScalar*rkVector.y;
  kProd.z = fScalar*rkVector.z;
  return kProd;
}

//----------------------------------------------------------------------------
Vector3& Vector3::operator+= (const Vector3& rkVector)
{
  x += rkVector.x;
  y += rkVector.y;
  z += rkVector.z;
  return *this;
}

//----------------------------------------------------------------------------
Vector3& Vector3::operator-= (const Vector3& rkVector)
{
  x -= rkVector.x;
  y -= rkVector.y;
  z -= rkVector.z;
  return *this;
}

//----------------------------------------------------------------------------
Vector3& Vector3::operator*= (double fScalar)
{
  x *= fScalar;
  y *= fScalar;
  z *= fScalar;
  return *this;
}

//----------------------------------------------------------------------------
double Vector3::SquaredLength () const
{
  return x*x + y*y + z*z;
}

//----------------------------------------------------------------------------
double Vector3::Dot (const Vector3& rkVector) const
{
  return x*rkVector.x + y*rkVector.y + z*rkVector.z;
}

//----------------------------------------------------------------------------
Vector3 Vector3::operator/ (double fScalar) const
{
  Vector3 kQuot;

  if ( fScalar != 0.0 ) {
    double fInvScalar = 1.0/fScalar;
    kQuot.x = fInvScalar*x;
    kQuot.y = fInvScalar*y;
    kQuot.z = fInvScalar*z;
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
    x *= fInvScalar;
    y *= fInvScalar;
    z *= fInvScalar;
  } else {
    cerr << "Vector3::operator/=: division by 0" << endl;
    exit(1);
  }

  return *this;
}

//----------------------------------------------------------------------------
double Vector3::Length () const
{
  return sqrt(x*x + y*y + z*z);
}

//----------------------------------------------------------------------------
double Vector3::Unitize (double fTolerance)
{
  double fLength = Length();

  if ( fLength > fTolerance ) {
    double fInvLength = 1.0/fLength;
    x *= fInvLength;
    y *= fInvLength;
    z *= fInvLength;
  } else {
    fLength = 0.0;
  }

  return fLength;
}

//----------------------------------------------------------------------------
Vector3 Vector3::Cross (const Vector3& rkVector) const
{
  Vector3 kCross;

  kCross.x = y*rkVector.z - z*rkVector.y;
  kCross.y = z*rkVector.x - x*rkVector.z;
  kCross.z = x*rkVector.y - y*rkVector.x;

  return kCross;
}

//----------------------------------------------------------------------------
Vector3 Vector3::UnitCross (const Vector3& rkVector) const
{
  Vector3 kCross;

  kCross.x = y*rkVector.z - z*rkVector.y;
  kCross.y = z*rkVector.x - x*rkVector.z;
  kCross.z = x*rkVector.y - y*rkVector.x;
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

  if (abs(rkW.x) >= abs(rkW.y) &&
      abs(rkW.x) >= abs(rkW.z)) {
    rkU.x = -rkW.y;
    rkU.y = +rkW.x;
    rkU.z = 0.0;
  } else {
    rkU.x = 0.0;
    rkU.y = +rkW.z;
    rkU.z = -rkW.y;
  }

  rkU.Unitize();
  rkV = rkW.Cross(rkU);
}


} // namespace mirtk
