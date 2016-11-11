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

#ifndef MIRTK_Vector3_H
#define MIRTK_Vector3_H

namespace mirtk {


class Point;


/**
 * 3D vector
 *
 * \todo Merge/replace uses of this class with Vector3D
 */
class Vector3
{
public:
  // construction
  Vector3 ();
  Vector3 (double fScalar);
  Vector3 (double fX, double fY, double fZ);
  Vector3 (const double afCoordinate[3]);
  Vector3 (const Vector3& rkVector);
  Vector3 (const Point &);

  // member access
  // (allows V._x or V(0) or V[0], V._y or V(1) or V[1], V._z or V(2) or or V[2])
  double _x, _y, _z;
  double& operator[] (int i);
  const double& operator[] (int i) const;
  operator double* ();
  operator const double* () const;

  // assignment and comparison
  Vector3& operator= (const Vector3& rkVector);
  Vector3& operator= (double fScalar);
  bool operator== (const Vector3& rkVector) const;
  bool operator!= (const Vector3& rkVector) const;

  // arithmetic operations
  Vector3 operator+ (const Vector3& rkVector) const;
  Vector3 operator- (const Vector3& rkVector) const;
  Vector3 operator* (double fScalar) const;
  Vector3 operator/ (double fScalar) const;
  Vector3 operator- () const;
  friend Vector3 operator* (double fScalar, const Vector3& rkVector);

  // arithmetic updates
  Vector3& operator+= (const Vector3& rkVector);
  Vector3& operator-= (const Vector3& rkVector);
  Vector3& operator*= (double fScalar);
  Vector3& operator/= (double fScalar);

  // vector operations
  double Length () const;
  double SquaredLength () const;
  double Dot (const Vector3& rkVector) const;
  double Normalize (double fTolerance = 1e-06);
  double Unitize (double fTolerance = 1e-06); ///< \deprecated Use Normalize instead
  Vector3 Cross (const Vector3& rkVector) const;
  Vector3 UnitCross (const Vector3& rkVector) const;

  // Gram-Schmidt orthonormalization.
  static void Orthonormalize (Vector3 akVector[3]);

  // Input W must be initialize to a nonzero vector, output is {U,V,W}
  // an orthonormal basis.  A hint is provided about whether or not W
  // is already unit length.
  static void GenerateOrthonormalBasis (Vector3& rkU, Vector3& rkV,
                                        Vector3& rkW, bool bUnitLengthW = true);

  // special points
  static const Vector3 ZERO;
  static const Vector3 UNIT_X;
  static const Vector3 UNIT_Y;
  static const Vector3 UNIT_Z;
};


} // namespace mirtk

#endif // MIRTK_Vector3_H
