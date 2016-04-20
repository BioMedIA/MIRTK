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

#ifndef MIRTK_Matrix3x3_H
#define MIRTK_Matrix3x3_H

#include "mirtk/NumericsExport.h"

#include "mirtk/Vector3.h"


namespace mirtk {


// NOTE.  The (x,y,z) coordinate system is assumed to be right-handed.
// Coordinate axis rotation matrices are of the form
//   RX =    1       0       0
//           0     cos(t) -sin(t)
//           0     sin(t)  cos(t)
// where t > 0 indicates a counterclockwise rotation in the yz-plane
//   RY =  cos(t)    0     sin(t)
//           0       1       0
//        -sin(t)    0     cos(t)
// where t > 0 indicates a counterclockwise rotation in the zx-plane
//   RZ =  cos(t) -sin(t)    0
//         sin(t)  cos(t)    0
//           0       0       1
// where t > 0 indicates a counterclockwise rotation in the xy-plane.

class Matrix3x3
{
public:
  // construction
  Matrix3x3 ();
  Matrix3x3 (int);
  Matrix3x3 (double);
  Matrix3x3 (const double aafEntry[3][3]);
  Matrix3x3 (const Matrix3x3& rkMatrix);
  Matrix3x3 (double fEntry00, double fEntry01, double fEntry02,
             double fEntry10, double fEntry11, double fEntry12,
             double fEntry20, double fEntry21, double fEntry22);

  // member access, allows use of construct mat[r][c]
  double*       operator[] (int iRow);
  const double* operator[] (int iRow) const;
  operator double* ();
  Vector3 GetColumn (int iCol) const;

  // assignment and comparison
  Matrix3x3& operator= (const Matrix3x3& rkMatrix);
  Matrix3x3& operator= (double);
  bool operator== (const Matrix3x3& rkMatrix) const;
  bool operator!= (const Matrix3x3& rkMatrix) const;

  // arithmetic operations
  Matrix3x3 &operator+=(double);
  Matrix3x3 &operator-=(double);
  Matrix3x3 &operator*=(double);
  Matrix3x3 &operator/=(double);

  Matrix3x3 operator+ (double) const;
  Matrix3x3 operator- (double) const;
  Matrix3x3 operator* (double) const;
  Matrix3x3 operator/ (double) const;

  Matrix3x3 &operator+=(const Matrix3x3& rkMatrix);
  Matrix3x3 &operator-=(const Matrix3x3& rkMatrix);
  Matrix3x3 &operator*=(const Matrix3x3& rkMatrix);
  Matrix3x3 &operator/=(const Matrix3x3& rkMatrix);

  Matrix3x3 operator+ (const Matrix3x3& rkMatrix) const;
  Matrix3x3 operator- (const Matrix3x3& rkMatrix) const;
  Matrix3x3 operator* (const Matrix3x3& rkMatrix) const;
  Matrix3x3 operator/ (const Matrix3x3& rkMatrix) const;
  Matrix3x3 operator- () const;

  // matrix * vector [3x3 * 3x1 = 3x1]
  Vector3 operator* (const Vector3& rkVector) const;

  // vector * matrix [1x3 * 3x3 = 1x3]
  friend Vector3 operator* (const Vector3& rkVector,
                                const Matrix3x3& rkMatrix);

  // scalar * matrix
  friend Matrix3x3 operator* (double fScalar, const Matrix3x3& rkMatrix);

  // utilities
  Matrix3x3 Transpose () const;
  bool Inverse (Matrix3x3& rkInverse, double fTolerance = 1e-06) const;
  Matrix3x3 Inverse (double fTolerance = 1e-06) const;
  double Determinant () const;
  Matrix3x3 Adjoint() const;
  double Trace() const;

  // singular value decomposition
  void SingularValueDecomposition (Matrix3x3& rkL, Vector3& rkS,
                                   Matrix3x3& rkR) const;
  void SingularValueComposition (const Matrix3x3& rkL,
                                 const Vector3& rkS, const Matrix3x3& rkR);

  // Gram-Schmidt orthonormalization (applied to columns of rotation matrix)
  void Orthonormalize ();

  // orthogonal Q, diagonal D, upper triangular U stored as (u01,u02,u12)
  void QDUDecomposition (Matrix3x3& rkQ, Vector3& rkD,
                         Vector3& rkU) const;

  double SpectralNorm () const;

  // matrix must be orthonormal
  void ToAxisAngle (Vector3& rkAxis, double& rfRadians) const;
  void FromAxisAngle (const Vector3& rkAxis, double fRadians);

  // The matrix must be orthonormal.  The decomposition is yaw*pitch*roll
  // where yaw is rotation about the Up vector, pitch is rotation about the
  // Right axis, and roll is rotation about the Direction axis.
  bool ToEulerAnglesXYZ (double& rfYAngle, double& rfPAngle, double& rfRAngle) const;
  bool ToEulerAnglesXZY (double& rfYAngle, double& rfPAngle, double& rfRAngle) const;
  bool ToEulerAnglesYXZ (double& rfYAngle, double& rfPAngle, double& rfRAngle) const;
  bool ToEulerAnglesYZX (double& rfYAngle, double& rfPAngle, double& rfRAngle) const;
  bool ToEulerAnglesZXY (double& rfYAngle, double& rfPAngle, double& rfRAngle) const;
  bool ToEulerAnglesZYX (double& rfYAngle, double& rfPAngle, double& rfRAngle) const;
  void FromEulerAnglesXYZ (double fYAngle, double fPAngle, double fRAngle);
  void FromEulerAnglesXZY (double fYAngle, double fPAngle, double fRAngle);
  void FromEulerAnglesYXZ (double fYAngle, double fPAngle, double fRAngle);
  void FromEulerAnglesYZX (double fYAngle, double fPAngle, double fRAngle);
  void FromEulerAnglesZXY (double fYAngle, double fPAngle, double fRAngle);
  void FromEulerAnglesZYX (double fYAngle, double fPAngle, double fRAngle);

  // eigensolver, matrix must be symmetric
  void EigenSolveSymmetric (double afEigenvalue[3],
                            Vector3 akEigenvector[3]) const;

  static void TensorProduct (const Vector3& rkU, const Vector3& rkV,
                             Matrix3x3& rkProduct);

  MIRTK_Numerics_EXPORT static const double EPSILON;
  MIRTK_Numerics_EXPORT static const Matrix3x3 ZERO;
  MIRTK_Numerics_EXPORT static const Matrix3x3 IDENTITY;

  // support for eigensolver
  void Tridiagonal(double afDiag[3], double afSubDiag[3]);
  bool QLAlgorithm(double afDiag[3], double afSubDiag[3]);

  // support for singular value decomposition
  MIRTK_Numerics_EXPORT static const double ms_fSvdEpsilon;
  MIRTK_Numerics_EXPORT static const int ms_iSvdMaxIterations;
  static void Bidiagonalize (Matrix3x3& kA, Matrix3x3& kL,
                             Matrix3x3& kR);
  static void GolubKahanStep (Matrix3x3& kA, Matrix3x3& kL,
                              Matrix3x3& kR);

  // support for spectral norm
  static double MaxCubicRoot (double afCoeff[3]);

protected:

  double m_aafEntry[3][3];
};


} // namespace mirtk

#endif // MIRTK_Matrix3x3_H
