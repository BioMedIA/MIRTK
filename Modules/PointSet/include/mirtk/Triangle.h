/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2016 Imperial College London
 * Copyright 2013-2016 Andreas Schuh
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

#ifndef MIRTK_Triangle_H
#define MIRTK_Triangle_H

#include "mirtk/Math.h"
#include "mirtk/Vector3.h"
#include "mirtk/Matrix3x3.h"
#include "mirtk/VtkMath.h"


namespace mirtk {


/**
 * Auxiliary class/static helper functions for triangulated surface meshes
 */
class Triangle
{
public:

  /// Compute center point of triangle
  ///
  /// \param[in]  a      Position of triangle vertex A.
  /// \param[in]  b      Position of triangle vertex B.
  /// \param[in]  c      Position of triangle vertex C.
  /// \param[out] center Position of triangle center.
  static void Center(const double a[3], const double b[3], const double c[3], double center[3]);

  /// Compute center point of triangle
  ///
  /// \param[in]  a      Position of triangle vertex A.
  /// \param[in]  b      Position of triangle vertex B.
  /// \param[in]  c      Position of triangle vertex C.
  /// \param[out] center Position of triangle center.
  ///
  /// \returns Radius of bounding sphere with optionally returned \p center point.
  static double BoundingSphereRadius(const double a[3], const double b[3], const double c[3], double *center = nullptr);

  /// Compute normal direction of triangle
  ///
  /// \param[in]  a Position of triangle vertex A.
  /// \param[in]  b Position of triangle vertex B.
  /// \param[in]  c Position of triangle vertex C.
  /// \param[out] n Non-normalized triangle normal vector.
  ///               The length of the resulting vector equals twice the triangle area.
  static void NormalDirection(const double a[3], const double b[3], const double c[3], double n[3]);

  /// Partial derivatives of normal direction w.r.t. position of vertex A
  ///
  /// \param[in] a Position of triangle vertex A.
  /// \param[in] b Position of triangle vertex B.
  /// \param[in] c Position of triangle vertex C.
  ///
  /// \returns Partial derivatives of normal direction components w.r.t. coordinates of vertex A.
  static Matrix3x3 NormalDirectionJacobian(const double a[3], const double b[3], const double c[3]);

  /// Compute normal of triangle
  ///
  /// \param[in]  a Position of triangle vertex A.
  /// \param[in]  b Position of triangle vertex B.
  /// \param[in]  c Position of triangle vertex C.
  /// \param[out] n Triangle normal vector.
  static void Normal(const double a[3], const double b[3], const double c[3], double n[3]);

  /// Partial derivatives of normal w.r.t. coordinates of vertex A
  ///
  /// \param[in]  a Position of triangle vertex A.
  /// \param[in]  b Position of triangle vertex B.
  /// \param[in]  c Position of triangle vertex C.
  ///
  /// \returns Partial derivatives of normal components w.r.t. coordinates of vertex A.
  static Matrix3x3 NormalJacobian(const double a[3], const double b[3], const double c[3]);

  /// Partial derivatives of normal w.r.t. coordinates of vertex A
  ///
  /// \param[in]  a Position of triangle vertex A.
  /// \param[in]  b Position of triangle vertex B.
  /// \param[in]  c Position of triangle vertex C.
  /// \param[in] dn Partial derivatives of normal direction vector.
  ///
  /// \returns Partial derivatives of normal components w.r.t. coordinates of vertex A.
  static Matrix3x3 NormalJacobian(const double a[3], const double b[3], const double c[3], const Matrix3x3 &dn);

  /// Partial derivatives of normal w.r.t. coordinates of vertex A
  ///
  /// \param[in]  n Triangle normal vector, i.e., normalized direction vector.
  /// \param[in] dn Partial derivatives of normal direction vector.
  ///
  /// \returns Partial derivatives of normal components w.r.t. coordinates of vertex A.
  static Matrix3x3 NormalJacobian(const double n[3], const Matrix3x3 &dn);

  /// Partial derivatives of normal w.r.t. coordinates of vertex A
  ///
  /// \param[in]  n Triangle normal vector, i.e., normalized direction vector.
  /// \param[in] dn Partial derivatives of normal direction vector.
  ///
  /// \returns Partial derivatives of normal components w.r.t. coordinates of vertex A.
  static Matrix3x3 NormalJacobian(const Vector3 &n, const Matrix3x3 &dn);

  /// Compute twice the area of triangle
  ///
  /// \param[in] a Position of triangle vertex A.
  /// \param[in] b Position of triangle vertex B.
  /// \param[in] c Position of triangle vertex C.
  ///
  /// \returns Twice the area of the triangle.
  static double DoubleArea(const double a[3], const double b[3], const double c[3]);

  /// Partial derivatives of twice the triangle area w.r.t. coordinates of vertex A
  ///
  /// \param[in] a Position of triangle vertex A.
  /// \param[in] b Position of triangle vertex B.
  /// \param[in] c Position of triangle vertex C.
  ///
  /// \returns Partial derivatives of twice the triangle area w.r.t. coordinates of vertex A.
  static Vector3 DoubleAreaGradient(const double a[3], const double b[3], const double c[3]);

  /// Compute area of triangle
  ///
  /// \param[in] a Position of triangle vertex A.
  /// \param[in] b Position of triangle vertex B.
  /// \param[in] c Position of triangle vertex C.
  ///
  /// \returns Area of triangle.
  static double Area(const double a[3], const double b[3], const double c[3]);

  /// Partial derivative of triangle area w.r.t. coordinates of vertex A
  ///
  /// \param[in] a Position of triangle vertex A.
  /// \param[in] b Position of triangle vertex B.
  /// \param[in] c Position of triangle vertex C.
  ///
  /// \returns Partial derivatives of triangle area w.r.t. coordinates of vertex A.
  static Vector3 AreaGradient(const double a[3], const double b[3], const double c[3]);

  /// Compute area of triangle in 2D
  ///
  /// \param[in] a Position of triangle vertex A.
  /// \param[in] b Position of triangle vertex B.
  /// \param[in] c Position of triangle vertex C.
  ///
  /// \returns Signed area of triangle.
  static double SignedArea2D(const double a[2], const double b[2], const double c[2]);

  /// Compute twice the area of triangle in 2D
  ///
  /// \param[in] a Position of triangle vertex A.
  /// \param[in] b Position of triangle vertex B.
  /// \param[in] c Position of triangle vertex C.
  ///
  /// \returns Twice the signed area of the triangle.
  static double DoubleSignedArea2D(const double a[2], const double b[2], const double c[2]);

  /// Compute area of triangle in 2D
  ///
  /// \param[in] a Position of triangle vertex A.
  /// \param[in] b Position of triangle vertex B.
  /// \param[in] c Position of triangle vertex C.
  ///
  /// \returns Area of triangle.
  static double Area2D(const double a[2], const double b[2], const double c[2]);

  /// Compute twice the area of triangle in 2D
  ///
  /// \param[in] a Position of triangle vertex A.
  /// \param[in] b Position of triangle vertex B.
  /// \param[in] c Position of triangle vertex C.
  ///
  /// \returns Twice the area of the triangle.
  static double DoubleArea2D(const double a[2], const double b[2], const double c[2]);

  /// Compute cotangent of angle ABC
  ///
  /// \param[in] a Position of triangle vertex A.
  /// \param[in] b Position of triangle vertex B.
  /// \param[in] c Position of triangle vertex C.
  ///
  /// \returns Cotangent of angle ABC (equals cotangent of angle CBA).
  static double Cotangent(double a[3], double b[3], double c[3]);

  /// Compute angle at triangle corners
  ///
  /// \param[in]  a     Position of triangle vertex A.
  /// \param[in]  b     Position of triangle vertex B.
  /// \param[in]  c     Position of triangle vertex C.
  /// \param[out] angle Angles at triangle vertices in radians.
  static void Angles(const double a[3], const double b[3], const double c[3], double angle[3]);

  /// Compute minimum angle at triangle vertices
  ///
  /// \return Minimum angle in radians.
  static double MinAngle(const double a[3], const double b[3], const double c[3]);

  /// Compute maximum angle at triangle corners
  ///
  /// \return Maximum angle in radians.
  static double MaxAngle(const double a[3], const double b[3], const double c[3]);

  /// Tests whether two triangles intersect each other
  static bool TriangleTriangleIntersection(const double a1[3], const double b1[3], const double c1[3],
                                           const double a2[3], const double b2[3], const double c2[3]);

  /// Tests whether two triangles intersect each other
  static bool TriangleTriangleIntersection(const double a1[3], const double b1[3], const double c1[3],
                                           const double a2[3], const double b2[3], const double c2[3],
                                           int &coplanar, double *p1, double *p2);

  /// Compute distance between triangle center points
  static double DistanceBetweenCenters(const double a1[3], const double b1[3], const double c1[3],
                                       const double a2[3], const double b2[3], const double c2[3],
                                       double *p1, double *p2);

  /// Compute distance between closest corner points of two triangles
  static double DistanceBetweenCorners(const double a1[3], const double b1[3], const double c1[3],
                                       const double a2[3], const double b2[3], const double c2[3],
                                       double *p1, double *p2);

  /// Compute distance between closest points of two triangles
  static double DistanceBetweenTriangles(const double a1[3], const double b1[3], const double c1[3], const double n1[3],
                                         const double a2[3], const double b2[3], const double c2[3], const double n2[3],
                                         double *p1, double *p2);

  /// Compute distance between closest points of two triangles
  static double DistanceBetweenTriangles(const double a1[3], const double b1[3], const double c1[3],
                                         const double a2[3], const double b2[3], const double c2[3],
                                         double *p1, double *p2);

};


////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Center
// =============================================================================

//----------------------------------------------------------------------------
inline void Triangle::Center(const double a[3], const double b[3], const double c[3], double center[3])
{
  center[0] = (a[0] + b[0] + c[0]) / 3.0;
  center[1] = (a[1] + b[1] + c[1]) / 3.0;
  center[2] = (a[2] + b[2] + c[2]) / 3.0;
}

//----------------------------------------------------------------------------
inline double Triangle::BoundingSphereRadius(const double a[3], const double b[3], const double c[3], double *center)
{
  double p[3];
  Center(a, b, c, p);
  if (center) memcpy(center, p, 3 * sizeof(double));
  return sqrt(max(max(vtkMath::Distance2BetweenPoints(a, p),
                      vtkMath::Distance2BetweenPoints(b, p)),
                      vtkMath::Distance2BetweenPoints(c, p)));
}

// =============================================================================
// Normals
// =============================================================================

// -----------------------------------------------------------------------------
inline void Triangle::NormalDirection(const double a[3], const double b[3], const double c[3], double n[3])
{
  double ab[3], ac[3];
  vtkMath::Subtract(b, a, ab);
  vtkMath::Subtract(c, a, ac);
  vtkMath::Cross(ab, ac, n);
}

// -----------------------------------------------------------------------------
inline void Triangle::Normal(const double a[3], const double b[3], const double c[3], double n[3])
{
  double l;
  NormalDirection(a, b, c, n);
  if ((l = vtkMath::Norm(n)) != 0.) {
    n[0] /= l, n[1] /= l, n[2] /= l;
  }
}

// =============================================================================
// Area
// =============================================================================

// -----------------------------------------------------------------------------
inline double Triangle::DoubleArea(const double a[3], const double b[3], const double c[3])
{
  double n[3];
  NormalDirection(a, b, c, n);
  return vtkMath::Norm(n);
}

// -----------------------------------------------------------------------------
inline double Triangle::Area(const double a[3], const double b[3], const double c[3])
{
  return .5 * DoubleArea(a, b, c);
}

// -----------------------------------------------------------------------------
inline double Triangle::DoubleSignedArea2D(const double a[2], const double b[2], const double c[2])
{
  return (a[0]*b[1] - a[1]*b[0]) + (b[0]*c[1] - b[1]*c[0]) + (c[0]*a[1] - c[1]*a[0]);
}

// -----------------------------------------------------------------------------
inline double Triangle::SignedArea2D(const double a[2], const double b[2], const double c[2])
{
  return .5 * DoubleSignedArea2D(a, b, c);
}

// -----------------------------------------------------------------------------
inline double Triangle::DoubleArea2D(const double a[2], const double b[2], const double c[2])
{
  return abs(DoubleSignedArea2D(a, b, c));
}

// -----------------------------------------------------------------------------
inline double Triangle::Area2D(const double a[2], const double b[2], const double c[2])
{
  return abs(SignedArea2D(a, b, c));
}

// =============================================================================
// Angles
// =============================================================================

// -----------------------------------------------------------------------------
inline void Triangle::Angles(const double a[3], const double b[3], const double c[3], double angle[3])
{
  const double ab2 = vtkMath::Distance2BetweenPoints(a, b), ab = sqrt(ab2);
  const double ac2 = vtkMath::Distance2BetweenPoints(a, c), ac = sqrt(ac2);
  const double bc2 = vtkMath::Distance2BetweenPoints(b, c), bc = sqrt(bc2);
  angle[0] = acos((ab2 + ac2 - bc2) / (2.0 * ab * ac));
  angle[1] = acos((ab2 + bc2 - ac2) / (2.0 * ab * bc));
  angle[2] = acos((ac2 + bc2 - ab2) / (2.0 * ac * bc));
}

// -----------------------------------------------------------------------------
inline double Triangle::MinAngle(const double a[3], const double b[3], const double c[3])
{
  double angle[3];
  Angles(a, b, c, angle);
  return min(min(angle[0], angle[1]), angle[2]);
}

// -----------------------------------------------------------------------------
inline double Triangle::MaxAngle(const double a[3], const double b[3], const double c[3])
{
  double angle[3];
  Angles(a, b, c, angle);
  return max(max(angle[0], angle[1]), angle[2]);
}

// -----------------------------------------------------------------------------
// Meyer et al. (2002). Generalized Barycentric Coordinates on Irregular Polygons.
inline double Triangle::Cotangent(double a[3], double b[3], double c[3])
{
  double ba[3], bc[3], n[3];
  vtkMath::Subtract(a, b, ba);
  vtkMath::Subtract(c, b, bc);
  vtkMath::Cross(ba, bc, n);
  return vtkMath::Dot(ba, bc) / vtkMath::Norm(n);
}

// =============================================================================
// Derivatives
// =============================================================================

// -----------------------------------------------------------------------------
inline Matrix3x3 Triangle
::NormalDirectionJacobian(const double a[3], const double b[3], const double c[3])
{
  // - Derivative of ab x ac w.r.t. ab:
  //   [     0.,     b[2] - a[2], a[1] - b[1];
  //    a[2] - b[2],      0.,     b[0] - a[0];
  //    b[1] - a[1], a[0] - b[0],      0.    ]
  // - Derivative of ab x ac w.r.t. ac:
  //   [         0., a[2] - c[2], c[1] - a[1];
  //    c[2] - a[2],          0., a[0] - c[0];
  //    a[1] - c[1], c[0] - a[0],          0.]
  // - Dervative of ab w.r.t. a: -I
  // - Dervative of ac w.r.t. a: -I
  //
  // Apply chain rule to obtain total derivative of ab x ac w.r.t. a.
  // This results in the following three subtractions and negations.
  const double m01 = c[2] - b[2];
  const double m02 = b[1] - c[1];
  const double m12 = c[0] - b[0];
  return Matrix3x3(0.,  m01, m02, -m01, 0., m12, -m02, -m12, 0.);
}

// -----------------------------------------------------------------------------
inline Matrix3x3 Triangle
::NormalJacobian(const Vector3 &n, const Matrix3x3 &dn)
{
  const double l2 = n.SquaredLength();
  if (l2 == 0.) return Matrix3x3(0.);

  const double l  = sqrt(l2);
  const double l3 = l2 * l;
  const double d  = 1. / l;

  const double m00 = n[0] * n[0] / l3;
  const double m01 = n[0] * n[1] / l3;
  const double m02 = n[0] * n[2] / l3;
  const double m11 = n[1] * n[1] / l3;
  const double m12 = n[1] * n[2] / l3;
  const double m22 = n[2] * n[2] / l3;

  return Matrix3x3(m00 + d, m01, m02,
                   m01, m11 + d, m12,
                   m02, m12, m22 + d) * dn;
}

// -----------------------------------------------------------------------------
inline Matrix3x3 Triangle
::NormalJacobian(const double n[3], const Matrix3x3 &dn)
{
  const double l2 = vtkMath::Dot(n, n);
  if (l2 == 0.) return Matrix3x3(0.);

  const double l  = sqrt(l2);
  const double l3 = l2 * l;
  const double d  = 1. / l;

  const double m00 = n[0] * n[0] / l3;
  const double m01 = n[0] * n[1] / l3;
  const double m02 = n[0] * n[2] / l3;
  const double m11 = n[1] * n[1] / l3;
  const double m12 = n[1] * n[2] / l3;
  const double m22 = n[2] * n[2] / l3;

  return Matrix3x3(m00 + d, m01, m02,
                   m01, m11 + d, m12,
                   m02, m12, m22 + d) * dn;
}

// -----------------------------------------------------------------------------
inline Matrix3x3 Triangle
::NormalJacobian(const double a[3], const double b[3], const double c[3], const Matrix3x3 &dn)
{
  double n[3];
  Normal(a, b, c, n);
  return NormalJacobian(n, dn);
}

// -----------------------------------------------------------------------------
inline Matrix3x3 Triangle
::NormalJacobian(const double a[3], const double b[3], const double c[3])
{
  return NormalJacobian(a, b, c, NormalDirectionJacobian(a, b, c));
}

// -----------------------------------------------------------------------------
inline Vector3 Triangle
::DoubleAreaGradient(const double a[3], const double b[3], const double c[3])
{
  double n[3];
  Normal(a, b, c, n);
  return - Vector3(n) * NormalDirectionJacobian(a, b, c);
}

// -----------------------------------------------------------------------------
inline Vector3 Triangle
::AreaGradient(const double a[3], const double b[3], const double c[3])
{
  return .5 * DoubleAreaGradient(a, b, c);
}

// =============================================================================
// Distance between two triangles
// =============================================================================

// -----------------------------------------------------------------------------
/// Compute distance between two triangles
inline double Triangle
::DistanceBetweenCenters(const double a1[3], const double b1[3], const double c1[3],
                         const double a2[3], const double b2[3], const double c2[3],
                         double *p1, double *p2)
{
  Center(a1, b1, c1, p1);
  Center(a2, b2, c2, p2);
  return sqrt(vtkMath::Distance2BetweenPoints(p1, p2));
}

// -----------------------------------------------------------------------------
/// Compute distance between two triangles
inline double Triangle
::DistanceBetweenCorners(const double a1[3], const double b1[3], const double c1[3],
                         const double a2[3], const double b2[3], const double c2[3],
                         double *p1, double *p2)
{
  double d, dmin = vtkMath::Distance2BetweenPoints(a1, a2);
  const double *x1 = a1, *x2 = a2;

  d = vtkMath::Distance2BetweenPoints(a1, b2);
  if (d < dmin) dmin = d, x1 = a1, x2 = b2;

  d = vtkMath::Distance2BetweenPoints(a1, c2);
  if (d < dmin) dmin = d, x1 = a1, x2 = c2;

  d = vtkMath::Distance2BetweenPoints(b1, a2);
  if (d < dmin) dmin = d, x1 = b1, x2 = a2;

  d = vtkMath::Distance2BetweenPoints(b1, b2);
  if (d < dmin) dmin = d, x1 = b1, x2 = b2;

  d = vtkMath::Distance2BetweenPoints(b1, c2);
  if (d < dmin) dmin = d, x1 = b1, x2 = c2;

  d = vtkMath::Distance2BetweenPoints(c1, a2);
  if (d < dmin) dmin = d, x1 = c1, x2 = a2;

  d = vtkMath::Distance2BetweenPoints(c1, b2);
  if (d < dmin) dmin = d, x1 = c1, x2 = b2;

  d = vtkMath::Distance2BetweenPoints(c1, c2);
  if (d < dmin) dmin = d, x1 = c1, x2 = c2;

  if (p1) memcpy(p1, x1, 3 * sizeof(double));
  if (p2) memcpy(p2, x2, 3 * sizeof(double));
  return sqrt(dmin);
}

// -----------------------------------------------------------------------------
/// Compute distance between two triangles
inline double Triangle
::DistanceBetweenTriangles(const double a1[3], const double b1[3], const double c1[3],
                           const double a2[3], const double b2[3], const double c2[3],
                           double *p1, double *p2)
{
  double n1[3], n2[3];
  Normal(a1, b1, c1, n1);
  Normal(a2, b2, c2, n2);
  return DistanceBetweenTriangles(a1, b1, c1, n1, a2, b2, c2, n2, p1, p2);
}


} // namespace mirtk

#endif // MIRTK_Triangle_H
