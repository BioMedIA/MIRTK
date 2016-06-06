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
#include "mirtk/VtkMath.h"


namespace mirtk {


/**
 * Auxiliary class/static helper functions for triangulated surface meshes
 */
class Triangle
{
public:

  /// Compute normal direction of triangle
  ///
  /// \param[in]  a Position of triangle vertex A.
  /// \param[in]  b Position of triangle vertex B.
  /// \param[in]  c Position of triangle vertex C.
  /// \param[out] n Non-normalized triangle normal vector. Length equals twice the triangle area.
  static void NormalDirection(const double a[3], const double b[3], const double c[3], double n[3]);

  /// Compute normal of triangle
  ///
  /// \param[in]  a Position of triangle vertex A.
  /// \param[in]  b Position of triangle vertex B.
  /// \param[in]  c Position of triangle vertex C.
  /// \param[out] n Triangle normal vector.
  static void Normal(const double a[3], const double b[3], const double c[3], double n[3]);

  /// Compute twice the area of triangle
  ///
  /// \param[in] a Position of triangle vertex A.
  /// \param[in] b Position of triangle vertex B.
  /// \param[in] c Position of triangle vertex C.
  ///
  /// \returns Twice the area of the triangle.
  static double DoubleArea(const double a[3], const double b[3], const double c[3]);

  /// Compute area of triangle
  ///
  /// \param[in] a Position of triangle vertex A.
  /// \param[in] b Position of triangle vertex B.
  /// \param[in] c Position of triangle vertex C.
  ///
  /// \returns Area of triangle.
  static double Area(const double a[3], const double b[3], const double c[3]);

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

  /// Tests whether two triangles intersect each other
  static bool TriangleTriangleIntersection(const double a1[3], const double b1[3], const double c1[3],
                                           const double a2[3], const double b2[3], const double c2[3]);

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
// Distance between two triangles
// =============================================================================

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
