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

#include "mirtk/Triangle.h"
#include "mirtk/Memory.h"

#include "mirtk/VtkMath.h"
#include "vtkVersion.h"
#include "vtkLine.h"
#include "vtkPlane.h"
#include "vtkTriangle.h"
#include "vtkIntersectionPolyDataFilter.h"


namespace mirtk {


// =============================================================================
// Distance between two triangles
// =============================================================================

// Copy 3D point coordinates stored in plain C array
#define MIRTK_RETURN_POINT(a, b) if (a) memcpy((a), (b), 3 * sizeof(double))

// -----------------------------------------------------------------------------
// Scaled vector subtraction
inline void VmsV(double *Vr, const double *V1, double s2, const double *V2)
{
  Vr[0] = V1[0] - s2 * V2[0];
  Vr[1] = V1[1] - s2 * V2[1];
  Vr[2] = V1[2] - s2 * V2[2];
}

// -----------------------------------------------------------------------------
inline double distance_between_triangles(double a1[3], double b1[3], double c1[3], double n1[3],
                                         double a2[3], double b2[3], double c2[3], double n2[3],
                                         double *p1 = nullptr, double *p2 = nullptr)
{
  const double TOL = 1e-6; // Tolerance for point coordinate equality check
  double q1[3], q2[3], t1, t2, d, min = inf;

  // Point in second triangle closest to any of the corners of the first triangle
  t2 = vtkPlane::Evaluate(n2, a2, a1);
  d  = abs(t2);
  if (d < min) {
    VmsV(q2, a1, t2, n2);
    if (vtkTriangle::PointInTriangle(q2, a2, b2, c2, TOL)) {
      MIRTK_RETURN_POINT(p1, a1);
      MIRTK_RETURN_POINT(p2, q2);
      min = d;
    }
  }
  t2 = vtkPlane::Evaluate(n2, a2, b1);
  d  = abs(t2);
  if (d < min) {
    VmsV(q2, b1, t2, n2);
    if (vtkTriangle::PointInTriangle(q2, a2, b2, c2, TOL)) {
      MIRTK_RETURN_POINT(p1, b1);
      MIRTK_RETURN_POINT(p2, q2);
      min = d;
    }
  }
  t2 = vtkPlane::Evaluate(n2, a2, c1);
  d  = abs(t2);
  if (d < min) {
    VmsV(q2, c1, t2, n2);
    if (vtkTriangle::PointInTriangle(q2, a2, b2, c2, TOL)) {
      MIRTK_RETURN_POINT(p1, c1);
      MIRTK_RETURN_POINT(p2, q2);
      min = d;
    }
  }

  // Point in first triangle closest to any of the corners of the second triangle
  t1 = vtkPlane::Evaluate(n1, a1, a2);
  d  = abs(t1);
  if (d < min) {
    VmsV(q1, a2, t1, n1);
    if (vtkTriangle::PointInTriangle(q1, a1, b1, c1, TOL)) {
      MIRTK_RETURN_POINT(p1, q1);
      MIRTK_RETURN_POINT(p2, a2);
      min = d;
    }
  }
  t1 = vtkPlane::Evaluate(n1, a1, b2);
  d  = abs(t1);
  if (d < min) {
    VmsV(q1, b2, t1, n1);
    if (vtkTriangle::PointInTriangle(q1, a1, b1, c1, TOL)) {
      MIRTK_RETURN_POINT(p1, q1);
      MIRTK_RETURN_POINT(p2, b2);
      min = d;
    }
  }
  t1 = vtkPlane::Evaluate(n1, a1, c2);
  d  = abs(t1);
  if (d < min) {
    VmsV(q1, c2, t1, n1);
    if (vtkTriangle::PointInTriangle(q1, a1, b1, c1, TOL)) {
      MIRTK_RETURN_POINT(p1, q1);
      MIRTK_RETURN_POINT(p2, c2);
      min = d;
    }
  }

  // Attention: vtkLine::DistanceBetweenLineSegments returns squared distance!
  // Take square root again before returning the final minimum distance value.
  min *= min;

  // Distance between each pair of edges
  d = vtkLine::DistanceBetweenLineSegments(a1, b1, a2, b2, q1, q2, t1, t2);
  if (d < min) {
    MIRTK_RETURN_POINT(p1, q1);
    MIRTK_RETURN_POINT(p2, q2);
    min = d;
  }
  d = vtkLine::DistanceBetweenLineSegments(a1, b1, b2, c2, q1, q2, t1, t2);
  if (d < min) {
    MIRTK_RETURN_POINT(p1, q1);
    MIRTK_RETURN_POINT(p2, q2);
    min = d;
  }
  d = vtkLine::DistanceBetweenLineSegments(a1, b1, c2, a2, q1, q2, t1, t2);
  if (d < min) {
    MIRTK_RETURN_POINT(p1, q1);
    MIRTK_RETURN_POINT(p2, q2);
    min = d;
  }
  d = vtkLine::DistanceBetweenLineSegments(b1, c1, a2, b2, q1, q2, t1, t2);
  if (d < min) {
    MIRTK_RETURN_POINT(p1, q1);
    MIRTK_RETURN_POINT(p2, q2);
    min = d;
  }
  d = vtkLine::DistanceBetweenLineSegments(b1, c1, b2, c2, q1, q2, t1, t2);
  if (d < min) {
    MIRTK_RETURN_POINT(p1, q1);
    MIRTK_RETURN_POINT(p2, q2);
    min = d;
  }
  d = vtkLine::DistanceBetweenLineSegments(b1, c1, c2, a2, q1, q2, t1, t2);
  if (d < min) {
    MIRTK_RETURN_POINT(p1, q1);
    MIRTK_RETURN_POINT(p2, q2);
    min = d;
  }
  d = vtkLine::DistanceBetweenLineSegments(c1, a1, a2, b2, q1, q2, t1, t2);
  if (d < min) {
    MIRTK_RETURN_POINT(p1, q1);
    MIRTK_RETURN_POINT(p2, q2);
    min = d;
  }
  d = vtkLine::DistanceBetweenLineSegments(c1, a1, b2, c2, q1, q2, t1, t2);
  if (d < min) {
    MIRTK_RETURN_POINT(p1, q1);
    MIRTK_RETURN_POINT(p2, q2);
    min = d;
  }
  d = vtkLine::DistanceBetweenLineSegments(c1, a1, c2, a2, q1, q2, t1, t2);
  if (d < min) {
    MIRTK_RETURN_POINT(p1, q1);
    MIRTK_RETURN_POINT(p2, q2);
    min = d;
  }

  return sqrt(min);
}

// -----------------------------------------------------------------------------
double Triangle
::DistanceBetweenTriangles(
    const double a1[3], const double b1[3], const double c1[3], const double n1[3],
    const double a2[3], const double b2[3], const double c2[3], const double n2[3],
    double *p1, double *p2
) {
  return distance_between_triangles(const_cast<double *>(a1),
                                    const_cast<double *>(b1),
                                    const_cast<double *>(c1),
                                    const_cast<double *>(n1),
                                    const_cast<double *>(a2),
                                    const_cast<double *>(b2),
                                    const_cast<double *>(c2),
                                    const_cast<double *>(n2), p1, p2);
}

// =============================================================================
// Triangle/triangle intersection test
// =============================================================================

#include "triangle_triangle_intersection.h"

// -----------------------------------------------------------------------------
bool Triangle
::TriangleTriangleOverlap(const double a1[2], const double b1[2], const double c1[2],
                          const double a2[2], const double b2[2], const double c2[2])
{
  return (0 != tri_tri_overlap_test_2d(const_cast<double *>(a1),
                                       const_cast<double *>(b1),
                                       const_cast<double *>(c1),
                                       const_cast<double *>(a2),
                                       const_cast<double *>(b2),
                                       const_cast<double *>(c2)));
}

// -----------------------------------------------------------------------------
bool Triangle
::TriangleTriangleIntersection(const double a1[3], const double b1[3], const double c1[3],
                               const double a2[3], const double b2[3], const double c2[3])
{
  return (0 != tri_tri_overlap_test_3d(const_cast<double *>(a1),
                                       const_cast<double *>(b1),
                                       const_cast<double *>(c1),
                                       const_cast<double *>(a2),
                                       const_cast<double *>(b2),
                                       const_cast<double *>(c2)));
}

// -----------------------------------------------------------------------------
// Copied from VTK 9.0.1 to fix segfault causing issue
//
// See: https://gitlab.kitware.com/vtk/vtk/-/issues/17722
bool Triangle
::TriangleTriangleIntersection(const double a1[3], const double b1[3], const double c1[3],
                               const double a2[3], const double b2[3], const double c2[3],
                               int &coplanar, double *pt1, double *pt2)
{
  double *p1 = const_cast<double *>(a1);
  double *q1 = const_cast<double *>(b1);
  double *r1 = const_cast<double *>(c1);

  double *p2 = const_cast<double *>(a2);
  double *q2 = const_cast<double *>(b2);
  double *r2 = const_cast<double *>(c2);

  double surfaceid[2], tolerance = 0;
  double n1[3], n2[3];

  // Compute supporting plane normals.
  vtkTriangle::ComputeNormal(p1, q1, r1, n1);
  vtkTriangle::ComputeNormal(p2, q2, r2, n2);
  double s1 = -vtkMath::Dot(n1, p1);
  double s2 = -vtkMath::Dot(n2, p2);

  // Compute signed distances of points p1, q1, r1 from supporting
  // plane of second triangle.
  double dist1[3];
  dist1[0] = vtkMath::Dot(n2, p1) + s2;
  dist1[1] = vtkMath::Dot(n2, q1) + s2;
  dist1[2] = vtkMath::Dot(n2, r1) + s2;

  // If signs of all points are the same, all the points lie on the
  // same side of the supporting plane, and we can exit early.
  if ((dist1[0] * dist1[1] > tolerance) && (dist1[0] * dist1[2] > tolerance))
  {
    // vtkDebugMacro(<<"Same side supporting plane 1!");
    return 0;
  }
  // Do the same for p2, q2, r2 and supporting plane of first
  // triangle.
  double dist2[3];
  dist2[0] = vtkMath::Dot(n1, p2) + s1;
  dist2[1] = vtkMath::Dot(n1, q2) + s1;
  dist2[2] = vtkMath::Dot(n1, r2) + s1;

  // If signs of all points are the same, all the points lie on the
  // same side of the supporting plane, and we can exit early.
  if ((dist2[0] * dist2[1] > tolerance) && (dist2[0] * dist2[2] > tolerance))
  {
    // vtkDebugMacro(<<"Same side supporting plane 2!");
    return 0;
  }
  // Check for coplanarity of the supporting planes.
  if (fabs(n1[0] - n2[0]) < 1e-9 && fabs(n1[1] - n2[1]) < 1e-9 && fabs(n1[2] - n2[2]) < 1e-9 &&
    fabs(s1 - s2) < 1e-9)
  {
    coplanar = 1;
    // vtkDebugMacro(<<"Coplanar!");
    return 0;
  }

  coplanar = 0;

  // There are more efficient ways to find the intersection line (if
  // it exists), but this is clear enough.
  double *pts1[3] = { p1, q1, r1 };
  double *pts2[3] = { p2, q2, r2 };

  // Find line of intersection (L = p + t*v) between two planes.
  double n1n2 = vtkMath::Dot(n1, n2);
  double a = (s1 - s2 * n1n2) / (n1n2 * n1n2 - 1.0);
  double b = (s2 - s1 * n1n2) / (n1n2 * n1n2 - 1.0);
  double p[3], v[3];
  p[0] = a * n1[0] + b * n2[0];
  p[1] = a * n1[1] + b * n2[1];
  p[2] = a * n1[2] + b * n2[2];
  vtkMath::Cross(n1, n2, v);
  vtkMath::Normalize(v);

  int index1 = 0, index2 = 0;
  double t1[3], t2[3];
  int ts1 = 50, ts2 = 50;
  for (int i = 0; i < 3; i++)
  {
    double t, x[3];
    int id1 = i, id2 = (i + 1) % 3;

    // Find t coordinate on line of intersection between two planes.
    double val1 = vtkPlane::IntersectWithLine(pts1[id1], pts1[id2], n2, p2, t, x);
    if (val1 == 1 || (t > (0 - tolerance) && t < (1 + tolerance)))
    {
      if (t < 1 + tolerance && t > 1 - tolerance)
      {
        ts1 = index1;
      }

      t1[index1++] = vtkMath::Dot(x, v) - vtkMath::Dot(p, v);
    }

    double val2 = vtkPlane::IntersectWithLine(pts2[id1], pts2[id2], n1, p1, t, x);
    if (val2 == 1 || (t > (0 - tolerance) && t < (1 + tolerance)))
    {
      if (t < 1 + tolerance && t > 1 - tolerance)
      {
        ts2 = index2;
      }

      t2[index2++] = vtkMath::Dot(x, v) - vtkMath::Dot(p, v);
    }
  }

  // If the value of the index is greater than 2, the intersecting point
  // actually is intersected by all three edges. In this case, set the two
  // edges to the two edges where the intersecting point is not the end point
  if (index1 > 2 && ts1 < 50)
  {
    index1--;
    std::swap(t1[ts1], t1[2]);
  }
  if (index2 > 2 && ts2 < 50)
  {
    index2--;
    std::swap(t2[ts2], t2[2]);
  }
  // Check if only one edge or all edges intersect the supporting
  // planes intersection.
  if (index1 != 2 || index2 != 2)
  {
    // vtkDebugMacro(<<"Only one edge intersecting!");
    return 0;
  }

  // Check for NaNs
  if (vtkMath::IsNan(t1[0]) || vtkMath::IsNan(t1[1]) || vtkMath::IsNan(t2[0]) || vtkMath::IsNan(t2[1]))
  {
    // vtkWarningMacro(<<"NaNs!");
    return 0;
  }

  if (t1[0] > t1[1])
  {
    std::swap(t1[0], t1[1]);
  }
  if (t2[0] > t2[1])
  {
    std::swap(t2[0], t2[1]);
  }
  // Handle the different interval configuration cases.
  double tt1, tt2;
  if (t1[1] < t2[0] || t2[1] < t1[0])
  {
    // vtkDebugMacro(<<"No Overlap!");
    return 0; // No overlap
  }
  else if (t1[0] < t2[0])
  {
    if (t1[1] < t2[1])
    {
      // First point on surface 2, second point on surface 1
      surfaceid[0] = 2;
      surfaceid[1] = 1;
      tt1 = t2[0];
      tt2 = t1[1];
    }
    else
    {
      // Both points belong to lines on surface 2
      surfaceid[0] = 2;
      surfaceid[1] = 2;
      tt1 = t2[0];
      tt2 = t2[1];
    }
  }
  else // t1[0] >= t2[0]
  {
    if (t1[1] < t2[1])
    {
      // Both points belong to lines on surface 1
      surfaceid[0] = 1;
      surfaceid[1] = 1;
      tt1 = t1[0];
      tt2 = t1[1];
    }
    else
    {
      // First point on surface 1, second point on surface 2
      surfaceid[0] = 1;
      surfaceid[1] = 2;
      tt1 = t1[0];
      tt2 = t2[1];
    }
  }

  // Create actual intersection points.
  pt1[0] = p[0] + tt1 * v[0];
  pt1[1] = p[1] + tt1 * v[1];
  pt1[2] = p[2] + tt1 * v[2];

  pt2[0] = p[0] + tt2 * v[0];
  pt2[1] = p[1] + tt2 * v[1];
  pt2[2] = p[2] + tt2 * v[2];

  return 1;
}


} // namespace mirtk
