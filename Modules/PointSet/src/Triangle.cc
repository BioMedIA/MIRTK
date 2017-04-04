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

  return min;
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
bool Triangle
::TriangleTriangleIntersection(const double a1[3], const double b1[3], const double c1[3],
                               const double a2[3], const double b2[3], const double c2[3],
                               int &coplanar, double *p1, double *p2)
{
  #if VTK_MAJOR_VERSION < 7 || (VTK_MAJOR_VERSION == 7 && VTK_MINOR_VERSION == 0)
    return vtkIntersectionPolyDataFilter
        ::TriangleTriangleIntersection(const_cast<double *>(a1),
                                       const_cast<double *>(b1),
                                       const_cast<double *>(c1),
                                       const_cast<double *>(a2),
                                       const_cast<double *>(b2),
                                       const_cast<double *>(c2),
                                       coplanar, p1, p2) != 0;
  #else
    // Tolerance value must not be close to zero
    // See https://gitlab.kitware.com/vtk/vtk/issues/17012
    double surfaceid[2], tol = 1e-6;
    return vtkIntersectionPolyDataFilter
        ::TriangleTriangleIntersection(const_cast<double *>(a1),
                                       const_cast<double *>(b1),
                                       const_cast<double *>(c1),
                                       const_cast<double *>(a2),
                                       const_cast<double *>(b2),
                                       const_cast<double *>(c2),
                                       coplanar, p1, p2, surfaceid, tol) != 0;
  #endif
}


} // namespace mirtk
