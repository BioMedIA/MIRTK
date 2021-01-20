/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2016 Imperial College London
 * Copyright 2013-2020 Andreas Schuh
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

#include "mirtk/SurfaceCollisions.h"

#include "mirtk/Assert.h"
#include "mirtk/Math.h"
#include "mirtk/Triangle.h"
#include "mirtk/PointSetUtils.h"
#include "mirtk/PointLocator.h"
#include "mirtk/Parallel.h"
#include "mirtk/Profiling.h"
#include "mirtk/Vtk.h"
#include "mirtk/VtkMath.h"

#include "vtkPlane.h"
#include "vtkTriangle.h"

#include "vtkPolyData.h"
#include "vtkCellData.h"
#include "vtkDataArray.h"
#include "vtkFloatArray.h"
#include "vtkUnsignedCharArray.h"

#include "vtkAbstractPointLocator.h"
#include "vtkOctreePointLocator.h"


namespace mirtk {


// =============================================================================
// Names of data arrays
// =============================================================================

MIRTK_PointSet_EXPORT const char * const SurfaceCollisions::BOUNDING_SPHERE_CENTER = "BoundingSphereCenter";
MIRTK_PointSet_EXPORT const char * const SurfaceCollisions::BOUNDING_SPHERE_RADIUS = "BoundingSphereRadius";
MIRTK_PointSet_EXPORT const char * const SurfaceCollisions::COLLISION_TYPE         = "CollisionType";

// =============================================================================
// Auxiliary functors
// =============================================================================

namespace SurfaceCollisionsUtils {


// -----------------------------------------------------------------------------
/// Compute center and radius of each bounding spheres
class ComputeBoundingSpheres
{
  vtkPolyData  *_Surface;
  vtkDataArray *_Center;
  vtkDataArray *_Radius;

  ComputeBoundingSpheres(vtkPolyData *surface, vtkDataArray *center, vtkDataArray *radius)
  :
    _Surface(surface), _Center(center), _Radius(radius)
  {}

public:

  /// Compute bounding spheres of specified range of cells
  void operator ()(const blocked_range<vtkIdType> &re) const
  {
    vtkNew<vtkIdList> ptIds;
    double a[3], b[3], c[3], origin[3], radius;

    for (vtkIdType cellId = re.begin(); cellId != re.end(); ++cellId) {
      // Get triangle vertices
      GetCellPoints(_Surface, cellId, ptIds.GetPointer());
      mirtkAssert(ptIds->GetNumberOfIds() == 3, "surface is triangular mesh");

      // Get triangle vertex positions
      _Surface->GetPoint(ptIds->GetId(0), a);
      _Surface->GetPoint(ptIds->GetId(1), b);
      _Surface->GetPoint(ptIds->GetId(2), c);

      // Get center of bounding sphere
      vtkTriangle::TriangleCenter(a, b, c, origin);
      _Center->SetTuple(cellId, origin);

      // Compute radius of bounding sphere
      radius = sqrt(max(max(vtkMath::Distance2BetweenPoints(a, origin),
                            vtkMath::Distance2BetweenPoints(b, origin)),
                            vtkMath::Distance2BetweenPoints(c, origin)));
      _Radius->SetTuple1(cellId, radius);
    }
  }

  /// Compute bounding spheres of surface faces
  static void Run(vtkPolyData *surface, vtkDataArray *center, vtkDataArray *radius)
  {
    center->SetNumberOfComponents(3);
    center->SetNumberOfTuples(surface->GetNumberOfCells());
    radius->SetNumberOfComponents(1);
    radius->SetNumberOfTuples(surface->GetNumberOfCells());
    if (surface->GetNumberOfCells() == 0) return;
    ComputeBoundingSpheres eval(surface, center, radius);
    parallel_for(blocked_range<vtkIdType>(0, surface->GetNumberOfCells()), eval);
  }
};

// -----------------------------------------------------------------------------
/// Check for collisions such as self-intersections and faces too close
class FindCollisions
{
  typedef SurfaceCollisions::IntersectionInfo   IntersectionInfo;
  typedef SurfaceCollisions::IntersectionsArray IntersectionsArray;
  typedef SurfaceCollisions::CollisionType      CollisionType;
  typedef SurfaceCollisions::CollisionsArray    CollisionsArray;
  typedef SurfaceCollisions::CollisionInfo      CollisionInfo;

  SurfaceCollisions       *_Filter;
  vtkAbstractPointLocator *_Locator;
  IntersectionsArray      *_Intersections;
  CollisionsArray         *_Collisions;
  double                   _MaxRadius;
  double                   _MinFrontfaceDistance;
  double                   _MinBackfaceDistance;
  double                   _MinAngleCos;

  /// Check whether point c is on the "left" of the line defined by a and b
  inline int IsLeft(const double a[2], const double b[2], const double c[2]) const
  {
     return sgn((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]));
  }

  /// Compute axes for projection of triangle points to 2D with origin at x
  inline int ProjectionAxes(const double x1[3], const double x2[3], const double n[3], double a[3], double b[3], double u1[2], double u2[2]) const
  {
    vtkMath::Subtract(x2, x1, a);
    u2[0] = vtkMath::Normalize(a);
    if (u2[0] <= .0) return 0;
    u1[0] = u1[1] = u2[1] = .0;
    vtkMath::Cross(n, a, b);
    return 1;
  }

  /// Project triangle point to 2D
  inline void ProjectTo2D(const double origin[3], const double a[3], const double b[3], const double x[3], double u[2]) const
  {
    double v[3];
    vtkMath::Subtract(x, origin, v);
    u[0] = vtkMath::Dot(v, a);
    u[1] = vtkMath::Dot(v, b);
  }

  /// Check if two adjacent triangles with a single common vertex are folding onto one another
  ///
  /// FIXME: A proper check requires an algorithm as outlined in CheckCollision().
  inline bool CheckCollisionOfAdjacentTriangles(vtkIdType cellId2, CollisionInfo &collision, double dp) const
  {
    if (dp < .0) {
      if (_Filter->BackfaceCollisionTest() && collision._Distance < _MinBackfaceDistance) {
        collision._Type = CollisionType::BackfaceCollision;
      } else {
        return false;
      }
    } else {
      if (_Filter->FrontfaceCollisionTest() && collision._Distance < _MinFrontfaceDistance) {
        collision._Type = CollisionType::FrontfaceCollision;
      } else {
        return false;
      }
    }
    collision._CellId = cellId2;
    return true;
  }

  /// Check whether the two faces nearly collide given the closest points recorded by CollisionInfo
  ///
  /// TODO: Use centroids of overlap of triangles after projection to 2D planes for collision test.
  ///
  /// Given a pair of tringles A and B, project A onto plane defined by B and clip it using the
  /// edges of triangle B. Compute centroid of clipped convex polygon and project it back to
  /// triangle A in 3D space (either using Barycentric coordinates or intersection with line
  /// along the normal of triangle B). Repeat these steps for B projected onto plane defined
  /// by triangle A. Use the distance between the overlap centroids and direction vector to
  /// determine whether there is a collision and on which side of triangle A.
  ///
  /// This algorithm should work for both adjacent and non-adjacent triangles. It can also
  /// implicitly determine if there is an overlap of the triangles after projection, i.e.,
  /// whether the clipped polygon has at least three non-colinear points.
  ///
  /// Algorithms:
  /// - Centroid of convex polygon: https://bell0bytes.eu/centroid-convex/
  /// - Clip polygon by another convex polygon (e.g., triangle):
  ///   https://en.wikipedia.org/wiki/Sutherland%E2%80%93Hodgman_algorithm
  inline bool CheckCollision(vtkIdType cellId1, const double n1[3],
                             vtkIdType cellId2, const double n2[3],
                             CollisionInfo &collision) const
  {
    double dp = vtkMath::Dot(n1, n2);
    if (dp >= -_MinAngleCos) {
      return false;
    }
    if (collision._Distance > 1e-9) {
      double v[3];
      vtkMath::Subtract(collision._Point2, collision._Point1, v);
      vtkMath::MultiplyScalar(v, 1 / collision._Distance);
      dp = vtkMath::Dot(v, n1);
    }
    if (abs(dp) <= _MinAngleCos) {
      return false;
    }
    return CheckCollisionOfAdjacentTriangles(cellId2, collision, dp);
  }

  /// Update cell collision type when a near-miss collision was detected
  inline void UpdateCollisionType(CollisionType &type, const CollisionInfo &collision) const
  {
    if (_Filter->IsIntersection(type)) {
      type = CollisionType::Ambiguous;
    } else if (type != CollisionType::NoCollision) {
      if (type != collision._Type) {
        type = CollisionType::Collision;
      }
    } else {
      type = collision._Type;
    }
  }

  /// Update cell collision type when an intersection was detected
  inline void UpdateCollisionType(CollisionType &type, const IntersectionInfo &intersection) const
  {
    if (_Filter->IsCollision(type)) {
      type = CollisionType::Ambiguous;
      return;
    }
    if (intersection._Adjacent) {
      if (type == CollisionType::SelfIntersection) {
        type = CollisionType::Intersection;
      } else {
        type = CollisionType::AdjacentIntersection;
      }
    } else {
      if (type == CollisionType::AdjacentIntersection) {
        type = CollisionType::Intersection;
      } else {
        type = CollisionType::SelfIntersection;
      }
    }
  }

  /// Store information about detected collision
  inline void SaveCollisionInfo(vtkIdType cellId, const CollisionInfo &collision) const
  {
    if (_Collisions == nullptr) return;
    (*_Collisions)[cellId].insert(collision);
  }

  /// Store information about detected collision
  inline void SaveIntersectionInfo(vtkIdType cellId1, const IntersectionInfo &intersection) const
  {
    if (_Intersections == nullptr) return;
    (*_Intersections)[cellId1].insert(intersection);
  }

public:

  void operator ()(const blocked_range<vtkIdType> &re) const
  {
    const double TOL = 1e-6; // Tolerance for point coordinate equality check

    vtkDataArray *mask      = _Filter->Mask();
    vtkPolyData  *surface   = _Filter->Output();
    vtkDataArray *center    = _Filter->GetCenterArray();
    vtkDataArray *radius    = _Filter->GetRadiusArray();
    vtkDataArray *coll_type = _Filter->GetCollisionTypeArray();

    const bool   coll_test     = (_MinFrontfaceDistance > .0 || _MinBackfaceDistance > .0);
    const bool   adj_coll_test = _Filter->AdjacentCollisionTest() && coll_test;
    const double min_distance  = max(_MinFrontfaceDistance, _MinBackfaceDistance);
    const double R             = _MaxRadius + 1.1 * min_distance;

    double tri1[3][3], tri2[3][3], tri1_2D[3][2], tri2_2D[3][2];
    double n1[3], n2[3], p1[3], p2[3], r1, c1[3], d[3], ax1[3], ax2[3], search_radius;
    int    tri12[3], i1, j1, k1, i2, j2, k2, shared_vertex1, shared_vertex2, coplanar;
    bool   intersect, same_side;

    vtkIdType cellId1, cellId2;
    UnorderedSet<vtkIdType> nearbyCellIds;
    vtkNew<vtkIdList> cellIds, ptIds, triPtIds1, triPtIds2;
  
    IntersectionInfo intersection;
    CollisionInfo collision;
    CollisionType type;

    for (cellId1 = re.begin(); cellId1 != re.end(); ++cellId1) {

      // No collision if none detected
      type = CollisionType::NoCollision;

      if (mask && mask->GetComponent(cellId1, 0) == .0) {
        continue; // Skip cell if (un)masked and keep collision type unmodified
      }

      // Get vertices and normal of this triangle
      GetCellPoints(surface, cellId1, triPtIds1.GetPointer());
      mirtkAssert(triPtIds1->GetNumberOfIds() == 3, "surface is triangular mesh");
      surface->GetPoint(triPtIds1->GetId(0), tri1[0]);
      surface->GetPoint(triPtIds1->GetId(1), tri1[1]);
      surface->GetPoint(triPtIds1->GetId(2), tri1[2]);
      vtkTriangle::ComputeNormal(tri1[0], tri1[1], tri1[2], n1);

      // Get bounding sphere
      center->GetTuple(cellId1, c1);
      r1 = radius->GetComponent(cellId1, 0);

      // Find other triangles within search radius
      nearbyCellIds.clear();
      search_radius = min(max(_Filter->MinSearchRadius(), r1 + R), _Filter->MaxSearchRadius());
      _Locator->FindPointsWithinRadius(search_radius, c1, ptIds.GetPointer());
      for (vtkIdType i = 0; i < ptIds->GetNumberOfIds(); ++i) {
        surface->GetPointCells(ptIds->GetId(i), cellIds.GetPointer());
        for (vtkIdType j = 0; j < cellIds->GetNumberOfIds(); ++j) {
          nearbyCellIds.insert(cellIds->GetId(j));
        }
      }
      nearbyCellIds.erase(cellId1);

      // Check for collisions between this triangle and the found nearby triangles
      for (const auto cellId2 : nearbyCellIds) {

        // Get vertex positions of nearby candidate triangle
        GetCellPoints(surface, cellId2, triPtIds2.GetPointer());
        mirtkAssert(triPtIds2->GetNumberOfIds() == 3, "surface is triangular mesh");
        surface->GetPoint(triPtIds2->GetId(0), tri2[0]);
        surface->GetPoint(triPtIds2->GetId(1), tri2[1]);
        surface->GetPoint(triPtIds2->GetId(2), tri2[2]);
        vtkTriangle::ComputeNormal(tri2[0], tri2[1], tri2[2], n2);

        // Get corresponding indices of shared vertices
        for (i1 = 0; i1 < 3; ++i1) {
          tri12[i1] = -1;
          for (i2 = 0; i2 < 3; ++i2) {
            if (triPtIds1->GetId(i1) == triPtIds2->GetId(i2)) {
              tri12[i1] = i2;
              break;
            }
          }
        }

        // Determine whether triangles share one vertex or even an edge
        shared_vertex1 = shared_vertex2 = -1;
        for (i1 = 0; i1 < 3; ++i1) {
          if (tri12[i1] != -1) {
            shared_vertex1 = i1;
            for (i2 = i1 + 1; i2 < 3; ++i2) {
              if (tri12[i2] != -1) {
                shared_vertex2 = i2;
                break;
              }
            }
            break;
          }
        }

        // Triangles with a shared edge can only overlap, but not intersect.
        // This includes near miss collision if a minimum face distance is set.
        if (shared_vertex2 != -1) {
          if (_Filter->AdjacentIntersectionTest() || adj_coll_test) {
            // Determine indices of non-shared vertices
            for (i1 = 0; i1 < 3; ++i1) {
              if (i1 != shared_vertex1 && i1 != shared_vertex2) break;
            }
            for (i2 = 0; i2 < 3; ++i2) {
              if (i2 != tri12[shared_vertex1] && i2 != tri12[shared_vertex2]) break;
            }
            double m[3], v1[3], v2[3];
            vtkMath::Add(tri1[shared_vertex1], tri1[shared_vertex2], m);
            vtkMath::MultiplyScalar(m, 0.5);
            vtkMath::Subtract(tri1[i1], m, v1);
            vtkMath::Normalize(v1);
            vtkMath::Subtract(tri2[i2], m, v2);
            vtkMath::Normalize(v2);
            double v_dot_v = vtkMath::Dot(v1, v2);
            // When angle is close to zero degrees, the triangles intersect
            if (1 - v_dot_v < TOL) {
              if (_Filter->AdjacentIntersectionTest()) {
                intersection._CellId = cellId2;
                intersection._Adjacent = true;
                UpdateCollisionType(type, intersection);
                SaveIntersectionInfo(cellId1, intersection);
              }
            // Otherwise, they may be folding onto one another
            } else if (adj_coll_test && v_dot_v > 0.75 * _MinAngleCos) {
              center->GetTuple(cellId2, collision._Point2);
              vtkPlane::ProjectPoint(collision._Point2, tri1[0], n1, collision._Point1);
              collision._Distance = sqrt(vtkMath::Distance2BetweenPoints(collision._Point1, collision._Point2));
              if (collision._Distance < min_distance) {
                double n_dot_v = vtkMath::Dot(n1, v2);
                if (CheckCollisionOfAdjacentTriangles(cellId2, collision, n_dot_v)) {
                  UpdateCollisionType(type, collision);
                  SaveCollisionInfo(cellId1, collision);
                }
              }
            }
          }
        }
        // If triangles share a single vertex, use different triangle/triangle
        // intersection check which gives us the points forming the line/point
        // of intersection, but does not handle the coplanar case itselfg
        else if (shared_vertex1 != -1) {
          intersection._CellId = -1;
          if (_Filter->AdjacentIntersectionTest()) {
            // Perform coplanarity test with our TOL instead of smaller tolerance used by vtkIntersectionPolyDataFilter
            coplanar = 0;
            if (fequal(n1[0], n2[0], TOL) && fequal(n1[1], n2[1], TOL) && fequal(n1[2], n2[2], TOL)) {
              const double b1 = -vtkMath::Dot(n1, tri1[0]);
              const double b2 = -vtkMath::Dot(n2, tri2[0]);
              if (fequal(b1, b2, TOL)) coplanar = 1;
            }
            if (coplanar) {
              // TODO: Either one triangle fully contained within the other or one edge of the first
              //       triangle intersects an edge of the second triangle.
            } else if (Triangle::TriangleTriangleIntersection(tri1[0], tri1[1], tri1[2],
                                                              tri2[0], tri2[1], tri2[2],
                                                              coplanar, p1, p2)) {
              // Ignore valid intersection of single shared vertex
              if (!fequal(p1[0], p2[0], TOL) || !fequal(p1[1], p2[1], TOL) || !fequal(p1[2], p2[2], TOL)) {
                intersection._CellId = cellId2;
              }
            }
          }
          if (intersection._CellId >= 0) {
            intersection._Adjacent = true;
            UpdateCollisionType(type, intersection);
            SaveIntersectionInfo(cellId1, intersection);
          } else if (adj_coll_test) {
            center->GetTuple(cellId1, collision._Point1);
            center->GetTuple(cellId2, collision._Point2);
            collision._Distance = min(vtkPlane::DistanceToPlane(collision._Point2, n1, tri1[0]),
                                      vtkPlane::DistanceToPlane(collision._Point1, n2, tri2[0]));
            if (collision._Distance < min_distance) {
              // Reorder vertex indices such that shared vertex is at first position
              i1 = shared_vertex1;
              j1 = i1 == 0 ? 1 : 0;
              k1 = i1 == 1 ? 2 : j1 + 1;
              i2 = tri12[shared_vertex1];
              j2 = i2 == 0 ? 1 : 0;
              k2 = i2 == 1 ? 2 : j2 + 1;
              // Vector from shared vertex to mid-point of opposite edge of first triangle
              double v1[3], p1[3];
              vtkMath::Add(tri1[j1], tri1[k1], p1);
              vtkMath::MultiplyScalar(p1, 0.5);
              vtkMath::Subtract(p1, tri1[i1], v1);
              vtkMath::Normalize(v1);
              // Vector from shared vertex to mid-point of opposite edge of second triangle
              double v2[3], p2[3];
              vtkMath::Add(tri2[j2], tri2[k2], p2);
              vtkMath::MultiplyScalar(p2, 0.5);
              vtkMath::Subtract(p2, tri2[i2], v2);
              vtkMath::Normalize(v2);
              // Check weak condition for triangles to fold onto one another
              double v_dot_v = vtkMath::Dot(v1, v2);
              double n_dot_n = vtkMath::Dot(n1, n2);
              if (n_dot_n < -0.4 && v_dot_v > 0.2) {
                double p3[3], v3[3];
                vtkPlane::ProjectPoint(p1, tri2[i2], n2, p3);
                vtkMath::Subtract(p3, tri2[i2], v3);
                vtkMath::Normalize(v3);
                v_dot_v = vtkMath::Dot(v1, v3);
                if (v_dot_v > _MinAngleCos) {
                  // Check if projection of triangles to common 2D plane intersects
                  if (ProjectionAxes(tri1[i1], tri1[j1], n1, ax1, ax2, tri1_2D[0], tri1_2D[1])) {
                    ProjectTo2D(tri1[i1], ax1, ax2, tri1[k1], tri1_2D[2]);
                    ProjectTo2D(tri1[i1], ax1, ax2, tri2[i2], tri2_2D[0]);
                    ProjectTo2D(tri1[i1], ax1, ax2, tri2[j2], tri2_2D[1]);
                    ProjectTo2D(tri1[i1], ax1, ax2, tri2[k2], tri2_2D[2]);
                    if (Triangle::TriangleTriangleOverlap(tri1_2D[0], tri1_2D[1], tri1_2D[2],
                                                          tri2_2D[0], tri2_2D[1], tri2_2D[2])) {
                      double n_dot_v = vtkMath::Dot(n1, v2);
                      if (CheckCollisionOfAdjacentTriangles(cellId2, collision, n_dot_v)) {
                        UpdateCollisionType(type, collision);
                        SaveCollisionInfo(cellId1, collision);
                      }
                    }
                  }
                }
              }
            }
          }
        }
        // In case of non-adjacent triangles, use fast intersection test which
        // also checks for overlap of coplanar triangles
        else if (_Filter->NonAdjacentIntersectionTest() &&
                 Triangle::TriangleTriangleIntersection(tri1[0], tri1[1], tri1[2],
                                                        tri2[0], tri2[1], tri2[2])) {
          intersection._CellId = cellId2;
          intersection._Adjacent = false;
          UpdateCollisionType(type, intersection);
          SaveIntersectionInfo(cellId1, intersection);
        }
        // If self-intersection check of non-adjacent triangles disabled or negative,
        // check for near miss collision if minimum distance set
        else if (coll_test) {
          if (_Filter->FastCollisionTest()) {
            center->GetTuple(cellId1, collision._Point1);
            center->GetTuple(cellId2, collision._Point2);
            collision._Distance = sqrt(vtkMath::Distance2BetweenPoints(collision._Point1, collision._Point2));
          } else {
            collision._Distance = Triangle::DistanceBetweenTriangles(tri1[0], tri1[1], tri1[2], n1,
                                                                     tri2[0], tri2[1], tri2[2], n2,
                                                                     collision._Point1, collision._Point2);
          }
          if (collision._Distance < min_distance) {
            if (CheckCollision(cellId1, n1, cellId2, n2, collision)) {
              UpdateCollisionType(type, collision);
              SaveCollisionInfo(cellId1, collision);
            }
          }
        }
      }

      // Set collision type of this cell
      coll_type->SetComponent(cellId1, 0, type);
    }
  }

  /// Find collision and self-intersections
  static void Run(SurfaceCollisions  *filter,
                  IntersectionsArray *intersections,
                  CollisionsArray    *collisions)
  {
    vtkDataArray *radius = filter->GetRadiusArray();
    vtkSmartPointer<vtkAbstractPointLocator> locator;
    locator = vtkSmartPointer<vtkOctreePointLocator>::New();
    locator->SetDataSet(filter->Output());
    locator->BuildLocator();
    FindCollisions body;
    body._Filter               = filter;
    body._Locator              = locator;
    body._Intersections        = intersections;
    body._Collisions           = collisions;
    body._MaxRadius            = radius->GetRange(0)[1];
    body._MinAngleCos          = cos(filter->MaxAngle() * rad_per_deg);
    body._MinFrontfaceDistance = filter->MinFrontfaceDistance();
    body._MinBackfaceDistance  = filter->MinBackfaceDistance();
    if (!filter->BackfaceCollisionTest())  body._MinFrontfaceDistance = .0;
    if (!filter->FrontfaceCollisionTest()) body._MinBackfaceDistance  = .0;
    blocked_range<vtkIdType> cellIds(0, filter->Output()->GetNumberOfCells());
    parallel_for(cellIds, body);
  }
};


} // namespace SurfaceCollisionsUtils
using namespace SurfaceCollisionsUtils;

// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
void SurfaceCollisions::CopyAttributes(const SurfaceCollisions &other)
{
  _Mask                        = other._Mask;
  _UseInputBoundingSpheres     = other._UseInputBoundingSpheres;
  _MinSearchRadius             = other._MinSearchRadius;
  _MaxSearchRadius             = other._MaxSearchRadius;
  _MinFrontfaceDistance        = other._MinFrontfaceDistance;
  _MinBackfaceDistance         = other._MinBackfaceDistance;
  _MaxAngle                    = other._MaxAngle;
  _AdjacentCollisionTest       = other._AdjacentCollisionTest;
  _AdjacentIntersectionTest    = other._AdjacentIntersectionTest;
  _NonAdjacentIntersectionTest = other._NonAdjacentIntersectionTest;
  _FrontfaceCollisionTest      = other._FrontfaceCollisionTest;
  _BackfaceCollisionTest       = other._BackfaceCollisionTest;
  _FastCollisionTest           = other._FastCollisionTest;
  _StoreIntersectionDetails    = other._StoreIntersectionDetails;
  _StoreCollisionDetails       = other._StoreCollisionDetails;
  _ResetCollisionType          = other._ResetCollisionType;
  _Intersections               = other._Intersections;
  _Collisions                  = other._Collisions;
}

// -----------------------------------------------------------------------------
SurfaceCollisions::SurfaceCollisions()
:
  _UseInputBoundingSpheres(false),
  _MinSearchRadius(.0),
  _MaxSearchRadius(inf),
  _MinFrontfaceDistance(.0),
  _MinBackfaceDistance(.0),
  _MaxAngle(90.0),
  _AdjacentCollisionTest(false),
  _AdjacentIntersectionTest(false),
  _NonAdjacentIntersectionTest(false),
  _FrontfaceCollisionTest(false),
  _BackfaceCollisionTest(false),
  _FastCollisionTest(false),
  _StoreIntersectionDetails(true),
  _StoreCollisionDetails(true),
  _ResetCollisionType(true)
{
}

// -----------------------------------------------------------------------------
SurfaceCollisions::SurfaceCollisions(const SurfaceCollisions &other)
:
  SurfaceFilter(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
SurfaceCollisions &SurfaceCollisions::operator =(const SurfaceCollisions &other)
{
  if (this != &other) {
    SurfaceFilter::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
SurfaceCollisions::~SurfaceCollisions()
{
}

// =============================================================================
// Output attributes
// =============================================================================

// -----------------------------------------------------------------------------
vtkDataArray *SurfaceCollisions::GetCenterArray() const
{
  return _Output->GetCellData()->GetArray(BOUNDING_SPHERE_CENTER);
}

// -----------------------------------------------------------------------------
vtkDataArray *SurfaceCollisions::GetRadiusArray() const
{
  return _Output->GetCellData()->GetArray(BOUNDING_SPHERE_RADIUS);
}

// -----------------------------------------------------------------------------
vtkDataArray *SurfaceCollisions::GetCollisionTypeArray() const
{
  return _Output->GetCellData()->GetArray(COLLISION_TYPE);
}

// -----------------------------------------------------------------------------
SurfaceCollisions::CollisionType SurfaceCollisions::GetCollisionType(int cellId) const
{
  return static_cast<CollisionType>(static_cast<int>(GetCollisionTypeArray()->GetComponent(cellId, 0)));
}

// -----------------------------------------------------------------------------
#ifdef VTK_USE_64BIT_IDS
SurfaceCollisions::CollisionType SurfaceCollisions::GetCollisionType(vtkIdType cellId) const
{
  return static_cast<CollisionType>(static_cast<int>(GetCollisionTypeArray()->GetComponent(cellId, 0)));
}
#endif // VTK_USE_64BIT_IDS

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void SurfaceCollisions::Initialize()
{
  // Initialize base class
  SurfaceFilter::Initialize();

  // Check input mesh
  if (!IsTriangularMesh(_Input)) {
    cerr << this->NameOfType() << "::Initialize: Input mesh must have triangular faces!" << endl;
    exit(1);
  }

  // Check mask
  if (_Mask != nullptr) {
    if (_Mask->GetNumberOfTuples() != _Input->GetNumberOfCells() || _Mask->GetNumberOfComponents() == 0) {
      cerr << this->NameOfType() << "::Initialize: Input mask must have one entry for each triangular face of input mesh" << endl;
      exit(1);
    }
  }

  // Reset stored collisions from previous run
  _Intersections.clear();
  _Collisions.clear();

  // Enable default collision tests if none selected
  if (!_AdjacentIntersectionTest && !_NonAdjacentIntersectionTest &&
      !_BackfaceCollisionTest    && !_FrontfaceCollisionTest) {
    _AdjacentIntersectionTest    = true;
    _NonAdjacentIntersectionTest = true;
    _FrontfaceCollisionTest      = (_MinFrontfaceDistance > .0);
    _BackfaceCollisionTest       = (_MinBackfaceDistance  > .0);
  }

  // Ensure that any collision tests are enabled
  if (!_AdjacentIntersectionTest && !_NonAdjacentIntersectionTest &&
     (_MinFrontfaceDistance <= .0 || !_FrontfaceCollisionTest) &&
     (_MinBackfaceDistance  <= .0 || !_BackfaceCollisionTest)) {
    cerr << this->NameOfType() << "::Initialize: No collision tests enabled or minimum distances not positive" << endl;
    exit(1);
  }

  // Ensure distance/radius parameters are non-negative
  if (_MinFrontfaceDistance < .0) _MinFrontfaceDistance = .0;
  if (_MinBackfaceDistance  < .0) _MinBackfaceDistance  = .0;
  if (_MinSearchRadius      < .0) _MinSearchRadius      = .0;

  // Precompute bounding spheres and prepare auxiliary data arrays
  bool compute_bounding_spheres = !_UseInputBoundingSpheres;
  bool reset_collision_types    = _ResetCollisionType;
  vtkSmartPointer<vtkDataArray> center, radius, types;
  center = _Output->GetCellData()->GetArray(BOUNDING_SPHERE_CENTER);
  radius = _Output->GetCellData()->GetArray(BOUNDING_SPHERE_RADIUS);
  types  = _Output->GetCellData()->GetArray(COLLISION_TYPE);

  if (!center || !radius) {
    if (!center) {
      center = vtkSmartPointer<vtkFloatArray>::New();
      center->SetName(BOUNDING_SPHERE_CENTER);
      _Output->GetCellData()->AddArray(center);
    }
    if (!radius) {
      radius = vtkSmartPointer<vtkFloatArray>::New();
      radius->SetName(BOUNDING_SPHERE_RADIUS);
      _Output->GetCellData()->AddArray(radius);
    }
    compute_bounding_spheres = true;
  }
  if (compute_bounding_spheres) {
    ComputeBoundingSpheres::Run(_Output, center, radius);
  }
  if (!types) {
    types = vtkSmartPointer<vtkUnsignedCharArray>::New();
    types->SetNumberOfComponents(1);
    types->SetNumberOfTuples(_Output->GetNumberOfCells());
    types->SetName(COLLISION_TYPE);
    _Output->GetCellData()->AddArray(types);
    reset_collision_types = true;
  }
  if (reset_collision_types) {
    types->FillComponent(0, static_cast<double>(NoCollision));
  }

  if (_StoreIntersectionDetails && (_AdjacentIntersectionTest || _NonAdjacentIntersectionTest)) {
    _Intersections.resize(_Output->GetNumberOfCells());
  }
  if (_StoreCollisionDetails && (_FrontfaceCollisionTest || _BackfaceCollisionTest)) {
    _Collisions.resize(_Output->GetNumberOfCells());
  }
}

// -----------------------------------------------------------------------------
void SurfaceCollisions::Execute()
{
  if (_Output->GetNumberOfCells() > 0) {
    FindCollisions::Run(this, _StoreIntersectionDetails ? &_Intersections : nullptr,
                              _StoreCollisionDetails    ? &_Collisions    : nullptr);
  }
}

// -----------------------------------------------------------------------------
bool SurfaceCollisions::FoundIntersections() const
{
  for (vtkIdType cellId = 0; cellId < _Output->GetNumberOfCells(); ++cellId) {
    if (IsIntersection(GetCollisionType(cellId))) return true;
  }
  return false;
}

// -----------------------------------------------------------------------------
bool SurfaceCollisions::FoundCollisions() const
{
  for (vtkIdType cellId = 0; cellId < _Output->GetNumberOfCells(); ++cellId) {
    if (IsCollision(GetCollisionType(cellId))) return true;
  }
  return false;
}


} // namespace mirtk
