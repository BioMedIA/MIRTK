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

#include "mirtk/SurfaceCollisions.h"

#include "mirtk/Assert.h"
#include "mirtk/Math.h"
#include "mirtk/Triangle.h"
#include "mirtk/PointSetUtils.h"
#include "mirtk/PointLocator.h"
#include "mirtk/Parallel.h"
#include "mirtk/Profiling.h"
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
    vtkIdType npts, *pts;
    double a[3], b[3], c[3], origin[3], radius;

    for (vtkIdType cellId = re.begin(); cellId != re.end(); ++cellId) {
      // Get triangle vertices
      _Surface->GetCellPoints(cellId, npts, pts);
      mirtkAssert(npts == 3, "surface is triangular mesh");

      // Get triangle vertex positions
      _Surface->GetPoint(pts[0], a);
      _Surface->GetPoint(pts[1], b);
      _Surface->GetPoint(pts[2], c);

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

public:

  void operator ()(const blocked_range<vtkIdType> &re) const
  {
    const double TOL = 1e-6; // Tolerance for point coordinate equality check

    vtkDataArray *mask      = _Filter->Mask();
    vtkPolyData  *surface   = _Filter->Output();
    vtkDataArray *center    = _Filter->GetCenterArray();
    vtkDataArray *radius    = _Filter->GetRadiusArray();
    vtkDataArray *coll_type = _Filter->GetCollisionTypeArray();

    const bool   coll_test    = (_MinFrontfaceDistance > .0 || _MinBackfaceDistance > .0);
    const double min_distance = max(_MinFrontfaceDistance, _MinBackfaceDistance);
    const double R            = _MaxRadius + 1.1 * min_distance;

    double         tri1[3][3], tri2[3][3], tri1_2D[3][2], tri2_2D[3][2];
    double         n1[3], n2[3], p1[3], p2[3], r1, c1[3], d[3], search_radius, dot;
    int            tri12[3], i1, i2, shared_vertex1, shared_vertex2, coplanar, s1, s2;
    vtkIdType      npts, *pts1, *pts2, *cells, cellId1, cellId2;
    unsigned short ncells;
    CollisionInfo  collision;
    CollisionType  type;

    vtkSmartPointer<vtkIdList> ptIds   = vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();

    for (cellId1 = re.begin(); cellId1 != re.end(); ++cellId1) {

      // No collision if none detected
      type = CollisionType::NoCollision;

      if (mask && mask->GetComponent(cellId1, 0) == .0) {
        continue; // Skip cell if (un)masked and keep collision type unmodified
      }

      // Get vertices and normal of this triangle
      surface->GetCellPoints(cellId1, npts, pts1);
      mirtkAssert(npts == 3, "surface is triangular mesh");
      surface->GetPoint(pts1[0], tri1[0]);
      surface->GetPoint(pts1[1], tri1[1]);
      surface->GetPoint(pts1[2], tri1[2]);
      vtkTriangle::ComputeNormal(tri1[0], tri1[1], tri1[2], n1);

      // Get bounding sphere
      center->GetTuple(cellId1, c1);
      r1 = radius->GetComponent(cellId1, 0);

      // Find other triangles within search radius
      search_radius = min(max(_Filter->MinSearchRadius(), r1 + R), _Filter->MaxSearchRadius());
      _Locator->FindPointsWithinRadius(search_radius, c1, ptIds);
      cellIds->Reset();
      for (vtkIdType i = 0; i < ptIds->GetNumberOfIds(); ++i) {
        surface->GetPointCells(ptIds->GetId(i), ncells, cells);
        for (unsigned short j = 0; j < ncells; ++j) {
          if (cells[j] != cellId1) cellIds->InsertUniqueId(cells[j]);
        }
      }

      // Check for collisions between this triangle and the found nearby triangles
      for (vtkIdType i = 0; i < cellIds->GetNumberOfIds(); ++i) {
        cellId2 = cellIds->GetId(i);

        // Get vertex positions of nearby candidate triangle
        surface->GetCellPoints(cellId2, npts, pts2);
        surface->GetPoint(pts2[0], tri2[0]);
        surface->GetPoint(pts2[1], tri2[1]);
        surface->GetPoint(pts2[2], tri2[2]);

        // Get corresponding indices of shared vertices
        for (i1 = 0; i1 < 3; ++i1) {
          tri12[i1] = -1;
          for (i2 = 0; i2 < 3; ++i2) {
            if (pts1[i1] == pts2[i2]) {
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

        // Triangles with a shared edge can only overlap, but not intersect
        if (shared_vertex2 != -1) {
          if (_Filter->AdjacentIntersectionTest()) {
            coplanar = 0;
            for (i2 = 0; i2 < 3; ++i2) {
              if (i2 != tri12[shared_vertex1] && i2 != tri12[shared_vertex2]) {
                coplanar = int(vtkPlane::DistanceToPlane(tri2[i2], tri1[0], n1) < TOL);
                break;
              }
            }
            if (coplanar) {
              // Project triangles into their common plane
              vtkTriangle::ProjectTo2D(tri1   [0], tri1   [1], tri1   [2],
                                       tri1_2D[0], tri1_2D[1], tri1_2D[2]);
              vtkTriangle::ProjectTo2D(tri2   [0], tri2   [1], tri2   [2],
                                       tri2_2D[0], tri2_2D[1], tri2_2D[2]);
              // Check on which side of the shared edge the non-shared vertices are
              s1 = 0;
              for (i1 = 0; i1 < 3; ++i1) {
                if (i1 != shared_vertex1 && i1 != shared_vertex2) {
                  s1 = IsLeft(tri1_2D[shared_vertex1], tri1_2D[shared_vertex2], tri1_2D[i1]);
                  break;
                }
              }
              // Note: re-uses i2 index found outside of this if-block
              s2 = IsLeft(tri1_2D[shared_vertex1], tri1_2D[shared_vertex2], tri2_2D[i2]);
              // If they are on the same side of the edge, the triangles must overlap
              if (s1 == s2) {
                //intersections.insert(IntersectionInfo(cellId2, true));
                if (_Filter->IsCollision(type)) {
                  type = CollisionType::Ambiguous;
                } else if (type == CollisionType::SelfIntersection) {
                  type = CollisionType::Intersection;
                } else {
                  type = CollisionType::AdjacentIntersection;
                }
              }
            }
          }
        }
        // If triangles share a single vertex, use different triangle/triangle
        // intersection check which gives us the points forming the line/point
        // of intersection, but does not handle the coplanar case itself
        else if (shared_vertex1 != -1) {
          if (_Filter->AdjacentIntersectionTest()) {
            // Perform coplanarity test with our TOL instead of the smaller
            // tolerance value used by vtkIntersectionPolyDataFilter
            coplanar = 0;
            Triangle::Normal(tri2[0], tri2[1], tri2[2], n2);
            if (fequal(n1[0], n2[0], TOL) &&
                fequal(n1[1], n2[1], TOL) &&
                fequal(n1[2], n2[2], TOL)) {
              const double b1 = -vtkMath::Dot(n1, tri1[0]);
              const double b2 = -vtkMath::Dot(n2, tri2[0]);
              if (fequal(b1, b2, TOL)) coplanar = 1;
            }
            if (coplanar) {
              // TODO: Either one triangle fully contained within the other
              //       or one edge of the first triangle intersects an edge of
              //       the second triangle.
            }
            else if (Triangle::TriangleTriangleIntersection(tri1[0], tri1[1], tri1[2],
                                                            tri2[0], tri2[1], tri2[2],
                                                            coplanar, p1, p2)) {
              // Ignore valid intersection of single shared vertex
              if (!fequal(p1[0], p2[0], TOL) ||
                  !fequal(p1[1], p2[1], TOL) ||
                  !fequal(p1[2], p2[2], TOL)) {
                if (_Intersections) {
                  (*_Intersections)[cellId1].insert(IntersectionInfo(cellId2, true));
                }
                if (_Filter->IsCollision(type)) {
                  type = CollisionType::Ambiguous;
                } else if (type == CollisionType::SelfIntersection) {
                  type = CollisionType::Intersection;
                } else {
                  type = CollisionType::AdjacentIntersection;
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
          if (_Intersections) {
            (*_Intersections)[cellId1].insert(IntersectionInfo(cellId2, false));
          }
          if (_Filter->IsCollision(type)) {
            type = CollisionType::Ambiguous;
          } else if (type == CollisionType::AdjacentIntersection) {
            type = CollisionType::Intersection;
          } else {
            type = CollisionType::SelfIntersection;
          }
        }
        // If self-intersection check of non-adjacent triangles disabled or negative,
        // check for near miss collision if minimum distance set
        else if (coll_test) {
          if (_Filter->FastCollisionTest()) {
            collision._Distance = Triangle::DistanceBetweenCenters(
                                      tri1[0], tri1[1], tri1[2],
                                      tri2[0], tri2[1], tri2[2],
                                      collision._Point1, collision._Point2);
          } else {
            Triangle::Normal(tri2[0], tri2[1], tri2[2], n2);
            collision._Distance = Triangle::DistanceBetweenTriangles(
                                      tri1[0], tri1[1], tri1[2], n1,
                                      tri2[0], tri2[1], tri2[2], n2,
                                      collision._Point1, collision._Point2);
          }
          if (collision._Distance < min_distance) {
            if (collision._Distance > 0.) {
              vtkMath::Subtract(collision._Point2, collision._Point1, d);
              dot = vtkMath::Dot(d, n1) / vtkMath::Norm(d);
            } else {
              dot = 0.;
            }
            if (abs(dot) >= _MinAngleCos) {
              if (dot < .0) {
                if (_Filter->BackfaceCollisionTest() && collision._Distance < _MinBackfaceDistance) {
                  collision._Type = CollisionType::BackfaceCollision;
                } else {
                  collision._Type = CollisionType::NoCollision;
                }
              } else {
                if (_Filter->FrontfaceCollisionTest() && collision._Distance < _MinFrontfaceDistance) {
                  collision._Type = CollisionType::FrontfaceCollision;
                } else {
                  collision._Type = CollisionType::NoCollision;
                }
              }
              if (collision._Type != CollisionType::NoCollision) {
                collision._CellId = cellId2;
                if (_Collisions) {
                  (*_Collisions)[cellId1].insert(collision);
                }
                if (_Filter->IsIntersection(type)) {
                  type = CollisionType::Ambiguous;
                } else if (type != CollisionType::NoCollision) {
                  if (type != collision._Type) type = CollisionType::Collision;
                } else {
                  type = collision._Type;
                }
              }
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
SurfaceCollisions::CollisionType SurfaceCollisions::GetCollisionType(vtkIdType cellId) const
{
  return static_cast<CollisionType>(static_cast<int>(GetCollisionTypeArray()->GetComponent(cellId, 0)));
}

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
