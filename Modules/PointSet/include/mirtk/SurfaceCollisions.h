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

#ifndef MIRTK_SurfaceCollisions_H
#define MIRTK_SurfaceCollisions_H

#include "mirtk/SurfaceFilter.h"

#include "mirtk/Math.h"
#include "mirtk/Array.h"
#include "mirtk/OrderedSet.h"
#include "mirtk/Memory.h"
#include "mirtk/PointSetExport.h"

#include "vtkSmartPointer.h"

class vtkPolyData;
class vtkDataArray;


namespace mirtk {


/**
 * Auxiliary class used to find self-collisions of a triangulated surface mesh
 *
 * Instances of this class detect different types of self-collisions such as
 * in particular self-intersections between non-adjacent as well as adjacent
 * triangular faces and a list of non-adjacent faces which are about to collide,
 * i.e., very close to each other. They are used to impose either hard or soft
 * non-self-intersection constraints on a deformable surface.
 */
class SurfaceCollisions : public SurfaceFilter
{
  mirtkObjectMacro(SurfaceCollisions);

  // ---------------------------------------------------------------------------
  // Types
public:

  /// Type of self-collision
  enum CollisionType
  {
    NoCollision,          ///< No collision with another triangle within search range
    Collision,            ///< Unknown or mixed collision types
    FrontfaceCollision,   ///< Second triangle is within minimum distance towards the front of first triangle
    BackfaceCollision,    ///< Second triangle is within minimum distance towards the back of first triangle
    Intersection,         ///< Unknown or mixed intersection types
    SelfIntersection,     ///< Non-adjacent triangles intersect each other
    AdjacentIntersection, ///< Adjacent triangles intersect each other
    Ambiguous             ///< Both collisions and self-intersections found
  };

  // Names of data arrays
  MIRTK_PointSet_EXPORT static const char * const BOUNDING_SPHERE_CENTER;
  MIRTK_PointSet_EXPORT static const char * const BOUNDING_SPHERE_RADIUS;
  MIRTK_PointSet_EXPORT static const char * const COLLISION_TYPE;

  /// Whether a given collision type indicates a near miss collision
  static bool IsCollision(CollisionType);

  /// Whether a given collision type indicates a self-intersection
  static bool IsIntersection(CollisionType);

  /// Structure storing information about detected self-intersection
  struct IntersectionInfo
  {
    int  _CellId;   ///< ID of other triangle
    bool _Adjacent; ///< Whether triangles are adjacent

    IntersectionInfo(int cellId = -1, bool adj = false)
    :
      _CellId(cellId), _Adjacent(adj)
    {}

    IntersectionInfo(const IntersectionInfo &other)
    :
      _CellId(other._CellId), _Adjacent(other._Adjacent)
    {}

    IntersectionInfo &operator =(const IntersectionInfo &other)
    {
      _CellId   = other._CellId;
      _Adjacent = other._Adjacent;
      return *this;
    }

    bool operator <(const IntersectionInfo &rhs) const
    {
      return _CellId < rhs._CellId;
    }
  };

  /// Structure storing information about detected collision
  struct CollisionInfo
  {
    int           _CellId;    ///< ID of other triangle
    double        _Point1[3]; ///< Point in this  triangle which is closest
    double        _Point2[3]; ///< Point in other triangle which is closest
    double        _Distance;  ///< Distance between closest points
    CollisionType _Type;      ///< Type of collision

    CollisionInfo(int cellId = -1)
    :
      _CellId(cellId),
      _Distance(numeric_limits<double>::quiet_NaN()),
      _Type(Collision)
    {
      memset(_Point1, 0, 3 * sizeof(double));
      memset(_Point2, 0, 3 * sizeof(double));
    }

    CollisionInfo(int cellId, const double p1[3], const double p2[3], CollisionType type = Collision)
    :
      _CellId(cellId),
      _Type(type)
    {
      memcpy(_Point1, p1, 3 * sizeof(double));
      memcpy(_Point2, p2, 3 * sizeof(double));
      const double a = _Point2[0] - _Point1[0];
      const double b = _Point2[1] - _Point1[1];
      const double c = _Point2[2] - _Point1[2];
      _Distance = sqrt(a * a + b * b + c * c);
    }

    CollisionInfo(int cellId, const double p1[3], const double p2[3], double dist, CollisionType type = Collision)
    :
      _CellId(cellId),
      _Distance(dist),
      _Type(type)
    {
      memcpy(_Point1, p1, 3 * sizeof(double));
      memcpy(_Point2, p2, 3 * sizeof(double));
    }

    CollisionInfo(const CollisionInfo &other)
    :
      _CellId(other._CellId),
      _Distance(other._Distance),
      _Type(other._Type)
    {
      memcpy(_Point1, other._Point1, 3 * sizeof(double));
      memcpy(_Point2, other._Point2, 3 * sizeof(double));
    }

    CollisionInfo &operator =(const CollisionInfo &other)
    {
      _CellId   = other._CellId;
      _Distance = other._Distance;
      _Type     = other._Type;
      memcpy(_Point1, other._Point1, 3 * sizeof(double));
      memcpy(_Point2, other._Point2, 3 * sizeof(double));
      return *this;
    }

    bool operator <(const CollisionInfo &rhs) const
    {
      return _CellId < rhs._CellId;
    }
  };

  typedef OrderedSet<IntersectionInfo>     IntersectionsSet;
  typedef IntersectionsSet::const_iterator IntersectionsIterator;
  typedef Array<IntersectionsSet>          IntersectionsArray;

  typedef OrderedSet<CollisionInfo>     CollisionsSet;
  typedef CollisionsSet::const_iterator CollisionsIterator;
  typedef Array<CollisionsSet>          CollisionsArray;

  // ---------------------------------------------------------------------------
  // Attributes
private:

  /// Optional mask of input triangles to check
  mirtkPublicAttributeMacro(vtkSmartPointer<vtkDataArray>, Mask);

  /// Use BoundingSphereCenter and BoundingSphereRadius cell data arrays of input
  mirtkPublicAttributeMacro(bool, UseInputBoundingSpheres);

  /// Minimum search radius around triangle center used to determine the set
  /// of nearby triangles to be tested for collisions with this triangle
  mirtkPublicAttributeMacro(double, MinSearchRadius);

  /// Maximum search radius around triangle center used to determine the set
  /// of nearby triangles to be tested for collisions with this triangle
  mirtkPublicAttributeMacro(double, MaxSearchRadius);

  /// Minimum required distance between non-adjacent triangular faces
  ///
  /// This distance threshold applies when the center point of the other
  /// face is in front of the current triangle whose collisions with other
  /// triangles is being determined (cf. CollisionType::FrontfaceCollision).
  mirtkPublicAttributeMacro(double, MinFrontfaceDistance);

  /// Minimum required distance between non-adjacent triangular faces
  ///
  /// This distance threshold applies when the center point of the other
  /// face is at the backside of the current triangle whose collisions with other
  /// triangles is being determined (cf. CollisionType::BackfaceCollision).
  mirtkPublicAttributeMacro(double, MinBackfaceDistance);

  /// Maximum angle between face normal and center to center point vector
  /// required for collision to be detected
  mirtkPublicAttributeMacro(double, MaxAngle);

  /// Whether to perform intersection tests between adjacent triangles
  mirtkPublicAttributeMacro(bool, AdjacentIntersectionTest);

  /// Whether to perform intersection tests between non-adjacent triangles
  mirtkPublicAttributeMacro(bool, NonAdjacentIntersectionTest);

  /// Whether to detect collisions with the front of a face
  mirtkPublicAttributeMacro(bool, FrontfaceCollisionTest);

  /// Whether to detect collisions with the back of a face
  mirtkPublicAttributeMacro(bool, BackfaceCollisionTest);

  /// Whether to use fast, approximate collision test
  mirtkPublicAttributeMacro(bool, FastCollisionTest);

  /// Whether to store detailed information about found self-intersections
  mirtkPublicAttributeMacro(bool, StoreIntersectionDetails);

  /// Whether to store detailed information about found collisions
  mirtkPublicAttributeMacro(bool, StoreCollisionDetails);

  /// Whether to clear collision type for unchecked triangles (cf. _Mask)
  mirtkPublicAttributeMacro(bool, ResetCollisionType);

  /// Found self-intersections per face
  /// \note Only non-empty after Run when _StoreIntersectionDetails is \c true.
  mirtkReadOnlyAttributeMacro(IntersectionsArray, Intersections);

  /// Found collisions per face
  /// \note Only non-empty after Run when _StoreCollisionDetails is \c true.
  mirtkReadOnlyAttributeMacro(CollisionsArray, Collisions);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const SurfaceCollisions &);

  // ---------------------------------------------------------------------------
  // Construction/destruction
public:

  /// Constructor
  SurfaceCollisions();

  /// Copy constructor
  SurfaceCollisions(const SurfaceCollisions &);

  /// Assignment operator
  SurfaceCollisions &operator =(const SurfaceCollisions &);

  /// Destructor
  virtual ~SurfaceCollisions();

  // ---------------------------------------------------------------------------
  // Execution
protected:

  /// Initialize filter after input and parameters are set
  virtual void Initialize();

  /// Execute filter
  virtual void Execute();

  // ---------------------------------------------------------------------------
  // Output
public:

  /// Get cell data array storing radii of bounding spheres
  vtkDataArray *GetRadiusArray() const;

  /// Get cell data array storing center points of bounding spheres
  vtkDataArray *GetCenterArray() const;

  /// Get cell data array storing type of cell collision
  vtkDataArray *GetCollisionTypeArray() const;

  /// Get type of cell collision
  CollisionType GetCollisionType(int) const;

  /// Get type of cell collision
  CollisionType GetCollisionType(vtkIdType) const;

  /// Set of intersections of other faces with the specified cell
  /// \note Use only when intersection tests enabled.
  const IntersectionsSet &Intersections(int) const;

  /// Set of collisions of other faces with the specified cell
  /// \note Use only when collision tests enabled.
  const CollisionsSet &Collisions(int) const;

  /// Whether at least one intersection was found
  bool FoundIntersections() const;

  /// Whether at least one collision was found
  bool FoundCollisions() const;

  /// Number of self-intersections with specified cell
  ///
  /// \note Use only when both intersection tests and storage of intersection
  ///       details are enabled.
  int NumberOfIntersections(int) const;

  /// Number of collisions with specified cell
  ///
  /// \note Use only when both collision tests and storage of collision details
  /// are enabled.
  int NumberOfCollisions(int) const;

  /// Enable/disable intersection test for adjacent triangles
  mirtkOnOffMacro(AdjacentIntersectionTest);

  /// Enable/disable intersection test for non-adjacent triangles
  mirtkOnOffMacro(NonAdjacentIntersectionTest);

  /// Enable/disable near miss collision test for front-facing triangles
  mirtkOnOffMacro(FrontfaceCollisionTest);

  /// Enable/disable near miss collision test for back-facing triangles
  mirtkOnOffMacro(BackfaceCollisionTest);

  /// Enable/disable fast, approximate collision test
  mirtkOnOffMacro(FastCollisionTest);

  /// Enable/disable storage of details about found intersections
  mirtkOnOffMacro(StoreIntersectionDetails);

  /// Enable/disable storage of details about found intersections
  mirtkOnOffMacro(StoreCollisionDetails);

  /// Enable/disable resetting of optional input collision type array
  mirtkOnOffMacro(ResetCollisionType);

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
inline bool SurfaceCollisions::IsCollision(CollisionType type)
{
  return (Collision <= type && type < Intersection) || (type == Ambiguous);
}

// -----------------------------------------------------------------------------
inline bool SurfaceCollisions::IsIntersection(CollisionType type)
{
  return Intersection <= type && type <= Ambiguous;
}

// -----------------------------------------------------------------------------
inline const SurfaceCollisions::IntersectionsSet &
SurfaceCollisions::Intersections(int cellId) const
{
  return _Intersections[cellId];
}

// -----------------------------------------------------------------------------
inline const SurfaceCollisions::CollisionsSet &
SurfaceCollisions::Collisions(int cellId) const
{
  return _Collisions[cellId];
}

// -----------------------------------------------------------------------------
inline int SurfaceCollisions::NumberOfIntersections(int cellId) const
{
  return static_cast<int>(Intersections(cellId).size());
}

// -----------------------------------------------------------------------------
inline int SurfaceCollisions::NumberOfCollisions(int cellId) const
{
  return static_cast<int>(Collisions(cellId).size());
}


} // namespace mirtk

#endif // MIRKT_SurfaceCollisions_H
