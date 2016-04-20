/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2016 Imperial College London
 * Copyright 2016 Andreas Schuh
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

#ifndef MIRTK_Plane_H
#define MIRTK_Plane_H

#include "mirtk/Object.h"

#include "mirtk/Vector3.h"
#include "mirtk/Point.h"
#include "mirtk/PointSet.h"


namespace mirtk {


/**
 * Plane in 3D space.
 *
 * The plane is stored in Hessian normal form. Additionally, the two orthonormal
 * tangent vectors to the plane normal are stored which are used for projecting
 * points onto the plane and expressing these as 2D (u, v) coordinates along these
 * in-plane axes. The vectors Plane::Tangent1, Plane::Tangent2, and Plane::Normal
 * form a right-handed coordinate system.
 */
class Plane : public Object
{
  mirtkObjectMacro(Plane);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Distance of plane from origin
  mirtkPublicAttributeMacro(double, Offset);

  /// First in-plane axis
  mirtkReadOnlyAttributeMacro(Vector3, Tangent1);

  /// Second in-plane axis
  mirtkReadOnlyAttributeMacro(Vector3, Tangent2);

  /// Plane normal vector
  mirtkReadOnlyAttributeMacro(Vector3, Normal);

  /// Origin of plane coordinate system
  ///
  /// When a plane is fit to a set of points, the origin is the centroid of the
  /// point cloud. Otherwise, the origin is the point on the plane that is
  /// closest to the world origin, i.e., the vector pointing from plane origin
  /// to the world origin is orthogonal to the plane normal.
  mirtkReadOnlyAttributeMacro(Point, Origin);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const Plane &);

protected:

  /// Update origin after change of plane normal
  void UpdateOrigin();

  /// Update tangent vectors after change of plane normal
  void UpdateTangents();

  // ---------------------------------------------------------------------------
  // Construction/Destruction
public:

  /// Default constructor
  Plane();

  /// Construct plane given its Hessian normal form parameters
  ///
  /// \param[in] n Normal vector.
  /// \param[in] b Distance of plane to origin.
  Plane(const Vector3 &n, double b);

  /// Construct plane given its Hessian normal form parameters
  ///
  /// \param[in] n Normal vector.
  /// \param[in] b Distance of plane to origin.
  Plane(double n[3], double b);

  /// Construct plane given its Hessian normal form parameters
  ///
  /// \param[in] nx First  component of normal vector.
  /// \param[in] ny Second component of normal vector.
  /// \param[in] nz Third  component of normal vector.
  /// \param[in] b  Distance of plane to origin.
  Plane(double nx, double ny, double nz, double b);

  /// Copy constructor
  Plane(const Plane &);

  /// Assignment operator
  Plane &operator =(const Plane &);

  /// Destructor
  virtual ~Plane();

  /// Set normal vector
  void Normal(const Vector3 &n);

  /// Set normal vector
  ///
  /// \param[in] n Normal vector.
  void Normal(double n[3]);

  /// Set normal vector
  ///
  /// \param[in] nx First  component of normal vector.
  /// \param[in] ny Second component of normal vector.
  /// \param[in] nz Third  component of normal vector.
  void Normal(double nx, double ny, double nz);

  // ---------------------------------------------------------------------------
  // Model fitting
public:

  /// Fit plane to a set of points
  ///
  /// \param[in] points Point cloud.
  ///
  /// \returns RMS error of the regression.
  double Fit(const PointSet &points);

  // ---------------------------------------------------------------------------
  // Evaluation

  /// Evaluates plane equation for a given point
  ///
  /// \param[in] p Point.
  ///
  /// \returns Signed distance of point to plane.
  double Evaluate(const Point &p) const;

  /// Evaluates plane equation for a given point
  ///
  /// \param[in] p Point.
  ///
  /// \returns Signed distance of point to plane.
  double Evaluate(const double p[3]) const;

  /// Evaluates plane equation for a given point
  ///
  /// \param[in] x First  coordinate of point.
  /// \param[in] y Second coordinate of point.
  /// \param[in] z Third  coordinate of point.
  ///
  /// \returns Signed distance of point to plane.
  double Evaluate(double x, double y, double z) const;

  /// Evaluates plane equation for a given point
  ///
  /// \param[in] p Point.
  ///
  /// \returns Signed distance of point to plane.
  double Distance(const Point &p) const;

  /// Evaluates plane equation for a given point
  ///
  /// \param[in] p Point.
  ///
  /// \returns Signed distance of point to plane.
  double Distance(const double p[3]) const;

  /// Evaluates plane equation for a given point
  ///
  /// \param[in] x First  coordinate of point.
  /// \param[in] y Second coordinate of point.
  /// \param[in] z Third  coordinate of point.
  ///
  /// \returns Signed distance of point to plane.
  double Distance(double x, double y, double z) const;

  // ---------------------------------------------------------------------------
  // Projection

  /// Project point onto plane
  ///
  /// \param[in]  p Point.
  /// \param[out] u Coordinate of projected point along Tangent1.
  /// \param[out] v Coordinate of projected point along Tangent2.
  void Project(const Point &p, double &u, double &v) const;

  /// Project point onto plane
  ///
  /// \param[in]  p Point.
  /// \param[out] u Coordinate of projected point along Tangent1.
  /// \param[out] v Coordinate of projected point along Tangent2.
  void Project(const double p[3], double &u, double &v) const;

  /// Project point onto plane
  ///
  /// \param[in]  x First  coordinate of point.
  /// \param[in]  y Second coordinate of point.
  /// \param[in]  z Third  coordinate of point.
  /// \param[out] u Coordinate of projected point along Tangent1.
  /// \param[out] v Coordinate of projected point along Tangent2.
  void Project(double x, double y, double z, double &u, double &v) const;

  // ---------------------------------------------------------------------------
  // Debugging

  /// Print plane equation
  ostream &Print(ostream &os, Indent = 0) const;

  /// Print plane equation to standard output stream
  void Print(Indent = 0) const;

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
inline Plane::Plane()
:
  _Offset(.0),
  _Tangent1(1.0, .0, .0),
  _Tangent2(.0, 1.0, .0),
  _Normal  (.0, .0, 1.0),
  _Origin  (.0, .0, .0)
{
}

// -----------------------------------------------------------------------------
inline Plane::Plane(const Vector3 &n, double b)
:
  _Offset(b),
  _Normal(n)
{
  UpdateOrigin();
  UpdateTangents();
}

// -----------------------------------------------------------------------------
inline Plane::Plane(double n[3], double b)
:
  _Offset(b),
  _Normal(n[0], n[1], n[2])
{
  UpdateOrigin();
  UpdateTangents();
}

// -----------------------------------------------------------------------------
inline Plane::Plane(double nx, double ny, double nz, double b)
:
  _Offset(b),
  _Normal(nx, ny, nz)
{
  UpdateOrigin();
  UpdateTangents();
}

// -----------------------------------------------------------------------------
inline void Plane::CopyAttributes(const Plane &other)
{
  _Offset   = other._Offset;
  _Normal   = other._Normal;
  _Tangent1 = other._Tangent1;
  _Tangent2 = other._Tangent2;
  _Origin   = other._Origin;
}

// -----------------------------------------------------------------------------
inline Plane::Plane(const Plane &other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
inline Plane &Plane::operator =(const Plane &other)
{
  if (this != &other) {
    Object::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
inline Plane::~Plane()
{
}

// -----------------------------------------------------------------------------
inline void Plane::Normal(const Vector3 &n)
{
  _Normal = n;
  UpdateOrigin();
  UpdateTangents();
}

// -----------------------------------------------------------------------------
inline void Plane::Normal(double n[3])
{
  _Normal[0] = n[0];
  _Normal[1] = n[1];
  _Normal[2] = n[2];
  UpdateOrigin();
  UpdateTangents();
}

// -----------------------------------------------------------------------------
inline void Plane::Normal(double nx, double ny, double nz)
{
  _Normal[0] = nx;
  _Normal[1] = ny;
  _Normal[2] = nz;
  UpdateOrigin();
  UpdateTangents();
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
inline double Plane::Evaluate(const Point &p) const
{
  return p._x * _Normal[0] + p._y * _Normal[1] + p._z * _Normal[2] - _Offset;
}

// -----------------------------------------------------------------------------
inline double Plane::Evaluate(const double p[3]) const
{
  return p[0] * _Normal[0] + p[1] * _Normal[1] + p[2] * _Normal[2] - _Offset;
}

// -----------------------------------------------------------------------------
inline double Plane::Evaluate(double x, double y, double z) const
{
  return x * _Normal[0] + y * _Normal[1] + z * _Normal[2] - _Offset;
}

// -----------------------------------------------------------------------------
inline double Plane::Distance(const Point &p) const
{
  return Evaluate(p);
}

// -----------------------------------------------------------------------------
inline double Plane::Distance(const double p[3]) const
{
  return Evaluate(p);
}

// -----------------------------------------------------------------------------
inline double Plane::Distance(double x, double y, double z) const
{
  return Evaluate(x, y, z);
}

// =============================================================================
// Projection
// =============================================================================

// -----------------------------------------------------------------------------
inline void Plane::Project(const Point &p, double &u, double &v) const
{
  Vector3 d(p._x - _Origin._x, p._y - _Origin._y, p._z - _Origin._z);
  u = _Tangent1.Dot(d);
  v = _Tangent2.Dot(d);
}

// -----------------------------------------------------------------------------
inline void Plane::Project(const double p[3], double &u, double &v) const
{
  Vector3 d(p[0] - _Origin._x, p[1] - _Origin._y, p[2] - _Origin._z);
  u = _Tangent1.Dot(d);
  v = _Tangent2.Dot(d);
}

// -----------------------------------------------------------------------------
inline void Plane::Project(double x, double y, double z, double &u, double &v) const
{
  Vector3 d(x - _Origin._x, y - _Origin._y, z - _Origin._z);
  u = _Tangent1.Dot(d);
  v = _Tangent2.Dot(d);
}

// =============================================================================
// Debugging
// =============================================================================

// -----------------------------------------------------------------------------
inline void Plane::Print(Indent indent) const
{
  Print(cout, indent);
}


} // namespace mirtk

#endif // MIRTK_Plane_H
