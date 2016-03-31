/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2015 Imperial College London
 * Copyright 2013-2015 Andreas Schuh
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

#ifndef MIRTK_Polyhedron_H
#define MIRTK_Polyhedron_H

#include "mirtk/Object.h"
#include "mirtk/Point.h"

#include "vtkSmartPointer.h"
#include "vtkPolyData.h"


namespace mirtk {


/**
 * Utility functions for dealing with VTK polydata representing a polyhedron.
 *
 * \note The point-in-polyhedron test based on the winding number method and
 * the polyhedron volume computation are a C++ implementation of the excellent
 * polyhedron.py Python module of Mark Dickinson found on GitHub at
 * https://github.com/mdickinson/polyhedron/blob/cd7361bcee8cbd9ef2e0aac58a7ce59bf9a52c4f/polyhedron.py
 */
class Polyhedron : public Object
{
  mirtkObjectMacro(Polyhedron);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Dataset defining the geometry and topology of this polyhedron
  mirtkPublicAttributeMacro(vtkSmartPointer<vtkPolyData>, DataSet);

public:

  // ---------------------------------------------------------------------------
  // Construction/destruction

  /// Constructor
  Polyhedron(vtkPolyData * = NULL);

  /// Copy constructor
  Polyhedron(const Polyhedron &);

  /// Assignment operator
  Polyhedron &operator =(const Polyhedron &);

  /// Destructor
  virtual ~Polyhedron();

  // ---------------------------------------------------------------------------
  // Geometry / Topology

  /// Number of vertices
  int NumberOfPoints() const;

  /// Get vertex position
  void GetPoint(int, double &, double &, double &) const;

  /// Get vertex position
  void GetPoint(int, double [3]) const;

  /// Get vertex position
  Point GetPoint(int) const;

  // ---------------------------------------------------------------------------
  // Properties

  /// Calculate volume enclosed by polyhedron
  double Volume() const;

  /// Calculate volume enclosed by polyhedron
  static double Volume(vtkPolyData *);

  /// Compute winding number of polyhedron around given point
  int WindingNumber(double, double, double) const;

  /// Compute winding number of polyhedron around given point
  int WindingNumber(double [3]) const;

  /// Compute winding number of polyhedron around given point
  int WindingNumber(const Point &) const;

  /// Compute winding number of polyhedron around given point
  static int WindingNumber(vtkPolyData *, double, double, double);

  /// Compute winding number of polyhedron around given point
  static int WindingNumber(vtkPolyData *, double [3]);

  /// Compute winding number of polyhedron around given point
  static int WindingNumber(vtkPolyData *, const Point &);

  // ---------------------------------------------------------------------------
  // Point-in-polyhedron test

  /// Test whether point is inside the polyhedron
  bool IsInside(double, double, double) const;

  /// Test whether point is inside the polyhedron
  bool IsInside(double [3]) const;

  /// Test whether point is inside the polyhedron
  bool IsInside(const Point &) const;

  /// Test whether point is inside the polyhedron
  static bool IsInside(vtkPolyData *, double, double, double);

  /// Test whether point is inside the polyhedron
  static bool IsInside(vtkPolyData *, double [3]);

  /// Test whether point is inside the polyhedron
  static bool IsInside(vtkPolyData *, const Point &);

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Geometry / Topology
// =============================================================================

// -----------------------------------------------------------------------------
inline int Polyhedron::NumberOfPoints() const
{
  return static_cast<int>(_DataSet->GetNumberOfPoints());
}

// -----------------------------------------------------------------------------
inline void Polyhedron::GetPoint(int i, double p[3]) const
{
  _DataSet->GetPoint(i, p);
}

// -----------------------------------------------------------------------------
inline void Polyhedron::GetPoint(int i, double &x, double &y, double &z) const
{
  double p[3];
  _DataSet->GetPoint(i, p);
  x = p[0], y = p[1], z = p[2];
}

// -----------------------------------------------------------------------------
inline Point Polyhedron::GetPoint(int i) const
{
  Point p;
  GetPoint(i, p._x, p._y, p._z);
  return p;
}

// =============================================================================
// Properties
// =============================================================================

// -----------------------------------------------------------------------------
inline double Polyhedron::Volume() const
{
  return Volume(_DataSet);
}

// -----------------------------------------------------------------------------
inline int Polyhedron::WindingNumber(double x, double y, double z) const
{
  return WindingNumber(_DataSet, x, y, z);
}

// -----------------------------------------------------------------------------
inline int Polyhedron::WindingNumber(double p[3]) const
{
  return WindingNumber(p[0], p[1], p[2]);
}

// -----------------------------------------------------------------------------
inline int Polyhedron::WindingNumber(const Point &p) const
{
  return WindingNumber(p._x, p._y, p._z);
}

// -----------------------------------------------------------------------------
inline int Polyhedron::WindingNumber(vtkPolyData *polydata, double p[3])
{
  return WindingNumber(polydata, p[0], p[1], p[2]);
}

// -----------------------------------------------------------------------------
inline int Polyhedron::WindingNumber(vtkPolyData *polydata, const Point &p)
{
  return WindingNumber(polydata, p._x, p._y, p._z);
}

// =============================================================================
// Point-in-polyhedron test
// =============================================================================

// -----------------------------------------------------------------------------
inline bool Polyhedron::IsInside(vtkPolyData *polydata, double x, double y, double z)
{
  return WindingNumber(polydata, x, y, z) != 0;
}

// -----------------------------------------------------------------------------
inline bool Polyhedron::IsInside(double x, double y, double z) const
{
  return IsInside(_DataSet, x, y, z);
}

// -----------------------------------------------------------------------------
inline bool Polyhedron::IsInside(double p[3]) const
{
  return IsInside(p[0], p[1], p[2]);
}

// -----------------------------------------------------------------------------
inline bool Polyhedron::IsInside(const Point &p) const
{
  return IsInside(p._x, p._y, p._z);
}

// -----------------------------------------------------------------------------
inline bool Polyhedron::IsInside(vtkPolyData *polydata, double p[3])
{
  return IsInside(polydata, p[0], p[1], p[2]);
}

// -----------------------------------------------------------------------------
inline bool Polyhedron::IsInside(vtkPolyData *polydata, const Point &p)
{
  return IsInside(polydata, p._x, p._y, p._z);
}


} // namespace mirtk

#endif // MIRTK_Polyhedron_H
