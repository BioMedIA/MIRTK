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

#ifndef MIRTK_RegisteredSurface_H
#define MIRTK_RegisteredSurface_H

#include "mirtk/RegisteredPointSet.h"

#include "vtkPolyData.h"
#include "vtkCellArray.h"


namespace mirtk {


/**
 * Registered boundary surface (vtkPolyData)
 */
class RegisteredSurface : public RegisteredPointSet
{
  mirtkObjectMacro(RegisteredSurface);

public:

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Constructor
  RegisteredSurface(vtkPolyData * = NULL, const class Transformation * = NULL);

  /// Copy constructor
  RegisteredSurface(const RegisteredSurface &);

  /// Assignment operator
  RegisteredSurface &operator =(const RegisteredSurface &);

  /// Destructor
  ~RegisteredSurface();

  /// Initialize dataset after input and parameters are set
  void Initialize();

  // ---------------------------------------------------------------------------
  // Polydata access

  /// Set input surface
  void InputSurface(vtkPolyData *);

  /// Get number of vertices
  int NumberOfVerts() const;

  /// Get number of lines
  int NumberOfLines() const;

  /// Get number of polygons
  int NumberOfPolys() const;

  /// Get number of triangle strips
  int NumberOfStrips() const;

  /// Get (transformed) polydata
  vtkPolyData *PolyData() const;

  /// Get the cell array defining vertices. If there are no vertices, an
  /// empty array will be returned (convenience to simplify traversal).
  vtkCellArray *Verts() const;

  /// Get the cell array defining lines. If there are no lines, an
  /// empty array will be returned (convenience to simplify traversal).
  vtkCellArray *Lines() const;

  /// Get the cell array defining polygons. If there are no polygons, an
  /// empty array will be returned (convenience to simplify traversal).
  vtkCellArray *Polys() const;

  /// Get the cell array defining triangle strips. If there are no triangle strips,
  /// an empty array will be returned (convenience to simplify traversal).
  vtkCellArray *Strips() const;

  /// Get pointer to IDs of points defining a cell
  void GetCellPoints(int, vtkIdType &, const vtkIdType *&) const;

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Polydata access
// =============================================================================

// -----------------------------------------------------------------------------
inline void RegisteredSurface::InputSurface(vtkPolyData *surface)
{
  _InputPointSet = surface;
}

// -----------------------------------------------------------------------------
inline vtkCellArray *RegisteredSurface::Verts() const
{
  return _OutputSurface->GetVerts();
}

// -----------------------------------------------------------------------------
inline vtkCellArray *RegisteredSurface::Lines() const
{
  return _OutputSurface->GetLines();
}

// -----------------------------------------------------------------------------
inline vtkCellArray *RegisteredSurface::Polys() const
{
  return _OutputSurface->GetPolys();
}

// -----------------------------------------------------------------------------
inline vtkCellArray *RegisteredSurface::Strips() const
{
  return _OutputSurface->GetStrips();
}

// -----------------------------------------------------------------------------
inline int RegisteredSurface::NumberOfVerts() const
{
  return static_cast<int>(Verts()->GetNumberOfCells());
}

// -----------------------------------------------------------------------------
inline int RegisteredSurface::NumberOfLines() const
{
  return static_cast<int>(Lines()->GetNumberOfCells());
}

// -----------------------------------------------------------------------------
inline int RegisteredSurface::NumberOfPolys() const
{
  return static_cast<int>(Polys()->GetNumberOfCells());
}

// -----------------------------------------------------------------------------
inline int RegisteredSurface::NumberOfStrips() const
{
  return static_cast<int>(Strips()->GetNumberOfCells());
}

// -----------------------------------------------------------------------------
inline void RegisteredSurface::GetCellPoints(int i, vtkIdType &npts, const vtkIdType *&pts) const
{
  _OutputSurface->GetCellPoints(i, npts, const_cast<vtkIdType *&>(pts));
}


} // namespace mirtk

#endif // MIRTK_RegisteredSurface_H
