/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2015 Imperial College London
 * Copyright 2008-2015 Daniel Rueckert, Julia Schnabel
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

#ifndef MIRTK_PointSet_H
#define MIRTK_PointSet_H

#include "mirtk/Object.h"

#include "mirtk/Memory.h"
#include "mirtk/Cfstream.h"
#include "mirtk/Point.h"
#include "mirtk/Array.h"
#include "mirtk/OrderedSet.h"

#include "mirtk/NumericsExport.h"


// Forward declaration of VTK type
class vtkAbstractArray;


namespace mirtk {


/**
 * Point set
 */
class PointSet : public Object
{
  mirtkObjectMacro(PointSet);

  // ---------------------------------------------------------------------------
  // Attributes

protected:

  /// Allocated size of PointSet
  int _m;

  /// Actual size of PointSet
  int _n;

  /// Pointer to Points
  Point *_data;

  MIRTK_Numerics_EXPORT static int POINTSET_SIZE;

public:

  // ---------------------------------------------------------------------------
  // Constructors and destructor

  /// Default constructor
  PointSet(int = 0);

  /// Copy constructor
  PointSet(const PointSet &);

  /// Copy subset constructor
  PointSet(const PointSet &, const Array<int> &);

  /// Copy subset constructor
  PointSet(const PointSet &, const OrderedSet<int> &);

  /// Destructor
  virtual ~PointSet();

  // ---------------------------------------------------------------------------
  // Size get acccessor and clearing of PointSet

  /// Resizes container so it contains n elements
  void Resize(int, const Point & = Point());

  /// Request a change in capacity
  void Reserve(int);

  /// Return size of allocated storage capacity
  int Capacity() const;

  /// Set container size (and capacity)
  ///
  /// The result of this function is equivalent to
  /// \code
  /// Resize(n);
  /// ShrinkToFit();
  /// \endcode
  /// but with only one reallocation operation. It sets both the size and
  /// the capacity of the point set container.
  void Size(int);

  /// Access function for size
  int Size() const;

  /// Requests the container to reduce its capacity to fit its size
  void ShrinkToFit();

  /// Clearing of PointSet
  ///
  /// \param[in] deallocate Whether to deallocate memory.
  void Clear(bool deallocate = true);

  // ---------------------------------------------------------------------------
  // Operators for access

  /// Get n-th point
  Point &GetPoint(int);

  /// Get n-th point
  const Point &GetPoint(int) const;

  /// Get n-th point
  void GetPoint(int, Point &) const;

  /// Get n-th point
  void GetPoint(int, double [3]) const;

  /// Set n-th point
  void SetPoint(int, const Point &);

  /// Set n-th point
  void SetPoint(int, const double [3]);

  /// Operator for Point put access
  Point &operator()(int);

  /// Operator for Point get access
  const Point &operator()(int) const;

  /// Operator
  PointSet operator()(int, int) const;

  // ---------------------------------------------------------------------------
  // Operators for PointSet arithmetic

  // PointSet operation for =
  PointSet &operator=(const PointSet &);

  // PointSet operator for Point adding
  PointSet& operator+=(const Point &);

  // PointSet operator for Point substraction
  PointSet& operator-=(const Point &);

  // PointSet operator for PointSet adding
  PointSet& operator+=(const PointSet &);

  // PointSet operator for PointSet substraction
  PointSet& operator-=(const PointSet &);

  /// Centre of gravity
  virtual Point CenterOfGravity() const;

  /// Closest point to given point
  virtual Point ClosestPoint(Point&);

  /// Point set distance to given point
  virtual double PointDistance(Point&);

  /// Bounding box
  virtual void BoundingBox(Point &, Point &) const;

  // ---------------------------------------------------------------------------
  // Operators for PointSet comparison

  /// Test for equality
  bool operator ==(const PointSet &) const;

  /// Test for inequality
  bool operator !=(const PointSet &) const;

  // ---------------------------------------------------------------------------
  // Explicit methods to add or delete points

  /// Adding of a point to point set
  void Add(const Point &);

  /// Deleting of a point from point set 
  void Del(const Point &);

  /// Adding of a PointSet to Pointset
  void Add(const PointSet &);

  /// Deleting of a PointSet from Pointset
  void Del(const PointSet &);

  /// Adding of a Point to Pointset
  void Add(double *);

  /// Deleting of a Point from Pointset
  void Del(double *);

  /// Delete all points without freeing already allocated memory
  void Del();

  // ---------------------------------------------------------------------------
  // I/O

  /// Read pointset from file
  void Read(const char *);

  /// Write pointset to file
  void Write(const char *) const;

  /// Read point set from file in VTK format
  ///
  /// \note Use this function only when MIRTK_Numerics_WITH_VTK is 1.
  void ReadVTK(const char *);

  /// Add point set from file in VTK format
  ///
  /// \note Use this function only when MIRTK_Numerics_WITH_VTK is 1.
  void AddVTK(const char *);

  /// Write point set to file in VTK format
  ///
  /// \note Use this function only when MIRTK_Numerics_WITH_VTK is 1.
  void WriteVTK(const char *, vtkAbstractArray * = NULL) const;

  // ---------------------------------------------------------------------------
  // Misc. functions

  /// Computes the standard deviation ellipsoid about the centroid of a point set.
  Point StandardDeviationEllipsoid() const;

  /// Tests if a point is inside the polygon defined by the point set
  int IsInside(double, double) const;

  /// Get (unweighted) centroid of point set
  Point Centroid() const;
};

////////////////////////////////////////////////////////////////////////////////
// Streaming operators
////////////////////////////////////////////////////////////////////////////////

/// Write point set to binary output file stream
Cofstream &operator <<(Cofstream &, const PointSet &);

/// Read point set from binary input file stream
Cifstream &operator >>(Cifstream &, PointSet &);

/// Write point set to text output stream
ostream &operator <<(ostream &, const PointSet &);

/// Read point set from text input stream
istream &operator >>(istream &, PointSet &);

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
inline PointSet::PointSet(int n)
:
  _m(n), _n(n), _data(Allocate<Point>(n))
{
}

// -----------------------------------------------------------------------------
inline PointSet::PointSet(const PointSet &pset)
:
  Object(pset),
  _m(0), _n(0), _data(NULL)
{
  (*this) = pset;
}

// -----------------------------------------------------------------------------
inline PointSet::PointSet(const PointSet &pset, const Array<int> &subset)
:
  Object(pset),
  _m(0), _n(0), _data(NULL)
{
  int i = 0;
  Size(static_cast<int>(subset.size()));
  for (Array<int>::const_iterator it = subset.begin(); it != subset.end(); ++it, ++i) {
    _data[i] = pset(*it);
  }
}

// -----------------------------------------------------------------------------
inline PointSet::PointSet(const PointSet &pset, const OrderedSet<int> &subset)
:
  Object(pset),
  _m(0), _n(0), _data(NULL)
{
  int i = 0;
  Size(static_cast<int>(subset.size()));
  for (OrderedSet<int>::const_iterator it = subset.begin(); it != subset.end(); ++it, ++i) {
    _data[i] = pset(*it);
  }
}

// -----------------------------------------------------------------------------
inline PointSet::~PointSet()
{
  Clear();
}

// -----------------------------------------------------------------------------
inline Point& PointSet::operator()(int j)
{
  return _data[j];
}

// -----------------------------------------------------------------------------
inline const Point &PointSet::operator()(int j) const
{
  return _data[j];
}

// -----------------------------------------------------------------------------
inline Point &PointSet::GetPoint(int i)
{
  return this->operator()(i);
}

// -----------------------------------------------------------------------------
inline const Point &PointSet::GetPoint(int i) const
{
  return this->operator()(i);
}

// -----------------------------------------------------------------------------
inline void PointSet::GetPoint(int i, Point &p) const
{
  p = this->operator()(i);
}

// -----------------------------------------------------------------------------
inline void PointSet::GetPoint(int i, double p[3]) const
{
  const Point &pt = this->operator()(i);
  p[0] = pt._x, p[1] = pt._y, p[2] = pt._z;
}

// -----------------------------------------------------------------------------
inline void PointSet::SetPoint(int i, const Point &p)
{
  this->operator()(i) = p;
}

// -----------------------------------------------------------------------------
inline void PointSet::SetPoint(int i, const double p[3])
{
  Point &pt = this->operator()(i);
  pt._x = p[0], pt._y = p[1], pt._z = p[2];
}

// -----------------------------------------------------------------------------
inline PointSet& PointSet::operator+=(const Point &p)
{
  this->Add(p);
  return *this;
}

// -----------------------------------------------------------------------------
inline PointSet& PointSet::operator-=(const Point &p)
{
  this->Del(p);
  return *this;
}

// -----------------------------------------------------------------------------
inline PointSet& PointSet::operator+=(const PointSet &pset)
{
  this->Add(pset);
  return *this;
}

// -----------------------------------------------------------------------------
inline PointSet& PointSet::operator-=(const PointSet &pset)
{
  this->Del(pset);
  return *this;
}

// -----------------------------------------------------------------------------
inline PointSet PointSet::operator()(int j, int k) const
{
  PointSet pset;
  for (int i = j; i < k; j++) {
    pset += _data[i];
  }
  return pset;
}

// -----------------------------------------------------------------------------
inline void PointSet::Resize(int n, const Point &p)
{
  int m = _m;
  while (n > m) m += POINTSET_SIZE;
  this->Reserve(m);
  _n = n;
  for (int i = _n; i < _m; ++i) _data[i] = p;
}

// -----------------------------------------------------------------------------
inline void PointSet::Reserve(int m)
{
  if (_m < m) {
    _m = m;
    Point *new_data = Allocate<Point>(m);
    for (int i = 0; i < _n; ++i) new_data[i] = _data[i];
    Deallocate(_data);
    _data = new_data;
  }
}

// -----------------------------------------------------------------------------
inline int PointSet::Capacity() const
{
  return _m;
}

// -----------------------------------------------------------------------------
inline void PointSet::Size(int n)
{
  if (_m != n) {
    Point *new_data = Allocate<Point>(n);
    for (int i = 0; i < _n; ++i) new_data[i] = _data[i];
    Deallocate(_data);
    _data = new_data;
    _m = _n = n;
  }
}

// -----------------------------------------------------------------------------
inline int PointSet::Size() const
{
  return _n;
}

// -----------------------------------------------------------------------------
inline void PointSet::ShrinkToFit()
{
  if (_n < _m) {
    Point *new_data = Allocate<Point>(_n);
    for (int i = 0; i < _n; ++i) new_data[i] = _data[i];
    Deallocate(_data);
    _data = new_data;
    _m    = _n;
  }
}

// -----------------------------------------------------------------------------
inline void PointSet::Add(double *p)
{
  this->Add(Point(p[0], p[1], p[2]));
}

// -----------------------------------------------------------------------------
inline void PointSet::Del(double *p)
{
  this->Del(Point(p[0], p[1], p[2]));
}

// -----------------------------------------------------------------------------
inline void PointSet::Del()
{
  this->Clear(false);
}


} // namespace mirtk

#endif // MIRTK_PointSet_H

