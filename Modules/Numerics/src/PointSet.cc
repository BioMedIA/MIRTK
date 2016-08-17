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

#include "mirtk/PointSet.h"

#include "mirtk/Math.h"
#include "mirtk/Path.h"

#include "mirtk/NumericsConfig.h"
#if MIRTK_Numerics_WITH_VTK
#  include "mirtk/Vtk.h"
#  include "vtkSmartPointer.h"
#  include "vtkPoints.h"
#  include "vtkCellArray.h"
#  include "vtkPointSet.h"
#  include "vtkPointData.h"
#  include "vtkPolyDataReader.h"
#  include "vtkXMLGenericDataObjectReader.h"
#  include "vtkOBJReader.h"
#  include "vtkPLYReader.h"
#  include "vtkSTLReader.h"
#  include "vtkPolyDataWriter.h"
#  include "vtkXMLPolyDataWriter.h"
#endif


namespace mirtk {


// -----------------------------------------------------------------------------
MIRTK_Numerics_EXPORT int PointSet::POINTSET_SIZE = 4096;

// -----------------------------------------------------------------------------
void PointSet::Clear(bool deallocate)
{
  if (deallocate) Deallocate(_data), _m = 0;
  _n = 0;
}

// -----------------------------------------------------------------------------
void PointSet::Add(const Point &p)
{
  if (_n + 1 <= _m) {
    // There is still enough memory left, so just add the point
    _data[_n] = p;
    _n++;
    return;
  }
  // There is not enough memory left, so allocate new point list and copy
  _m += POINTSET_SIZE;
  Point *new_data = Allocate<Point>(_m);
  for (int i = 0; i < _n; i++) {
    new_data[i] = _data[i];
  }
  new_data[_n] = p;
  Deallocate(_data);
  _data = new_data;
  _n++;
}

// -----------------------------------------------------------------------------
void PointSet::Del(const Point &p)
{
  Point *new_data = Allocate<Point>(_m);

  int new_n = 0;
  for (int i = 0; i < _n; ++i) {
    if (_data[i] == p) continue;
    new_data[new_n] = _data[i];
    ++new_n;
  }

  Deallocate(_data);
  _data = new_data;
  _n    = new_n;

  if ((new_n % POINTSET_SIZE) == 0) {
    _m = (new_n / POINTSET_SIZE) * POINTSET_SIZE;
    new_data = Allocate<Point>(_m);
    for (int i = 0; i < _n; ++i) {
      new_data[i] = _data[i];
    }
    Deallocate(_data);
    _data = new_data;
  }
}

// -----------------------------------------------------------------------------
void PointSet::Add(const PointSet &pset)
{
  for (int i = 0; i < pset.Size(); i++) this->Add(pset(i));
}

// -----------------------------------------------------------------------------
void PointSet::Del(const PointSet &pset)
{
  for (int i = 0; i < pset.Size(); i++) this->Del(pset(i));
}

// -----------------------------------------------------------------------------
PointSet& PointSet::operator=(const PointSet &pset)
{
  Reserve(pset._n);
  _n = pset._n;
  for (int i = 0; i < _n; ++i) {
    _data[i] = pset._data[i];
  }
  return *this;
}

// -----------------------------------------------------------------------------
Point PointSet::CenterOfGravity() const
{
  int i;
  Point p;

  if (this->Size() == 0) {
    cerr << "PointSet::CenterOfGravity(): No points in point set" << endl;
    return p;
  }
  for (i = 0; i < this->Size(); i++) {
    p += (*this)(i);
  }
  return p / (double)this->Size();
}

// -----------------------------------------------------------------------------
void PointSet::BoundingBox(Point &p1, Point &p2) const
{
  int i;

  if (this->Size() == 0) {
    p1 = Point();
    p2 = Point();
    return;
  } else {
    p1 = (*this)(0);
    p2 = (*this)(0);
  }
  for (i = 1; i < this->Size(); i++) {
    p1._x = p1._x < (*this)(i)._x ? p1._x : (*this)(i)._x;
    p1._y = p1._y < (*this)(i)._y ? p1._y : (*this)(i)._y;
    p1._z = p1._z < (*this)(i)._z ? p1._z : (*this)(i)._z;
    p2._x = p2._x > (*this)(i)._x ? p2._x : (*this)(i)._x;
    p2._y = p2._y > (*this)(i)._y ? p2._y : (*this)(i)._y;
    p2._z = p2._z > (*this)(i)._z ? p2._z : (*this)(i)._z;
  }
}

// -----------------------------------------------------------------------------
bool PointSet::operator ==(const PointSet &rhs) const
{
  if (_n != rhs._n) return false;
  for (int i = 0; i < _n; ++i) {
    if (_data[i] != rhs._data[i]) return false;
  }
  return true;
}

// -----------------------------------------------------------------------------
bool PointSet::operator !=(const PointSet &rhs) const
{
  return !(*this == rhs);
}

// -----------------------------------------------------------------------------
Cofstream &operator <<(Cofstream &os, const PointSet &pset)
{
  os.WriteAsChar("PointSet", 13);
  os.WriteAsInt(pset.Size());
  for (int i = 0; i < pset.Size(); ++i) {
    const Point &p = pset(i);
    os.WriteAsDouble(p._x);
    os.WriteAsDouble(p._y);
    os.WriteAsDouble(p._z);
  }
  return os;
}

// -----------------------------------------------------------------------------
Cifstream &operator >>(Cifstream &is, PointSet &pset)
{
  char buffer[13];
  is.ReadAsChar(buffer, 13);
  if (strncmp(buffer, "PointSet", 13) != 0) {
    cerr << "Failed to read point set from binary file stream" << endl;
    exit(1);
  }

  int n;
  is.ReadAsInt(&n, 1);
  pset.Reserve(pset.Size() + n);

  Point p;
  for (int i = 0; i < n; ++i) {
    is.ReadAsDouble(&p._x, 1);
    is.ReadAsDouble(&p._y, 1);
    is.ReadAsDouble(&p._z, 1);
    pset.Add(p);
  }

  return is;
}

// -----------------------------------------------------------------------------
ostream &operator <<(ostream &os, const PointSet &pset)
{
  os << "PointSet " << pset.Size() << endl;
  os.setf(ios::right);
  os.setf(ios::fixed);
  os.precision(10);
  for (int i = 0; i < pset.Size(); ++i) {
    os << setw(15) << pset(i) << endl;
  }
  os.precision(6);
  os.unsetf(ios::right);
  os.unsetf(ios::fixed);
  return os;
}

// -----------------------------------------------------------------------------
istream &operator >>(istream &is, PointSet &pset)
{
  // Read header
  char buffer[256];
  is >> buffer;
  if (strcmp(buffer, "PointSet") != 0) {
    cerr << "Can't read PointSet file: " << buffer << endl;
    exit(1);
  }

  // Read size
  int n;
  is >> n;
  pset.Reserve(pset.Size() + n);

  // Read PointSet
  Point p;
  for (int i = 0; i < n; ++i) {
    is >> p;
    pset.Add(p);
  }
  return is;
}

// -----------------------------------------------------------------------------
void PointSet::Read(const char *filename)
{
  // Read VTK file if extension matches known format supported by VTK
  const string ext = Extension(filename);
  if (ext == ".vtk" ||
     (ext.length() == 4 && ext.substr(0, 3) == ".vt") ||
      ext == ".obj" || ext == ".ply" || ext == ".stl") {
#if MIRTK_Numerics_WITH_VTK
    ReadVTK(filename);
#else
    cerr << "PointSet::Read: Cannot read input file " << filename << " when Numerics module not built WITH_VTK" << endl;
    exit(1);
#endif // MIRTK_Numerics_WITH_VTK
    return;
  }

  // Clear pointset first
  this->Clear();

  // Open file stream
  ifstream from(filename);

  // Check whether file opened ok
  if (!from) {
    cerr << "PointSet::Read: Can't open file " << filename << endl;
    exit(1);
  }

  // Read PointSet
  from >> *this;
}

// -----------------------------------------------------------------------------
void PointSet::Write(const char *filename) const
{
  // Open file stream
  ofstream to(filename);

  // Check whether file opened ok
  if (!to) {
    cerr << "PointSet::Write: Can't open file " << filename << endl;
    exit(1);
  }

  // Write PointSet
  to << *this;
}

#if MIRTK_Numerics_WITH_VTK

// -----------------------------------------------------------------------------
void PointSet::AddVTK(const char *fname)
{
  vtkSmartPointer<vtkPointSet> pointset;
  const string ext = Extension(fname);
  if (ext == ".obj") {
    vtkSmartPointer<vtkOBJReader> reader;
    reader = vtkSmartPointer<vtkOBJReader>::New();
    reader->SetFileName(fname);
    reader->Update();
    pointset = reader->GetOutput();
  } else if (ext == ".ply") {
    vtkSmartPointer<vtkPLYReader> reader;
    reader = vtkSmartPointer<vtkPLYReader>::New();
    reader->SetFileName(fname);
    reader->Update();
    pointset = reader->GetOutput();
  } else if (ext == ".stl") {
    vtkSmartPointer<vtkSTLReader> reader;
    reader = vtkSmartPointer<vtkSTLReader>::New();
    reader->SetFileName(fname);
    reader->Update();
    pointset = reader->GetOutput();
  } else if (ext.length() == 4 && ext.substr(0, 3) == ".vt" && ext[3] != 'k') {
    vtkSmartPointer<vtkXMLGenericDataObjectReader> reader;
    reader = vtkSmartPointer<vtkXMLGenericDataObjectReader>::New();
    reader->SetFileName(fname);
    reader->Update();
    pointset = vtkPointSet::SafeDownCast(reader->GetOutput());
  } else {
    vtkSmartPointer<vtkPolyDataReader> reader;
    reader = vtkSmartPointer<vtkPolyDataReader>::New();
    reader->SetFileName(fname);
    reader->Update();
    pointset = reader->GetOutput();
  }

  double p[3];
  vtkPoints * const points = pointset->GetPoints();
  this->Reserve(this->Size() + points->GetNumberOfPoints());
  for (vtkIdType i = 0; i < points->GetNumberOfPoints(); i++) {
    points->GetPoint(i, p);
    this->Add(Point(p));
  }
}

void PointSet::ReadVTK(const char *fname)
{
  this->Clear();
  AddVTK(fname);
}

void PointSet::WriteVTK(const char *fname, vtkAbstractArray *data) const
{
  const vtkIdType npoints = static_cast<vtkIdType>(this->Size());

  vtkSmartPointer<vtkPoints>    points;
  vtkSmartPointer<vtkCellArray> vertices;
  vtkSmartPointer<vtkPolyData>  output;

  points = vtkSmartPointer<vtkPoints>::New();
  points->SetNumberOfPoints(npoints);

  vertices = vtkSmartPointer<vtkCellArray>::New();
  vertices->Allocate(npoints);

  double p[3];
  for (vtkIdType i = 0; i < npoints; ++i) {
    this->GetPoint(i, p);
    points->InsertPoint(i, p);
    vertices->InsertNextCell(1, &i);
  }

  output = vtkSmartPointer<vtkPolyData>::New();
  output->SetPoints(points);
  output->SetVerts(vertices);
  if (data) output->GetPointData()->AddArray(data);

  if (Extension(fname) == ".vtp") {
    vtkSmartPointer<vtkXMLPolyDataWriter> writer;
    writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(fname);
    SetVTKInput(writer, output);
    writer->Update();
  } else {
    vtkSmartPointer<vtkPolyDataWriter> writer;
    writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetFileName(fname);
    SetVTKInput(writer, output);
    writer->Update();
  }
}

#endif // MIRTK_Numerics_WITH_VTK

// -----------------------------------------------------------------------------
Point PointSet::ClosestPoint(Point &p_input)
{
  int i,j;
  double mindistance = 100000;
  double tmpdistance;
  Point e;

  if (this->Size() == 0) {
    cerr << "PointSet::ClosestPoint(): No points in pointset" << endl;
    return e;
  }
  j = 0;
  for (i = 0; i < this->Size(); i++) {
    tmpdistance = _data[i].Distance(p_input);
    if (tmpdistance < mindistance) {
      mindistance = tmpdistance;
      j = i;
    }
  }
  e = _data[j];
  return e;
}

// -----------------------------------------------------------------------------
double PointSet::PointDistance(Point &p_input)
{
  int i;
  double mindistance = 100000;
  double tmpdistance;
  if (this->Size() == 0) {
    cerr << "PointSet::PointDistant(): No points in pointset" << endl;
    return 0;
  }
  for (i = 0; i < this->Size(); i++) {
    tmpdistance = _data[i].Distance(p_input);
    if(tmpdistance < mindistance){
      mindistance = tmpdistance;
    }
  }
  return mindistance ;
}

// -----------------------------------------------------------------------------
Point PointSet::StandardDeviationEllipsoid() const
{
  int i;
  Point p;
  Point e;

  if (this->Size() == 0) {
    cerr << "PointSet::StandardDeviationEllipsoid(): No points in pointset" << endl;
    return e;
  }
  p = this->CenterOfGravity();

  for (i = 0; i < this->Size(); i++) {
    e._x += ((*this)(i)._x-p._x)*((*this)(i)._x-p._x);
    e._y += ((*this)(i)._y-p._y)*((*this)(i)._y-p._y);
    e._z += ((*this)(i)._z-p._z)*((*this)(i)._z-p._z);
  }
  e._x = sqrt(e._x / ((double)this->Size()-1));
  e._y = sqrt(e._y / ((double)this->Size()-1));
  e._z = sqrt(e._z / ((double)this->Size()-1));
  return e ;
}

// -----------------------------------------------------------------------------
static inline int intersection(double x1, double y1, double x2, double y2,
                               double x3, double y3, double x4, double y4)
{
  double a = (x4 - x3)*(y1 - y3) - (y4 - y3)*(x1 - x3);
  double b = (x2 - x1)*(y1 - y3) - (y2 - y1)*(x1 - x3);
  if ((y4 - y3)*(x2 - x1) - (x4 - x3)*(y2 - y1) != 0) {
    a /= (y4 - y3)*(x2 - x1) - (x4 - x3)*(y2 - y1);
    b /= (y4 - y3)*(x2 - x1) - (x4 - x3)*(y2 - y1);
    if ((a >= 0) && (a < 1) && (b >= 0) && (b < 1)) {
      return 1;
    } else {
      return 0;
    }
  } else {
    return 0;
  }
}

// -----------------------------------------------------------------------------
int PointSet::IsInside(double x, double y) const
{
  int i;

  // compute no of points
  if (_n == 0) return false;

  // compute centre
  double cx = 0;
  double cy = 0;
  for (i = 0; i < _n; ++i) {
    cx += _data[i]._x;
    cy += _data[i]._y;
  }
  cx /= _n;
  cy /= _n;

  // compute a point outside the polygon.
  double ox = 0, oy = 0;
  for (i = 0; i < _n; ++i) {
    double tmp;

    tmp = _data[i]._x-cx;
    if (tmp<0) tmp = -tmp;
    if (tmp>ox) ox = tmp;

    tmp = _data[i]._y-cy;
    if (tmp<0) tmp = -tmp;
    if (tmp>oy) oy = tmp;
  }
  ox = cx + ox + oy + 1;
  oy = cy + ox + oy + 1;

  // count crossings.
  int crossings = 0;
  for (i = 0; i < _n; ++i) {
    crossings += intersection(_data[i]._x, _data[i]._y, _data[(i+1)%_n]._x,
                              _data[(i+1)%_n]._y, ox, oy, x, y);
  }

  // inside iff there was an odd number of crossings.
  return crossings % 2 != 0;
}

// -----------------------------------------------------------------------------
Point PointSet::Centroid() const
{
  Point c;
  if (_n > 0) {
    for (int i = 0; i < _n; ++i) {
      c += _data[i];
    }
    c /= _n;
  }
  return c;
}


} // namespace mirtk
