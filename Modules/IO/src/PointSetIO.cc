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

#include "mirtk/PointSetIO.h"

#include "mirtk/Path.h"
#include "mirtk/Stream.h"
#include "mirtk/Vtk.h"

#include "vtkPoints.h"
#include "vtkPointData.h"
#include "vtkDataArray.h"
#include "vtkCellArray.h"
#include "vtkIdList.h"

#include "vtkUnsignedShortArray.h"
#include "vtkFloatArray.h"

#include "vtkUnstructuredGrid.h"
#include "vtkStructuredGrid.h"

#include "vtkGenericDataObjectReader.h"
#include "vtkDataSetWriter.h"

#include "vtkPolyDataReader.h"
#include "vtkPolyDataWriter.h"

#include "vtkXMLGenericDataObjectReader.h"
#include "vtkXMLDataSetWriter.h"

#include "vtkXMLPolyDataReader.h"
#include "vtkXMLPolyDataWriter.h"

#include "vtkOBJReader.h"
#include "vtkPLYReader.h"
#include "vtkSTLReader.h"

#include "vtkPLYWriter.h"
#include "vtkSTLWriter.h"

#include "brainsuite/dfsurface.h"


namespace mirtk {


// =============================================================================
// File name extension
// =============================================================================

// -----------------------------------------------------------------------------
const char *DefaultExtension(vtkDataSet *dataset)
{
  if      (vtkPolyData        ::SafeDownCast(dataset)) return ".vtp";
  else if (vtkUnstructuredGrid::SafeDownCast(dataset)) return ".vtu";
  else if (vtkStructuredGrid  ::SafeDownCast(dataset)) return ".vts";
  return ".vtk";
}

// =============================================================================
// Generic I/O functions
// =============================================================================

// -----------------------------------------------------------------------------
vtkSmartPointer<vtkPointSet> ReadPointSet(const char *fname, int *ftype, bool exit_on_failure)
{
  vtkSmartPointer<vtkPointSet> pointset;
  const string ext = Extension(fname);
  if (ext == ".vtp" || ext == ".stl" || ext == ".ply" || ext == ".obj" || ext == ".dfs" || ext == ".off") {
    pointset = ReadPolyData(fname);
  } else if (ext.length() == 4  && ext.substr(0, 3) == ".vt" && ext != ".vtk") {
    vtkSmartPointer<vtkXMLGenericDataObjectReader> reader;
    reader = vtkSmartPointer<vtkXMLGenericDataObjectReader>::New();
    reader->SetFileName(fname);
    reader->Update();
    pointset = vtkPointSet::SafeDownCast(reader->GetOutput());
  } else {
    vtkSmartPointer<vtkGenericDataObjectReader> reader;
    reader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
    reader->SetFileName(fname);
    reader->Update();
    if (ftype) *ftype = reader->GetFileType();
    pointset = vtkPointSet::SafeDownCast(reader->GetOutput());
  }
  if (exit_on_failure && (!pointset || pointset->GetNumberOfPoints() == 0)) {
    cerr << "File " << fname << " either contains no points or could not be read" << endl;
    exit(1);
  }
  return pointset;
}

// -----------------------------------------------------------------------------
bool WritePointSet(const char *fname, vtkPointSet *pointset, bool compress, bool ascii)
{
  vtkPolyData *polydata = vtkPolyData::SafeDownCast(pointset);
  if (polydata) return WritePolyData(fname, polydata, compress, ascii);
  const string ext = Extension(fname);
  int success = 0;
  if (ext.length() == 4 && ext.substr(0, 3) == ".vt" && ext != ".vtk") {
    vtkSmartPointer<vtkXMLDataSetWriter> writer;
    writer = vtkSmartPointer<vtkXMLDataSetWriter>::New();
    SetVTKInput(writer, pointset);
    writer->SetFileName(fname);
    if (compress) writer->SetCompressorTypeToZLib();
    else          writer->SetCompressorTypeToNone();
    success = writer->Write();
  } else {
    vtkSmartPointer<vtkDataSetWriter> writer;
    writer = vtkSmartPointer<vtkDataSetWriter>::New();
    SetVTKInput(writer, pointset);
    writer->SetFileName(fname);
    if (ascii) writer->SetFileTypeToASCII();
    else       writer->SetFileTypeToBinary();
    success = writer->Write();
  }
  return (success == 1);
}

// -----------------------------------------------------------------------------
vtkSmartPointer<vtkPolyData> ReadPolyData(const char *fname, int *ftype, bool exit_on_failure)
{
  vtkSmartPointer<vtkPolyData> polydata;
  const string ext = Extension(fname);
  if (ext == ".vtp") {
    vtkSmartPointer<vtkXMLPolyDataReader> reader;
    reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
    reader->SetFileName(fname);
    reader->Update();
    polydata = reader->GetOutput();
  } else if (ext == ".stl") {
    vtkSmartPointer<vtkSTLReader> reader;
    reader = vtkSmartPointer<vtkSTLReader>::New();
    reader->SetFileName(fname);
    reader->Update();
    polydata = reader->GetOutput();
  } else if (ext == ".ply") {
    vtkSmartPointer<vtkPLYReader> reader;
    reader = vtkSmartPointer<vtkPLYReader>::New();
    reader->SetFileName(fname);
    reader->Update();
    polydata = reader->GetOutput();
  } else if (ext == ".obj") {
    vtkSmartPointer<vtkOBJReader> reader;
    reader = vtkSmartPointer<vtkOBJReader>::New();
    reader->SetFileName(fname);
    reader->Update();
    polydata = reader->GetOutput();
  } else if (ext == ".dfs") {
    polydata = ReadDFS(fname);
  } else if (ext == ".off") {
    polydata = ReadOFF(fname);
  } else {
    vtkSmartPointer<vtkPolyDataReader> reader;
    reader = vtkSmartPointer<vtkPolyDataReader>::New();
    reader->SetFileName(fname);
    reader->Update();
    if (ftype) *ftype = reader->GetFileType();
    polydata = reader->GetOutput();
  }
  if (exit_on_failure && polydata->GetNumberOfPoints() == 0) {
    cerr << "File " << fname << " either contains no points or could not be read" << endl;
    exit(1);
  }
  return polydata;
}

// -----------------------------------------------------------------------------
bool WritePolyData(const char *fname, vtkPolyData *polydata, bool compress, bool ascii)
{
  const string ext = Extension(fname);
  int success = 0;
  if (ext == ".vtp") {
    vtkSmartPointer<vtkXMLPolyDataWriter> writer;
    writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    SetVTKInput(writer, polydata);
    writer->SetFileName(fname);
    if (compress) writer->SetCompressorTypeToZLib();
    else          writer->SetCompressorTypeToNone();
    success = writer->Write();
  } else if (ext == ".stl") {
    vtkSmartPointer<vtkSTLWriter> writer;
    writer = vtkSmartPointer<vtkSTLWriter>::New();
    SetVTKInput(writer, polydata);
    if (ascii) writer->SetFileTypeToASCII();
    else       writer->SetFileTypeToBinary();
    writer->SetFileName(fname);
    success = writer->Write();
  } else if (ext == ".ply") {
    vtkSmartPointer<vtkPLYWriter> writer;
    writer = vtkSmartPointer<vtkPLYWriter>::New();
    SetVTKInput(writer, polydata);
    if (ascii) writer->SetFileTypeToASCII();
    else       writer->SetFileTypeToBinary();
    writer->SetFileName(fname);
    success = writer->Write();
  } else if (ext == ".node") {
    success = WriteTetGenNode(fname, polydata);
  } else if (ext == ".poly") {
    success = WriteTetGenPoly(fname, polydata);
  } else if (ext == ".smesh") {
    success = WriteTetGenSMesh(fname, polydata);
  } else if (ext == ".dfs") {
    success = WriteDFS(fname, polydata);
  } else if (ext == ".off") {
    success = WriteOFF(fname, polydata);
  } else {
    vtkSmartPointer<vtkPolyDataWriter> writer;
    writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    SetVTKInput(writer, polydata);
    writer->SetFileName(fname);
    if (ascii) writer->SetFileTypeToASCII();
    else       writer->SetFileTypeToBinary();
    success = writer->Write();
  }
  return (success == 1);
}

// =============================================================================
// BrainSuite I/O functions
// =============================================================================

// -----------------------------------------------------------------------------
vtkSmartPointer<vtkPolyData> ReadDFS(const char *fname)
{
  vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
  // Read .dfs file
  SILT::DFSurface surface;
  if (!surface.readDFS(fname)) return polydata;
  const vtkIdType npoints = static_cast<vtkIdType>(surface.vertices .size());
  const vtkIdType ncells  = static_cast<vtkIdType>(surface.triangles.size());
  // Copy vertex coordinates
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  points->SetNumberOfPoints(npoints);
  double p[3];
  for (vtkIdType i = 0; i < npoints; ++i) {
    p[0] = static_cast<double>(surface.vertices[i].x);
    p[1] = static_cast<double>(surface.vertices[i].y);
    p[2] = static_cast<double>(surface.vertices[i].z);
    points->SetPoint(i, p);
  }
  polydata->SetPoints(points);
  // Copy triangle face list
  vtkIdType pts[3];
  vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
  cells->Allocate(cells->EstimateSize(ncells, 3));
  for (vtkIdType i = 0; i < ncells; ++i) {
    pts[0] = static_cast<vtkIdType>(surface.triangles[i].a);
    pts[1] = static_cast<vtkIdType>(surface.triangles[i].b);
    pts[2] = static_cast<vtkIdType>(surface.triangles[i].c);
    cells->InsertNextCell(3, pts);
  }
  polydata->SetPolys(cells);
  // Copy vertex normals
  if (!surface.vertexNormals.empty()) {
    double n[3];
    vtkSmartPointer<vtkFloatArray> normals = vtkSmartPointer<vtkFloatArray>::New();
    normals->SetName("Normals");
    normals->SetNumberOfComponents(3);
    normals->SetNumberOfTuples(npoints);
    for (vtkIdType i = 0; i < npoints; ++i) {
      n[0] = static_cast<double>(surface.vertexNormals[i].x);
      n[1] = static_cast<double>(surface.vertexNormals[i].y);
      n[2] = static_cast<double>(surface.vertexNormals[i].z);
      normals->SetTuple(i, n);
    }
    polydata->GetPointData()->SetNormals(normals);
  }
  // Copy vertex colors
  if (!surface.vertexColors.empty()) {
    double rgb[3];
    vtkSmartPointer<vtkFloatArray> colors = vtkSmartPointer<vtkFloatArray>::New();
    colors->SetName("Colors");
    colors->SetNumberOfComponents(3);
    colors->SetNumberOfTuples(npoints);
    for (vtkIdType i = 0; i < npoints; ++i) {
      rgb[0] = static_cast<double>(surface.vertexColors[i].x);
      rgb[1] = static_cast<double>(surface.vertexColors[i].y);
      rgb[2] = static_cast<double>(surface.vertexColors[i].z);
      colors->SetTuple(i, rgb);
    }
    polydata->GetPointData()->AddArray(colors);
  }
  // Copy vertex UV coordinates
  if (!surface.vertexUV.empty()) {
    double uv[3] = {.0};
    vtkSmartPointer<vtkFloatArray> coords = vtkSmartPointer<vtkFloatArray>::New();
    coords->SetName("UV");
    coords->SetNumberOfComponents(3);
    coords->SetNumberOfTuples(npoints);
    for (vtkIdType i = 0; i < npoints; ++i) {
      uv[0] = static_cast<double>(surface.vertexUV[i].u);
      uv[1] = static_cast<double>(surface.vertexUV[i].v);
      coords->SetTuple(i, uv);
    }
    polydata->GetPointData()->SetTCoords(coords);
  }
  // Copy vertex labels
  if (!surface.vertexLabels.empty()) {
    vtkSmartPointer<vtkUnsignedShortArray> labels = vtkSmartPointer<vtkUnsignedShortArray>::New();
    labels->SetName("Labels");
    labels->SetNumberOfComponents(1);
    labels->SetNumberOfTuples(npoints);
    for (vtkIdType i = 0; i < npoints; ++i) {
      labels->SetValue(i, surface.vertexLabels[i]);
    }
    polydata->GetPointData()->AddArray(labels);
  }
  // Copy vertex attributes
  if (!surface.vertexAttributes.empty()) {
    vtkSmartPointer<vtkFloatArray> scalars = vtkSmartPointer<vtkFloatArray>::New();
    scalars->SetName("Attributes");
    scalars->SetNumberOfComponents(1);
    scalars->SetNumberOfTuples(npoints);
    for (vtkIdType i = 0; i < npoints; ++i) {
      scalars->SetValue(i, surface.vertexAttributes[i]);
    }
    polydata->GetPointData()->SetScalars(scalars);
  }
  return polydata;
}

// -----------------------------------------------------------------------------
bool WriteDFS(const char *fname, vtkPolyData *polydata)
{
  SILT::DFSurface surface;
  // Copy vertex coordinates
  double p[3];
  surface.vertices.resize(polydata->GetNumberOfPoints());
  for (vtkIdType i = 0; i < polydata->GetNumberOfPoints(); ++i) {
    polydata->GetPoint(i, p);
    surface.vertices[i].x = static_cast<float>(p[0]);
    surface.vertices[i].y = static_cast<float>(p[1]);
    surface.vertices[i].z = static_cast<float>(p[2]);
  }
  // Copy triangular faces
  vtkIdType npts, *pts;
  surface.triangles.resize(polydata->GetNumberOfCells());
  for (vtkIdType i = 0; i < polydata->GetNumberOfCells(); ++i) {
    polydata->GetCellPoints(i, npts, pts);
    if (npts != 3) return false;
    surface.triangles[i].a = static_cast<int>(pts[0]);
    surface.triangles[i].b = static_cast<int>(pts[1]);
    surface.triangles[i].c = static_cast<int>(pts[2]);
  }
  // Copy vertex normals
  vtkDataArray *normals = polydata->GetPointData()->GetNormals();
  if (!normals) normals = polydata->GetPointData()->GetArray("Normals");
  if (normals) {
    double n[3];
    surface.vertexNormals.resize(polydata->GetNumberOfPoints());
    for (vtkIdType i = 0; i < polydata->GetNumberOfPoints(); ++i) {
      normals->GetTuple(i, n);
      surface.vertexNormals[i].x = static_cast<float>(n[0]);
      surface.vertexNormals[i].y = static_cast<float>(n[1]);
      surface.vertexNormals[i].z = static_cast<float>(n[2]);
    }
  }
  // Copy vertex colors
  vtkDataArray *colors = polydata->GetPointData()->GetArray("Colors");
  if (colors) {
    double rgb[3];
    surface.vertexColors.resize(polydata->GetNumberOfPoints());
    for (vtkIdType i = 0; i < polydata->GetNumberOfPoints(); ++i) {
      colors->GetTuple(i, rgb);
      surface.vertexColors[i].x = static_cast<float>(rgb[0]);
      surface.vertexColors[i].y = static_cast<float>(rgb[1]);
      surface.vertexColors[i].z = static_cast<float>(rgb[2]);
    }
  }
  // Copy vertex UV coordinates
  vtkDataArray *coords = polydata->GetPointData()->GetTCoords();
  if (!coords)  coords = polydata->GetPointData()->GetArray("UV");
  if (coords && (coords->GetNumberOfComponents() == 2 || coords->GetNumberOfComponents() == 3)) {
    surface.vertexUV.resize(polydata->GetNumberOfPoints());
    for (vtkIdType i = 0; i < polydata->GetNumberOfPoints(); ++i) {
      surface.vertexUV[i].u = static_cast<float>(coords->GetComponent(i, 0));
      surface.vertexUV[i].v = static_cast<float>(coords->GetComponent(i, 1));
    }
  }
  // Copy vertex labels
  vtkDataArray *labels = polydata->GetPointData()->GetArray("Labels");
  if (labels && labels->GetNumberOfComponents() == 1) {
    surface.vertexLabels.resize(polydata->GetNumberOfPoints());
    for (vtkIdType i = 0; i < polydata->GetNumberOfPoints(); ++i) {
      surface.vertexLabels[i] = static_cast<unsigned short>(labels->GetComponent(i, 0));
    }
  }
  // Copy vertex attributes
  vtkDataArray *scalars = polydata->GetPointData()->GetScalars();
  if (!scalars) scalars = polydata->GetPointData()->GetArray("Attributes");
  if (scalars && scalars->GetNumberOfComponents() == 1) {
    surface.vertexAttributes.resize(polydata->GetNumberOfPoints());
    for (vtkIdType i = 0; i < polydata->GetNumberOfPoints(); ++i) {
      surface.vertexAttributes[i] = static_cast<float>(scalars->GetComponent(i, 0));
    }
  }
  // Write .dfs file
  return surface.writeDFS(fname);
}

// =============================================================================
// Object File Format I/O functions
// =============================================================================

// -----------------------------------------------------------------------------
vtkSmartPointer<vtkPolyData> ReadOFF(const char *fname)
{
  vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();

  ifstream ifs(fname);
  if (!ifs) return polydata;

  string keyword;
  ifs >> keyword;
  if (keyword != "OFF" && keyword != "off") return polydata;

  int numVertices = -1, numFaces = -1, numEdges = -1;
  ifs >> numVertices >> numFaces >> numEdges;
  if (ifs.fail()) return polydata;

  if (numVertices < 0) numVertices = 0;
  if (numFaces    < 0) numFaces    = 0;

  vtkSmartPointer<vtkPoints> points;
  points = vtkSmartPointer<vtkPoints>::New();
  points->SetNumberOfPoints(numVertices);

  double p[3];
  for (int i = 0; i < numVertices; ++i) {
    ifs >> p[0] >> p[1] >> p[2];
    if (ifs.fail()) break;
    points->SetPoint(i, p);
  }
  if (ifs.fail()) return polydata;

  vtkSmartPointer<vtkCellArray> verts, lines, polys;
  vtkSmartPointer<vtkIdList> cell = vtkSmartPointer<vtkIdList>::New();
  verts = vtkSmartPointer<vtkCellArray>::New();
  lines = vtkSmartPointer<vtkCellArray>::New();
  polys = vtkSmartPointer<vtkCellArray>::New();
  for (int i = 0, ptId, n; i < numFaces; ++i) {
    ifs >> n;
    if (ifs.fail()) break;
    if (n > 0) {
      cell->Reset();
      for (int j = 0; j < n; ++j) {
        ifs >> ptId;
        cell->InsertNextId(ptId);
      }
      if      (n == 1) verts->InsertNextCell(cell);
      else if (n == 2) lines->InsertNextCell(cell);
      else             polys->InsertNextCell(cell);
      if (!ifs.good()) break;
    }
  }
  if (ifs.fail()) return polydata;

  verts->Squeeze();
  lines->Squeeze();
  polys->Squeeze();

  polydata->SetPoints(points);
  if (verts->GetNumberOfCells() > 0) polydata->SetVerts(verts);
  if (lines->GetNumberOfCells() > 0) polydata->SetLines(lines);
  if (polys->GetNumberOfCells() > 0) polydata->SetPolys(polys);

  return polydata;
}

// -----------------------------------------------------------------------------
bool WriteOFF(const char *fname, vtkPolyData *polydata)
{
  ofstream ofs(fname);
  ofs.precision(12);

  ofs << "OFF\n";
  ofs << polydata->GetNumberOfPoints() << " ";
  ofs << polydata->GetNumberOfCells()  << " 0\n";

  double p[3];
  for (vtkIdType ptId = 0; ptId < polydata->GetNumberOfPoints(); ++ptId) {
    polydata->GetPoint(ptId, p);
    ofs << p[0] << " " << p[1] << " " << p[2] << "\n";
  }

  vtkIdType numPts, *ptIds;
  polydata->BuildCells();
  for (vtkIdType cellId = 0; cellId < polydata->GetNumberOfCells(); ++cellId) {
    polydata->GetCellPoints(cellId, numPts, ptIds);
    ofs << numPts;
    for (vtkIdType i = 0; i < numPts; ++i) {
      ofs << " " << ptIds[i];
    }
    ofs << "\n";
  }

  return !ofs.fail();
}

// =============================================================================
// TetGen I/O functions
// =============================================================================

// -----------------------------------------------------------------------------
/// Write point set to TetGen .node output stream
///
/// @param[in, out] os       Output stream.
/// @param[in]      pointset Point set.
///
/// @return Whether point set was written successfully to the given output stream.
bool WriteTetGenNode(ostream &os, vtkPointSet *pointset)
{
  vtkPointData *pd = pointset->GetPointData();
  vtkDataArray *ar;
  int nattributes = 0;
  for (int i = 0; i < pd->GetNumberOfArrays(); ++i) {
    nattributes += pd->GetArray(i)->GetNumberOfComponents();
  }
  os << pointset->GetNumberOfPoints() << " 3 " << nattributes << " 0\n";
  double p[3];
  streamsize precision = os.precision();
  for (vtkIdType ptId = 0; ptId < pointset->GetNumberOfPoints(); ++ptId) {
    os << (ptId + 1) << " ";
    os.precision(8); // default TetGen tolerance is 1e-8
    pointset->GetPoint(ptId, p);
    os << " " << p[0] << " " << p[1] << " " << p[2];
    os.precision(5);
    for (int i = 0; i < pd->GetNumberOfArrays(); ++i) {
      ar = pd->GetArray(i);
      for (int j = 0; j < ar->GetNumberOfComponents(); ++j) {
        os << " " << ar->GetComponent(ptId, j);
      }
    }
    os << "\n";
  }
  os.precision(precision);
  return !os.fail();
}

// -----------------------------------------------------------------------------
bool WriteTetGenNode(const char *fname, vtkPointSet *pointset)
{
  ofstream os(fname);
  if (!os.is_open()) return false;
  return WriteTetGenNode(os, pointset);
}

// -----------------------------------------------------------------------------
bool WriteTetGenPoly(const char *fname, vtkPolyData *polydata, const PointSet *holes)
{
  ofstream os(fname);
  if (!os.is_open()) return false;
  os << "# part 1: nodes\n";
  WriteTetGenNode(os, polydata);
  vtkIdType npts, *pts;
  vtkCellArray *verts  = polydata->GetVerts();
  vtkCellArray *lines  = polydata->GetLines();
  vtkCellArray *polys  = polydata->GetPolys();
  vtkCellArray *strips = polydata->GetStrips();
  int nfacets = 0;
  if (verts ->GetNumberOfCells() > 0) nfacets += 1;
  if (lines ->GetNumberOfCells() > 0) nfacets += 1;
  if (polys ->GetNumberOfCells() > 0) nfacets += 1;
  if (strips->GetNumberOfCells() > 0) nfacets += 1;
  os << "\n# part 2: facets\n";
  os << nfacets << " 0\n";
  if (verts->GetNumberOfCells() > 0) {
    os << "# verts\n";
    os << verts->GetNumberOfCells() << "\n";
    verts->InitTraversal();
    while (verts->GetNextCell(npts, pts)) {
      os << npts << " ";
      for (vtkIdType i = 0; i < npts; ++i) os << " " << (pts[i] + 1);
      os << "\n";
    }
  }
  if (lines->GetNumberOfCells() > 0) {
    os << "# lines\n";
    os << lines->GetNumberOfCells() << "\n";
    lines->InitTraversal();
    while (lines->GetNextCell(npts, pts)) {
      os << npts << " ";
      for (vtkIdType i = 0; i < npts; ++i) os << " " << (pts[i] + 1);
      os << "\n";
    }
  }
  if (polys->GetNumberOfCells() > 0) {
    os << "# polys\n";
    os << polys->GetNumberOfCells() << "\n";
    polys->InitTraversal();
    while (polys->GetNextCell(npts, pts)) {
      os << npts << " ";
      for (vtkIdType i = 0; i < npts; ++i) os << " " << (pts[i] + 1);
      os << "\n";
    }
  }
  if (strips->GetNumberOfCells() > 0) {
    os << "# strips\n";
    os << strips->GetNumberOfCells() << "\n";
    strips->InitTraversal();
    while (strips->GetNextCell(npts, pts)) {
      os << npts << " ";
      for (vtkIdType i = 0; i < npts; ++i) os << " " << (pts[i] + 1);
      os << "\n";
    }
  }
  os << "\n# part 3: hole list\n";
  if (holes) {
    os << holes->Size() << "\n";
    for (int i = 0; i < holes->Size(); ++i) {
      const Point &p = holes->GetPoint(i);
      os << (i+1) << " " << p._x << " " << p._y << " " << p._z << "\n";
    }
  } else {
    os << "0\n";
  }
  os << "\n# part 4: region list\n";
  os << "0\n";
  return !os.fail();
}

// -----------------------------------------------------------------------------
bool WriteTetGenSMesh(const char *fname, vtkPolyData *polydata, const PointSet *holes)
{
  ofstream os(fname);
  if (!os.is_open()) return false;
  os << "# part 1: nodes\n";
  WriteTetGenNode(os, polydata);
  vtkIdType npts, *pts;
  vtkCellArray *verts  = polydata->GetVerts();
  vtkCellArray *lines  = polydata->GetLines();
  vtkCellArray *polys  = polydata->GetPolys();
  vtkCellArray *strips = polydata->GetStrips();
  int nfacets = 0;
  if (verts ->GetNumberOfCells() > 0) nfacets += verts ->GetNumberOfCells();
  if (lines ->GetNumberOfCells() > 0) nfacets += lines ->GetNumberOfCells();
  if (polys ->GetNumberOfCells() > 0) nfacets += polys ->GetNumberOfCells();
  if (strips->GetNumberOfCells() > 0) nfacets += strips->GetNumberOfCells();
  os << "\n# part 2: facets\n";
  os << nfacets << " 0\n";
  if (verts->GetNumberOfCells() > 0) {
    verts->InitTraversal();
    while (verts->GetNextCell(npts, pts)) {
      os << npts << " ";
      for (vtkIdType i = 0; i < npts; ++i) os << " " << (pts[i] + 1);
      os << "\n";
    }
  }
  if (lines->GetNumberOfCells() > 0) {
    lines->InitTraversal();
    while (lines->GetNextCell(npts, pts)) {
      os << npts << " ";
      for (vtkIdType i = 0; i < npts; ++i) os << " " << (pts[i] + 1);
      os << "\n";
    }
  }
  if (polys->GetNumberOfCells() > 0) {
    polys->InitTraversal();
    while (polys->GetNextCell(npts, pts)) {
      os << npts << " ";
      for (vtkIdType i = 0; i < npts; ++i) os << " " << (pts[i] + 1);
      os << "\n";
    }
  }
  if (strips->GetNumberOfCells() > 0) {
    strips->InitTraversal();
    while (strips->GetNextCell(npts, pts)) {
      os << npts;
      for (vtkIdType i = 0; i < npts; ++i) os << " " << (pts[i] + 1);
      os << "\n";
    }
  }
  os << "\n# part 3: hole list\n";
  if (holes) {
    os << holes->Size() << "\n";
    for (int i = 0; i < holes->Size(); ++i) {
      const Point &p = holes->GetPoint(i);
      os << (i+1) << " " << p._x << " " << p._y << " " << p._z << "\n";
    }
  } else {
    os << "0\n";
  }
  os << "\n# part 4: region list\n";
  os << "0\n";
  return !os.fail();
}


} // namespace mirtk
