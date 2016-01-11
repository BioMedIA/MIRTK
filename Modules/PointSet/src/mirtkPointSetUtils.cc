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

#include <mirtkPointSetUtils.h>

#include <mirtkVtk.h>
#include <mirtkUtils.h>
#include <mirtkMath.h>
#include <mirtkMemory.h>
#include <mirtkStream.h>
#include <mirtkPath.h>
#include <mirtkParallel.h>
#include <mirtkEdgeTable.h>
#include <mirtkDataStatistics.h>
#include <mirtkPointSet.h>
#include <mirtkVoxel.h>
#include <mirtkBaseImage.h>
#include <mirtkVector3D.h>
#include <mirtkMatrix3x3.h>
#include <mirtkVector3.h>
#include <mirtkAlgorithm.h>

#include <vtkSmartPointer.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkStructuredGrid.h>
#include <vtkDataArray.h>
#include <vtkUnsignedShortArray.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkDataSetAttributes.h>
#include <vtkCell.h>
#include <vtkTriangle.h>
#include <vtkTetra.h>
#include <vtkGenericCell.h>
#include <vtkImageData.h>
#include <vtkImageStencilData.h>

#include <vtkOBJReader.h>
#include <vtkSTLReader.h>
#include <vtkSTLWriter.h>
#include <vtkPLYReader.h>
#include <vtkPLYWriter.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkGenericDataObjectReader.h>
#include <vtkXMLGenericDataObjectReader.h>
#include <vtkXMLDataSetWriter.h>
#include <vtkDataSetWriter.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>

#include <vtkImageStencilIterator.h>
#include <vtkPolyDataToImageStencil.h>
#include <vtkCleanPolyData.h>
#include <vtkTriangleFilter.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkHull.h>
#include <vtkMaskPoints.h>
#include <vtkDelaunay3D.h>
#include <vtkMassProperties.h>

#include "dfsurface.h" // BrainSuite surface (.dfs) I/O functions

#include <algorithm>


namespace mirtk {


// =============================================================================
// I/O
// =============================================================================

// -----------------------------------------------------------------------------
const char *DefaultExtension(vtkDataSet *dataset)
{
  if      (vtkPolyData        ::SafeDownCast(dataset)) return ".vtp";
  else if (vtkUnstructuredGrid::SafeDownCast(dataset)) return ".vtu";
  else if (vtkStructuredGrid  ::SafeDownCast(dataset)) return ".vts";
  return ".vtk";
}

// -----------------------------------------------------------------------------
vtkSmartPointer<vtkPointSet> ReadPointSet(const char *fname, int *ftype, bool exit_on_failure)
{
  vtkSmartPointer<vtkPointSet> pointset;
  const string ext = Extension(fname);
  if (ext == ".vtp" || ext == ".stl" || ext == ".ply" || ext == ".obj" || ext == ".dfs") {
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
/// Read BrainSuite surface from .dfs file
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
/// Write surface in BrainSuite .dfs format
int WriteDFS(const char *fname, vtkPolyData *polydata)
{
  // Ensure surface is triangular mesh
  vtkSmartPointer<vtkPolyData> mesh = Triangulate(polydata);
  SILT::DFSurface surface;
  // Copy vertex coordinates
  double p[3];
  surface.vertices.resize(mesh->GetNumberOfPoints());
  for (vtkIdType i = 0; i < mesh->GetNumberOfPoints(); ++i) {
    mesh->GetPoint(i, p);
    surface.vertices[i].x = static_cast<float>(p[0]);
    surface.vertices[i].y = static_cast<float>(p[1]);
    surface.vertices[i].z = static_cast<float>(p[2]);
  }
  // Copy triangular faces
  vtkIdType npts, *pts;
  surface.triangles.resize(mesh->GetNumberOfCells());
  for (vtkIdType i = 0; i < mesh->GetNumberOfCells(); ++i) {
    mesh->GetCellPoints(i, npts, pts);
    if (npts != 3) return false;
    surface.triangles[i].a = static_cast<int>(pts[0]);
    surface.triangles[i].b = static_cast<int>(pts[1]);
    surface.triangles[i].c = static_cast<int>(pts[2]);
  }
  // Copy vertex normals
  vtkDataArray *normals = mesh->GetPointData()->GetNormals();
  if (!normals) normals = mesh->GetPointData()->GetArray("Normals");
  if (normals) {
    double n[3];
    surface.vertexNormals.resize(mesh->GetNumberOfPoints());
    for (vtkIdType i = 0; i < mesh->GetNumberOfPoints(); ++i) {
      normals->GetTuple(i, n);
      surface.vertexNormals[i].x = static_cast<float>(n[0]);
      surface.vertexNormals[i].y = static_cast<float>(n[1]);
      surface.vertexNormals[i].z = static_cast<float>(n[2]);
    }
  }
  // Copy vertex colors
  vtkDataArray *colors = mesh->GetPointData()->GetArray("Colors");
  if (colors) {
    double rgb[3];
    surface.vertexColors.resize(mesh->GetNumberOfPoints());
    for (vtkIdType i = 0; i < mesh->GetNumberOfPoints(); ++i) {
      colors->GetTuple(i, rgb);
      surface.vertexColors[i].x = static_cast<float>(rgb[0]);
      surface.vertexColors[i].y = static_cast<float>(rgb[1]);
      surface.vertexColors[i].z = static_cast<float>(rgb[2]);
    }
  }
  // Copy vertex UV coordinates
  vtkDataArray *coords = mesh->GetPointData()->GetTCoords();
  if (!coords)  coords = mesh->GetPointData()->GetArray("UV");
  if (coords && (coords->GetNumberOfComponents() == 2 || coords->GetNumberOfComponents() == 3)) {
    surface.vertexUV.resize(mesh->GetNumberOfPoints());
    for (vtkIdType i = 0; i < mesh->GetNumberOfPoints(); ++i) {
      surface.vertexUV[i].u = static_cast<float>(coords->GetComponent(i, 0));
      surface.vertexUV[i].v = static_cast<float>(coords->GetComponent(i, 1));
    }
  }
  // Copy vertex labels
  vtkDataArray *labels = mesh->GetPointData()->GetArray("Labels");
  if (labels && labels->GetNumberOfComponents() == 1) {
    surface.vertexLabels.resize(mesh->GetNumberOfPoints());
    for (vtkIdType i = 0; i < mesh->GetNumberOfPoints(); ++i) {
      surface.vertexLabels[i] = static_cast<unsigned short>(labels->GetComponent(i, 0));
    }
  }
  // Copy vertex attributes
  vtkDataArray *scalars = mesh->GetPointData()->GetScalars();
  if (!scalars) scalars = mesh->GetPointData()->GetArray("Attributes");
  if (scalars && scalars->GetNumberOfComponents() == 1) {
    surface.vertexAttributes.resize(mesh->GetNumberOfPoints());
    for (vtkIdType i = 0; i < mesh->GetNumberOfPoints(); ++i) {
      surface.vertexAttributes[i] = static_cast<float>(scalars->GetComponent(i, 0));
    }
  }
  // Write .dfs file
  return static_cast<int>(surface.writeDFS(fname));
}

// -----------------------------------------------------------------------------
/// Write dataset points to TetGen .node file
int WriteTetGenNode(ostream &os, vtkPolyData *polydata)
{
  vtkPointData *pd = polydata->GetPointData();
  vtkDataArray *ar;
  int nattributes = 0;
  for (int i = 0; i < pd->GetNumberOfArrays(); ++i) {
    nattributes += pd->GetArray(i)->GetNumberOfComponents();
  }
  os << polydata->GetNumberOfPoints() << " 3 " << nattributes << " 0\n";
  double p[3];
  streamsize precision = os.precision();
  for (vtkIdType ptId = 0; ptId < polydata->GetNumberOfPoints(); ++ptId) {
    os << (ptId + 1) << " ";
    os.precision(8); // default TetGen tolerance is 1e-8
    polydata->GetPoint(ptId, p);
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
  return os.fail() ? 0 : 1;
}

// -----------------------------------------------------------------------------
/// Write dataset points to TetGen .node file
int WriteTetGenNode(const char *fname, vtkPolyData *polydata)
{
  ofstream os(fname);
  if (!os.is_open()) return false;
  return WriteTetGenNode(os, polydata);
}

// -----------------------------------------------------------------------------
/// Write polygonal dataset to TetGen .poly file
int WriteTetGenPoly(const char *fname, vtkPolyData *polydata, const PointSet *holes)
{
  ofstream os(fname);
  if (!os.is_open()) return 0;
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
  return os.fail() ? 0 : 1;
}

// -----------------------------------------------------------------------------
/// Write polygonal dataset to TetGen .smesh file
int WriteTetGenSMesh(const char *fname, vtkPolyData *polydata, const PointSet *holes)
{
  ofstream os(fname);
  if (!os.is_open()) return 0;
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
  return os.fail() ? 0 : 1;
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
// Point set domain
// =============================================================================

// -----------------------------------------------------------------------------
ImageAttributes PointSetDomain(vtkPointSet *data, double dx, double dy, double dz)
{
  ImageAttributes attr;
  attr._dx = dx;
  attr._dy = dy >= .0 ? dy : dx;
  attr._dz = dz >= .0 ? dz : dx;
  attr._dt = .0;
  // Compute eigenvectors of covariance matrix
  double c[3], p[3];
  data->GetCenter(c);
  Matrix3x3 covar(.0);
  for (vtkIdType i = 0; i < data->GetNumberOfPoints(); ++i) {
    data->GetPoint(i, p);
    for (int d = 0; d < 3; ++d) p[d] -= c[d];
    for (int r = 0; r < 3; ++r) {
      for (int c = 0; c < 3; ++c) {
        covar[r][c] += p[r] * p[c];
      }
    }
  }
  double  eigen[3];
  Vector3 axis [3];
  covar.EigenSolveSymmetric(eigen, axis);
  Vector3::GenerateOrthonormalBasis(axis[0], axis[1], axis[2]);
  // Set output origin and orientation
  attr._xorigin = c[0];
  attr._yorigin = c[1];
  attr._zorigin = c[2];
  for (int d = 0; d < 3; ++d) {
    attr._xaxis[d] = axis[0][d];
    attr._yaxis[d] = axis[1][d];
    attr._zaxis[d] = axis[2][d];
  }
  // Determine bounds of reoriented data set
  double x, y, z;
  double xmin = numeric_limits<double>::max(), ymin = xmin, zmin = xmin;
  double xmax = -xmin, ymax = -ymin, zmax = -zmin;
  OrderedSet<double> xs, ys, zs;
  for (vtkIdType i = 0; i < data->GetNumberOfPoints(); ++i) {
    data->GetPoint(i, p);
    x = attr._xaxis[0] * p[0] + attr._xaxis[1] * p[1] + attr._xaxis[2] * p[2];
    y = attr._yaxis[0] * p[0] + attr._yaxis[1] * p[1] + attr._yaxis[2] * p[2];
    z = attr._zaxis[0] * p[0] + attr._zaxis[1] * p[1] + attr._zaxis[2] * p[2];
    if (x < xmin) xmin = x;
    if (x > xmax) xmax = x;
    if (y < ymin) ymin = y;
    if (y > ymax) ymax = y;
    if (z < zmin) zmin = z;
    if (z > zmax) zmax = z;
    if (attr._dx <= .0) xs.insert(x);
    if (attr._dy <= .0) ys.insert(y);
    if (attr._dz <= .0) zs.insert(z);
  }
  // Set output resolution and size
  const double extent[3]  = {xmax - xmin, ymax - ymin, zmax - zmin};
  const double avg_extent = (extent[0] + extent[1] + extent[2]) / 3.0;
  if (attr._dx <= .0) attr._dx = ((extent[0] / avg_extent > 1e-3) ? AverageInterval(xs) : .0);
  if (attr._dy <= .0) attr._dy = ((extent[1] / avg_extent > 1e-3) ? AverageInterval(ys) : .0);
  if (attr._dz <= .0) attr._dz = ((extent[2] / avg_extent > 1e-3) ? AverageInterval(zs) : .0);
  attr._x  = (attr._dx > .0 ? iround(extent[0] / attr._dx) : 0) + 1;
  attr._y  = (attr._dy > .0 ? iround(extent[1] / attr._dy) : 0) + 1;
  attr._z  = (attr._dz > .0 ? iround(extent[2] / attr._dz) : 0) + 1;
  attr._dx = (attr._x  >  1 ? extent[0] / (attr._x - 1) : .0);
  attr._dy = (attr._y  >  1 ? extent[1] / (attr._y - 1) : .0);
  attr._dz = (attr._z  >  1 ? extent[2] / (attr._z - 1) : .0);
  return attr;
}

// -----------------------------------------------------------------------------
ImageAttributes PointSetDomain(vtkPointSet *data, const Vector3D<double> &ds)
{
  return PointSetDomain(data, ds._x, ds._y, ds._z);
}

// -----------------------------------------------------------------------------
ImageAttributes PolyDataDomain(vtkPolyData *data, double dx, double dy, double dz)
{
  return PointSetDomain(data, dx, dy, dz);
}

// -----------------------------------------------------------------------------
ImageAttributes PolyDataDomain(vtkPolyData *data, const Vector3D<double> &ds)
{
  return PointSetDomain(data, ds);
}

// =============================================================================
// Point/cell data
// =============================================================================

// -----------------------------------------------------------------------------
int PolyDataAttributeType(const char *type)
{
  string ltype(type);
  transform(ltype.begin(), ltype.end(), ltype.begin(), ::tolower);
  if (ltype == "scalars")     return vtkDataSetAttributes::SCALARS;
  if (ltype == "vectors")     return vtkDataSetAttributes::VECTORS;
  if (ltype == "normals")     return vtkDataSetAttributes::NORMALS;
  if (ltype == "tcoords")     return vtkDataSetAttributes::TCOORDS;
  if (ltype == "tensors")     return vtkDataSetAttributes::TENSORS;
  if (ltype == "globalids")   return vtkDataSetAttributes::GLOBALIDS;
  if (ltype == "pedigreeids") return vtkDataSetAttributes::PEDIGREEIDS;
  if (ltype == "edgeflag")    return vtkDataSetAttributes::EDGEFLAG;
  return -1;
}

// -----------------------------------------------------------------------------
vtkDataArray *GetArrayByCaseInsensitiveName(vtkDataSetAttributes *data, const char *name, int *loc)
{
  string lname = ToLower(name), lower_name;
  for (int i = 0; i < data->GetNumberOfArrays(); ++i) {
    const char *array_name = data->GetArrayName(i);
    if (array_name) {
      lower_name = ToLower(array_name);
      if (lower_name == lname) {
        if (loc) *loc = i;
        return data->GetArray(i);
      }
    }
  }
  if (loc) *loc = -1;
  return NULL;
}

// -----------------------------------------------------------------------------
int DeepCopyArrayUsingCaseInsensitiveName(vtkDataSetAttributes *dst, vtkDataSetAttributes *src, const char *name, bool deep)
{
  int loc = -1;
  vtkDataArray *src_array = GetArrayByCaseInsensitiveName(src, name);
  if (src_array) {
    vtkDataArray *dst_array = GetArrayByCaseInsensitiveName(dst, name, &loc);
    if (dst_array) {
      dst_array->DeepCopy(src_array);
    } else {
      vtkSmartPointer<vtkDataArray> copy;
      copy = vtkSmartPointer<vtkDataArray>::NewInstance(src_array);
      copy->DeepCopy(src_array);
      loc = dst->AddArray(copy);
    }
  }
  return loc;
}

// =============================================================================
// Cells
// =============================================================================

// -----------------------------------------------------------------------------
// Attention: vtkGenericCell is not derived from the respective cell types
double ComputeArea(vtkCell *cell)
{
  if (cell->GetCellDimension() < 2) return .0;
  vtkPoints *points = cell->GetPoints();
  switch (cell->GetCellType()) {
    case VTK_TRIANGLE: {
      double p1[3], p2[3], p3[3];
      points->GetPoint(0, p1);
      points->GetPoint(1, p2);
      points->GetPoint(2, p3);
      return vtkTriangle::TriangleArea(p1, p2, p3);
    }
    default: return numeric_limits<double>::quiet_NaN();
  }
}

// -----------------------------------------------------------------------------
// Attention: vtkGenericCell is not derived from the respective cell types
double ComputeVolume(vtkCell *cell)
{
  if (cell->GetCellDimension() < 3) return .0;
  vtkPoints *points = cell->GetPoints();
  switch (cell->GetCellType()) {
    case VTK_TETRA: {
      double p1[3], p2[3], p3[3], p4[3];
      points->GetPoint(0, p1);
      points->GetPoint(1, p2);
      points->GetPoint(2, p3);
      points->GetPoint(3, p4);
      return vtkTetra::ComputeVolume(p1, p2, p3, p4);
    }
    default: return numeric_limits<double>::quiet_NaN();
  }
}

// =============================================================================
// Basic point set manipulation
// =============================================================================

// -----------------------------------------------------------------------------
vtkSmartPointer<vtkPolyData> DataSetSurface(vtkSmartPointer<vtkDataSet> dataset, bool passPtIds, bool passCellIds)
{
  vtkSmartPointer<vtkDataSetSurfaceFilter> surface;
  surface = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
  SetVTKInput(surface, dataset);
  surface->SetPassThroughPointIds(passPtIds);
  surface->SetPassThroughCellIds(passCellIds);
  surface->Update();
  return surface->GetOutput();
}

// -----------------------------------------------------------------------------
void Center(vtkSmartPointer<vtkPointSet> pointset)
{
  double c[3], p[3];
  pointset->GetCenter(c);
  vtkPoints *points = pointset->GetPoints();
  for (vtkIdType ptId = 0; ptId < pointset->GetNumberOfPoints(); ++ptId) {
    points->GetPoint(ptId, p);
    p[0] -= c[0];
    p[1] -= c[1];
    p[2] -= c[2];
    points->SetPoint(ptId, p);
  }
  pointset->Modified();
}

// -----------------------------------------------------------------------------
void Scale(vtkSmartPointer<vtkPointSet> pointset, double scale)
{
  double c[3], p[3];
  pointset->GetCenter(c);
  vtkPoints *points = pointset->GetPoints();
  for (vtkIdType ptId = 0; ptId < pointset->GetNumberOfPoints(); ++ptId) {
    points->GetPoint(ptId, p);
    p[0] = c[0] + scale * (p[0] - c[0]);
    p[1] = c[1] + scale * (p[1] - c[1]);
    p[2] = c[2] + scale * (p[2] - c[2]);
    points->SetPoint(ptId, p);
  }
  pointset->Modified();
}

// =============================================================================
// Surface meshes
// =============================================================================

// -----------------------------------------------------------------------------
namespace AreaUtils {

struct ComputeSurfaceArea
{
  vtkPolyData  *_Surface;
  vtkDataArray *_CellArea;
  double        _SurfaceArea;

  ComputeSurfaceArea(vtkPolyData *surface)
  :
    _Surface(surface),
    _CellArea(surface->GetCellData()->GetArray("Area")),
    _SurfaceArea(.0)
  {}

  ComputeSurfaceArea(const ComputeSurfaceArea &other, split)
  :
    _Surface(other._Surface), _CellArea(other._CellArea), _SurfaceArea(.0)
  {}

  void join(const ComputeSurfaceArea &other)
  {
    _SurfaceArea += other._SurfaceArea;
  }

  void operator ()(const blocked_range<vtkIdType> cellIds)
  {
    vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();
    vtkIdType npts, *pts;
    double p1[3], p2[3], p3[3], area;
    for (vtkIdType cellId = cellIds.begin(); cellId != cellIds.end(); ++cellId)
    {
      _Surface->GetCellPoints(cellId, npts, pts);
      if (npts == 3) {
        _Surface->GetPoint(pts[0], p1);
        _Surface->GetPoint(pts[1], p2);
        _Surface->GetPoint(pts[2], p3);
        area = vtkTriangle::TriangleArea(p1, p2, p3);
      } else {
        _Surface->GetCell(cellId, cell);
        area = ComputeArea(cell);
      }
      if (_CellArea) _CellArea->SetComponent(cellId, 0, area);
      _SurfaceArea += area;
    }
  }
};

} // namespace AreaUtils

// -----------------------------------------------------------------------------
double Area(vtkSmartPointer<vtkPolyData> polydata, bool per_cell)
{
  vtkSmartPointer<vtkDataArray> area;
  if (per_cell) {
    area = vtkSmartPointer<vtkFloatArray>::New();
    area->SetName("Area");
    area->SetNumberOfComponents(1);
    area->SetNumberOfTuples(polydata->GetNumberOfCells());
    polydata->GetCellData()->AddArray(area);
  }
  AreaUtils::ComputeSurfaceArea eval(polydata);
  blocked_range<vtkIdType> cellIds(0, polydata->GetNumberOfCells());
  parallel_reduce(cellIds, eval);
  return eval._SurfaceArea;
}

// -----------------------------------------------------------------------------
double Area(vtkSmartPointer<vtkPointSet> pointset)
{
  return Area(DataSetSurface(pointset), false);
}


namespace EdgeLengthUtils {

// -----------------------------------------------------------------------------
/// Compute lengths of all edges
struct ComputeEdgeLength
{
  vtkPoints       *_Points;
  const EdgeTable *_EdgeTable;
  double          *_EdgeLength;

  void operator ()(const blocked_range<int> &re) const
  {
    int    ptId1, ptId2, edgeId;
    double p1[3], p2[3];

    EdgeIterator it(*_EdgeTable);
    for (it.InitTraversal(re); (edgeId = it.GetNextEdge(ptId1, ptId2)) != -1;) {
      _Points->GetPoint(ptId1, p1);
      _Points->GetPoint(ptId2, p2);
      _EdgeLength[edgeId] = sqrt(vtkMath::Distance2BetweenPoints(p1, p2));
    }
  }
};

// -----------------------------------------------------------------------------
/// Determine minimum/maximum edge length
struct MinMaxEdgeLength
{
  vtkPoints       *_Points;
  const EdgeTable *_EdgeTable;
  double           _Min;
  double           _Max;

  MinMaxEdgeLength() : _Min(numeric_limits<double>::infinity()), _Max(-_Min) {}

  MinMaxEdgeLength(const MinMaxEdgeLength &other, split)
  :
    _Points(other._Points),
    _EdgeTable(other._EdgeTable),
    _Min(other._Min),
    _Max(other._Max)
  {}

  void join(const MinMaxEdgeLength &other)
  {
    if (other._Min < _Min) _Min = other._Min;
    if (other._Max < _Max) _Max = other._Max;
  }

  void operator ()(const blocked_range<int> &re)
  {
    int    ptId1, ptId2;
    double p1[3], p2[3], d;

    EdgeIterator it(*_EdgeTable);
    for (it.InitTraversal(re); it.GetNextEdge(ptId1, ptId2) != -1;) {
      _Points->GetPoint(ptId1, p1);
      _Points->GetPoint(ptId2, p2);
      d = sqrt(vtkMath::Distance2BetweenPoints(p1, p2));
      if (d < _Min) _Min = d;
      if (d > _Max) _Max = d;
    }
  }
};

// -----------------------------------------------------------------------------
/// Calculate sum of edge lengths
struct SumEdgeLengths
{
  vtkPoints       *_Points;
  const EdgeTable *_EdgeTable;
  double           _Sum;

  SumEdgeLengths() : _Sum(.0) {}

  SumEdgeLengths(const SumEdgeLengths &other, split)
  :
    _Points(other._Points),
    _EdgeTable(other._EdgeTable),
    _Sum(.0)
  {}

  void join(const SumEdgeLengths &other)
  {
    _Sum += other._Sum;
  }

  void operator ()(const blocked_range<int> &re)
  {
    int    ptId1, ptId2;
    double p1[3], p2[3];

    EdgeIterator it(*_EdgeTable);
    for (it.InitTraversal(re); it.GetNextEdge(ptId1, ptId2) != -1;) {
      _Points->GetPoint(ptId1, p1);
      _Points->GetPoint(ptId2, p2);
      _Sum += sqrt(vtkMath::Distance2BetweenPoints(p1, p2));
    }
  }
};

// -----------------------------------------------------------------------------
/// Calculate sum of edge lengths
struct SumEdgeLengthsWithDuplicates
{
  vtkPoints       *_Points;
  const EdgeTable *_EdgeTable;
  double           _Sum;

  SumEdgeLengthsWithDuplicates() : _Sum(.0) {}

  SumEdgeLengthsWithDuplicates(const SumEdgeLengthsWithDuplicates &other, split)
  :
    _Points(other._Points),
    _EdgeTable(other._EdgeTable),
    _Sum(.0)
  {}

  void join(const SumEdgeLengthsWithDuplicates &other)
  {
    _Sum += other._Sum;
  }

  void operator ()(const blocked_range<int> &ptIds)
  {
    const int *adjPtIds;
    int        numAdjPts;
    double     p1[3], p2[3];

    for (int ptId = ptIds.begin(); ptId != ptIds.end(); ++ptId) {
      _Points->GetPoint(ptId, p1);
      _EdgeTable->GetAdjacentPoints(ptId, numAdjPts, adjPtIds);
      for (int i = 0; i < numAdjPts; ++i) {
        _Points->GetPoint(adjPtIds[i], p2);
        _Sum += sqrt(vtkMath::Distance2BetweenPoints(p1, p2));
      }
    }
  }
};


} // namespace EdgeLengthUtils

// -----------------------------------------------------------------------------
double AverageEdgeLength(vtkSmartPointer<vtkPoints> points, const EdgeTable &edgeTable)
{
  if (edgeTable.NumberOfEdges() == 0) return .0;
  // FIXME: Parallel reduction of sum of edge lengths yields different result
  //        from the code below or a single-threaded execution. This
  //        could indicate a bug in irtkEdgeIterator::InitTraversal.
//  EdgeLengthUtils::SumEdgeLengths eval;
//  eval._Points    = points;
//  eval._EdgeTable = &edgeTable;
//  parallel_reduce(blocked_range<int>(0, edgeTable.NumberOfEdges()), eval);
//  eval(blocked_range<int>(0, edgeTable.NumberOfEdges()));
//  return eval._Sum / edgeTable.NumberOfEdges();
  EdgeLengthUtils::SumEdgeLengthsWithDuplicates eval;
  eval._Points    = points;
  eval._EdgeTable = &edgeTable;
  parallel_reduce(blocked_range<int>(0, static_cast<int>(points->GetNumberOfPoints())), eval);
  return eval._Sum / (2 * edgeTable.NumberOfEdges());
}

// -----------------------------------------------------------------------------
double AverageEdgeLength(vtkSmartPointer<vtkPointSet> pointset)
{
  EdgeTable edgeTable(pointset);
  return AverageEdgeLength(pointset->GetPoints(), edgeTable);
}

// -----------------------------------------------------------------------------
double RobustAverageEdgeLength(vtkSmartPointer<vtkPoints> points, const EdgeTable &edgeTable)
{
  const int n = edgeTable.NumberOfEdges();
  if (n == 0) return .0;
  double *edgeLength = new double[n];
  EdgeLengthUtils::ComputeEdgeLength eval;
  eval._Points     = points;
  eval._EdgeTable  = &edgeTable;
  eval._EdgeLength = edgeLength;
  parallel_for(blocked_range<int>(0, n), eval);
  double mean = data::statistic::RobustMean::Calculate(5, n, edgeLength);
  delete[] edgeLength;
  return mean;
}

// -----------------------------------------------------------------------------
double RobustAverageEdgeLength(vtkSmartPointer<vtkPointSet> pointset)
{
  EdgeTable edgeTable(pointset);
  return RobustAverageEdgeLength(pointset->GetPoints(), edgeTable);
}

// -----------------------------------------------------------------------------
double MedianEdgeLength(vtkSmartPointer<vtkPoints> points, const EdgeTable &edgeTable)
{
  const int n = edgeTable.NumberOfEdges();
  if (n == 0) return .0;
  double *edgeLength = new double[n];
  EdgeLengthUtils::ComputeEdgeLength eval;
  eval._Points     = points;
  eval._EdgeTable  = &edgeTable;
  eval._EdgeLength = edgeLength;
  parallel_for(blocked_range<int>(0, n), eval);
  sort(edgeLength, edgeLength + n);
  double median = edgeLength[n / 2];
  delete[] edgeLength;
  return median;
}

// -----------------------------------------------------------------------------
double MedianEdgeLength(vtkSmartPointer<vtkPointSet> pointset)
{
  EdgeTable edgeTable(pointset);
  return MedianEdgeLength(pointset->GetPoints(), edgeTable);
}

// -----------------------------------------------------------------------------
void GetMinMaxEdgeLength(vtkSmartPointer<vtkPoints> points, const EdgeTable &edgeTable, double &min, double &max)
{
  min = max = .0;
  if (edgeTable.NumberOfEdges() == 0) return;
  EdgeLengthUtils::MinMaxEdgeLength eval;
  eval._Points    = points;
  eval._EdgeTable = &edgeTable;
  parallel_reduce(blocked_range<int>(0, edgeTable.NumberOfEdges()), eval);
  min = eval._Min;
  max = eval._Max;
}

// -----------------------------------------------------------------------------
void GetMinMaxEdgeLength(vtkSmartPointer<vtkPointSet> pointset, double &min, double &max)
{
  EdgeTable edgeTable(pointset);
  GetMinMaxEdgeLength(pointset->GetPoints(), edgeTable, min, max);
}

// -----------------------------------------------------------------------------
double MinEdgeLength(vtkSmartPointer<vtkPoints> points, const EdgeTable &edgeTable)
{
  double min, max;
  GetMinMaxEdgeLength(points, edgeTable, min, max);
  return min;
}

// -----------------------------------------------------------------------------
double MinEdgeLength(vtkSmartPointer<vtkPointSet> pointset)
{
  EdgeTable edgeTable(pointset);
  return MinEdgeLength(pointset->GetPoints(), edgeTable);
}

// -----------------------------------------------------------------------------
double MaxEdgeLength(vtkSmartPointer<vtkPoints> points, const EdgeTable &edgeTable)
{
  double min, max;
  GetMinMaxEdgeLength(points, edgeTable, min, max);
  return max;
}

// -----------------------------------------------------------------------------
double MaxEdgeLength(vtkSmartPointer<vtkPointSet> pointset)
{
  EdgeTable edgeTable(pointset);
  return MaxEdgeLength(pointset->GetPoints(), edgeTable);
}

// -----------------------------------------------------------------------------
void EdgeLengthNormalDistribution(vtkSmartPointer<vtkPoints> points, const EdgeTable &edgeTable, double &mean, double &sigma)
{
  mean = sigma = .0;
  if (edgeTable.NumberOfEdges() == 0) return;

  int    ptId1, ptId2, n = 0;
  double p1[3], p2[3], d, delta;

  EdgeIterator it(edgeTable);
  for (it.InitTraversal(); it.GetNextEdge(ptId1, ptId2) != -1;) {
    points->GetPoint(ptId1, p1);
    points->GetPoint(ptId2, p2);
    d = sqrt(vtkMath::Distance2BetweenPoints(p1, p2));
    ++n;
    delta = d - mean;
    mean  += delta / n;
    sigma += delta * (d - mean);
  }

  if (n > 1) sigma /= n - 1;
  sigma = sqrt(sigma);
}

// -----------------------------------------------------------------------------
void EdgeLengthNormalDistribution(vtkSmartPointer<vtkPointSet> pointset, double &mean, double &sigma)
{
  EdgeTable edgeTable(pointset);
  EdgeLengthNormalDistribution(pointset->GetPoints(), edgeTable, mean, sigma);
}

// -----------------------------------------------------------------------------
double GetVolume(vtkSmartPointer<vtkPolyData> surface)
{
  vtkSmartPointer<vtkMassProperties> mp = vtkMassProperties::New();
  SetVTKInput(mp, surface);
  return mp->GetVolume();
}

// -----------------------------------------------------------------------------
bool IsSurfaceMesh(vtkDataSet *dataset)
{
  vtkPolyData *poly = vtkPolyData::SafeDownCast(dataset);
  return poly && poly->GetPolys() && poly->GetPolys()->GetNumberOfCells() > 0 &&
         (!poly->GetLines() || poly->GetLines()->GetNumberOfCells() == 0) &&
         (!poly->GetVerts() || poly->GetVerts()->GetNumberOfCells() == 0);
}

// -----------------------------------------------------------------------------
vtkSmartPointer<vtkPolyData> Triangulate(vtkSmartPointer<vtkPolyData> mesh)
{
  vtkSmartPointer<vtkPolyData> output;
  if (IsTriangularMesh(mesh)) {
    output = vtkSmartPointer<vtkPolyData>::NewInstance(mesh);
    output->ShallowCopy(mesh);
  } else {
    vtkSmartPointer<vtkTriangleFilter> filter;
    filter = vtkSmartPointer<vtkTriangleFilter>::New();
    SetVTKInput(filter, mesh);
    filter->PassVertsOff();
    filter->PassLinesOff();
    vtkSmartPointer<vtkCleanPolyData> cleaner;
    cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
    cleaner->PointMergingOff();
    cleaner->ConvertLinesToPointsOn();
    cleaner->ConvertPolysToLinesOn();
    cleaner->ConvertStripsToPolysOn();
    SetVTKConnection(cleaner, filter);
    cleaner->Update();
    cleaner->GetOutput()->SetLines(NULL);
    cleaner->GetOutput()->SetVerts(NULL);
    output = cleaner->GetOutput();
  }
  return output;
}

// -----------------------------------------------------------------------------
// Implements two different filter pipelines to extract the convex hull
vtkSmartPointer<vtkPolyData> ConvexHull(vtkSmartPointer<vtkPointSet> pointset, int levels)
{
  // Spatially stratify points to prevent "Unable to factor linear system"
  // warning of vtkDelaunay3D filter due to numerical imprecisions
  vtkNew<vtkMaskPoints> stratify;
  stratify->RandomModeOn();
  stratify->SetRandomModeType(2);
  stratify->SetMaximumNumberOfPoints(.75 * pointset->GetNumberOfPoints());
  SetVTKInput(stratify, pointset);

  // Get convex hull of largest component
  vtkNew<vtkHull> hull;
  hull->AddRecursiveSpherePlanes(levels);
  //SetVTKConnection(hull, stratify);
  SetVTKInput(hull, pointset);

  // Compute Delaunay triangulation
  vtkNew<vtkDelaunay3D> delaunay;
  SetVTKConnection(delaunay, hull);

  // Construct surface mesh
  vtkNew<vtkDataSetSurfaceFilter> mesher;
  SetVTKConnection(mesher, delaunay);

  mesher->Update();
  return mesher->GetOutput();
}

// -----------------------------------------------------------------------------
bool IsTriangularMesh(vtkDataSet *input)
{
  if (vtkPointSet::SafeDownCast(input) == NULL) return false;
  for (vtkIdType cellId = 0; cellId < input->GetNumberOfCells(); ++cellId) {
    int type = input->GetCellType(cellId);
    if (type != VTK_EMPTY_CELL && type != VTK_TRIANGLE) return false;
  }
  return true;
}

// -----------------------------------------------------------------------------
bool IsTetrahedralMesh(vtkDataSet *input)
{
  if (vtkPointSet::SafeDownCast(input) == NULL) return false;
  for (vtkIdType cellId = 0; cellId < input->GetNumberOfCells(); ++cellId) {
    int type = input->GetCellType(cellId);
    if (type != VTK_EMPTY_CELL && type != VTK_TETRA) return false;
  }
  return true;
}

// -----------------------------------------------------------------------------
vtkSmartPointer<vtkPointSet> Tetrahedralize(vtkSmartPointer<vtkPointSet> input)
{
  vtkSmartPointer<vtkPointSet> mesh;
  if (IsTetrahedralMesh(input)) {
    mesh = vtkSmartPointer<vtkPointSet>::NewInstance(input);
    mesh->ShallowCopy(input);
  } else {
    // TODO: Use TetGen library to tetrahedralize interior of input PLC
    cerr << "irtkPolyDataUtils::Tetrahedralize: Not implemented, use TetGen command instead" << endl;
    exit(1);
  }
  return mesh;
}

// -----------------------------------------------------------------------------
vtkSmartPointer<vtkImageData> NewVtkMask(int nx, int ny, int nz)
{
  vtkSmartPointer<vtkImageData> imagedata = vtkSmartPointer<vtkImageData>::New();
  imagedata->SetOrigin(.0, .0, .0);
  imagedata->SetDimensions(nx, ny, nz);
  imagedata->SetSpacing(1.0, 1.0, 1.0);
#if VTK_MAJOR_VERSION >= 6
  imagedata->AllocateScalars(ToVTKDataType(MIRTK_VOXEL_BINARY), 1);
#else
  imagedata->SetScalarType(ToVTKDataType(MIRTK_VOXEL_BINARY));
  imagedata->AllocateScalars();
#endif
  return imagedata;
}

// -----------------------------------------------------------------------------
vtkSmartPointer<vtkPointSet> WorldToImage(vtkSmartPointer<vtkPointSet> pointset,
                                          const BaseImage             *image)
{
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  points->SetNumberOfPoints(pointset->GetNumberOfPoints());
  double p[3];
  for (vtkIdType ptId = 0; ptId < points->GetNumberOfPoints(); ++ptId) {
    pointset->GetPoint(ptId, p);
    image->WorldToImage(p[0], p[1], p[2]);
    points->SetPoint(ptId, p);
  }
  vtkSmartPointer<vtkPointSet> output;
  output = vtkSmartPointer<vtkPointSet>::NewInstance(pointset);
  output->ShallowCopy(pointset);
  output->SetPoints(points);
  return output;
}

// -----------------------------------------------------------------------------
vtkSmartPointer<vtkImageStencilData> ImageStencil(vtkSmartPointer<vtkImageData> image,
                                                  vtkSmartPointer<vtkPointSet>  pointset)
{
  vtkSmartPointer<vtkPolyData> surface = DataSetSurface(pointset);
  vtkSmartPointer<vtkPolyDataToImageStencil> filter;
  filter = vtkSmartPointer<vtkPolyDataToImageStencil>::New();
  SetVTKInput(filter, surface);
  filter->SetOutputOrigin(image->GetOrigin());
  filter->SetOutputSpacing(image->GetSpacing());
  filter->SetOutputWholeExtent(image->GetExtent());
  filter->Update();
  return filter->GetOutput();
}

// -----------------------------------------------------------------------------
void ImageStencilToMask(vtkSmartPointer<vtkImageStencilData> stencil,
                        vtkSmartPointer<vtkImageData>        image)
{
  if (image->GetScalarType() != ToVTKDataType(MIRTK_VOXEL_BINARY)) {
    cerr << "ImageStencilToMask: vtkImageData must have scalar type MIRTK_VOXEL_BINARY" << endl;
    exit(1);
  }
  const vtkIdType nvox = image->GetNumberOfPoints();
  memset(image->GetScalarPointer(), 0, nvox * sizeof(BinaryPixel));
  vtkImageStencilIterator<BinaryPixel> it;
  it.Initialize(image, stencil, image->GetExtent());
  while (!it.IsAtEnd()) {
    if (it.IsInStencil()) {
      for (BinaryPixel *cur = it.BeginSpan(); cur != it.EndSpan(); ++cur) {
        *cur = static_cast<BinaryPixel>(true);
      }
    }
    it.NextSpan();
  }
}


} // namespace mirtk
