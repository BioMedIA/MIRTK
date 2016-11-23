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
#include "mirtk/System.h" // GetUser, GetDateTime
#include "mirtk/UnorderedMap.h"
#include "mirtk/Vtk.h"

#include "vtkPoints.h"
#include "vtkPointData.h"
#include "vtkDataArray.h"
#include "vtkCellArray.h"
#include "vtkIdTypeArray.h"
#include "vtkIdList.h"
#include "vtkInformation.h"
#include "vtkInformationIterator.h"

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

#if MIRTK_IO_WITH_GIFTI
  #include "mirtk/NiftiImageInfo.h"
  #include "gifti/gifti_io.h"
#endif


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
vtkSmartPointer<vtkPointSet> ReadPointSet(const char *fname, bool exit_on_failure)
{
  FileOption fopt;
  return ReadPointSet(fname, fopt, exit_on_failure);
}

// -----------------------------------------------------------------------------
vtkSmartPointer<vtkPointSet> ReadPointSet(const char *fname, FileOption &fopt, bool exit_on_failure)
{
  const string ext = Extension(fname);
  if (ext == ".vtp" || ext == ".stl" || ext == ".ply" || ext == ".obj" || ext == ".dfs" || ext == ".off" || ext == ".gii") {
    return ReadPolyData(fname, fopt, exit_on_failure);
  }
  fopt = FO_Default;
  vtkSmartPointer<vtkPointSet> pointset;
  if (ext == ".txt" || ext == ".csv" || ext == ".tsv") {
    char sep    = ',';
    if (ext == ".tsv") sep = '\t';
    pointset = ReadPointSetTable(fname, sep);
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
    pointset = vtkPointSet::SafeDownCast(reader->GetOutput());
    fopt = (reader->GetFileType() == VTK_ASCII ? FO_ASCII : FO_Binary);
  }
  if (exit_on_failure && (!pointset || pointset->GetNumberOfPoints() == 0)) {
    cerr << "File " << fname << " either contains no points or could not be read" << endl;
    exit(1);
  }
  return pointset;
}

// -----------------------------------------------------------------------------
bool WritePointSet(const char *fname, vtkPointSet *pointset, FileOption fopt)
{
  vtkPolyData *polydata = vtkPolyData::SafeDownCast(pointset);
  if (polydata) return WritePolyData(fname, polydata, fopt);

  int success = 0;
  const string ext = Extension(fname);
  if (ext == ".txt" || ext == ".csv" || ext == ".tsv") {
    string type = fname;
    type = type.substr(0, type.length() - ext.length());
    type = Extension(type, EXT_Last);
    bool coords = (type != ".attr");
    char sep    = ',';
    if (ext == ".tsv") sep = '\t';
    WritePointSetTable(fname, pointset, sep, true, coords);
  } else if (ext.length() == 4 && ext.substr(0, 3) == ".vt" && ext != ".vtk") {
    vtkSmartPointer<vtkXMLDataSetWriter> writer;
    writer = vtkSmartPointer<vtkXMLDataSetWriter>::New();
    SetVTKInput(writer, pointset);
    writer->SetFileName(fname);
    if (fopt == FO_NoCompress) writer->SetCompressorTypeToNone();
    else                       writer->SetCompressorTypeToZLib();
    success = writer->Write();
  } else {
    vtkSmartPointer<vtkDataSetWriter> writer;
    writer = vtkSmartPointer<vtkDataSetWriter>::New();
    SetVTKInput(writer, pointset);
    writer->SetFileName(fname);
    if (fopt == FO_ASCII) writer->SetFileTypeToASCII();
    else                  writer->SetFileTypeToBinary();
    success = writer->Write();
  }
  return (success == 1);
}

// -----------------------------------------------------------------------------
vtkSmartPointer<vtkPolyData> ReadPolyData(const char *fname, bool exit_on_failure)
{
  FileOption fopt;
  return ReadPolyData(fname, fopt, exit_on_failure);
}

// -----------------------------------------------------------------------------
vtkSmartPointer<vtkPolyData> ReadPolyData(const char *fname, FileOption &fopt, bool exit_on_failure)
{
  fopt = FO_Default;
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
  } else if (ext == ".gii") {
    #if MIRTK_IO_WITH_GIFTI
      polydata = ReadGIFTI(fname, nullptr, exit_on_failure);
    #else
      if (exit_on_failure) {
        cerr << "Error: File '" << fname << "' cannot be read because MIRTK I/O library was built without GIFTI support!" << endl;
        exit(1);
      }
    #endif
  } else {
    vtkSmartPointer<vtkPolyDataReader> reader;
    reader = vtkSmartPointer<vtkPolyDataReader>::New();
    reader->SetFileName(fname);
    reader->Update();
    polydata = reader->GetOutput();
    fopt = (reader->GetFileType() == VTK_ASCII ? FO_ASCII : FO_Binary);
  }
  if (exit_on_failure && polydata->GetNumberOfPoints() == 0) {
    cerr << "Error: File '" << fname << "' either contains no points or could not be read!" << endl;
    exit(1);
  }
  return polydata;
}

// -----------------------------------------------------------------------------
bool WritePolyData(const char *fname, vtkPolyData *polydata, FileOption fopt)
{
  const string ext = Extension(fname);
  int success = 0;
  if (ext == ".vtp") {
    vtkSmartPointer<vtkXMLPolyDataWriter> writer;
    writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    SetVTKInput(writer, polydata);
    writer->SetFileName(fname);
    if (fopt == FO_NoCompress) writer->SetCompressorTypeToNone();
    else                       writer->SetCompressorTypeToZLib();
    success = writer->Write();
  } else if (ext == ".stl") {
    vtkSmartPointer<vtkSTLWriter> writer;
    writer = vtkSmartPointer<vtkSTLWriter>::New();
    SetVTKInput(writer, polydata);
    if (fopt == FO_ASCII) writer->SetFileTypeToASCII();
    else                  writer->SetFileTypeToBinary();
    writer->SetFileName(fname);
    success = writer->Write();
  } else if (ext == ".ply") {
    vtkSmartPointer<vtkPLYWriter> writer;
    writer = vtkSmartPointer<vtkPLYWriter>::New();
    SetVTKInput(writer, polydata);
    if (fopt == FO_ASCII) writer->SetFileTypeToASCII();
    else                  writer->SetFileTypeToBinary();
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
  } else if (ext == ".gii") {
    #if MIRTK_IO_WITH_GIFTI
      success = WriteGIFTI(fname, polydata, fopt);
    #else
      cerr << "Error: Cannot write surface to GIFTI file because MIRTK I/O library was built without GIFTI support!" << endl;
    #endif
  } else if (ext == ".txt" || ext == ".csv" || ext == ".tsv") {
    string type = fname;
    type = type.substr(0, type.length() - ext.length());
    type = Extension(type, EXT_Last);
    bool coords = (type != ".attr");
    char sep    = ',';
    if (ext == ".tsv") sep = '\t';
    WritePointSetTable(fname, polydata, sep, true, coords);
  } else {
    vtkSmartPointer<vtkPolyDataWriter> writer;
    writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    SetVTKInput(writer, polydata);
    writer->SetFileName(fname);
    if (fopt == FO_ASCII) writer->SetFileTypeToASCII();
    else                  writer->SetFileTypeToBinary();
    success = writer->Write();
  }
  return (success == 1);
}

// =============================================================================
// CSV I/O functions
// =============================================================================

// -----------------------------------------------------------------------------
vtkSmartPointer<vtkPointSet> ReadPointSetTable(const char *fname, char sep, vtkPointSet *pointset)
{
  const bool errmsg = true;

  ifstream ifs(fname);
  if (!ifs.is_open()) {
    if (errmsg) cerr << "Failed to open file " << fname << endl;
    return vtkSmartPointer<vtkPolyData>::New();
  }

  string line;
  string lstr;
  Array<string> cols;
  double value;

  // header / column names
  if (!getline(ifs, line) || line.empty()) {
    if (errmsg) cerr << "IOError: Failed to read column header" << endl;
    return vtkSmartPointer<vtkPolyData>::New();
  }
  if (line.back() == '\r') line.pop_back();
  Array<string> header = Split(line, sep);
  if (header.empty()) {
    if (errmsg) cerr << "Header is an empty line" << endl;
    return vtkSmartPointer<vtkPolyData>::New();
  }
  const size_t ncols = header.size();
  for (size_t i = 0; i < ncols; ++i) {
    if (header[i].length() >= 2 && header[i].front() == '"' && header[i].back() == '"') {
      header[i] = header[i].substr(1, header[i].length()-2);
    }
  }

  // get index of point ID column
  int id_col = -1;
  for (size_t i = 0; i < ncols; ++i) {
    if (header[i].empty()) continue;
    lstr = ToLower(header[i]);
    if (lstr == "id" || lstr == "point.id") {
      if (id_col == -1) {
        id_col = static_cast<int>(i);
      } else {
        if (errmsg) cerr << "More than one point ID column found" << endl;
        return vtkSmartPointer<vtkPolyData>::New();
      }
    }
  }

  // get indices of point coordinate columns
  int x_col = -1, y_col = -1, z_col = -1;
  for (size_t i = 0; i < ncols; ++i) {
    if (header[i].empty()) continue;
    lstr = ToLower(header[i]);
    if (lstr.back() == 'x') {
      if (lstr == "x"       ||
          lstr == "p.x"     ||
          lstr == "point.x" ||
          lstr == "coord.x") {
        if (x_col == -1) {
          x_col = static_cast<int>(i);
        } else {
          if (errmsg) cerr << "More than one point x coordinate column found" << endl;
          return vtkSmartPointer<vtkPolyData>::New();
        }
      }
    } else if (lstr.back() == 'y') {
      if (lstr == "y"       ||
          lstr == "p.y"     ||
          lstr == "point.y" ||
          lstr == "coord.y") {
        if (y_col == -1) {
          y_col = static_cast<int>(i);
        } else {
          if (errmsg) cerr << "More than one point y coordinate column found" << endl;
          return vtkSmartPointer<vtkPolyData>::New();
        }
      }
    } else if (lstr.back() == 'z') {
      if (lstr == "z"       ||
          lstr == "p.z"     ||
          lstr == "point.z" ||
          lstr == "coord.z") {
        if (z_col == -1) {
          z_col = static_cast<int>(i);
        } else {
          if (errmsg) cerr << "More than one point z coordinate column found" << endl;
          return vtkSmartPointer<vtkPolyData>::New();
        }
      }
    }
  }

  // read points
  vtkSmartPointer<vtkPoints> points;
  if (x_col != -1 && y_col != -1) { // z coord optional for 2D point sets
    const auto pos = ifs.tellg();
    Array<Point> coords;
    double p[3] = {0.};
    vtkIdType ptId = 0;
    size_t irow = 0;
    while (getline(ifs, line)) {
      ++irow;
      if (line.back() == '\r') line.pop_back();
      if (line.empty()) continue;
      cols = Split(line, sep);
      if (cols.size() != ncols) {
        if (errmsg) cerr << "Row " << irow << " has different number of columsn" << endl;
        return vtkSmartPointer<vtkPolyData>::New();
      }
      if (!FromString(cols[x_col], p[0])) {
        if (errmsg) cerr << "Invalid point x coordinate in row " << irow << endl;
        return vtkSmartPointer<vtkPolyData>::New();
      }
      if (!FromString(cols[y_col], p[1])) {
        if (errmsg) cerr << "Invalid point y coordinate in row " << irow << endl;
        return vtkSmartPointer<vtkPolyData>::New();
      }
      if (z_col != -1) {
        if (!FromString(cols[z_col], p[2])) {
          if (errmsg) cerr << "Invalid point z coordinate in row " << irow << endl;
          return vtkSmartPointer<vtkPolyData>::New();
        }
      }
      if (id_col != -1) {
        if (!FromString(cols[id_col], value) || static_cast<double>(static_cast<vtkIdType>(value)) != value) {
          if (errmsg) cerr << "Invalid point ID in row " << irow << endl;
          return vtkSmartPointer<vtkPolyData>::New();
        }
        ptId = static_cast<vtkIdType>(value);
        if (ptId >= static_cast<vtkIdType>(coords.size())) {
          coords.resize(ptId + 1);
        }
        coords[ptId] = Point(p);
      } else {
        coords.push_back(Point(p));
        ++ptId;
      }
    }
    if (pointset) points.TakeReference(pointset->GetPoints()->NewInstance());
    else          points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(static_cast<vtkIdType>(coords.size()));
    for (ptId = 0; ptId < points->GetNumberOfPoints(); ++ptId) {
      points->SetPoint(ptId, coords[ptId]);
    }
    ifs.seekg(pos);
  } else {
    points = pointset->GetPoints();
  }
  const size_t npoints = static_cast<size_t>(points->GetNumberOfPoints());
  if (npoints == 0) {
    return vtkSmartPointer<vtkPolyData>::New();
  }

  // determine which columns to extract and group them
  UnorderedMap<string, OrderedSet<Pair<int, string>>> comps;
  for (size_t i = 0; i < ncols; ++i) {
    if (header[i].empty() || i == id_col || i == x_col || i == y_col || i == z_col) continue;
    const auto pos = header[i].find('.');
    if (pos != string::npos && pos != header[i].length()-1) {
      const string name = header[i].substr(0, pos);
      const string comp = header[i].substr(pos + 1);
      comps[name].insert(MakePair(static_cast<int>(i), move(comp)));
    } else {
      comps[header[i]].insert(MakePair(static_cast<int>(i), ""));
    }
  }

  // allocate data arrays
  Array<vtkSmartPointer<vtkDataArray>> data(comps.size());
  Array<Pair<int, int>> col2data(ncols, MakePair(-1, -1));

  int i = 0, j;
  for (const auto &kv : comps) {
    data[i] = NewVTKDataArray(VTK_FLOAT);
    data[i]->SetName(kv.first.c_str());
    data[i]->SetNumberOfComponents(static_cast<int>(kv.second.size()));
    data[i]->SetNumberOfTuples(npoints);
    j = 0;
    for (const auto comp : kv.second) {
      if (!comp.second.empty()) data[i]->SetComponentName(j, comp.second.c_str());
      col2data[comp.first] = MakePair(i, j);
      ++j;
    }
    ++i;
  }

  // read
  size_t ptId = 0;
  size_t irow = 0;
  while (getline(ifs, line)) {
    ++irow;
    if (line.back() == '\r') line.pop_back();
    if (line.empty()) continue;
    cols = Split(line, sep);
    if (cols.size() != ncols) {
      if (errmsg) cerr << "Row " << irow << " has different number of columsn" << endl;
      return vtkSmartPointer<vtkPolyData>::New();
    }
    if (id_col != -1) {
      if (!FromString(cols[id_col], value) || static_cast<double>(static_cast<size_t>(value)) != value) {
        if (errmsg) cerr << "Invalid point ID in row " << irow << endl;
        return vtkSmartPointer<vtkPolyData>::New();
      }
      ptId = static_cast<size_t>(value);
      if (ptId >= npoints) continue;
    } else {
      if (ptId >= npoints) break;
    }
    for (size_t c = 0; c < ncols; ++c) {
      i = col2data[c].first;
      j = col2data[c].second;
      if (i == -1 || j == -1) continue;
      if (!FromString(cols[c], value)) {
        if (errmsg) cerr << "Failed to parse value in row " << irow << ", column " << c+1 << endl;
        return vtkSmartPointer<vtkPolyData>::New();
      }
      data[i]->SetComponent(ptId, j, value);
    }
    ++ptId;
  }

  // TODO: convert to integral types when possible

  vtkSmartPointer<vtkPointSet> output;
  if (pointset) {
    output.TakeReference(pointset->NewInstance());
    output->ShallowCopy(pointset);
  } else {
    output = vtkSmartPointer<vtkPolyData>::New();
  }
  output->SetPoints(points);
  for (auto &arr : data) {
    output->GetPointData()->AddArray(arr);
  }
  return output;
}

// -----------------------------------------------------------------------------
bool WritePointSetTable(const char *fname, vtkPointSet *pointset, char sep, bool ids, bool coords)
{
  ofstream ofs(fname);
  if (!ofs.is_open()) return false;

  // output arrays and precision
  vtkPointData * const pd = pointset->GetPointData();

  Array<vtkDataArray *> arrays;
  Array<streamsize>     digits;

  arrays.reserve(pd->GetNumberOfArrays());
  digits.reserve(arrays.size() + 1u);

  for (int i = 0; i < pd->GetNumberOfArrays(); ++i) {
    vtkDataArray *arr = pd->GetArray(i);
    if (arr->GetName() != nullptr) {
      arrays.push_back(arr);
      switch (arr->GetDataType()) {
        case VTK_FLOAT:  { digits.push_back(numeric_limits<float >::max_digits10); } break;
        case VTK_DOUBLE: { digits.push_back(numeric_limits<double>::max_digits10); } break;
        default:         { digits.push_back(0); } break;
      }
    }
  }
  if (coords) {
    if (pointset->GetPoints()->GetDataType() == VTK_FLOAT) {
      digits.push_back(numeric_limits<float>::max_digits10);
    } else {
      digits.push_back(numeric_limits<double>::max_digits10);
    }
  } else {
    if (arrays.empty()) return false;
  }

  // header
  int col = 0;
  if (ids) {
    ofs << "ID";
    col += 1;
  }
  if (coords) {
    if (col > 0) ofs << sep;
    ofs << "X" << sep << "Y" << sep << "Z";
    col += 3;
  }
  for (size_t i = 0; i < arrays.size(); ++i) {
    vtkDataArray *arr = arrays[i];
    for (int j = 0; j < arr->GetNumberOfComponents(); ++j) {
      if (++col > 1) ofs << sep;
      ofs << arr->GetName();
      const char *comp = arr->GetComponentName(j);
      if (comp != nullptr) {
        ofs << '.' << comp;
      } else if (arr->GetNumberOfComponents() > 1) {
        ofs << '.' << (j+1);
      }
    }
  }
  ofs << "\n";

  // data rows
  double p[3];
  for (vtkIdType ptId = 0; ptId < pointset->GetNumberOfPoints(); ++ptId) {
    int col = 0;
    if (ids) {
      ofs.precision(0);
      ofs << ptId;
      col = 1;
    }
    if (coords) {
      pointset->GetPoint(ptId, p);
      ofs.precision(digits.back());
      if (col > 0) ofs << sep;
      ofs << p[0] << sep << p[1] << sep << p[2];
      col += 3;
    }
    for (size_t i = 0; i < arrays.size(); ++i) {
      vtkDataArray *arr = arrays[i];
      ofs.precision(digits[i]);
      for (int j = 0; j < arr->GetNumberOfComponents(); ++j) {
        if (++col > 1) ofs << sep;
        ofs << arr->GetComponent(ptId, j);
      }
    }
    ofs << "\n";
  }

  ofs.close();
  return !ofs.fail();
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

// =============================================================================
// GIFTI I/O functions
// =============================================================================
#if MIRTK_IO_WITH_GIFTI

// -----------------------------------------------------------------------------
/// Get VTK data type enumeration value corresponding to given GIFTI datatype
static int GiftiDataTypeToVtk(int datatype)
{
  switch (datatype) {
    case NIFTI_TYPE_INT8:    return VTK_CHAR;
    case NIFTI_TYPE_INT16:   return VTK_SHORT;
    case NIFTI_TYPE_INT32:   return VTK_INT;
    case NIFTI_TYPE_INT64:   return VTK_LONG_LONG;
    case NIFTI_TYPE_UINT8:   return VTK_UNSIGNED_CHAR;
    case NIFTI_TYPE_UINT16:  return VTK_UNSIGNED_SHORT;
    case NIFTI_TYPE_UINT32:  return VTK_UNSIGNED_INT;
    case NIFTI_TYPE_UINT64:  return VTK_UNSIGNED_LONG_LONG;
    case NIFTI_TYPE_FLOAT32: return VTK_FLOAT;
    case NIFTI_TYPE_FLOAT64: return VTK_DOUBLE;
    default:                 return VTK_VOID;
  }
}

// -----------------------------------------------------------------------------
/// Get GIFTI data type suitable for given VTK data array
///
/// @note GIFTI only supports NIFTI_TYPE_UINT8, NIFTI_TYPE_INT32, and NIFTI_TYPE_FLOAT32.
static int GiftiDataType(vtkDataArray *data)
{
  switch (data->GetDataType()) {
    case VTK_UNSIGNED_CHAR:
      return NIFTI_TYPE_UINT8;
    case VTK_CHAR:
    case VTK_SHORT:
    case VTK_INT:
    case VTK_UNSIGNED_SHORT:
    case VTK_UNSIGNED_INT:
      return NIFTI_TYPE_INT32;
    default:
      return NIFTI_TYPE_FLOAT32;
  }
}

// -----------------------------------------------------------------------------
/// Get GIFTI intent code for given VTK point data array
///
/// @param[in] data Point data array.
/// @param[in] attr VTK attribute type or -1 if it is no VTK attribute.
static int GiftiIntentCode(vtkDataArray *data, int attr = -1)
{
  if (data->GetName()) {
    NiftiIntent intent = NIFTI_INTENT_NONE;
    if (FromString(data->GetName(), intent) && intent != NIFTI_INTENT_NONE) {
      return intent;
    }
  }
  if (attr == vtkDataSetAttributes::NORMALS ||
      attr == vtkDataSetAttributes::VECTORS) {
    return NIFTI_INTENT_VECTOR;
  }
  if (data->GetNumberOfComponents() > 1) {
    return NIFTI_INTENT_VECTOR;
  } else if (data->GetName()) {
    const string lname = ToLower(data->GetName());
    if (lname.find("curvature") != string::npos ||
        lname == "curv" || lname == "sulc" || // Sulcal depth
        lname == "sulcal depth" ||
        lname == "sulcaldepth"  ||
        lname == "sulcal_depth" ||
        lname == "k1" || lname == "k2" || // Principle curvature
        lname == "h" || // Mean curvature
        lname == "k" || // Gaussian curvature
        lname == "c" || lname == "curvedness" || // Curvedness
        lname.find("shape") != string::npos) { 
      return NIFTI_INTENT_SHAPE;
    }
  }
  return NIFTI_INTENT_NONE;
}

// -----------------------------------------------------------------------------
/// Copy GIFTI data array of data type matching the template argument to vtkDataArray
template <class T>
void CopyDataArrayWithType(vtkDataArray *dst, const giiDataArray *src,
                           vtkIdTypeArray *indices = nullptr)
{
  const int m = src->dims[0];
  const int n = static_cast<int>(src->nvals / static_cast<long long>(m));
  const T *v = reinterpret_cast<const T *>(src->data);
  if (indices) {
    vtkIdType index;
    for (int i = 0; i < m; ++i)
    for (int j = 0; j < n; ++j) {
      dst->SetComponent(i, j, .0);
    }
    if (src->ind_ord == GIFTI_IND_ORD_COL_MAJOR) {
      for (int j = 0; j < n; ++j)
      for (int i = 0; i < m; ++i, ++v) {
        index = indices->GetComponent(i, 0);
        dst->SetComponent(index, j, static_cast<double>(*v));
      }
    } else {
      for (int i = 0; i < m; ++i) {
        index = indices->GetComponent(i, 0);
        for (int j = 0; j < n; ++j, ++v) {
          dst->SetComponent(index, j, static_cast<double>(*v));
        }
      }
    }
  } else {
    if (src->ind_ord == GIFTI_IND_ORD_COL_MAJOR) {
      for (int j = 0; j < n; ++j)
      for (int i = 0; i < m; ++i, ++v) {
        dst->SetComponent(i, j, static_cast<double>(*v));
      }
    } else {
      for (int i = 0; i < m; ++i)
      for (int j = 0; j < n; ++j, ++v) {
        dst->SetComponent(i, j, static_cast<double>(*v));
      }
    }
  }
}

// -----------------------------------------------------------------------------
/// Copy vtkDataArray to GIFTI data array of data type matching the template argument
template <class T>
void CopyDataArrayWithType(giiDataArray *dst, vtkDataArray *src)
{
  const int m = static_cast<int>(src->GetNumberOfTuples());
  const int n = static_cast<int>(src->GetNumberOfComponents());
  T *v = reinterpret_cast<T *>(dst->data);
  if (dst->ind_ord == GIFTI_IND_ORD_COL_MAJOR) {
    for (int j = 0; j < n; ++j)
    for (int i = 0; i < m; ++i, ++v) {
      (*v) = static_cast<T>(src->GetComponent(i, j));
    }
  } else {
    for (int i = 0; i < m; ++i)
    for (int j = 0; j < n; ++j, ++v) {
      (*v) = static_cast<T>(src->GetComponent(i, j));
    }
  }
}

// -----------------------------------------------------------------------------
/// Copy GIFTI data array to vtkDataArray
static void CopyDataArray(vtkDataArray *dst, const giiDataArray *src,
                          vtkIdTypeArray *indices = nullptr)
{
  switch (src->datatype) {
    case NIFTI_TYPE_INT8:    CopyDataArrayWithType<int8_t  >(dst, src, indices); break;
    case NIFTI_TYPE_INT16:   CopyDataArrayWithType<int16_t >(dst, src, indices); break;
    case NIFTI_TYPE_INT32:   CopyDataArrayWithType<int32_t >(dst, src, indices); break;
    case NIFTI_TYPE_INT64:   CopyDataArrayWithType<int64_t >(dst, src, indices); break;
    case NIFTI_TYPE_UINT8:   CopyDataArrayWithType<uint8_t >(dst, src, indices); break;
    case NIFTI_TYPE_UINT16:  CopyDataArrayWithType<uint16_t>(dst, src, indices); break;
    case NIFTI_TYPE_UINT32:  CopyDataArrayWithType<uint32_t>(dst, src, indices); break;
    case NIFTI_TYPE_UINT64:  CopyDataArrayWithType<uint64_t>(dst, src, indices); break;
    case NIFTI_TYPE_FLOAT32: CopyDataArrayWithType<float   >(dst, src, indices); break;
    case NIFTI_TYPE_FLOAT64: CopyDataArrayWithType<double  >(dst, src, indices); break;
    default:
      cerr << "GIFTI data array has unknown/invalid data type: " << src->datatype << endl;
      exit(1);
  }
}

// -----------------------------------------------------------------------------
/// Copy vtkDataArray to GIFTI data array
static void CopyDataArray(giiDataArray *dst, vtkDataArray *src)
{
  switch (dst->datatype) {
    case NIFTI_TYPE_UINT8:   CopyDataArrayWithType<uint8_t >(dst, src); break;
    case NIFTI_TYPE_INT32:   CopyDataArrayWithType<int32_t >(dst, src); break;
    case NIFTI_TYPE_FLOAT32: CopyDataArrayWithType<float   >(dst, src); break;
    default:
      cerr << "GIFTI data array has unknown/invalid data type: " << dst->datatype << endl;
      exit(1);
  }
}

// -----------------------------------------------------------------------------
// vtkInformation keys of standard GIFTI meta data entries
#define GiftiMetaDataKeyMacro(getter, name, type) \
  vtkInformation##type##Key *GiftiMetaData::getter() \
  { \
    static vtkSmartPointer<vtkInformation##type##Key> key \
     = new vtkInformation##type##Key(name, "GiftiMetaData"); \
    return key; \
  }

GiftiMetaDataKeyMacro(DATE, "Date", String);
GiftiMetaDataKeyMacro(USER_NAME, "UserName", String);
GiftiMetaDataKeyMacro(SUBJECT_ID, "SubjectID", String);
GiftiMetaDataKeyMacro(SURFACE_ID, "SurfaceID", String);
GiftiMetaDataKeyMacro(UNIQUE_ID,  "UniqueID", String);
GiftiMetaDataKeyMacro(NAME, "Name", String);
GiftiMetaDataKeyMacro(DESCRIPTION, "Description", String);
GiftiMetaDataKeyMacro(TIME_STEP, "TimeStep", Double);
GiftiMetaDataKeyMacro(DATA_SPACE, "DataSpace", String);
GiftiMetaDataKeyMacro(ANATOMICAL_STRUCTURE_PRIMARY, "AnatomicalStructurePrimary", String);
GiftiMetaDataKeyMacro(ANATOMICAL_STRUCTURE_SECONDARY, "AnatomicalStructureSecondary", String);
GiftiMetaDataKeyMacro(GEOMETRIC_TYPE, "GeometricType", String);
GiftiMetaDataKeyMacro(TOPOLOGICAL_TYPE, "TopologicalType", String);
GiftiMetaDataKeyMacro(INTENT_CODE, "Intent_code", Integer);
GiftiMetaDataKeyMacro(INTENT_P1, "intent_p1", Double);
GiftiMetaDataKeyMacro(INTENT_P2, "intent_p2", Double);
GiftiMetaDataKeyMacro(INTENT_P3, "intent_p3", Double);

// HCP Workbench
GiftiMetaDataKeyMacro(PROGRAM_PROVENANCE, "ProgramProvenance", String);
GiftiMetaDataKeyMacro(PROVENANCE, "Provenance", String);
GiftiMetaDataKeyMacro(PARENT_PROVENANCE, "ParentProvenance", String);
GiftiMetaDataKeyMacro(WORKING_DIRECTORY, "WorkingDirectory", String);

// -----------------------------------------------------------------------------
Array<vtkInformationKey *> GiftiMetaData::KeysForFile()
{
  Array<vtkInformationKey *> keys;
  keys.push_back(DATE());
  keys.push_back(USER_NAME());
  keys.push_back(DESCRIPTION());
  keys.push_back(SUBJECT_ID());
  keys.push_back(UNIQUE_ID());
  keys.push_back(TIME_STEP());
  keys.push_back(PROGRAM_PROVENANCE());
  keys.push_back(PROVENANCE());
  keys.push_back(PARENT_PROVENANCE());
  keys.push_back(WORKING_DIRECTORY());
  keys.push_back(ANATOMICAL_STRUCTURE_PRIMARY()); // HCP Workbench .func.gii
  keys.push_back(ANATOMICAL_STRUCTURE_SECONDARY());
  return keys;
}

// -----------------------------------------------------------------------------
Array<vtkInformationKey *> GiftiMetaData::KeysForDataArray(int intent)
{
  Array<vtkInformationKey *> keys;
  keys.push_back(NAME());
  keys.push_back(DESCRIPTION());
  keys.push_back(UNIQUE_ID());
  keys.push_back(SUBJECT_ID());
  keys.push_back(SURFACE_ID());
  if (intent < 0) {
    keys.push_back(ANATOMICAL_STRUCTURE_PRIMARY());
    keys.push_back(ANATOMICAL_STRUCTURE_SECONDARY());
    keys.push_back(GEOMETRIC_TYPE());
    keys.push_back(TOPOLOGICAL_TYPE());
    keys.push_back(INTENT_CODE());
    keys.push_back(INTENT_P1());
    keys.push_back(INTENT_P2());
    keys.push_back(INTENT_P3());
  } else if (intent == NIFTI_INTENT_POINTSET) {
    keys.push_back(ANATOMICAL_STRUCTURE_PRIMARY());
    keys.push_back(ANATOMICAL_STRUCTURE_SECONDARY());
    keys.push_back(GEOMETRIC_TYPE());
  } else if (intent == NIFTI_INTENT_TRIANGLE) {
    keys.push_back(TOPOLOGICAL_TYPE());
  } else if (NIFTI_FIRST_STATCODE <= intent && intent <= NIFTI_LAST_STATCODE) {
    keys.push_back(INTENT_CODE());
    keys.push_back(INTENT_P1());
    keys.push_back(INTENT_P2());
    keys.push_back(INTENT_P3());
  }
  return keys;
}

// -----------------------------------------------------------------------------
/// Get GIFTI meta data from vtkInformation given a vtkInformationKey as string
string GiftiMetaData::Get(vtkInformation *info, vtkInformationKey *key)
{
  vtkInformationStringKey *skey = vtkInformationStringKey::SafeDownCast(key);
  if (skey) return info->Get(skey);
  vtkInformationDoubleKey *dkey = vtkInformationDoubleKey::SafeDownCast(key);
  if (dkey) return ToString(info->Get(dkey));
  vtkInformationIntegerKey *ikey = vtkInformationIntegerKey::SafeDownCast(key);
  if (dkey) return ToString(info->Get(ikey));
  return string();
}

// -----------------------------------------------------------------------------
/// Copy standard GIFTI meta data to vtkInformation
static void CopyMetaData(vtkInformation *info, const giiMetaData &meta)
{
  for (int i = 0; i < meta.length; ++i) {
    if (strcmp(meta.name[i], GiftiMetaData::DATE()->GetName()) == 0) {
      info->Set(GiftiMetaData::DATE(), meta.value[i]);
    } else if (strcmp(meta.name[i], GiftiMetaData::USER_NAME()->GetName()) == 0) {
      info->Set(GiftiMetaData::USER_NAME(), meta.value[i]);
    } else if (strcmp(meta.name[i], GiftiMetaData::SUBJECT_ID()->GetName()) == 0) {
      info->Set(GiftiMetaData::SUBJECT_ID(), meta.value[i]);
    } else if (strcmp(meta.name[i], GiftiMetaData::SURFACE_ID()->GetName()) == 0) {
      info->Set(GiftiMetaData::SURFACE_ID(), meta.value[i]);
    } else if (strcmp(meta.name[i], GiftiMetaData::UNIQUE_ID()->GetName()) == 0) {
      info->Set(GiftiMetaData::UNIQUE_ID(), meta.value[i]);
    } else if (strcmp(meta.name[i], GiftiMetaData::NAME()->GetName()) == 0) {
      info->Set(GiftiMetaData::NAME(), meta.value[i]);
    } else if (strcmp(meta.name[i], GiftiMetaData::DESCRIPTION()->GetName()) == 0) {
      info->Set(GiftiMetaData::DESCRIPTION(), meta.value[i]);
    } else if (strcmp(meta.name[i], GiftiMetaData::TIME_STEP()->GetName()) == 0) {
      double tr;
      if (FromString(meta.value[i], tr)) {
        info->Set(GiftiMetaData::TIME_STEP(), tr);
      }
    } else if (strcmp(meta.name[i], GiftiMetaData::ANATOMICAL_STRUCTURE_PRIMARY()->GetName()) == 0) {
      info->Set(GiftiMetaData::ANATOMICAL_STRUCTURE_PRIMARY(), meta.value[i]);
    } else if (strcmp(meta.name[i], GiftiMetaData::ANATOMICAL_STRUCTURE_SECONDARY()->GetName()) == 0) {
      info->Set(GiftiMetaData::ANATOMICAL_STRUCTURE_SECONDARY(), meta.value[i]);
    } else if (strcmp(meta.name[i], GiftiMetaData::GEOMETRIC_TYPE()->GetName()) == 0) {
      info->Set(GiftiMetaData::GEOMETRIC_TYPE(), meta.value[i]);
    } else if (strcmp(meta.name[i], GiftiMetaData::TOPOLOGICAL_TYPE()->GetName()) == 0) {
      info->Set(GiftiMetaData::TOPOLOGICAL_TYPE(), meta.value[i]);
    } else if (strcmp(meta.name[i], GiftiMetaData::INTENT_CODE()->GetName()) == 0 ||
               strcmp(meta.name[i], "Intent")     == 0 ||
               strcmp(meta.name[i], "IntentCode") == 0) {
      int intent_code;
      if (FromString(meta.value[i], intent_code)) {
        info->Set(GiftiMetaData::INTENT_CODE(), intent_code);
      }
    } else if (strcmp(meta.name[i], GiftiMetaData::INTENT_P1()->GetName())  == 0 ||
               strcmp(meta.name[i], "Intent_p1") == 0 ||
               strcmp(meta.name[i], "IntentP1")  == 0) {
      double intent_p1;
      if (FromString(meta.value[i], intent_p1)) {
        info->Set(GiftiMetaData::INTENT_P1(), intent_p1);
      }
    } else if (strcmp(meta.name[i], GiftiMetaData::INTENT_P2()->GetName())  == 0 ||
               strcmp(meta.name[i], "Intent_p2") == 0 ||
               strcmp(meta.name[i], "IntentP2")  == 0) {
      double intent_p2;
      if (FromString(meta.value[i], intent_p2)) {
        info->Set(GiftiMetaData::INTENT_P2(), intent_p2);
      }
    } else if (strcmp(meta.name[i], GiftiMetaData::INTENT_P3()->GetName())  == 0 ||
               strcmp(meta.name[i], "Intent_p3") == 0 ||
               strcmp(meta.name[i], "IntentP3")  == 0) {
      double intent_p3;
      if (FromString(meta.value[i], intent_p3)) {
        info->Set(GiftiMetaData::INTENT_P3(), intent_p3);
      }
    }
  }
}

// -----------------------------------------------------------------------------
/// Copy standard GIFTI meta data from vtkInformation if present
#if 0 // unused
static void CopyMetaData(giiMetaData &meta, vtkInformation *info)
{
  vtkNew<vtkInformationIterator> it;
  it->SetInformation(info);
  for (it->InitTraversal(); !it->IsDoneWithTraversal(); it->GoToNextItem()) {
    vtkInformationKey * const key = it->GetCurrentKey();
    if (info->Has(key)) {
      const string value = GiftiMetaData::Get(info, key);
      if (!value.empty()) {
        gifti_add_to_meta(&meta, key->GetName(), value.c_str(), 1);
      }
    }
  }
}
#endif

// -----------------------------------------------------------------------------
/// Copy specified standard GIFTI meta data from vtkInformation if present
static void CopyMetaData(giiMetaData &meta, vtkInformation *info,
                         const Array<vtkInformationKey *> &keys)
{
  for (auto it = keys.begin(); it != keys.end(); ++it) {
    vtkInformationKey * const key = *it;
    if (info->Has(key)) {
      const string value = GiftiMetaData::Get(info, key);
      if (!value.empty()) {
        gifti_add_to_meta(&meta, key->GetName(), value.c_str(), 1);
      }
    }
  }
}

// -----------------------------------------------------------------------------
/// Copy GIFTI point set to vtkPoints
static vtkSmartPointer<vtkPoints>
GetPoints(const gifti_image *gim, vtkInformation *info = nullptr, bool errmsg = false)
{
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  for (int i = 0; i < gim->numDA; ++i) {
    giiDataArray *da = gim->darray[i];
    if (da->intent == NIFTI_INTENT_POINTSET) {
      if (da->datatype != NIFTI_TYPE_FLOAT32) {
        if (errmsg) {
          cerr << "Error: GIFTI coordinates array must have datatype NIFTI_TYPE_FLOAT32!" << endl;
        }
        break;
      }
      if (da->num_dim != 2) {
        if (errmsg) {
          cerr << "Error: GIFTI coordinates array must have 2 dimensions!" << endl;
        }
        break;
      }
      if (da->dims[1] != 3) {
        if (errmsg) {
          cerr << "Error: Second dimension of GIFTI coordinates array must have size 3!" << endl;
        }
        break;
      }
      const int n = da->dims[0];
      points->SetNumberOfPoints(n);
      const float *x = reinterpret_cast<float *>(da->data);
      if (da->ind_ord == GIFTI_IND_ORD_COL_MAJOR) {
        const float *y = x + n, *z = y + n;
        for (int j = 0; j < n; ++j, ++x, ++y, ++z) {
          points->SetPoint(j, static_cast<double>(*x), static_cast<double>(*y), static_cast<double>(*z));
        }
      } else {
        const float *y = x + 1, *z = x + 2;
        for (int j = 0; j < n; ++j, x += 3, y += 3, z += 3) {
          points->SetPoint(j, static_cast<double>(*x), static_cast<double>(*y), static_cast<double>(*z));
        }
      }
      if (info) {
        CopyMetaData(info, da->meta);
        const char *dataspace;
        if (da->numCS > 0) {
          dataspace = da->coordsys[0]->dataspace;
          for (int c = 1; c < da->numCS; ++c) {
            if (strcmp(dataspace, da->coordsys[c]->dataspace) != 0) {
              dataspace = "NIFTI_XFORM_UNKNOWN";
              break;
            }
          }
        } else {
          dataspace = "NIFTI_XFORM_UNKNOWN";
        }
        info->Set(GiftiMetaData::DATA_SPACE(), dataspace);
      }
      break;
    }
  }
  return points;
}

// -----------------------------------------------------------------------------
/// Copy GIFTI topology information to vtkCellArray
static vtkSmartPointer<vtkCellArray>
GetTriangles(const gifti_image *gim, vtkInformation *info = nullptr, bool errmsg = false)
{
  vtkSmartPointer<vtkCellArray> triangles;
  for (int i = 0; i < gim->numDA; ++i) {
    giiDataArray *da = gim->darray[i];
    if (da->intent == NIFTI_INTENT_TRIANGLE) {
      if (da->datatype != NIFTI_TYPE_INT32) {
        if (errmsg) {
          cerr << "Error: GIFTI topology array must have datatype NIFTI_TYPE_INT32!" << endl;
        }
        break;
      }
      if (da->num_dim != 2) {
        if (errmsg) {
          cerr << "Error: GIFTI topology array must have 2 dimensions!" << endl;
        }
        break;
      }
      if (da->dims[1] != 3) {
        if (errmsg) {
          cerr << "Error: Second dimension of GIFTI topology array must have size 3!" << endl;
        }
        break;
      }
      vtkIdType pts[3];
      const int n = da->dims[0];
      triangles = vtkSmartPointer<vtkCellArray>::New();
      triangles->Allocate(3 * n);
      const int *a = reinterpret_cast<int *>(da->data);
      if (da->ind_ord == GIFTI_IND_ORD_COL_MAJOR) {
        const int *b = a + n, *c = b + n;
        for (int j = 0; j < n; ++j, ++a, ++b, ++c) {
          pts[0] = static_cast<vtkIdType>(*a);
          pts[1] = static_cast<vtkIdType>(*b);
          pts[2] = static_cast<vtkIdType>(*c);
          triangles->InsertNextCell(3, pts);
        }
      } else {
        const int *b = a + 1, *c = a + 2;
        for (int j = 0; j < n; ++j, a += 3, b += 3, c += 3) {
          pts[0] = static_cast<vtkIdType>(*a);
          pts[1] = static_cast<vtkIdType>(*b);
          pts[2] = static_cast<vtkIdType>(*c);
          triangles->InsertNextCell(3, pts);
        }
      }
      if (info) CopyMetaData(info, da->meta);
      break;
    }
  }
  return triangles;
}

// -----------------------------------------------------------------------------
/// Convert GIFTI node indices array to vtkDataArray
static vtkSmartPointer<vtkIdTypeArray>
GetNodeIndices(const gifti_image *gim, bool errmsg = false)
{
  vtkSmartPointer<vtkIdTypeArray> indices;
  for (int i = 0; i < gim->numDA; ++i) {
    giiDataArray *da = gim->darray[i];
    if (da->intent == NIFTI_INTENT_NODE_INDEX) {
      if (da->num_dim != 1) {
        if (errmsg) {
          cerr << "Error: GIFTI node indices array must have 1 dimension!" << endl;
        }
        break;
      }
      if (da->dims[0] <= 0) {
        if (errmsg) {
          cerr << "Error: GIFTI node indices array must contain at least one index!" << endl;
        }
        break;
      }
      indices = vtkSmartPointer<vtkIdTypeArray>::New();
      indices->SetNumberOfComponents(1);
      indices->SetNumberOfTuples(da->dims[0]);
      CopyDataArray(indices, da);
    }
  }
  return indices;
}

// -----------------------------------------------------------------------------
/// Convert GIFTI data arrays to vtkDataArray instances of a vtkPointData
static vtkSmartPointer<vtkPointData>
GetPointData(const gifti_image *gim, vtkIdType npoints = 0, vtkIdTypeArray *indices = nullptr, bool errmsg = false)
{
  vtkIdType nindices = 0;
  if (indices) {
    nindices = indices->GetNumberOfTuples();
    if (npoints == 0) {
      cerr << "Error: Number of points cannot be zero when reading sparse GIFTI point data arrays!" << endl;
      exit(1);
    }
    if (nindices > npoints) {
      cerr << "Error: Number of points cannot be less then number of GIFTI node indices!" << endl;
      exit(1);
    }
  }
  bool ok = true;
  vtkSmartPointer<vtkPointData> pd = vtkSmartPointer<vtkPointData>::New();
  for (int i = 0; i < gim->numDA; ++i) {
    giiDataArray *da = gim->darray[i];
    if (da->intent != NIFTI_INTENT_POINTSET &&
        da->intent != NIFTI_INTENT_TRIANGLE &&
        da->intent != NIFTI_INTENT_NODE_INDEX &&
        da->num_dim > 0 && da->dims[0] > 0 && da->nvals > 0) {
      const int ncomp = static_cast<int>(da->nvals / static_cast<long long>(da->dims[0]));
      vtkSmartPointer<vtkDataArray> data;
      data = NewVTKDataArray(GiftiDataTypeToVtk(da->datatype));
      data->SetNumberOfComponents(ncomp);
      if (npoints) {
        if (( indices && static_cast<vtkIdType>(da->dims[0]) != nindices) ||
            (!indices && static_cast<vtkIdType>(da->dims[0]) != npoints)) {
          if (errmsg) {
            cerr << "Error: GIFTI array size does not match point set or node indices array size!" << endl;
          }
          ok = false;
          break;
        }
        data->SetNumberOfTuples(npoints);
      } else {
        data->SetNumberOfTuples(da->dims[0]);
      }
      CopyDataArray(data, da, indices);
      vtkInformation * const info = data->GetInformation();
      CopyMetaData(info, da->meta);
      if (info->Has(GiftiMetaData::NAME())) {
        data->SetName(info->Get(GiftiMetaData::NAME()));
      } else {
        data->SetName(ToString(da->intent).c_str());
      }
      int idx = -1;
      switch (da->intent) {
        case NIFTI_INTENT_SHAPE: {
          if (!pd->GetScalars()) idx = pd->SetScalars(data);
        } break;
        case NIFTI_INTENT_VECTOR: {
          if (ncomp == 3) {
            const string lname = ToLower(data->GetName());
            if (lname == "normals" || lname == "normal") {
              if (!pd->GetNormals()) idx = pd->SetNormals(data);
            } else {
              if (!pd->GetVectors()) idx = pd->SetVectors(data);
            }
          }
        } break;
      }
      if (idx == -1) idx = pd->AddArray(data);
    }
  }
  if (!ok) pd->Initialize();
  return pd;
}

// -----------------------------------------------------------------------------
vtkSmartPointer<vtkPoints> ReadGIFTICoordinates(const char *fname, vtkInformation *info, bool errmsg)
{
  gifti_image *gim = gifti_read_image(fname, 1);
  if (gim == nullptr) {
    if (errmsg) {
      cerr << "Error: Could not read GIFTI file: " << fname << endl;
    }
    return vtkSmartPointer<vtkPoints>::New();
  }
  return GetPoints(gim, info, errmsg);
}

// -----------------------------------------------------------------------------
vtkSmartPointer<vtkCellArray> ReadGIFTITopology(const char *fname, vtkInformation *info, bool errmsg)
{
  gifti_image *gim = gifti_read_image(fname, 1);
  if (gim == nullptr) {
    if (errmsg) {
      cerr << "Error: Could not read GIFTI file: " << fname << endl;
    }
    return nullptr;
  }
  return GetTriangles(gim, info, errmsg);
}

// -----------------------------------------------------------------------------
vtkSmartPointer<vtkPointData> ReadGIFTIPointData(const char *fname, bool errmsg)
{
  gifti_image *gim = gifti_read_image(fname, 1);
  if (gim == nullptr) {
    if (errmsg) {
      cerr << "Error: Could not read GIFTI file: " << fname << endl;
    }
    return nullptr;
  }
  return GetPointData(gim, 0, nullptr, errmsg);
}

// -----------------------------------------------------------------------------
vtkSmartPointer<vtkPolyData> ReadGIFTI(const char *fname, vtkPolyData *surface, bool errmsg)
{
  vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();

  // Read GIFTI
  const int read_data = 1;
  gifti_image *gim = gifti_read_image(fname, read_data);
  if (gim == nullptr) return polydata;

  // Convert geometry and topology arrays including their meta data
  vtkSmartPointer<vtkInformation> geom_info = vtkSmartPointer<vtkInformation>::New();
  vtkSmartPointer<vtkInformation> topo_info = vtkSmartPointer<vtkInformation>::New();
  vtkSmartPointer<vtkPoints>    points = GetPoints   (gim, geom_info,  errmsg);
  vtkSmartPointer<vtkCellArray> polys  = GetTriangles(gim, topo_info,  errmsg);

  // Polygonal dataset requires a point set
  if (points->GetNumberOfPoints() == 0) {
    if (surface && surface->GetNumberOfPoints() > 0) {
      points = surface->GetPoints();
    } else {
      if (errmsg) {
        cerr << "Error: Cannot read GIFTI point data without input point set (e.g., from .coords.gii or .surf.gii file)!" << endl;
      }
      return polydata;
    }
  }
  const vtkIdType npoints = points->GetNumberOfPoints();

  // Check topology information
  if (polys) {
    polys->InitTraversal();
    vtkIdType npts, *pts;
    while (polys->GetNextCell(npts, pts) != 0) {
      if (npts != 3 || pts[0] >= npoints || pts[1] >= npoints || pts[2] >= npoints) {
        if (errmsg) {
          cerr << "Error: GIFTI topology array has invalid point index!" << endl;
        }
        return polydata;
      }
    }
  } else if (surface) {
    polys = surface->GetPolys();
  }

  // Get node indices array
  vtkSmartPointer<vtkIdTypeArray> indices = GetNodeIndices(gim, errmsg);
  if (indices) {
    for (vtkIdType i = 0; i < indices->GetNumberOfTuples(); ++i) {
      const vtkIdType index = indices->GetComponent(i, 0);
      if (index >= npoints) {
        if (errmsg) {
          cerr << "Error: Index of GIFTI node indices array element is out of range!" << endl;
          cerr << "       - Number of points = " << npoints << endl;
          cerr << "       - Node index       = " << index << endl;
        }
        return polydata;
      }
    }
  }

  // Convert possibly sparse point data arrays
  vtkSmartPointer<vtkPointData> pd = GetPointData(gim, npoints, indices, errmsg);

  // Copy file meta data to vtkPolyData information
  vtkInformation * const info = polydata->GetInformation();
  CopyMetaData(info, gim->meta);

  // Free gifti_image instance
  gifti_free_image(gim);
  gim = nullptr;

  // Check number of tuples of point data arrays
  bool ok = true;
  if (pd) {
    for (int i = 0; i < pd->GetNumberOfArrays(); ++i) {
      vtkDataArray *array = pd->GetArray(i);
      if (array->GetNumberOfTuples() != npoints) {
        cerr << "Error: GIFTI array '" << array->GetName()
             << "' at index " << i << " has mismatching size!" << endl;
        cerr << "       - Number of points = " << npoints << endl;
        cerr << "       - Number of tuples = " << array->GetNumberOfTuples() << endl;
        ok = false;
      }
    }
  }

  // Finalize polygonal dataset
  if (ok) {
    // Set geometry, topology, and point data
    polydata->SetPoints(points);
    polydata->SetPolys(polys);
    polydata->GetPointData()->ShallowCopy(pd);
    // In the GIFTI files released by the Human Connectome Project (HCP),
    // the name of point set and triangle list is identical and contains the
    // path of the FreeSurfer surface file. We store this data array name in
    // the vtkPolyData information and in WriteGIFTI we then set the Name meta
    // data of point set and triangle list to the name stored in this map.
    if (geom_info->Has(GiftiMetaData::NAME()) && topo_info->Has(GiftiMetaData::NAME())) {
      if (strcmp(geom_info->Get(GiftiMetaData::NAME()),
                 topo_info->Get(GiftiMetaData::NAME())) == 0) {
        info->CopyEntry(geom_info, GiftiMetaData::NAME());
      }
    } else if (geom_info->Has(GiftiMetaData::NAME())) {
      info->CopyEntry(geom_info, GiftiMetaData::NAME());
    } else if (topo_info->Has(GiftiMetaData::NAME())) {
      info->CopyEntry(topo_info, GiftiMetaData::NAME());
    }
    // Copy meta data of geometry and topology data arrays to vtkPolyData information
    if (geom_info->Has(GiftiMetaData::SUBJECT_ID())) {
      info->CopyEntry(geom_info, GiftiMetaData::SUBJECT_ID());
    }
    if (geom_info->Has(GiftiMetaData::SURFACE_ID())) {
      info->CopyEntry(geom_info, GiftiMetaData::SURFACE_ID());
    }
    if (geom_info->Has(GiftiMetaData::UNIQUE_ID())) {
      info->CopyEntry(geom_info, GiftiMetaData::UNIQUE_ID());
    }
    if (geom_info->Has(GiftiMetaData::DESCRIPTION())) {
      info->CopyEntry(geom_info, GiftiMetaData::DESCRIPTION());
    }
    if (geom_info->Has(GiftiMetaData::DATA_SPACE())) {
      info->CopyEntry(geom_info, GiftiMetaData::DATA_SPACE());
    }
    if (geom_info->Has(GiftiMetaData::ANATOMICAL_STRUCTURE_PRIMARY())) {
      info->CopyEntry(geom_info, GiftiMetaData::ANATOMICAL_STRUCTURE_PRIMARY());
    }
    if (geom_info->Has(GiftiMetaData::ANATOMICAL_STRUCTURE_SECONDARY())) {
      info->CopyEntry(geom_info, GiftiMetaData::ANATOMICAL_STRUCTURE_SECONDARY());
    }
    if (geom_info->Has(GiftiMetaData::GEOMETRIC_TYPE())) {
      info->CopyEntry(geom_info, GiftiMetaData::GEOMETRIC_TYPE());
    }
    if (topo_info->Has(GiftiMetaData::TOPOLOGICAL_TYPE())) {
      info->CopyEntry(topo_info, GiftiMetaData::TOPOLOGICAL_TYPE());
    }
  } else {
    info->Clear();
  }

  return polydata;
}

// -----------------------------------------------------------------------------
static bool AddPoints(gifti_image *gim, vtkPoints *points, vtkInformation *info = nullptr)
{
  if (gifti_add_empty_darray(gim, 1) != 0) return false;
  giiDataArray *da = gim->darray[gim->numDA-1];

  // Set data array attributes
  da->intent     = NIFTI_INTENT_POINTSET;
  da->datatype   = NIFTI_TYPE_FLOAT32;
  da->ind_ord    = GIFTI_IND_ORD_ROW_MAJOR;
  da->num_dim    = 2;
  da->dims[0]    = static_cast<int>(points->GetNumberOfPoints());
  da->dims[1]    = 3;
  #ifdef HAVE_ZLIB
    da->encoding = GIFTI_ENCODING_B64GZ;
  #else
    da->encoding = GIFTI_ENCODING_B64BIN;
  #endif
  da->endian     = gifti_get_this_endian();
  da->ext_fname  = nullptr;
  da->ext_offset = 0;
  da->nvals      = gifti_darray_nvals(da);
  gifti_datatype_sizes(da->datatype, &da->nbyper, nullptr);

  // Allocate memory for point set coordinates
  da->data = calloc(da->nvals * da->nbyper, sizeof(char));
  if (da->data == nullptr) {
    gifti_free_DataArray(da);
    gim->darray[--gim->numDA] = nullptr;
    return false;
  }

  // Copy point set coordinates
  double p[3];
  float *pdata = reinterpret_cast<float *>(da->data);
  for (int i = 0; i < da->dims[0]; ++i) {
    points->GetPoint(i, p);
    (*pdata) = static_cast<float>(p[0]), ++pdata;
    (*pdata) = static_cast<float>(p[1]), ++pdata;
    (*pdata) = static_cast<float>(p[2]), ++pdata;
  }

  // Add coordinate system with identity matrix
  if (gifti_add_empty_CS(da) != 0) {
    gifti_free_DataArray(da);
    gim->darray[--gim->numDA] = nullptr;
    return false;
  }
  giiCoordSystem *cs = da->coordsys[da->numCS-1];
  const char *dataspace;
  if (info && info->Has(GiftiMetaData::DATA_SPACE())) {
    dataspace = info->Get(GiftiMetaData::DATA_SPACE());
  } else {
    dataspace = "NIFTI_XFORM_UNKNOWN";
  }
  cs->dataspace  = gifti_strdup(dataspace);
  cs->xformspace = gifti_strdup(dataspace);
  cs->xform[0][0] = cs->xform[1][1] = cs->xform[2][2] = cs->xform[3][3] = 1.0;

  // Copy meta data from vtkPolyData information
  if (info) {
    CopyMetaData(da->meta, info, GiftiMetaData::KeysForDataArray(da->intent));
  }

  return true;
}

// -----------------------------------------------------------------------------
static bool AddTriangles(gifti_image *gim, vtkCellArray *triangles, vtkInformation *info = nullptr)
{
  if (triangles->GetMaxCellSize() != 3) return false;

  if (gifti_add_empty_darray(gim, 1) != 0) return false;
  giiDataArray *da = gim->darray[gim->numDA-1];

  // Set data array attributes
  da->intent     = NIFTI_INTENT_TRIANGLE;
  da->datatype   = NIFTI_TYPE_INT32;
  da->ind_ord    = GIFTI_IND_ORD_ROW_MAJOR;
  da->num_dim    = 2;
  da->dims[0]    = static_cast<int>(triangles->GetNumberOfCells());
  da->dims[1]    = 3;
  #ifdef HAVE_ZLIB
    da->encoding = GIFTI_ENCODING_B64GZ;
  #else
    da->encoding = GIFTI_ENCODING_B64BIN;
  #endif
  da->endian     = gifti_get_this_endian();
  da->ext_fname  = nullptr;
  da->ext_offset = 0;
  da->nvals      = gifti_darray_nvals(da);
  gifti_datatype_sizes(da->datatype, &da->nbyper, nullptr);

  // Allocate memory for point indices
  da->data = calloc(da->nvals * da->nbyper, sizeof(char));
  if (da->data == nullptr) {
    gifti_free_DataArray(da);
    gim->darray[--gim->numDA] = nullptr;
    return false;
  }

  // Copy triangles
  vtkIdType npts, *pts;
  int *pdata = reinterpret_cast<int *>(da->data);
  triangles->InitTraversal();
  for (int i = 0; i < da->dims[0]; ++i) {
    triangles->GetNextCell(npts, pts);
    if (npts != 3) {
      gifti_free_DataArray(da);
      gim->darray[--gim->numDA] = nullptr;
      return false;
    }
    (*pdata) = static_cast<int>(pts[0]), ++pdata;
    (*pdata) = static_cast<int>(pts[1]), ++pdata;
    (*pdata) = static_cast<int>(pts[2]), ++pdata;
  }

  // Copy meta data from vtkPolyData information
  if (info) {
    CopyMetaData(da->meta, info, GiftiMetaData::KeysForDataArray(da->intent));
  }

  return true;
}

// -----------------------------------------------------------------------------
static bool AddDataArray(gifti_image *gim, vtkDataArray *data, int intent)
{
  if (gifti_add_empty_darray(gim, 1) != 0) return false;
  giiDataArray *da = gim->darray[gim->numDA-1];

  // Set data array attributes
  da->intent     = intent;
  da->datatype   = GiftiDataType(data);
  da->ind_ord    = GIFTI_IND_ORD_ROW_MAJOR;
  da->num_dim    = 1;
  da->dims[0]    = static_cast<int>(data->GetNumberOfTuples());
  if (data->GetNumberOfComponents() > 1) {
    da->num_dim  = 2;
    da->dims[1]  = static_cast<int>(data->GetNumberOfComponents());
  }
  #ifdef HAVE_ZLIB
    da->encoding = GIFTI_ENCODING_B64GZ;
  #else
    da->encoding = GIFTI_ENCODING_B64BIN;
  #endif
  da->endian     = gifti_get_this_endian();
  da->ext_fname  = nullptr;
  da->ext_offset = 0;
  da->nvals      = gifti_darray_nvals(da);
  gifti_datatype_sizes(da->datatype, &da->nbyper, nullptr);

  // Allocate memory
  da->data = malloc(da->nvals * da->nbyper);
  if (da->data == nullptr) {
    gifti_free_DataArray(da);
    gim->darray[--gim->numDA] = nullptr;
    return false;
  }

  // Convert and copy data
  CopyDataArray(da, data);

  // Copy meta data from vtkDataArray information
  CopyMetaData(da->meta, data->GetInformation(), GiftiMetaData::KeysForDataArray(da->intent));
  if (data->GetName()) {
    gifti_add_to_meta(&da->meta, GiftiMetaData::NAME()->GetName(), data->GetName(), 1);
  }

  return true;
}

// -----------------------------------------------------------------------------
static bool AddPointData(gifti_image *gim, vtkPointData *pd, const string &type)
{
  const int numDA = gim->numDA;
  for (int i = 0; i < pd->GetNumberOfArrays(); ++i) {
    vtkSmartPointer<vtkDataArray> array = pd->GetArray(i);
    const int intent = GiftiIntentCode(array, pd->IsArrayAnAttribute(i));

    if (type == ".shape"){
      if (intent != NIFTI_INTENT_SHAPE || array->GetNumberOfComponents() > 1) continue;

      if (GiftiDataType(array) != NIFTI_TYPE_FLOAT32) {
        const vtkIdType n = array->GetNumberOfTuples();
        vtkDataArray * const orig = array;
        array = NewVTKDataArray(VTK_FLOAT);
        array->SetNumberOfComponents(1);
        array->SetNumberOfTuples(n);
        for (vtkIdType i = 0; i < n; ++i) {
          array->SetComponent(i, 0, static_cast<float>(orig->GetComponent(i, 0)));
        }
      }
    }
    
    if (!AddDataArray(gim, array, intent)) {
      for (int j = numDA; j < gim->numDA; ++j) {
        gifti_free_DataArray(gim->darray[j]);
        gim->darray[j] = nullptr;
      }
      gim->numDA = numDA;
      return false;
    }
  }
  return true;
}

// -----------------------------------------------------------------------------
bool WriteGIFTI(const char *fname, vtkPolyData *polydata, FileOption fopt)
{
  // Determine type of GIFTI file from file name extensions
  const string ext  = Extension(fname, EXT_Last);

  string type;
  if (ext == ".gii") {
    const string name = fname;
    type = Extension(name.substr(0, name.length() - ext.length()), EXT_Last);
    if (type != ".coord" &&
        type != ".func" &&
        type != ".label" &&
        type != ".rgba" &&
        type != ".shape" &&
        type != ".surf" &&
        type != ".tensor" &&
        type != ".time" &&
        type != ".topo" &&
        type != ".vector") {
      type.clear();
    } else if (type != ".coord" && type != ".topo" && type != ".surf" && type != ".shape") {
      cerr << "WriteGIFTI: Output file type " << type << ext << " not supported!\n"
              "            Can only write .coord.gii, .topo.gii, .shape.gii, .surf.gii, or generic .gii file.\n"
              "            To write a generic GIFTI file, remove the " << type << " infix before\n"
              "            the .gii file name extension." << endl;
      return false;
    }
  }

  // Allocate new GIFTI structure
  gifti_image *gim = gifti_create_image(0, 0, 0, 0, nullptr, 0);
  if (gim == nullptr) return false;

  // Set extra attributes for XML validation
  gifti_add_to_nvpairs(&gim->ex_atrs, "xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance");
  gifti_add_to_nvpairs(&gim->ex_atrs, "xsi:noNamespaceSchemaLocation", "http://brainvis.wustl.edu/caret6/xml_schemas/GIFTI_Caret.xsd");

  // Remove file meta data that is no longer valid
  vtkSmartPointer<vtkInformation> info = vtkSmartPointer<vtkInformation>::New();
  info->Copy(polydata->GetInformation());
  info->Remove(GiftiMetaData::DATE());
  info->Remove(GiftiMetaData::USER_NAME());
  info->Remove(GiftiMetaData::WORKING_DIRECTORY());
  info->Remove(GiftiMetaData::PROGRAM_PROVENANCE());
  info->Remove(GiftiMetaData::PROVENANCE());
  info->Remove(GiftiMetaData::PARENT_PROVENANCE());

  // Copy file level meta data from vtkPolyData information
  CopyMetaData(gim->meta, info, GiftiMetaData::KeysForFile());

  // Set UserName and Date
  gifti_add_to_meta(&gim->meta, GiftiMetaData::DATE()     ->GetName(), GetDateTime().c_str(), 1);
  gifti_add_to_meta(&gim->meta, GiftiMetaData::USER_NAME()->GetName(), GetUser().c_str(), 1);

  // Add point coordinates
  if (polydata->GetNumberOfPoints() > 0) {
    if (type.empty() || type == ".coord" || type == ".surf") {
      if (!AddPoints(gim, polydata->GetPoints(), info)) {
        gifti_free_image(gim);
        return false;
      }
    }
  }

  // Add triangles
  if (polydata->GetPolys() && polydata->GetPolys()->GetNumberOfCells() > 0) {
    if (type.empty() || type == ".topo" || type == ".surf") {
      if (!AddTriangles(gim, polydata->GetPolys(), info)) {
        gifti_free_image(gim);
        return false;
      }
    }
  }

  // Add point data arrays
  if (type.empty() || (type != ".coord" && type != ".topo" && type != ".surf")) {
    if (!AddPointData(gim, polydata->GetPointData(), type)) {
      gifti_free_image(gim);
      return false;
    }
  }

  // Set encoding of all data arrays
  #ifndef HAVE_ZLIB
    compress = false;
  #endif
  int encoding = (fopt == FO_ASCII      ? GIFTI_ENCODING_ASCII  :
                 (fopt == FO_NoCompress ? GIFTI_ENCODING_B64BIN :
                                          GIFTI_ENCODING_B64GZ));
  for (int i = 0; i < gim->numDA; ++i) {
    gim->darray[i]->encoding = encoding;
  }

  // Write GIFTI file
  const int write_data = 1;
  bool success = (gifti_write_image(gim, fname, write_data) == 0);
  gifti_free_image(gim);

  return success;
}

#endif // MIRTK_IO_WITH_GIFTI


} // namespace mirtk
