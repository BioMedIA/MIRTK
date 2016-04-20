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

#include "mirtk/Common.h"
#include "mirtk/Options.h"

#include "mirtk/Matrix.h"
#include "mirtk/Matlab.h"
#include "mirtk/PointSetIO.h"
#include "mirtk/PointSetUtils.h"

#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkCellArray.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkTriangleFilter.h"

using namespace mirtk;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
/// Print help screen
void PrintHelp(const char *name)
{
  cout << endl;
  cout << "Usage: " << name << " <input> <mat> [options]" << endl;
  cout << endl;
  cout << "Description:" << endl;
  cout << "  Saves the points and faces of a vtkPolyData to a MAT file." << endl;
  cout << endl;
  cout << "Arguments:" << endl;
  cout << "  input   File name of input dataset." << endl;
  cout << "  mat     File name of output MAT file." << endl;
  cout << endl;
  cout << "Optional arguments:" << endl;
  cout << "  -points [<var>]             Name of points variable. (default: X)" << endl;
  cout << "  -nopoints                   Exclude points from MAT file." << endl;
  cout << "  -faces [<var>]              Name of faces variable. (default: F)" << endl;
  cout << "  -nofaces                    Exclude faces from MAT file." << endl;
  cout << "  -pointdata <name> [<var>]   Include named point data. (default: none)" << endl;
  cout << "  -celldata  <name> [<var>]   Include named cell  data. (default: none)" << endl;
  cout << endl;
}

// =============================================================================
// Auxiliaries
// =============================================================================

// -----------------------------------------------------------------------------
/// Make valid MATLAB variable name from given string
string ValidVariableName(const char *str)
{
  string name(str);
  for (size_t i = 0; i < name.length(); ++i) {
    if (!isalnum(name[i])) name[i] = '_';
  }
  if (isdigit(name.front())) name = string("_") + name;
  return name;
}

// -----------------------------------------------------------------------------
inline mxArray *MatrixToMxArray(const Matrix &m)
{
  mxArray *mx = mxCreateDoubleMatrix(m.Rows(), m.Cols(), mxREAL);
  memcpy(mxGetPr(mx), m.RawPointer(), m.NumberOfElements() * sizeof(double));
  return mx;
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  // Positional arguments
  EXPECTS_POSARGS(2);

  const char *vtk_name = POSARG(1);
  const char *mat_name = POSARG(2);

  // Optional arguments
  const char *x_name = "X";
  const char *f_name = "F";
  Array<Pair<const char *, const char *> > pd_names, cd_names;

  for (ALL_OPTIONS) {
    if (OPTION("-points")) {
      if (HAS_ARGUMENT) x_name = ARGUMENT;
      else              x_name = "X";
    }
    else if (OPTION("-nopoints")) x_name = NULL;
    else if (OPTION("-faces")) {
      if (HAS_ARGUMENT) f_name = ARGUMENT;
      else              f_name = "F";
    }
    else if (OPTION("-nofaces"))  f_name = NULL;
    else if (OPTION("-pointdata")) {
      const char *name = ARGUMENT;
      const char *var  = (HAS_ARGUMENT ? ARGUMENT : name);
      pd_names.push_back(MakePair(name, var));
    }
    else if (OPTION("-celldata")) {
      const char *name = ARGUMENT;
      const char *var  = (HAS_ARGUMENT ? ARGUMENT : name);
      cd_names.push_back(MakePair(name, var));
    }
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  // Initialize MCR
  Matlab::Initialize();

  // Read input dataset
  if (verbose > 1) cout << "Reading dataset from " << vtk_name;
  vtkSmartPointer<vtkPolyData> polydata = ReadPolyData(vtk_name);
  if (verbose > 1) cout << " done" << endl;

  // Open output MAT file
  MATFile *fp = matOpen(mat_name, "w");
  if (fp == NULL) FatalError("Failed to open MAT file " << mat_name);

  // Write points
  if (x_name) {
    double x[3];
    Matrix X(polydata->GetNumberOfPoints(), 3);
    for (vtkIdType i = 0; i < polydata->GetNumberOfPoints(); ++i) {
      polydata->GetPoint(i, x);
      X(i, 0) = x[0], X(i, 1) = x[1], X(i, 2) = x[2];
    }
    mxArray *var = MatrixToMxArray(X);
    if (matPutVariable(fp, x_name, var) != 0) {
      mxDestroyArray(var);
      matClose(fp);
      FatalError("Failed to write points to MAT file!");
    }
    mxDestroyArray(var);
  }

  // Write polygons as triangle faces
  if (f_name) {
    vtkSmartPointer<vtkTriangleFilter> triangulate;
    triangulate = vtkSmartPointer<vtkTriangleFilter>::New();
    SetVTKInput(triangulate, polydata);
    triangulate->Update();
    polydata = triangulate->GetOutput();
    if (polydata->GetPolys()->GetMaxCellSize() != 3) {
      FatalError("Input dataset must be triangulated surface mesh!");
    }
    Matrix F(polydata->GetNumberOfCells(),  3);
    vtkIdType i = 0;
    vtkCellArray *cells = polydata->GetPolys();
    cells->InitTraversal();
    vtkIdType npts, *pts;
    while (cells->GetNextCell(npts, pts)) {
      F(i, 0) = pts[0]+1, F(i, 1) = pts[1]+1, F(i, 2) = pts[2]+1;
      ++i;
    }
    mxArray *var = MatrixToMxArray(F);
    if (matPutVariable(fp, x_name, var) != 0) {
      mxDestroyArray(var);
      matClose(fp);
      FatalError("Failed to write faces to MAT file!");
    }
    mxDestroyArray(var);
  }

  // Write point data
  for (size_t i = 0; i < pd_names.size(); ++i) {
    string var_name = ValidVariableName(pd_names[i].second);
    vtkDataArray *array = GetArrayByCaseInsensitiveName(polydata->GetPointData(), pd_names[i].first);
    if (!array) {
      matClose(fp);
      FatalError("Input polydata has no point data array named " << pd_names[i].first);
    }
    Matrix S(array->GetNumberOfTuples(), array->GetNumberOfComponents());
    for (vtkIdType i = 0; i < array->GetNumberOfTuples(); ++i) {
      for (int j = 0; j < array->GetNumberOfComponents(); ++j) {
        S(i, j) = array->GetComponent(i, j);
      }
    }
    mxArray *var = MatrixToMxArray(S);
    if (matPutVariable(fp, var_name.c_str(), var) != 0) {
      mxDestroyArray(var);
      matClose(fp);
      FatalError("Failed to write point data array " << pd_names[i].first << " to MAT file!");
    }
    mxDestroyArray(var);
  }

  // Write cell data
  for (size_t i = 0; i < cd_names.size(); ++i) {
    string var_name = ValidVariableName(cd_names[i].second);
    vtkDataArray *array = GetArrayByCaseInsensitiveName(polydata->GetCellData(), cd_names[i].first);
    if (!array) {
      matClose(fp);
      FatalError("Input polydata has no cell data array named " << cd_names[i].first);
    }
    Matrix S(array->GetNumberOfTuples(), array->GetNumberOfComponents());
    for (vtkIdType i = 0; i < array->GetNumberOfTuples(); ++i) {
      for (int j = 0; j < array->GetNumberOfComponents(); ++j) {
        S(i, j) = array->GetComponent(i, j);
      }
    }
    mxArray *var = MatrixToMxArray(S);
    if (matPutVariable(fp, var_name.c_str(), var) != 0) {
      mxDestroyArray(var);
      matClose(fp);
      FatalError("Failed to write cell data array " << cd_names[i].first << " to MAT file!");
    }
    mxDestroyArray(var);
  }

  // Close MAT file
  return (matClose(fp) == 0 ? 0 : 1);
}
