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

#include "mirtk/Common.h"
#include "mirtk/Options.h"

#include "mirtk/EdgeConnectivity.h"
#include "mirtk/PointSetIO.h"

#include "vtkSmartPointer.h"
#include "vtkPointSet.h"
#include "vtkPointData.h"
#include "vtkDataArray.h"

using namespace mirtk;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  cout << "\n";
  cout << "Usage: " << name << " <input> <output> [options]\n";
  cout << "\n";
  cout << "Description:\n";
  cout << "  Erodes scalar point data of an input point set by replacing a point's\n";
  cout << "  value by the minimum of the values of its neighboring points. When the\n";
  cout << "  input data array has more than one component, each component is processed\n";
  cout << "  separately.\n";
  cout << "\n";
  cout << "Arguments:\n";
  cout << "  input    Input  point set.\n";
  cout << "  output   Output point set with dilated scalars.\n";
  cout << "\n";
  cout << "Optional arguments:\n";
  cout << "  -a, -array, -scalars <name>\n";
  cout << "      Name of point data array. (default: active scalars)\n";
  cout << "  -c, -connectivity, -neighborhood <int>\n";
  cout << "      Maximum node edge-connectivity of neighbors. (default: 1)\n";
  cout << "  -r, -radius <float>\n";
  cout << "      Maximum point distance of neighboring points. Used instead of :option:`-c` when positive. (default: 0)\n";
  cout << "  -n, -iterations, -iter <n>\n";
  cout << "      Number of iterations. (default: 1)\n";
  cout << "\n";
  cout << "Output format options:\n";
  cout << "  -ascii, -nobinary" << endl;
  cout << "      Write legacy VTK in ASCII format. (default: input)" << endl;
  cout << "  -binary, -noascii" << endl;
  cout << "      Write legacy VTK in binary format. (default: input)" << endl;
  cout << "  -[no]compress" << endl;
  cout << "      Write XML VTK file with or without compression. (default: on)" << endl;
  PrintStandardOptions(cout);
  cout << "\n";
  cout.flush();
}

// =============================================================================
// Auxiliaries
// =============================================================================

// -----------------------------------------------------------------------------
struct ErodeScalars
{
  vtkDataArray           *_Input;
  vtkDataArray           *_Output;
  const EdgeConnectivity *_Neighbors;

  void operator ()(const blocked_range<int> &ptIds) const
  {
    int        nbrPts;
    const int *nbrIds;
    double     value;

    for (int ptId = ptIds.begin(); ptId != ptIds.end(); ++ptId) {
      for (int j = 0; j < _Input->GetNumberOfComponents(); ++j) {
        value = _Input->GetComponent(ptId, j);
        _Neighbors->GetConnectedPoints(ptId, nbrPts, nbrIds);
        for (int i = 0; i < nbrPts; ++i) {
          value = min(value, _Input->GetComponent(nbrIds[i], j));
        }
        _Output->SetComponent(ptId, j, value);
      }
    }
  }
};

// -----------------------------------------------------------------------------
/// Erode point data
///
/// \todo Implement this as MeshFilter as part of the PointSet module.
vtkSmartPointer<vtkPointSet> Erode(vtkSmartPointer<vtkPointSet> input,
                                   const EdgeConnectivity &neighbors,
                                   vtkSmartPointer<vtkDataArray> arr, int niter = 1)
{
  vtkSmartPointer<vtkPointSet> output;
  output.TakeReference(input->NewInstance());
  output->ShallowCopy(input);

  int attr = -1;
  vtkPointData * const pd = output->GetPointData();
  for (int i = 0; i < pd->GetNumberOfArrays(); ++i) {
    if (pd->GetArray(i) == arr) {
      attr = pd->IsArrayAnAttribute(i);
      pd->RemoveArray(i);
      break;
    }
  }

  vtkSmartPointer<vtkDataArray> res;
  res.TakeReference(arr->NewInstance());
  res->DeepCopy(arr); // in particular size, name, and component names if any
  const int idx = pd->AddArray(res);
  if (attr == -1) pd->SetActiveAttribute(idx, attr);

  for (int iter = 0; iter < niter; ++iter) {
    if (iter > 0) {
      if (iter == 1) arr.TakeReference(res->NewInstance());
      arr->DeepCopy(res);
    }
    ErodeScalars body;
    body._Input     = arr;
    body._Output    = res;
    body._Neighbors = &neighbors;
    parallel_for(blocked_range<int>(0, static_cast<int>(input->GetNumberOfPoints())), body);
  }

  return output;
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  EXPECTS_POSARGS(2);

  const char *input_name  = POSARG(1);
  const char *output_name = POSARG(2);
  const char *array_name  = nullptr;

  int iterations    = 1;
  int connectivity  = 1;
  double     radius = 0.;
  FileOption fopt   = FO_Default;

  for (ALL_OPTIONS) {
    if (OPTION("-a") || OPTION("-array") || OPTION("-scalars")) {
      array_name = ARGUMENT;
    }
    else if (OPTION("-n") || OPTION("-iterations") || OPTION("-iter")) {
      PARSE_ARGUMENT(iterations);
    }
    else if (OPTION("-c") || OPTION("-connectivity") || OPTION("-neighborhood") || OPTION("-neighbourhood")) {
      PARSE_ARGUMENT(connectivity);
    }
    else if (OPTION("-r") || OPTION("-radius")) {
      PARSE_ARGUMENT(radius);
    }
    else HANDLE_POINTSETIO_OPTION(fopt);
    else HANDLE_STANDARD_OR_UNKNOWN_OPTION();
  }

  if (verbose) cout << "Reading " << input_name << "...", cout.flush();
  vtkSmartPointer<vtkPointSet> pset = ReadPointSet(input_name);
  if (verbose) cout << " done" << endl;

  vtkDataArray *arr = nullptr;
  if (array_name) {
    arr = pset->GetPointData()->GetArray(array_name);
    if (arr == nullptr) {
      FatalError("Input point set has not point data array named " << array_name << "!");
    }
  } else {
    arr = pset->GetPointData()->GetScalars();
    if (arr == nullptr) {
      FatalError("Input point set has no active scalars, use -array, -scalars option!");
    }
  }

  EdgeConnectivity neighbors;
  if (radius > 0.) {
    neighbors.Initialize(pset, radius);
  } else {
    neighbors.Initialize(pset, connectivity);
  }

  if (verbose) cout << "Eroding with c=" << neighbors.Maximum() << "...", cout.flush();
  pset = Erode(pset, neighbors, arr, iterations);
  arr  = nullptr;
  if (verbose) cout << " done" << endl;

  if (!WritePointSet(output_name, pset, fopt)) {
    FatalError("Failed to write result to " << output_name << "!");
  }

  return 0;
}
